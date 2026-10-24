# Prepare the bioassay data for modelling: read the five databases, put them on
# a common format, remove copies of the same bioassay, and save the An. gambiae
# complex records to data/clean/all_gambiae_complex_data.RDS. The supplement
# (doc/supp.qmd, "Bioassay data processing") describes the method;
# outputs/review/bioassay_attrition.csv counts the records each step removes,
# outputs/review/bioassay_links.csv lists every copy removed, and
# outputs/review/bioassay_review_pairs.csv the pairs flagged for review.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_prep_functions.R")

# ---- read the databases ----

# The 20261005 MTM discriminating concentration download fills missing values
# of MORTALITY_ADJUSTED with 0: it has no missing values, and all 306 records
# missing in the 20251031 download are exactly 0, 87 of them labelled
# susceptible. The papers checked report nonzero mortality for these (#27). So
# we drop them, by ID against the 20251031 download; the 232 zeros already in
# that download are kept.
mtm_filled_mortality_ids <- read_xlsx(
  "data/raw/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20251031.xlsx",
  sheet = "Data", col_types = "text") %>%
  filter(is.na(MORTALITY_ADJUSTED)) %>%
  pull(ID)

mtm_columns <- c(record_id = "ID", citation = "CITATION",
                 link = "CITATION_URL", site = "SITE_NAME",
                 country_name = "COUNTRY_NAME", latitude = "LATITUDE",
                 longitude = "LONGITUDE", year_start = "YEAR_START",
                 generation = "STAGE_ORIGIN",
                 insecticide_type = "INSECTICIDE_TYPE",
                 concentration = "INSECTICIDE_CONC", test_type = "TEST_TYPE",
                 mosquito_number = "MOSQUITO_NUMBER",
                 mortality_adjusted = "MORTALITY_ADJUSTED")

# An IR Mapper analysis-ready extract, as a source (below), given the column
# naming its species. The extracts are in Windows-1252, with "NR" and "NF" for
# missing values. Records pooled over several sites ("multipoint") are
# dropped.
ir_mapper_extract <- function(source, file, species_column) {
  force(file)
  list(
    source = source,
    database = "IR Mapper extracts",
    read = function() {
      read_csv(file, col_types = cols(.default = "c"),
               na = c("", "NA", "NR", "NF"),
               locale = locale(encoding = "windows-1252"))
    },
    keep = function(d) filter(d, `Site type` == "point"),
    columns = c(ir_mapper_extract_columns, species = species_column)
  )
}
ir_mapper_extract_columns <- c(record_id = "Bioassay ID",
                               citation = "Source 1 citation",
                               site = "Site name", country_name = "Country",
                               latitude = "Latitude", longitude = "Longitude",
                               year_start = "Start year",
                               year_end = "End year",
                               month = "Start month",
                               generation = "Generation",
                               insecticide_type = "Insecticide tested",
                               concentration = "Concentration (%)",
                               test_type = "Test protocol",
                               mosquito_number = "No. mosquitoes tested",
                               mortality_adjusted = "Percent mortality")

# Vector Atlas cites a study by DOI, so we build a citation from its author and
# publication year columns. These sometimes vary within a study (a fill-down
# error gives one study a year per record, 2014 to 2031), so we take the most
# common, and the earliest of equally common years.
read_vector_atlas <- function() {
  va <- read_csv("data/raw/VA_FULL_DATA_20250721.csv",
                 col_types = cols(.default = "c")) %>%
    filter(insecticide_resistance_data == "phenotypic")
  study_citations <- va %>%
    count(source_id, author, publication_year, article_title,
          journal_title) %>%
    arrange(source_id, desc(n), publication_year) %>%
    distinct(source_id, .keep_all = TRUE) %>%
    transmute(source_id,
              citation = paste0(author, " (", publication_year, ") ",
                                article_title, ". ", journal_title, "."))
  va %>%
    left_join(study_citations, by = "source_id")
}

# each source: its database, how to read it, its own filters, and its columns.
# Sources are listed in order of preference, for choosing which copy of a
# bioassay to keep.
sources <- list(
  list(
    source = "va_new",
    database = "Vector Atlas",
    read = read_vector_atlas,
    # records at a point, with a positive number of mosquitoes if recorded
    keep = function(d) {
      d %>%
        filter(area_type == "point",
               is.na(mosquitoes_tested_n) |
                 as.numeric(mosquitoes_tested_n) != 0)
    },
    columns = c(citation = "citation", link = "citation_doi", site = "site",
                country_name = "country", latitude = "latitude_1",
                longitude = "longitude_1", year_start = "year_start",
                year_end = "year_end",
                month = "month_start", species = "species",
                generation = "generation",
                insecticide_type = "insecticide_tested",
                concentration = "concentration_percent",
                test_type = "test_protocol",
                mosquito_number = "mosquitoes_tested_n",
                mortality_adjusted = "percent_mortality")
  ),
  list(
    source = "mtm_dis",
    database = "MTM",
    read = function() {
      read_xlsx(
        "data/raw/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20261005.xlsx",
        sheet = "Data", col_types = "text")
    },
    # drop the filled mortalities (above) and a record in the Americas
    keep = function(d) {
      d %>%
        filter(!ID %in% mtm_filled_mortality_ids,
               as.numeric(LONGITUDE) > -50)
    },
    columns = c(mtm_columns, species = "SPECIES")
  ),
  list(
    source = "mtm_int",
    database = "MTM",
    read = function() {
      read_xlsx(
        "data/raw/MTM_INTENSITY_CONCENTRATION_BIOASSAY_20261005.xlsx",
        sheet = "Data", col_types = "text")
    },
    keep = function(d) filter(d, as.numeric(LONGITUDE) > -50),
    columns = c(mtm_columns, species = "VECTOR_SPECIES")
  ),
  ir_mapper_extract("mapper_dis_species",
                    "data/raw/2_standard-WHO-susc-test_species.csv",
                    "Species"),
  ir_mapper_extract("mapper_int_complex",
                    "data/raw/5_intensity-bioassays_complex-subgroup.csv",
                    "Complex/Subgroup"),
  ir_mapper_extract("mapper_dis_complex",
                    "data/raw/1_standard-WHO-susc-test_complex-subgroup.csv",
                    "Complex/Subgroup"),
  list(
    source = "mapper_source_2022",
    database = "IR Mapper 2022",
    read = function() {
      read_xlsx("data/raw/IR Mapper data August 2022.xlsx", col_types = "text")
    },
    keep = identity,
    columns = c(citation = "Reference_Name", link = "URL", site = "Locality",
                country_name = "Country", latitude = "Latitude",
                longitude = "Longitude",
                year_start = "Mosquito_collection_Year_start",
                year_end = "Mosquito_collection_year_end",
                species = "Vector_species", insecticide_type = "Chemical_type",
                concentration = "Insecticide_dosage",
                test_type = "IR_Test_Method",
                mosquito_number = "IR_Test_NumExposed",
                mortality_adjusted = "IR_Test_mortality")
  ),
  list(
    source = "mapper_source_2024",
    database = "IR Mapper 2024",
    read = function() {
      read_xlsx("data/raw/IR_Mapper_Anopheles_data_2021_Sept_2024.xlsx",
                col_types = "text")
    },
    keep = identity,
    columns = c(citation = "reference_Name", link = "url", site = "locality",
                country_name = "country", latitude = "latitude",
                longitude = "longitude", year_start = "collection_Year_Start",
                year_end = "collection_Year_End",
                species = "vector_Species", insecticide_type = "chemical_Type",
                concentration = "insecticide_Dosage",
                test_type = "iR_Test_Method",
                mosquito_number = "iR_Test_NumExposed/Tested",
                mortality_adjusted = "iR_Test_Mortality")
  )
)

bioassays <- lapply(sources, read_bioassay_source) %>%
  bind_rows() %>%
  mutate(
    record = row_number(),
    source_rank = match(source, map_chr(sources, "source")),
    country_name = standardise_country(country_name),
    insecticide_type = standardise_insecticide(insecticide_type),
    insecticide_class = insecticide_class_of(insecticide_type),
    species = standardise_species(species),
    # the paper a record cites, as its first author and year of publication
    # (see the study labels, below)
    paper = coalesce(citation_author_year(citation), citation),
    # the years a test ran: one year where no end year is given (as in MTM),
    # or where the end year is before the start (3 records)
    year_end = if_else(is.na(year_end) | year_end < year_start, year_start,
                       year_end),
    mosquito_number = round(mosquito_number),
    n_recorded = !is.na(mosquito_number)
  )

# ---- select records ----

# Each step keeps the records meeting a condition, and counts by source the
# records it removes
attrition <- count(bioassays, source, name = "read")
keep_records <- function(records, condition, step) {
  kept <- filter(records, {{ condition }})
  removed <- anti_join(records, kept, by = "record") %>%
    count(source, name = step)
  attrition <<- attrition %>%
    left_join(removed, by = "source") %>%
    mutate(across(all_of(step), ~ coalesce(.x, 0L)))
  kept
}

bioassays <- bioassays %>%
  keep_records(!is.na(country_name), "not in Africa") %>%
  keep_records(grepl("who", tolower(test_type)) &
                 !grepl("bottle|larval|synergist", tolower(test_type)),
               "not a WHO tube test") %>%
  keep_records(!is.na(latitude) & !is.na(longitude) & !is.na(year_start) &
                 !is.na(concentration) & !is.na(mortality_adjusted) &
                 !is.na(insecticide_type),
               "missing a field") %>%
  keep_records(species %in% gambiae_complex_species,
               "not An. gambiae complex")

# map pixels, which the sums of species-level records (below) are matched on
mask <- rast("data/clean/raster_mask.tif")
bioassays <- bioassays %>%
  mutate(cell = terra::cellFromXY(mask, cbind(longitude, latitude))) %>%
  keep_records(!is.na(cell), "outside the map")

# ---- the number of mosquitoes, and proportions ----

# Where the number of mosquitoes isn't recorded we take the number closest to
# 100 that gives a whole number dead, or 100 if none between 50 and 200 does
# (infer_sample_size()). Copies are found from the recorded numbers only
# (n_recorded), below.
bioassays <- bioassays %>%
  mutate(
    mosquito_number = replace(
      mosquito_number, !n_recorded,
      infer_sample_size(mortality_adjusted[!n_recorded])),
    died = round(mosquito_number * mortality_adjusted / 100)
  ) %>%
  # a positive mortality that gives no dead mosquito is an error: a proportion
  # read as a percentage, e.g. 1.0 for 100% (#60), or a number dead, e.g. 3
  # of 3 or 4
  keep_records(!(died == 0 & mortality_adjusted > 0),
               "positive mortality below one mosquito")

# ---- copies of the same bioassay ----

# A download sometimes holds a paper more than once (imported twice), or
# repeats its rows (one row copied down for many tests): at least half of the
# paper's records in the download are repeated, identical in every column but
# their ID (row_hash). One of each is kept. Records of one paper in one
# download that are repeated otherwise are taken as the paper's separate tests:
# a paper reporting many tests at one site can give two the same result. The
# supplement gives the shares of repeated records and the checks by hand.
bioassays <- bioassays %>%
  group_by(source, citation) %>%
  add_count(row_hash, name = "row_copies") %>%
  mutate(repeated = mean(row_copies > 1) >= 0.5 & duplicated(row_hash)) %>%
  ungroup() %>%
  select(-row_copies) %>%
  keep_records(!repeated, "paper repeated in a download")

# The databases hold many of the same bioassays, and some hold a bioassay more
# than once. Copies have the same country, years, insecticide, concentration
# and result (agreeing_pairs()), but not always the same coordinates: the
# databases georeference sites differently, and copies in different databases
# are often tens of km apart. So copies are looked for within max_km of each
# other in different databases, and within max_km_within in the same
# database.
#
# The distances are a choice. Matching on the percentage alone, where a number
# of mosquitoes is missing, is less specific, so its distance is shorter.
# Within a database, copies are mostly at the same coordinates.
max_km <- c(count = 100, percent = 50, misscaled = 100)
max_km_within <- 10

# Databases sometimes give one test different years. Where two databases'
# records of one paper match exactly in years that don't overlap
# (cross_year_pairs()), they take the years of the better record
# (cross_year_dates()), and so are found as copies below. The rules come from
# checks by hand in October 2026: the year of publication was given where a
# paper gives none, a database lumped a paper's tests into one year, and
# otherwise the IR Mapper extracts (Moyes et al. 2019) were most often right.
year_changes <- cross_year_dates(
  bioassays,
  database_order = c("IR Mapper extracts", "Vector Atlas", "MTM",
                     "IR Mapper 2022", "IR Mapper 2024"),
  max_km = max_km[["count"]])
year_corrections <- bioassays %>%
  inner_join(year_changes, by = "record") %>%
  left_join(select(bioassays, from = record, from_source = source,
                   from_record_id = record_id),
            by = "from") %>%
  select(source, record_id, citation, site, insecticide_type,
         mosquito_number, died, recorded_start = year_start,
         recorded_end = year_end, from_start, from_end, from_source,
         from_record_id, rule)
bioassays <- bioassays %>%
  left_join(select(year_changes, record, from_start, from_end),
            by = "record") %>%
  mutate(year_start = coalesce(from_start, year_start),
         year_end = coalesce(from_end, year_end)) %>%
  select(-from_start, -from_end)

# Which copy to keep: one with a recorded number of mosquitoes, then from the
# first source in the list above
bioassays <- bioassays %>%
  arrange(!n_recorded, source_rank, record) %>%
  mutate(preference = row_number(),
         record_key = paste(source, record_id)) %>%
  arrange(record)

# A reviewer records whether a pair of records are copies ("copy") or
# separate tests ("different") in review/bioassay_duplicates.csv, by pair_id.
# A pair decided "different" is never taken as copies; any other decision
# (blank, "unsure") leaves the rules below to decide. The pairs flagged for
# review are written to outputs/review/bioassay_review_pairs.csv (below).
decisions <- read_csv("review/bioassay_duplicates.csv",
                      col_types = cols(.default = "c")) %>%
  filter(decision %in% c("copy", "different")) %>%
  select(pair_id, decision)

record_info <- bioassays %>%
  select(record, record_key, database, source, citation, site, month,
         generation,
         preference)
pairs <- agreeing_pairs(bioassays, max_km, max_km_within) %>%
  left_join(record_info, by = c("a" = "record")) %>%
  left_join(record_info, by = c("b" = "record"), suffix = c("_a", "_b")) %>%
  mutate(pair_id = pair_id(record_key_a, record_key_b),
         name_match = same_site_name(site_a, site_b),
         # the two citations in a fixed order, a missing one as ""
         first = pmin(coalesce(citation_a, ""), coalesce(citation_b, "")),
         second = pmax(coalesce(citation_a, ""), coalesce(citation_b, ""))) %>%
  left_join(select(decisions, pair_id, decision), by = "pair_id")

# Within a database. Two records with the same citation are one test in two
# downloads of the database (e.g. the IR Mapper extracts' species and complex
# files) if they are at the same coordinates in the same month and generation,
# and else separate tests. Within one download, a paper's records are separate
# tests (bar a paper repeated in a download, above): we trust the authors and
# the databases' extractors. Records with different citations are copies if
# the two citations share at least two results between 0% and 100% (the same
# data published twice, e.g. in a paper and a report), and else are flagged
# for review.
same_text <- function(x, y) {
  coalesce(as.character(x), "") == coalesce(as.character(y), "")
}
shared_results <- pairs %>%
  filter(within, !extreme, first != second) %>%
  count(database_a, first, second, name = "shared")
within_pairs <- pairs %>%
  filter(within) %>%
  left_join(shared_results, by = c("database_a", "first", "second")) %>%
  mutate(
    shared = coalesce(shared, 0L),
    rule = case_when(
      first == second & source_a != source_b & km < 0.01 &
        same_text(month_a, month_b) &
        same_text(generation_a, generation_b) ~ "in two downloads",
      first == second ~ "separate tests",
      shared >= 2 ~ "published twice",
      .default = "review: different citations"
    )
  )

# The records to drop from pairs of copies, each with the preferred record
# kept in its place. Copies are paired one-to-one (one_to_one()) between the
# two citations, or downloads (`between`), as one paper can hold separate tests
# with one result, such as 100%.
drop_copies <- function(pairs, between = c("citation_a", "citation_b")) {
  if (nrow(pairs) == 0) {
    return(tibble(dropped = integer(), kept = integer(), rule = character(),
                  decision = character(), km = numeric()))
  }
  pairs %>%
    one_to_one(between = between) %>%
    mutate(dropped = if_else(preference_a > preference_b, a, b),
           kept = if_else(dropped == a, b, a)) %>%
    select(dropped, kept, rule, decision, km)
}
# pairs under these rules, unless a reviewer decided they are different tests
is_copy <- function(pairs, rules) {
  pairs$rule %in% rules & coalesce(pairs$decision == "copy", TRUE)
}
not_dropped <- function(pairs, dropped) {
  filter(pairs, !a %in% dropped$dropped, !b %in% dropped$dropped)
}
dropped_entered <- drop_copies(
  filter(within_pairs, is_copy(within_pairs, "in two downloads")),
  between = c("source_a", "source_b"))
within_left <- not_dropped(within_pairs, dropped_entered)
dropped_published <- drop_copies(
  filter(within_left, is_copy(within_left, "published twice")))
within_left <- not_dropped(within_left, dropped_published)

# Then, among the records left: records with different citations flagged for
# review, which are copies unless a reviewer found them separate tests. Every
# flagged pair needs a decision, so the prep warns of any without one.
dropped_reviewed <- drop_copies(
  filter(within_left, is_copy(within_left, "review: different citations")))
undecided <- sum(startsWith(within_left$rule, "review") &
                   is.na(within_left$decision))
if (undecided > 0) {
  warning(undecided, " pairs flagged for review have no decision; see ",
          "outputs/review/bioassay_review_pairs.csv")
}
dropped_within <- c(dropped_entered$dropped, dropped_published$dropped,
                    dropped_reviewed$dropped)

# Between databases, pairs are chosen one-to-one (one_to_one()). In a pair on
# different scales, the record on the wrong scale is dropped, unless it also
# has a copy on its own scale in another database (then the match was a
# coincidence); else the preferred record is kept.
has_copy <- pairs %>%
  filter(!within, agreement != "misscaled") %>%
  select(a, b) %>%
  unlist() %>%
  unique()
cross_pairs <- pairs %>%
  filter(!within, !a %in% dropped_within, !b %in% dropped_within,
         agreement != "misscaled" | !a %in% has_copy) %>%
  one_to_one() %>%
  mutate(dropped = case_when(agreement == "misscaled" ~ a,
                             preference_a > preference_b ~ a,
                             .default = b),
         kept = if_else(dropped == a, b, a))

dropped_sums <- aggregated_copies(
  filter(bioassays, !record %in% dropped_within))
dropped_pooled <- pooled_copies(
  filter(bioassays, !record %in% dropped_within))

# every copy dropped, with the record kept in its place
links <- bind_rows(
  bind_rows(dropped_entered, dropped_published, dropped_reviewed) %>%
    mutate(rule = paste("within a database:", coalesce(decision, rule))),
  cross_pairs %>%
    mutate(rule = if_else(agreement == "misscaled",
                          "between databases: wrong scale",
                          "between databases: copy")),
  tibble(dropped = dropped_sums,
         rule = "sum of species-level records of one study"),
  tibble(dropped = dropped_pooled,
         rule = "separate tests of one study pooled into one record")
) %>%
  select(rule, km, dropped, kept)

dropped_scale <- filter(cross_pairs, agreement == "misscaled")$dropped
dropped_cross <- setdiff(cross_pairs$dropped, dropped_scale)
candidates <- bioassays
bioassays <- bioassays %>%
  keep_records(!record %in% dropped_within, "copy within a database") %>%
  keep_records(!record %in% dropped_cross, "copy in another database") %>%
  keep_records(!record %in% dropped_scale, "number dead given as mortality") %>%
  keep_records(!record %in% dropped_sums, "sum of species records") %>%
  keep_records(!record %in% dropped_pooled, "pooled separate tests")

# Each record kept takes the best location and the finest years of its copies
# (locate_from_copies(), date_from_copies()). The IR Mapper extracts are the
# analysis-ready datasets of Moyes et al. (2019), abstracted and checked to one
# standard, so their coordinates are preferred among equally precise copies.
# The models use one year for each test: the middle of its years, rounded
# down, so a test over two years takes the first.
members <- copy_members(bioassays, links)
bioassays <- bioassays %>%
  locate_from_copies(candidates, members,
                     analysis_ready = "IR Mapper extracts") %>%
  date_from_copies(candidates, members) %>%
  mutate(year = (year_start + year_end) %/% 2)

# ---- country labels ----

# the country labels must agree with the coordinates; see check_country_labels()
country_check <- check_country_labels(bioassays)
dir.create("outputs/review", showWarnings = FALSE, recursive = TRUE)
bioassays %>%
  select(longitude, latitude, year_start, insecticide_type, source,
         citation) %>%
  bind_cols(country_check) %>%
  filter(action != "keep") %>%
  write.csv("outputs/review/bioassay_country_relabels.csv", row.names = FALSE)
bioassays <- bioassays %>%
  mutate(
    country_name = if_else(country_check$action == "relabel",
                           country_check$country_geometry, country_name),
    longitude = if_else(country_check$action == "move",
                        country_check$new_longitude, longitude),
    latitude = if_else(country_check$action == "move",
                       country_check$new_latitude, latitude)
  ) %>%
  keep_records(country_check$action != "drop", "far from its country") %>%
  left_join(country_region_lookup(), by = "country_name")

# ---- study labels ----

# A label for each study, to model noise shared by the bioassays of a study
# (the survey is the study and year). The databases cite the same paper
# differently, e.g. "Dabire KR, Diabate A, ... 2009;103(11):1113-20.", "Dabire
# et al. 2009" and "Dabire et al. (2009) Distribution of ...", so we label a
# study by the first author's surname and the publication year parsed from the
# citation, or by the citation itself where none parses (institutional or
# missing citations such as "Ministry of Health"). Labels are split by country,
# which stops institutional labels and personal communications of
# multi-country datasets pooling surveys across countries. See #53.
#
# Some studies' data are cited differently again (in a report and a paper, or
# by an institution and an author). Labels that share copies are listed in
# outputs/review/bioassay_study_links.csv, but not merged: merging them chains
# through reports that compile many studies (e.g. the PMI and WHO reports) into
# one label for most of the data.
labels <- candidates %>%
  select(record, label = paper)
links %>%
  filter(!is.na(kept)) %>%
  left_join(labels, by = c("dropped" = "record")) %>%
  left_join(select(labels, record, label), by = c("kept" = "record"),
            suffix = c("_dropped", "_kept")) %>%
  filter(label_dropped != label_kept) %>%
  count(label_kept, label_dropped, name = "copies") %>%
  arrange(desc(copies)) %>%
  write.csv("outputs/review/bioassay_study_links.csv", row.names = FALSE)

bioassays <- bioassays %>%
  left_join(select(labels, record, label), by = "record") %>%
  mutate(study = paste(label, country_name, sep = ", "),
         species_complex = "gambiae complex",
         n_imputed = !n_recorded,
         link = paper_link(link, citation))

# ---- records for review, and reports ----

# The pairs flagged for review, with what the reviewer needs: each record's
# citation, a link to the paper (its DOI, the database's URL, or a Google
# Scholar search), the fields that might distinguish two tests, and any
# decision already recorded. A reviewer copies a pair's row into
# review/bioassay_duplicates.csv and fills in `decision` ("copy" or
# "different") and `notes`. Undecided pairs are taken as copies.
review_fields <- c("source", "record_id", "citation", "link", "year_start",
                   "year_end", "site", "latitude", "longitude", "month",
                   "generation", "species", "mosquito_number", "n_recorded",
                   "died", "mortality_adjusted")
record_fields <- candidates %>%
  mutate(link = paper_link(link, citation)) %>%
  select(record, insecticide_type, concentration,
         all_of(review_fields))
side_by_side <- function(pairs) {
  a <- record_fields %>%
    rename_with(~ paste0(.x, "_a"), all_of(review_fields))
  b <- record_fields %>%
    select(record, all_of(review_fields)) %>%
    rename_with(~ paste0(.x, "_b"), all_of(review_fields))
  pairs %>%
    select(any_of(c("pair_id", "reason", "database_a", "database_b", "km",
                    "shared", "a", "b"))) %>%
    left_join(a, by = c("a" = "record")) %>%
    left_join(b, by = c("b" = "record")) %>%
    select(any_of(c("pair_id", "reason", "database_a", "database_b", "km",
                    "shared", "insecticide_type",
                    "concentration")),
           all_of(as.vector(rbind(paste0(review_fields, "_a"),
                                  paste0(review_fields, "_b")))))
}
filter(within_left, startsWith(rule, "review")) %>%
  mutate(reason = rule, km = round(km, 2)) %>%
  side_by_side() %>%
  left_join(decisions, by = "pair_id") %>%
  relocate(pair_id, decision) %>%
  arrange(!is.na(decision), reason, citation_a) %>%
  write_excel_csv("outputs/review/bioassay_review_pairs.csv", na = "")

write_excel_csv(year_corrections,
                "outputs/review/bioassay_year_corrections.csv", na = "")

links %>%
  rename(a = dropped, b = kept) %>%
  mutate(reason = rule, km = round(km, 2)) %>%
  side_by_side() %>%
  rename_with(~ sub("_a$", "_dropped", sub("_b$", "_kept", .x))) %>%
  write_excel_csv("outputs/review/bioassay_links.csv", na = "")

attrition <- attrition %>%
  left_join(count(bioassays, source, name = "kept"), by = "source") %>%
  arrange(match(source, map_chr(sources, "source")))
write.csv(attrition, "outputs/review/bioassay_attrition.csv", row.names = FALSE)
print(as.data.frame(attrition))

bioassays %>%
  select(mortality_adjusted, mosquito_number, died, n_imputed, year,
         year_start, year_end, month, latitude, longitude, site, species,
         species_complex, insecticide_type, insecticide_class, concentration,
         test_type, country_name, region, citation, link, study, source,
         database, record_id, coordinates_from) %>%
  saveRDS("data/clean/all_gambiae_complex_data.RDS")
