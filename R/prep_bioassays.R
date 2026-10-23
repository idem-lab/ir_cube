# Prepare the bioassay data for modelling: read the five databases, put them on
# a common format, remove copies of the same bioassay, and save the An. gambiae
# complex records to data/clean/all_gambiae_complex_data.RDS. The supplement
# (doc/supp.qmd, "Bioassay data processing") describes the method;
# outputs/review/bioassay_attrition.csv counts the records each step removes,
# and outputs/review/bioassay_links.csv lists every copy removed.

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

# The IR Mapper analysis-ready extracts are in Windows-1252, with "NR" and "NF"
# for missing values. Records pooled over several sites ("multipoint") are
# dropped.
read_ir_mapper_extract <- function(file) {
  read_csv(file, col_types = cols(.default = "c"),
           na = c("", "NA", "NR", "NF"),
           locale = locale(encoding = "windows-1252"))
}
ir_mapper_extract_columns <- c(record_id = "Bioassay ID",
                               citation = "Source 1 citation",
                               site = "Site name", country_name = "Country",
                               latitude = "Latitude", longitude = "Longitude",
                               year_start = "Start year",
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
  list(
    source = "mapper_dis_species",
    database = "IR Mapper extracts",
    read = function() {
      read_ir_mapper_extract("data/raw/2_standard-WHO-susc-test_species.csv")
    },
    keep = function(d) filter(d, `Site type` == "point"),
    columns = c(ir_mapper_extract_columns, species = "Species")
  ),
  list(
    source = "mapper_int_complex",
    database = "IR Mapper extracts",
    read = function() {
      read_ir_mapper_extract(
        "data/raw/5_intensity-bioassays_complex-subgroup.csv")
    },
    keep = function(d) filter(d, `Site type` == "point"),
    columns = c(ir_mapper_extract_columns, species = "Complex/Subgroup")
  ),
  list(
    source = "mapper_dis_complex",
    database = "IR Mapper extracts",
    read = function() {
      read_ir_mapper_extract(
        "data/raw/1_standard-WHO-susc-test_complex-subgroup.csv")
    },
    keep = function(d) filter(d, `Site type` == "point"),
    columns = c(ir_mapper_extract_columns, species = "Complex/Subgroup")
  ),
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
  keep_records(grepl("who", tolower(test_type)), "not a WHO tube test") %>%
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
  # a positive mortality that gives no dead mosquito is a proportion read as a
  # percentage, e.g. 1.0 for 100% (#60)
  keep_records(!(died == 0 & mortality_adjusted > 0),
               "mortality a proportion")

# ---- copies of the same bioassay ----

# The databases hold many of the same bioassays, and some hold a bioassay more
# than once. Copies have the same country, year, insecticide and result
# (agreeing_pairs()), but not always the same coordinates: the databases
# georeference sites differently, and copies in different databases are often
# tens of km apart. So copies are looked for within max_km of each other in
# different databases, and within max_km_within in the same database, where the
# concentration must match too.
#
# Each distance is where looking further would find more coincidences than
# copies. Coincidences are counted by pairs whose results differ by 3-5 dead
# mosquitoes (or 3-5 percentage points), which copies don't, so this is for
# results between 0% and 100%. Between databases, records with the same number
# of mosquitoes and number dead are copies up to 100 km apart: 98% of copies
# are within 100 km, against 62% within a pixel, and about 3% of the matches
# are coincidences. Matching on the percentage alone, where a number of
# mosquitoes is missing, is less specific: 10% of matches within 50 km are
# coincidences, and beyond 50 km coincidences outnumber copies. Within a
# database, copies are mostly at the same coordinates, and coincidences grow
# faster with distance, as one study's sites tend to have similar results. At
# 0% or 100% (where coincidences are counted from deltamethrin records matched
# to permethrin ones) up to about a fifth of matches are coincidences at any
# distance. See #63.
max_km <- c(count = 100, percent = 50, misscaled = 100)
max_km_within <- 10

# Which copy to keep: one with a recorded number of mosquitoes, then with the
# most common concentration of its insecticide (so that it stays in the
# modelled subset; the other is probably mislabelled), then with a
# species-level identification, then from the first source in the list above
bioassays <- bioassays %>%
  group_by(insecticide_type) %>%
  mutate(modal_concentration = concentration == sample_mode(concentration)) %>%
  ungroup() %>%
  arrange(!n_recorded, !modal_concentration, species == "gambiae complex",
          source_rank, record) %>%
  mutate(preference = row_number(),
         record_key = paste(source, record_id)) %>%
  arrange(record)

# A reviewer can record whether a pair of records are copies ("copy") or
# separate tests ("different") in review/bioassay_duplicates.csv, and for
# copies which of the two is located more finely (`finer`: "a" or "b"; see
# merge_copies()). A pair decided "different" is never taken as copies; any
# other decision (blank, "unsure") leaves the rules below to decide. The pairs
# flagged for review are written there (below).
review_file <- "review/bioassay_duplicates.csv"
reviewed <- tibble(pair_id = character(), decision = character(),
                   notes = character(), finer = character())
if (file.exists(review_file)) {
  reviewed <- read_csv(review_file, col_types = cols(.default = "c")) %>%
    bind_rows(reviewed)
}
decisions <- filter(reviewed, decision %in% c("copy", "different"))
finer_pairs <- decisions %>%
  filter(decision == "copy", finer %in% c("a", "b")) %>%
  transmute(key_a = paste(source_a, record_id_a),
            key_b = paste(source_b, record_id_b),
            finer = if_else(finer == "a", key_a, key_b),
            coarser = if_else(finer == key_a, key_b, key_a))

record_info <- bioassays %>%
  select(record, record_key, database, citation, site, month, generation,
         preference)
pairs <- agreeing_pairs(bioassays, max_km, max_km_within) %>%
  left_join(record_info, by = c("a" = "record")) %>%
  left_join(record_info, by = c("b" = "record"), suffix = c("_a", "_b")) %>%
  mutate(pair_id = pair_id(record_key_a, record_key_b),
         name_match = same_site_name(site_a, site_b)) %>%
  left_join(select(decisions, pair_id, decision), by = "pair_id")

# Within a database. Two records with the same citation are one entered twice
# if they are at the same coordinates in the same month and generation, and
# else separate tests. Records with different citations are copies if the two
# citations share at least two results between 0% and 100% (the same data
# published twice, e.g. in a paper and a report), and else are flagged for
# review.
same_text <- function(x, y) {
  coalesce(as.character(x), "") == coalesce(as.character(y), "")
}
shared_results <- pairs %>%
  filter(within, agreement != "misscaled", !extreme,
         !same_text(citation_a, citation_b)) %>%
  count(database_a, first = pmin(citation_a, citation_b),
        second = pmax(citation_a, citation_b), name = "shared")
within_pairs <- pairs %>%
  filter(within, agreement != "misscaled") %>%
  mutate(first = pmin(citation_a, citation_b),
         second = pmax(citation_a, citation_b)) %>%
  left_join(shared_results, by = c("database_a", "first", "second")) %>%
  mutate(
    shared = coalesce(shared, 0L),
    rule = case_when(
      same_text(citation_a, citation_b) & km < 0.01 &
        same_text(month_a, month_b) &
        same_text(generation_a, generation_b) ~ "entered twice",
      same_text(citation_a, citation_b) ~ "separate tests",
      shared >= 2 ~ "published twice",
      .default = "review: different citations"
    )
  )

# The records to drop from pairs of copies, each with the preferred record
# kept in its place. A record entered several times is grouped (a copy of a
# copy is a copy). Copies under different citations are paired one-to-one
# (one_to_one()), as two citations can hold the same tests at several sites
# with one result, such as 100%.
drop_copies <- function(pairs, grouped) {
  if (nrow(pairs) == 0) {
    return(tibble(dropped = integer(), kept = integer(), rule = character(),
                  decision = character(), km = numeric()))
  }
  if (grouped) {
    records <- unique(c(pairs$a, pairs$b))
    rules <- bind_rows(select(pairs, record = a, rule, decision, km),
                       select(pairs, record = b, rule, decision, km)) %>%
      distinct(record, .keep_all = TRUE)
    tibble(record = records,
           group = connected_groups(records, pairs$a, pairs$b)) %>%
      left_join(select(record_info, record, preference), by = "record") %>%
      group_by(group) %>%
      mutate(kept = record[which.min(preference)]) %>%
      ungroup() %>%
      filter(record != kept) %>%
      left_join(rules, by = "record") %>%
      select(dropped = record, kept, rule, decision, km)
  } else {
    pairs %>%
      one_to_one(between = c("citation_a", "citation_b")) %>%
      mutate(dropped = if_else(preference_a > preference_b, a, b),
             kept = if_else(dropped == a, b, a)) %>%
      select(dropped, kept, rule, decision, km)
  }
}
# pairs under these rules, unless a reviewer decided they are different tests
is_copy <- function(pairs, rules) {
  pairs$rule %in% rules & coalesce(pairs$decision == "copy", TRUE)
}
not_dropped <- function(pairs, dropped) {
  filter(pairs, !a %in% dropped$dropped, !b %in% dropped$dropped)
}
dropped_entered <- drop_copies(
  filter(within_pairs, is_copy(within_pairs, "entered twice")),
  grouped = TRUE)
within_left <- not_dropped(within_pairs, dropped_entered)
dropped_published <- drop_copies(
  filter(within_left, is_copy(within_left, "published twice")),
  grouped = FALSE)
within_left <- not_dropped(within_left, dropped_published)

# Then, among the records left: a result between 0% and 100% that one citation
# gives more than twice, at different sites or in different months or
# generations, is implausible as separate tests (e.g. one result entered 22
# times in 3 months), so these are flagged for review too. Pairs flagged for
# review are taken as copies until reviewed.
repeated <- within_left %>%
  filter(rule == "separate tests", !extreme) %>%
  select(a, b) %>%
  pivot_longer(everything(), values_to = "record") %>%
  count(record) %>%
  filter(n >= 2) %>%
  pull(record)
within_left <- within_left %>%
  mutate(rule = if_else(rule == "separate tests" &
                          (a %in% repeated | b %in% repeated),
                        "review: one result repeated", rule))
dropped_reviewed <- bind_rows(
  drop_copies(filter(within_left,
                     is_copy(within_left, "review: one result repeated")),
              grouped = TRUE),
  drop_copies(filter(within_left,
                     is_copy(within_left, "review: different citations")),
              grouped = FALSE)
)
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
         is.na(decision) | decision == "copy",
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

# Each record kept takes the best location and species of its copies (see
# merge_copies()). The IR Mapper extracts are the analysis-ready datasets of
# Moyes et al. (2019). In October 2026, of the tests where their copy was more
# than 5 km from IR Mapper 2022's, MTM's copy agreed with the extracts in 392
# and with IR Mapper 2022 in 192; against the Vector Atlas, 94 and 11.
bioassays <- merge_copies(bioassays, candidates, links, finer_pairs,
                          analysis_ready = "IR Mapper extracts")

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
  transmute(record,
            label = coalesce(citation_author_year(citation), citation))
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

# The pairs flagged for review, and any other pair a reviewer has decided, with
# what the reviewer needs: each record's citation, a link to the paper (its
# DOI, the database's URL, or a Google Scholar search), and the fields that
# might distinguish two tests. The reviewer fills in `decision` ("copy" or
# "different"), `finer` and `notes`, which are kept when the file is
# rewritten. Undecided pairs are taken as copies.
review_fields <- c("source", "record_id", "citation", "link", "site",
                   "latitude", "longitude", "month", "generation", "species",
                   "mosquito_number", "n_recorded", "died",
                   "mortality_adjusted")
record_fields <- candidates %>%
  mutate(link = paper_link(link, citation)) %>%
  select(record, year_start, insecticide_type, concentration,
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
                    "shared", "year_start", "insecticide_type",
                    "concentration")),
           all_of(as.vector(rbind(paste0(review_fields, "_a"),
                                  paste0(review_fields, "_b")))))
}
review <- bind_rows(
  filter(within_left, startsWith(rule, "review")),
  filter(pairs, !is.na(decision))
) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  mutate(reason = coalesce(rule, "between databases"),
         km = round(km, 2)) %>%
  side_by_side()
review <- review %>%
  mutate(across(everything(), as.character)) %>%
  left_join(select(reviewed, pair_id, decision, notes, finer),
            by = "pair_id") %>%
  bind_rows(reviewed %>%
              filter(!pair_id %in% review$pair_id, !is.na(decision)) %>%
              mutate(reason = "no longer matched")) %>%
  relocate(pair_id, decision, finer, notes, reason) %>%
  arrange(!is.na(decision), reason, citation_a)
dir.create(dirname(review_file), showWarnings = FALSE)
write_excel_csv(review, review_file, na = "")

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
  select(mortality_adjusted, mosquito_number, died, n_imputed, year_start,
         month, latitude, longitude, site, species, species_complex,
         insecticide_type, insecticide_class, concentration, test_type,
         country_name, region, citation, link, study, source, database,
         record_id, coordinates_from) %>%
  saveRDS("data/clean/all_gambiae_complex_data.RDS")
