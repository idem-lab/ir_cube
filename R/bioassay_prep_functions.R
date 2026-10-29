# Functions used by R/prep_bioassays.R to read the bioassay databases, put them
# on a common format, find copies of the same bioassay, and label studies.

# ---- reading ----

# the columns every source is mapped to; any a source lacks are NA
bioassay_columns <- c(
  "record_id", "citation", "link", "site", "country_name", "latitude",
  "longitude", "year_start", "year_end", "month", "species", "generation",
  "insecticide_type", "concentration", "test_type", "mosquito_number",
  "mortality_adjusted", "synergist", "intensity"
)
numeric_bioassay_columns <- c("latitude", "longitude", "year_start",
                              "year_end", "month",
                              "concentration", "mosquito_number",
                              "mortality_adjusted")

# Read one source, given its specification: `read()` returns the raw table with
# every column as text, `keep()` applies the source's own filters to it, and
# `columns` maps our column names to the source's. Numbers are parsed after
# dropping "%" and spaces. Sources without their own record IDs get one from a
# hash of the raw row, numbered within identical rows, so that a record keeps
# its ID when the data are re-read.
read_bioassay_source <- function(spec) {
  raw <- spec$keep(spec$read())
  if (!"record_id" %in% names(spec$columns)) {
    row_hash <- apply(raw, 1, function(row) {
      substr(rlang::hash(paste(row, collapse = "|")), 1, 12)
    })
    raw$record_id <- paste0(row_hash, "-", ave(seq_along(row_hash), row_hash,
                                                FUN = seq_along))
    spec$columns <- c(spec$columns, record_id = "record_id")
  }
  out <- raw %>%
    select(all_of(spec$columns))
  out[setdiff(bioassay_columns, names(out))] <- NA_character_
  out %>%
    mutate(
      across(all_of(numeric_bioassay_columns),
             ~ suppressWarnings(as.numeric(str_remove_all(.x, "[%\\s]")))),
      citation = str_squish(citation),
      source = spec$source,
      database = spec$database
    )
}

# ---- standardising ----

# A name in lower-case ASCII letters and digits only, for matching names that
# differ in case, accents or punctuation
name_key <- function(x) {
  gsub("[^a-z0-9]", "", tolower(stringi::stri_trans_general(x, "Latin-ASCII")))
}

# Our name for each country: the UNSD name as printed in this project (see
# country_region_lookup()), matched ignoring case, accents and punctuation, or
# one of the databases' other spellings
standardise_country <- function(country) {
  names <- country_region_lookup()$country_name
  lookup <- c(
    setNames(names, name_key(names)),
    unitedrepublicoftanzania = "Tanzania",
    tanzaniaunitedrepublicof = "Tanzania",
    thegambia = "Gambia",
    democraticrepublicofthecongo = "DR Congo",
    drcongo = "DR Congo",
    centralafricanrepublic = "CAR",
    saotomeandprincipe = "Sao Tome & Principe",
    cotedivoire = "Côte d’Ivoire",
    swaziland = "Eswatini"
  )
  unname(lookup[name_key(country)])
}

# Our name for each insecticide, ignoring case, punctuation and the multiplier
# of an intensity assay ("deltamethrin 5x"), or NA for those not listed
standardise_insecticide <- function(insecticide) {
  names <- c("Alpha-cypermethrin", "Bendiocarb", "Carbosulfan", "Chlorfenapyr",
             "Chlorpyrifos-methyl", "Clothianidin", "Cyfluthrin", "DDT",
             "Deltamethrin", "Dieldrin", "Etofenprox", "Fenitrothion",
             "Lambda-cyhalothrin", "Malathion", "Permethrin",
             "Pirimiphos-methyl", "Propoxur")
  key <- function(x) gsub("[^a-z]", "", sub("\\s*[0-9]+x$", "", tolower(x)))
  unname(setNames(names, key(names))[key(insecticide)])
}

# the class of each insecticide
insecticide_class_of <- function(insecticide) {
  classes <- c(
    "Alpha-cypermethrin" = "Pyrethroids", "Cyfluthrin" = "Pyrethroids",
    "Deltamethrin" = "Pyrethroids", "Etofenprox" = "Pyrethroids",
    "Lambda-cyhalothrin" = "Pyrethroids", "Permethrin" = "Pyrethroids",
    "DDT" = "Organochlorines", "Dieldrin" = "Organochlorines",
    "Bendiocarb" = "Carbamates", "Carbosulfan" = "Carbamates",
    "Propoxur" = "Carbamates",
    "Chlorpyrifos-methyl" = "Organophosphates",
    "Fenitrothion" = "Organophosphates", "Malathion" = "Organophosphates",
    "Pirimiphos-methyl" = "Organophosphates",
    "Clothianidin" = "Neonicotinoids", "Chlorfenapyr" = "Pyrroles"
  )
  unname(classes[insecticide])
}

# Our name for each species: genus spelled out, S and M molecular forms as An.
# gambiae s.s. and An. coluzzii (Coetzee et al. 2013), and mixed or unresolved
# identifications as the complex
standardise_species <- function(species) {
  species <- gsub("An. ", "Anopheles ", species, fixed = TRUE)
  case_when(
    species == "arabiensis" ~ "Anopheles arabiensis",
    species %in% c("stephensi", "Anopheles stephensi sl") ~
      "Anopheles stephensi",
    species %in% c("coluzzii (gambiae m)", "coluzzii_gambiae_m form") ~
      "Anopheles coluzzii",
    species %in% c("gambiae (s)", "gambiae_s form") ~ "Anopheles gambiae s.s.",
    species %in% c("Funestus Subgroup", "funestus", "Anopheles funestus s.l.",
                   "Anopheles funestus sl", "FUNESTUS COMPLEX") ~
      "funestus complex",
    species %in% c("Gambiae Complex", "GAMBIAE COMPLEX", "gambiae", "coluzzi",
                   "Anopheles gambiae s.l.", "Anopheles gambiae sl",
                   "gambiae (S_M)", "gambiae_s form_m form",
                   "Anopheles coluzzii/gambiae") ~ "gambiae complex",
    .default = species
  )
}

gambiae_complex_species <- c(
  "gambiae complex", "Anopheles gambiae s.s.", "Anopheles arabiensis",
  "Anopheles coluzzii", "Anopheles merus", "Anopheles gambiae",
  "Anopheles quadriannulatus", "Anopheles gambiae/Anopheles coluzzii",
  "Anopheles gambiae ss"
)

# A link to each record's paper: its DOI, else the database's URL, else a
# Google Scholar search for the citation. DOIs have no spaces, but some are
# given with spaces for hyphens ("10.1186/s12936 017 1960 3") or after a
# slash or dot ("10.1371/journal. pone.0249440"); in a link that is a DOI, a
# space after a slash or dot is dropped and any other is read as a hyphen.
paper_link <- function(link, citation) {
  is_doi <- grepl("^\\s*10\\.\\d{4,}/", link)
  link[is_doi] <- gsub("\\s+", "-",
                       gsub("([./])\\s+", "\\1", str_squish(link[is_doi])))
  doi <- str_extract(link, "10\\.\\d{4,}/\\S+")
  url <- if_else(grepl("^https?://", link), link, NA_character_)
  search <- if_else(is.na(citation), NA_character_,
                    paste0("https://scholar.google.com/scholar?q=",
                           vapply(coalesce(citation, ""), URLencode,
                                  character(1), reserved = TRUE)))
  coalesce(if_else(is.na(doi), NA_character_, paste0("https://doi.org/", doi)),
           url, search)
}

# the great-circle distance in km between points (longitude, latitude) and
# (longitude_b, latitude_b)
km_between <- function(longitude, latitude, longitude_b, latitude_b) {
  fields::rdist.earth.vec(cbind(longitude, latitude),
                          cbind(longitude_b, latitude_b), miles = FALSE)
}

# ---- finding copies ----

# Pairs of records whose results agree, as a tibble with the two records (`a`,
# `b`), how they agree, and the distance between them in km. Records must be
# in the same country, with overlapping years (from year_start to year_end),
# and use the same insecticide at the same concentration, and be within
# `max_km` (between databases, a distance for each way of agreeing) or
# `max_km_within` (within a database). Results agree in one of three ways:
#
# "count": both record the number of mosquitoes, and it and the number dead
# are equal.
#
# "percent": one or both don't record the number of mosquitoes, and the
# percentage mortalities are within half a percentage point.
#
# "misscaled": record a's mortality, read as a number dead or as a proportion,
# gives the number dead of record b, in another database, with the same number
# of mosquitoes, where read as a percentage it gives a number dead more than
# one away. Record a is then on the wrong scale.
#
# Where the agreed result is 0% or 100%, which genuinely separate tests often
# share, the species must also be compatible (the same, or one recorded only to
# the complex).
agreeing_pairs <- function(records, max_km, max_km_within) {

  # a row for each year of each record, so that records whose years overlap
  # share a year
  key <- c("year", "country_name", "insecticide_type")
  r <- records %>%
    mutate(year = map2(year_start, year_end, seq)) %>%
    unnest(year) %>%
    select(record, database, all_of(key), concentration, species,
           mosquito_number, died, mortality_adjusted, latitude, longitude,
           n_recorded)
  counted <- filter(r, n_recorded)

  count <- counted %>%
    inner_join(counted, by = c(key, "mosquito_number", "died"),
               suffix = c("", "_b"), relationship = "many-to-many") %>%
    filter(record < record_b) %>%
    mutate(agreement = "count",
           extreme = died == 0 | died == mosquito_number)

  # percentages within half a point are at most one apart when rounded
  percent <- r %>%
    filter(!n_recorded) %>%
    mutate(percent_key = round(mortality_adjusted)) %>%
    inner_join(r %>%
                 cross_join(tibble(offset = -1:1)) %>%
                 mutate(percent_key = round(mortality_adjusted) + offset) %>%
                 select(-offset),
               by = c(key, "percent_key"),
               suffix = c("", "_b"), relationship = "many-to-many") %>%
    filter(record != record_b,
           abs(mortality_adjusted - mortality_adjusted_b) <= 0.5,
           n_recorded_b | record < record_b) %>%
    mutate(agreement = "percent",
           extreme = mortality_adjusted %in% c(0, 100) |
             mortality_adjusted_b %in% c(0, 100))

  misscaled <- bind_rows(
    counted %>%
      mutate(reading = mortality_adjusted),
    counted %>%
      filter(mortality_adjusted <= 1) %>%
      mutate(reading = round(mortality_adjusted * mosquito_number))
  ) %>%
    filter(reading == round(reading), reading <= mosquito_number,
           abs(reading - died) > 1) %>%
    inner_join(counted, by = c(key, "mosquito_number", "reading" = "died"),
               suffix = c("", "_b"), relationship = "many-to-many") %>%
    filter(database != database_b) %>%
    mutate(agreement = "misscaled",
           extreme = reading == 0 | reading == mosquito_number) %>%
    select(-reading)

  complexes <- c("gambiae complex", "funestus complex")
  bind_rows(count, percent, misscaled) %>%
    mutate(
      km = km_between(longitude, latitude, longitude_b, latitude_b),
      within = database == database_b
    ) %>%
    filter(
      km <= if_else(within, max_km_within, max_km[agreement]),
      concentration == concentration_b,
      !extreme | species == species_b | species %in% complexes |
        species_b %in% complexes
    ) %>%
    distinct(a = record, b = record_b, agreement, extreme, km, within)
}

# Whether two site names are the same, ignoring case, accents and punctuation.
# Similar names are not taken as the same, as some are different places (e.g.
# Kawama and Kasama in Zambia).
same_site_name <- function(x, y) {
  kx <- coalesce(name_key(x), "")
  ky <- coalesce(name_key(y), "")
  kx != "" & kx == ky
}

# Choose pairs of copies one-to-one: pairs of records with the same site name
# (`name_match`) first, then in order of distance, accept a pair unless either
# record is already in an accepted pair with a record of the other's database
# (or citation, given `between`). So one record can't remove several separate
# tests of another database that share its result, and where the databases
# place sites inconsistently, records are paired by site name rather than by
# the nearest coordinates.
one_to_one <- function(pairs, between = c("database_a", "database_b")) {
  pairs <- arrange(pairs, !name_match, km)
  key_a <- paste(pairs$a, pairs[[between[2]]])
  key_b <- paste(pairs$b, pairs[[between[1]]])
  used <- new.env(hash = TRUE)
  accept <- logical(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    if (is.null(used[[key_a[i]]]) && is.null(used[[key_b[i]]])) {
      accept[i] <- TRUE
      used[[key_a[i]]] <- TRUE
      used[[key_b[i]]] <- TRUE
    }
  }
  pairs[accept, ]
}

# Some databases pool a study's separate tests at a site (in different months,
# say) into one record, summing the mosquitoes tested, where another database
# holds the tests. The pooled records: those whose number tested equals the sum
# over two or more records of one study at one site (to 0.01 degrees), year
# and insecticide in another database, within 25 km, and whose mortality is
# within 5 percentage points of theirs pooled (supplement step 8).
pooled_copies <- function(records) {
  tests <- records %>%
    filter(n_recorded) %>%
    mutate(site_key = paste(round(latitude, 2), round(longitude, 2)))
  sets <- tests %>%
    group_by(database, paper, site_key, year_start, insecticide_type) %>%
    filter(n() >= 2) %>%
    summarise(n_sum = sum(mosquito_number),
              pooled_mortality = 100 * sum(died) / n_sum,
              latitude = mean(latitude), longitude = mean(longitude),
              .groups = "drop")
  sets %>%
    inner_join(tests %>%
                 select(record, database_pooled = database, year_start,
                        insecticide_type, n_sum = mosquito_number,
                        mortality_adjusted, lat_p = latitude,
                        lon_p = longitude),
               by = c("year_start", "insecticide_type", "n_sum"),
               relationship = "many-to-many") %>%
    filter(database_pooled != database,
           abs(mortality_adjusted - pooled_mortality) <= 5,
           km_between(longitude, latitude, lon_p, lat_p) <= 25) %>%
    pull(record) %>%
    unique()
}

# Some versions (one database's records of one paper or citation,
# `version_key`) hold a bioassay both split by species (e.g. An. arabiensis and
# An. gambiae s.s.) and summed to complex level, e.g. in the IR Mapper
# extracts' species and complex files. The
# complex-level records whose died and mosquito_number equal exactly the sums
# over the version's species-level records of one source at the same pixel,
# year, insecticide and concentration. See issue #53.
aggregated_copies <- function(records) {
  species_sums <- records %>%
    filter(species != "gambiae complex") %>%
    group_by(cell, year_start, insecticide_type, concentration, version_key,
             source) %>%
    filter(n() > 1) %>%
    summarise(
      died_sum = sum(died),
      mosquito_number_sum = sum(mosquito_number),
      .groups = "drop"
    )
  records %>%
    filter(species == "gambiae complex") %>%
    inner_join(species_sums,
               by = c("cell", "year_start", "insecticide_type",
                      "concentration", "version_key"),
               relationship = "many-to-many") %>%
    filter(died == died_sum,
           mosquito_number == mosquito_number_sum) %>%
    pull(record) %>%
    unique()
}

# A stable ID for each pair of records, from their sources and record IDs
pair_id <- function(key_a, key_b) {
  first <- pmin(key_a, key_b)
  second <- pmax(key_a, key_b)
  vapply(paste(first, second, sep = " & "),
         function(x) substr(rlang::hash(x), 1, 12),
         character(1), USE.NAMES = FALSE)
}

# ---- study labels ----

# The first author's surname and the publication year of each citation, e.g.
# "dabire 2009", or NA where either can't be found, the citation is a DOI, or
# its first word is generic ("Unpublished", "Personal communication", "pers.
# comm."). The surname is the citation's first word, ASCII-transliterated, in
# lower case and without apostrophes or hyphens (the databases write
# "Wat'senga" and "Watsenga"). The year is the first of: "et al. YYYY",
# "(YYYY)" near the start, the journal year in Vancouver style (a year, optionally with month
# and day, followed by ";volume"), or else any 19xx or 20xx. A stray "Â" (UTF-8
# text decoded twice, as in some Vector Atlas citations, e.g. "Â Riveron
# (2017)") is dropped first.
citation_author_year <- function(citation) {
  x <- str_squish(stringi::stri_trans_general(gsub("\u00c2", "", citation),
                                              "Latin-ASCII"))
  is_doi <- grepl("^(doi:\\s*|https?://(dx\\.)?doi\\.org/)?10\\.\\d{4,}/", x,
                  ignore.case = TRUE)
  is_generic <- grepl("^(unpublished|personal communication|pers\\.? comm)",
                      x, ignore.case = TRUE)
  surname <- str_match(x, "^([A-Za-z][A-Za-z'\\-]+)")[, 2]
  months <- "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*"
  year <- coalesce(
    str_match(x, "et al\\.*,?\\s*\\(?((19|20)\\d{2})")[, 2],
    str_match(x, "^[^()]{0,80}\\(((19|20)\\d{2})[a-z]?\\)")[, 2],
    str_match(x, paste0("\\b((19|20)\\d{2})[a-z]?(\\s+", months,
                        "(\\s+\\d{1,2})?)?\\s*;"))[, 2],
    str_match(x, "\\b((19|20)\\d{2})\\b")[, 2]
  )
  if_else(is_doi | is_generic | is.na(surname) | is.na(year),
          NA_character_,
          paste(gsub("['-]", "", tolower(surname)), year))
}

# ---- one version of each paper ----

is_personal_communication <- function(citation) {
  grepl("personal communication|pers\\.? comm", citation, ignore.case = TRUE)
}

# Whether each citation is a paper: it gives a first author and a year
# (citation_author_year()), and isn't a personal communication or a report: of
# the PMI, WHO (OMS in French) or a national malaria control programme (NMCP),
# which compile many studies' data, or one that opens with its title rather
# than its authors (eight or more words, with no comma, before its year in
# brackets, e.g. "Rapport Interimaire de la Troisieme Annee ... (2010)").
is_paper_citation <- function(citation) {
  before_year <- str_match(citation, "^([^,(]*?)\\(\\s*(?:\\S+\\s+)?(?:19|20)\\d{2}\\)")[, 2]
  title_first <- str_count(coalesce(before_year, ""), "\\S+") >= 8
  !is.na(citation_author_year(citation)) &
    !is_personal_communication(citation) &
    !grepl("^(PMI|WHO|OMS)([^a-z]|$)|^President\\S* Malaria Initiative|\\bNMCP\\b",
           citation, ignore.case = TRUE) &
    !title_first
}

# Whether each of a paper's versions is dominated: another version is at least
# as good on every measure (the columns of `measures`, larger being better)
# and better on at least one.
dominated <- function(measures) {
  m <- as.matrix(measures)
  apply(m, 1, function(x) {
    any(apply(m, 1, function(y) all(y >= x) && any(y > x)))
  })
}

# The versions of the papers held by two or more databases, ranked. A version
# is the records of one paper (`paper`, in one country) in one database. Each
# version is measured on the records the model uses (`modelled`, record
# numbers), so that tests it doesn't use (e.g. intensity tests, or a
# mis-entered concentration) can't decide the choice:
# - whether its recorded numbers tested vary (a database giving every test one
#   number, e.g. 100, often gives a placeholder);
# - its places: the number of model pixels (`cell`) its records are in;
# - its years: the number of years its records give, as the model reads them
#   (the middle of a record's years, rounded down), so a version putting every
#   test in one year, or giving every test one range of years, has one.
# The best version is at least as good as each other on all three and better
# on one. Where several versions tie on all three, the one with more records
# the model uses is kept, then the one first in `database_order`. Where no
# version is best on all three (e.g. one gives varying numbers tested and
# another more years), the versions conflict: the paper needs a decision by
# hand, and the version ranked first meanwhile is the first of the versions
# not bettered, by records the model uses and then `database_order`.
# Returns one row per version, with its rank, whether the paper's versions
# conflict, and, for the first, the step that chose it (`reason`).
paper_versions <- function(records, modelled, database_order) {
  records %>%
    filter(is_paper) %>%
    mutate(modelled = record %in% modelled,
           year = (year_start + year_end) %/% 2) %>%
    group_by(paper, database) %>%
    summarise(records = n(),
              modelled_records = sum(modelled),
              numbers_vary = n_distinct(
                mosquito_number[modelled & n_recorded]) > 1,
              places = n_distinct(cell[modelled]),
              years = n_distinct(year[modelled]),
              .groups = "drop") %>%
    group_by(paper) %>%
    mutate(best = !dominated(tibble(numbers_vary, places, years))) %>%
    arrange(!best, desc(modelled_records), match(database, database_order),
            .by_group = TRUE) %>%
    mutate(rank = row_number(),
           conflict = n_distinct(numbers_vary[best], places[best],
                                 years[best]) > 1,
           reason = case_when(
             rank > 1 ~ NA_character_,
             n() == 1 ~ "the only version",
             conflict ~ "conflict: decide by hand",
             sum(best) == 1 ~ "best on numbers tested, places and years",
             modelled_records > lead(modelled_records) ~
               "most modelled records",
             .default = "database order")) %>%
    ungroup()
}

# The records of the versions not kept that the kept version doesn't match.
# A record is matched if the kept version has a record of the same test: the
# same insecticide and concentration, a site within `max_km` or of the same
# name within `name_km` (similar names far apart are different places, e.g.
# Kawama and Kasama in Zambia, 450 km apart), overlapping years, and a
# mortality within `tolerance` percentage points; or a record that agrees with
# it exactly (the same number tested and number dead, or mortality within half
# a point where a number tested is missing), within the distances used to find
# copies between databases, which allow for the databases placing one site
# differently. `copies` gives these
# exact agreements, as pairs of record numbers (`a`, `b`; agreeing_pairs()).
# `kept` gives each paper's kept database (`kept_database`).
unmatched_in_kept <- function(records, kept, copies, tolerance = 10,
                              max_km = 10, name_km = 100) {
  records <- inner_join(records, kept, by = "paper")
  kept_records <- filter(records, database == kept_database)
  kept_tests <- kept_records %>%
    select(paper, insecticide_type, concentration,
           site_k = site, latitude_k = latitude, longitude_k = longitude,
           year_start_k = year_start, year_end_k = year_end,
           mortality_k = mortality_adjusted)
  others <- filter(records, database != kept_database)
  same_test <- others %>%
    inner_join(kept_tests,
               by = c("paper", "insecticide_type", "concentration"),
               relationship = "many-to-many") %>%
    filter(pmax(year_start, year_start_k) <= pmin(year_end, year_end_k),
           abs(mortality_adjusted - mortality_k) <= tolerance) %>%
    mutate(km = km_between(longitude, latitude, longitude_k, latitude_k)) %>%
    filter(km <= max_km | (km <= name_km & same_site_name(site, site_k))) %>%
    pull(record)
  exact <- bind_rows(select(copies, record = a, other = b),
                     select(copies, record = b, other = a)) %>%
    filter(record %in% others$record, other %in% kept_records$record) %>%
    pull(record)
  others %>%
    filter(!record %in% c(same_test, exact)) %>%
    select(record, paper, database)
}

# ---- country labels ----

# check the country labels against the coordinates (supplement step 3). The
# model takes each pixel's country (and so its initial condition) from this
# label, so the label and the pixel have to agree. See issue #13. A label is
# kept if the point is within country_tolerance_km of the labelled country, and
# else relabelled to the country it lies in. Records outside every polygon take
# their nearest country if within the tolerance and are dropped otherwise.
# Records labelled with a country that has no mask cells, and so no border
# polygon (Cabo Verde), are left alone.
#
# Returns, for each of `records`, the country label, the country of its
# coordinates and the action to take.
check_country_labels <- function(records, country_tolerance_km = 5) {
  country_borders <- readRDS("data/clean/country_borders.RDS")

  record_points <- records %>%
    st_as_sf(coords = c("longitude", "latitude"),
             crs = st_crs(country_borders))
  polygon_index <- st_intersects(record_points, country_borders) %>%
    vapply(function(i) i[1], integer(1))
  nearest_index <- st_nearest_feature(record_points, country_borders)
  label_index <- match(records$country_name,
                       country_borders$country_name)
  km_from_label <- as.numeric(st_distance(record_points,
                                          country_borders[label_index, ],
                                          by_element = TRUE)) / 1000
  km_from_nearest <- as.numeric(st_distance(record_points,
                                            country_borders[nearest_index, ],
                                            by_element = TRUE)) / 1000

  tibble(
    country_label = records$country_name,
    country_geometry = country_borders$country_name[
      ifelse(is.na(polygon_index), nearest_index, polygon_index)],
    km_from_label = km_from_label,
    km_outside_polygons = ifelse(is.na(polygon_index), km_from_nearest, 0),
    action = case_when(
      is.na(label_index) ~ "keep",
      km_from_label <= country_tolerance_km ~ "keep",
      km_outside_polygons <= country_tolerance_km ~ "relabel",
      .default = "drop"
    )
  )
}
