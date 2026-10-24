# Functions used by R/prep_bioassays.R to read the bioassay databases, put them
# on a common format, find copies of the same bioassay, and label studies.

# ---- reading ----

# the columns every source is mapped to; any a source lacks are NA
bioassay_columns <- c(
  "record_id", "citation", "link", "site", "country_name", "latitude",
  "longitude", "year_start", "year_end", "month", "species", "generation",
  "insecticide_type", "concentration", "test_type", "mosquito_number",
  "mortality_adjusted"
)
numeric_bioassay_columns <- c("latitude", "longitude", "year_start",
                              "year_end", "month",
                              "concentration", "mosquito_number",
                              "mortality_adjusted")

# Read one source, given its specification: `read()` returns the raw table with
# every column as text, `keep()` applies the source's own filters to it, and
# `columns` maps our column names to the source's. Numbers are parsed after
# dropping "%" and spaces. Each record gets `row_hash`, a hash of its raw row
# without the source's record ID, which is the same for rows identical in
# every column. Sources without their own record IDs get one from this hash,
# numbered within identical rows, so that a record keeps its ID when the data
# are re-read.
read_bioassay_source <- function(spec) {
  raw <- spec$keep(spec$read())
  row_hash <- apply(raw[setdiff(names(raw), spec$columns["record_id"])], 1,
                    function(row) {
                      substr(rlang::hash(paste(row, collapse = "|")), 1, 12)
                    })
  if (!"record_id" %in% names(spec$columns)) {
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
      database = spec$database,
      row_hash = row_hash
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

# the number of decimal places a coordinate is given to, up to 6
decimal_places <- function(x) {
  places <- rep(6L, length(x))
  for (d in 5:0) places[abs(x * 10 ^ d - round(x * 10 ^ d)) < 1e-7] <- d
  places
}

# the great-circle distance in km between points (longitude, latitude) and
# (longitude_b, latitude_b)
km_between <- function(longitude, latitude, longitude_b, latitude_b) {
  fields::rdist.earth.vec(cbind(longitude, latitude),
                          cbind(longitude_b, latitude_b), miles = FALSE)
}

# ---- finding copies ----

# The groups of records joined by a set of pairs (`a`, `b` record numbers), as
# a group number for each of `records`
connected_groups <- function(records, a, b) {
  parent <- setNames(records, records)
  find <- function(x) {
    while (parent[[as.character(x)]] != x) x <- parent[[as.character(x)]]
    x
  }
  for (i in seq_along(a)) {
    ra <- find(a[i])
    rb <- find(b[i])
    if (ra != rb) parent[[as.character(max(ra, rb))]] <- min(ra, rb)
  }
  vapply(records, function(x) as.numeric(find(x)), numeric(1))
}

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

# Whether two site names name the same place, ignoring case, accents and
# punctuation: the same, one within the other (for names of four or more
# characters), or differing in at most a fifth of their characters (spelling),
# provided any numbers in them are the same ("Houeyiho (1)" and "Houeyiho (2)"
# are two sites)
same_site_name <- function(x, y) {
  kx <- coalesce(name_key(x), "")
  ky <- coalesce(name_key(y), "")
  ok <- nchar(kx) >= 4 & nchar(ky) >= 4 &
    gsub("[^0-9]", "", kx) == gsub("[^0-9]", "", ky)
  same <- rep(FALSE, length(kx))
  i <- which(ok)
  edits <- mapply(function(a, b) utils::adist(a, b)[1, 1], kx[i], ky[i])
  same[i] <- kx[i] == ky[i] | str_detect(kx[i], fixed(ky[i])) |
    str_detect(ky[i], fixed(kx[i])) |
    edits <= 0.2 * pmax(nchar(kx[i]), nchar(ky[i]))
  same
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

# Whether each record is imprecisely located: its site is named as a centre
# point or centroid, or a coordinate has at most one decimal place (about 5-11
# km)
imprecise_location <- function(site, latitude, longitude) {
  grepl("cent(re|er) ?point|centroid", site, ignore.case = TRUE) |
    pmin(decimal_places(latitude), decimal_places(longitude)) <= 1
}

# The copies of each record kept (`kept`), itself included (`member`), given
# the copies dropped (`links`). A copy of a copy is a copy (e.g. a record
# dropped within a database, whose kept record was then dropped for a copy in
# another database).
copy_members <- function(records, links) {
  copies <- links %>%
    filter(!is.na(kept)) %>%
    distinct(dropped, kept)
  members <- transmute(records, kept = record, member = record)
  repeat {
    found <- members %>%
      inner_join(copies, by = c("member" = "kept"),
                 relationship = "many-to-many") %>%
      distinct(kept, member = dropped) %>%
      anti_join(members, by = c("kept", "member"))
    if (nrow(found) == 0) break
    members <- bind_rows(members, found)
  }
  members
}

# The records kept, each with its coordinates and site name from the best
# located of its copies (copy_members()), since copies are one test. Precisely
# located copies are preferred (imprecise_location()), then a copy in an
# `analysis_ready` database, then the record's own coordinates, then a copy
# with the same site name, then the nearest. `coordinates_from` names the copy
# the coordinates came from.
locate_from_copies <- function(records, all_records, members,
                               analysis_ready = character()) {
  own <- select(records, kept = record, own_lat = latitude,
                own_lon = longitude, own_site = site)
  location <- members %>%
    left_join(all_records %>%
                select(record, record_key, database, latitude, longitude,
                       site),
              by = c("member" = "record")) %>%
    left_join(own, by = "kept") %>%
    mutate(imprecise = imprecise_location(site, latitude, longitude),
           same_name = same_site_name(site, own_site),
           km = km_between(longitude, latitude, own_lon, own_lat)) %>%
    arrange(kept, imprecise, !database %in% analysis_ready, member != kept,
            !same_name, km) %>%
    distinct(kept, .keep_all = TRUE) %>%
    select(record = kept, latitude, longitude, site,
           coordinates_from = record_key)
  records %>%
    select(-latitude, -longitude, -site) %>%
    left_join(location, by = "record")
}

# The records kept, each with the years (year_start to year_end) of the copy
# (copy_members()) that dates the test most finely: the shortest range, then
# the record's own years, then the preferred record's (`preference`).
date_from_copies <- function(records, all_records, members) {
  years <- members %>%
    left_join(select(all_records, record, year_start, year_end, preference),
              by = c("member" = "record")) %>%
    arrange(kept, year_end - year_start, member != kept, preference) %>%
    distinct(kept, .keep_all = TRUE) %>%
    select(record = kept, year_start, year_end)
  records %>%
    select(-year_start, -year_end) %>%
    left_join(years, by = "record")
}

# Some databases pool a study's separate tests at a site (in different months,
# say) into one record, summing the mosquitoes tested, where another database
# holds the tests. The pooled records: those whose number tested equals the sum
# over two or more records of one study at one site (to 0.01 degrees), year
# and insecticide in another database, within 25 km, and whose mortality is
# within 5 percentage points of theirs pooled (see the supplement for how
# these were chosen).
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

# Copies recorded in different years: pairs of records (`a`, `b`) of two
# databases citing the same paper (`paper`, the study label without its
# country), with the same country, insecticide, concentration,
# site name, number tested and number dead, from 5% to 95% mortality, within
# `max_km` (`km`), whose years don't overlap. One database sometimes gives a
# test another year, but a site is also visited again in later years, and below
# 5% or above 95% repeat visits often give the same result (e.g. 1 or 99 of
# 100). Matches between different papers are taken as separate tests, though
# some are one test published twice.
cross_year_pairs <- function(records, max_km = 100) {
  tests <- records %>%
    filter(n_recorded, died >= 0.05 * mosquito_number,
           died <= 0.95 * mosquito_number) %>%
    filter(!is.na(paper)) %>%
    select(record, database, paper, country_name, insecticide_type,
           concentration, mosquito_number, died, year_start, year_end, site,
           latitude, longitude)
  tests %>%
    inner_join(tests,
               by = c("paper", "country_name", "insecticide_type",
                      "concentration", "mosquito_number", "died"),
               suffix = c("", "_b"), relationship = "many-to-many") %>%
    filter(record < record_b, database != database_b,
           year_end < year_start_b | year_end_b < year_start) %>%
    mutate(km = km_between(longitude, latitude, longitude_b, latitude_b)) %>%
    filter(km <= max_km, same_site_name(site, site_b)) %>%
    select(a = record, b = record_b, km)
}

# The years to give copies recorded in different years (cross_year_pairs()).
# Records joined by these pairs are one test, and take the years of its best
# record, chosen by, in turn:
# 1. not being in the paper's publication year where another record is earlier
#    (a test is rarely run in the year it is published, and databases give the
#    publication year where a paper gives none);
# 2. its database giving the paper's records in these pairs more different
#    years (another lumps them into one);
# 3. its database coming first in `database_order`.
# A test whose best records differ in their years keeps its years. Returns the
# records whose years change, with the best record (`from`), its years
# (`from_start`, `from_end`) and the step that chose it (`rule`).
cross_year_dates <- function(records, database_order, max_km = 100) {
  pairs <- cross_year_pairs(records, max_km)
  records %>%
    filter(record %in% c(pairs$a, pairs$b)) %>%
    mutate(published = if_else(paper == citation, NA,
                               as.numeric(str_extract(paper, "\\d{4}$"))),
           test = connected_groups(record, pairs$a, pairs$b)) %>%
    group_by(paper, country_name, database) %>%
    mutate(years = n_distinct(year_start, year_end)) %>%
    group_by(test) %>%
    mutate(publication_year = coalesce(year_start == published, FALSE) &
             year_start > min(year_start),
           rank = match(database, database_order)) %>%
    arrange(publication_year, desc(years), rank, .by_group = TRUE) %>%
    mutate(best = publication_year == first(publication_year) &
             years == first(years) & rank == first(rank)) %>%
    filter(n_distinct(year_start[best], year_end[best]) == 1) %>%
    mutate(from = first(record),
           from_start = first(year_start),
           from_end = first(year_end),
           rule = case_when(
             publication_year != first(publication_year) ~ "publication year",
             years != first(years) ~ "more years",
             .default = "database order")) %>%
    ungroup() %>%
    filter(year_start != from_start | year_end != from_end) %>%
    select(record, from, from_start, from_end, rule)
}

# Some sources hold a bioassay split by species (e.g. An. arabiensis and An.
# gambiae s.s.) and another, or the same, source holds it summed to complex
# level. The complex-level records whose died and mosquito_number equal exactly
# the sums over the species-level records of one source and citation at the
# same pixel, year, insecticide and concentration. In October 2026, comparing
# against species-level records a year earlier or later found no such matches,
# so these are not coincidences. See issue #53.
aggregated_copies <- function(records) {
  species_sums <- records %>%
    filter(species != "gambiae complex") %>%
    group_by(cell, year_start, insecticide_type, concentration, source,
             citation) %>%
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
                      "concentration"),
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
# "(YYYY)", the journal year in Vancouver style (a year, optionally with month
# and day, followed by ";volume"), or else any 19xx or 20xx.
citation_author_year <- function(citation) {
  x <- str_squish(stringi::stri_trans_general(citation, "Latin-ASCII"))
  is_doi <- grepl("^(doi:\\s*|https?://(dx\\.)?doi\\.org/)?10\\.\\d{4,}/", x,
                  ignore.case = TRUE)
  is_generic <- grepl("^(unpublished|personal communication|pers\\.? comm)",
                      x, ignore.case = TRUE)
  surname <- str_match(x, "^([A-Za-z][A-Za-z'\\-]+)")[, 2]
  months <- "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*"
  year <- coalesce(
    str_match(x, "et al\\.*,?\\s*\\(?((19|20)\\d{2})")[, 2],
    str_match(x, "\\(((19|20)\\d{2})[a-z]?\\)")[, 2],
    str_match(x, paste0("\\b((19|20)\\d{2})[a-z]?(\\s+", months,
                        "(\\s+\\d{1,2})?)?\\s*;"))[, 2],
    str_match(x, "\\b((19|20)\\d{2})\\b")[, 2]
  )
  if_else(is_doi | is_generic | is.na(surname) | is.na(year),
          NA_character_,
          paste(gsub("['-]", "", tolower(surname)), year))
}

# ---- country labels ----

# check the country labels against the coordinates. The model takes each
# pixel's country (and so its initial condition) from this label, so the label
# and the pixel have to agree. See issue #13. Three rules, in order:
#
# 1. A record at rounded coordinates that falls outside its labelled country,
# but whose rounding box (the area the rounding implies) overlaps that country,
# is moved into it: to the pixel of the labelled country nearest the recorded
# point within the box, or failing that the nearest pixel of that country
# anywhere. Its coordinates become that pixel's centre, so the model puts it
# there. "Rounded" is at most one decimal place in both coordinates, a box of
# +/- 0.5 or 0.05 degrees. In practice only whole degrees trigger it: 76 records
# (6 locations) are at whole degrees and 256 at one decimal place, and the only
# rounded records outside their country are 47 in Rwanda at (29, -2), (29, -3)
# and (30, -1), 16-36 km into DR Congo and Uganda. Two decimals is +/- 0.5 km,
# well inside the tolerance of rule 2.
#
# 2. Otherwise a label is kept if the point is within country_tolerance_km of
# the labelled country, and else relabelled to the country it lies in. Of the
# 148 other records in a neighbouring country's polygon, 139 are within 3.5 km
# of their labelled country's border (GADM generalisation and coordinate
# rounding) and 9 are at 7.8 km or more, with nothing in between.
#
# 3. Records outside every polygon take their nearest country if within the
# tolerance (106 of 126, all within 4 km of their label) and are dropped
# otherwise (20 records 40 km offshore of Madagascar, outside the mask).
#
# Records labelled with a country that has no mask cells, and so no border
# polygon (8 in Cabo Verde), are left alone.
#
# Returns, for each of `records`, the country label, the country of its
# coordinates, the action to take and any new coordinates.
check_country_labels <- function(records, country_tolerance_km = 5) {
  country_borders <- readRDS("data/clean/country_borders.RDS")
  country_raster <- rast("data/clean/country_raster.tif")

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
  country_geometry <- country_borders$country_name[
    ifelse(is.na(polygon_index), nearest_index, polygon_index)]

  # half-width of the rounding box in degrees, from the number of decimal places
  # of the less precise coordinate; NA if either has more than one
  latitude_places <- decimal_places(records$latitude)
  longitude_places <- decimal_places(records$longitude)
  box_half_width <- if_else(
    pmax(latitude_places, longitude_places) <= 1,
    0.5 * 10 ^ -pmin(latitude_places, longitude_places), NA_real_)

  box_overlaps_label <- rep(FALSE, nrow(records))
  candidates <- which(!is.na(box_half_width) & !is.na(label_index) &
                        km_from_label > 0)
  for (i in candidates) {
    box <- st_as_sfc(st_bbox(c(
      xmin = records$longitude[i] - box_half_width[i],
      xmax = records$longitude[i] + box_half_width[i],
      ymin = records$latitude[i] - box_half_width[i],
      ymax = records$latitude[i] + box_half_width[i]),
      crs = st_crs(country_borders)))
    box_overlaps_label[i] <- lengths(st_intersects(
      box, country_borders[label_index[i], ])) > 0
  }

  # the pixel of the labelled country to move each rounded record to, per
  # location
  country_cells <- data.frame(cell = cells(country_raster)) %>%
    mutate(country_name = as.character(
             terra::extract(country_raster, cell)$country_name)) %>%
    bind_cols(as.data.frame(xyFromCell(country_raster, .$cell)))
  move_to <- records %>%
    mutate(half_width = box_half_width) %>%
    filter(box_overlaps_label) %>%
    distinct(longitude, latitude, country_name, half_width)
  move_to$new_longitude <- NA_real_
  move_to$new_latitude <- NA_real_
  move_to$within_box <- NA
  for (j in seq_len(nrow(move_to))) {
    options <- country_cells %>%
      filter(country_name == move_to$country_name[j])
    in_box <- abs(options$x - move_to$longitude[j]) <= move_to$half_width[j] &
      abs(options$y - move_to$latitude[j]) <= move_to$half_width[j]
    move_to$within_box[j] <- any(in_box)
    if (any(in_box)) {
      options <- options[in_box, ]
    }
    km <- as.numeric(st_distance(
      st_as_sf(move_to[j, ], coords = c("longitude", "latitude"),
               crs = st_crs(country_borders)),
      st_as_sf(options, coords = c("x", "y"),
               crs = st_crs(country_borders)))) / 1000
    move_to$new_longitude[j] <- options$x[which.min(km)]
    move_to$new_latitude[j] <- options$y[which.min(km)]
  }
  move_index <- match(paste(records$longitude,
                            records$latitude,
                            records$country_name),
                      paste(move_to$longitude, move_to$latitude,
                            move_to$country_name))
  move_index[!box_overlaps_label] <- NA

  country_check <- tibble(
    country_label = records$country_name,
    country_geometry = country_geometry,
    km_from_label = km_from_label,
    km_outside_polygons = ifelse(is.na(polygon_index), km_from_nearest, 0),
    new_longitude = move_to$new_longitude[move_index],
    new_latitude = move_to$new_latitude[move_index],
    new_pixel_within_box = move_to$within_box[move_index],
    action = case_when(
      is.na(label_index) ~ "keep",
      box_overlaps_label ~ "move",
      km_from_label <= country_tolerance_km ~ "keep",
      km_outside_polygons <= country_tolerance_km ~ "relabel",
      .default = "drop"
    )
  )

  country_check
}
