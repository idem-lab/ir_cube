# Functions behind R/bioassay_vs_map.R (#30): survey identifiers, the #31
# duplicate rule applied to held-out records, replicate pairs within held-out
# pixel-years, the clusters for the survey-structured bootstrap, and the n*
# summaries. Functions only: source from the repo root, after
# R/validation_scoring.R and R/two_stage_correction.R (for project_km()).

suppressMessages({
  library(dplyr)
  library(sf)
})


# surveys -----------------------------------------------------------------------

# Survey identifier: citation x country x year. Restored from tag
# two-stage-options-6edcf0d (R/two_stage_correction.R), where it defined the
# survey effect removed in #29. The citation alone does not identify a study:
# aggregated sources ("Ministry of Health", "PMI 2016", "VectorBase", personal
# communications) span countries, so the country splits them (110 of 1434
# citation-years span more than one country, with 9,902 assays). Assays with no
# citation (12 in the cleaned data, all Madagascar 2023) fall back to spatial
# clusters within a year: pixels within cluster_km of each other (single
# linkage) in the same country-year.
survey_id <- function(citation, country, year, coords_km, cluster_km = 25) {
  id <- paste(citation, country, year, sep = " | ")
  missing <- is.na(citation) | citation == ""
  if (any(missing)) {
    group <- paste(country, year)[missing]
    xy <- coords_km[missing, , drop = FALSE]
    cluster <- integer(sum(missing))
    for (g in unique(group)) {
      rows <- which(group == g)
      cluster[rows] <- if (length(rows) == 1) 1L else
        stats::cutree(stats::hclust(stats::dist(xy[rows, , drop = FALSE]),
                                    method = "single"), h = cluster_km)
    }
    id[missing] <- paste("no citation", group, "cluster", cluster, sep = " | ")
  }
  id
}


# held-out records --------------------------------------------------------------

# Attach the bioassay records (`data`: the modelled subset, with cell) to the
# held-out records of one fold (`held_out`), to recover the fields the per-record
# scores don't carry: citation, source, species and coordinates. Records are
# matched on pixel, year, insecticide and result. Records sharing all of these
# are interchangeable as far as the scores go, so the k-th held-out copy is given
# the k-th data row with that key; a held-out record with no data row left gets
# NA in `data_row`.
attach_bioassay_records <- function(held_out, data) {
  key <- c("cell", "year_start", "insecticide_type", "died", "mosquito_number")
  data <- data %>%
    mutate(data_row = row_number()) %>%
    group_by(across(all_of(key))) %>%
    mutate(copy = row_number()) %>%
    ungroup()
  held_out %>%
    group_by(across(all_of(key))) %>%
    mutate(copy = row_number()) %>%
    ungroup() %>%
    left_join(data %>% select(all_of(key), copy, data_row, longitude, latitude,
                              citation, source, species, species_complex,
                              concentration),
              by = c(key, "copy")) %>%
    select(-copy)
}

# the bioassay database each source table comes from (as R/prep_bioassays.R,
# commit 73b0963, #31)
bioassay_database <- function(source) {
  case_when(
    grepl("^mtm", source) ~ "MTM",
    grepl("^mapper", source) ~ "IR Mapper",
    grepl("^va", source) ~ "Vector Atlas"
  )
}

# The data rows of `records` that #31's rule (R/prep_bioassays.R, commit
# 73b0963, find_cross_database_duplicates()) drops as copies of one bioassay
# held by more than one database. The September folds predate the rule, so it
# is applied here to the held-out records; once #31 is merged and the folds are
# rerun (#35) the records are already deduplicated and this drops nothing.
#
# The rule as there: records sharing a pixel, year, insecticide and result
# (died and mosquito_number), from more than one database, are copies, except
# at 0% or 100% mortality, where species and concentration must also match and
# the sites be within max_extreme_km. Within a match set, the k-th preferred
# record of each database is matched to the k-th of each other database, and
# the most preferred of each matched set is kept. Preference: the modal
# concentration (every held-out record has it), an assigned species complex, a
# species-level identification, then data row order. Held-out records are
# matched only to each other, which is all the rule can act on here: the copies
# of one assay share a pixel-year, so they are held out together.
cross_database_duplicates <- function(records, max_extreme_km = 5) {
  matched <- records %>%
    mutate(
      database = bioassay_database(source),
      extreme = died == 0 | died == mosquito_number,
      species_key = if_else(extreme, species, NA),
      concentration_key = if_else(extreme, concentration, NA)
    ) %>%
    group_by(cell, year_start, insecticide_type, died, mosquito_number,
             species_key, concentration_key) %>%
    filter(n_distinct(database) > 1) %>%
    arrange(is.na(species_complex),
            species %in% c("gambiae complex", "funestus complex"),
            data_row,
            .by_group = TRUE) %>%
    group_by(database, .add = TRUE) %>%
    mutate(k = row_number()) %>%
    group_by(cell, year_start, insecticide_type, died, mosquito_number,
             species_key, concentration_key, k) %>%
    mutate(kept_row = first(data_row),
           kept_longitude = first(longitude),
           kept_latitude = first(latitude)) %>%
    ungroup() %>%
    filter(data_row != kept_row)
  if (nrow(matched) == 0) {
    return(integer())
  }
  km <- as.numeric(st_distance(
    st_as_sf(matched, coords = c("longitude", "latitude"), crs = 4326),
    st_as_sf(matched, coords = c("kept_longitude", "kept_latitude"),
             crs = 4326),
    by_element = TRUE)) / 1000
  matched$data_row[!matched$extreme | km <= max_extreme_km]
}


# replicate pairs ---------------------------------------------------------------

# Per held-out assay in a pixel-year (cell x year x insecticide) with at least
# one other held-out assay: the mean squared difference to its partners, over
# all of them (`pair_any`) and over those from a different survey only
# (`pair_other_survey`, NA without one), and the number of each.
#
# Averaging over each assay's partners, rather than over all pairs, weights
# every held-out assay once, as the map's squared error does. Averaging over
# pairs would weight a pixel-year of n assays by n^2 in D_b but by n in D_m.
replicate_partners <- function(y, survey) {
  difference <- outer(y, y, "-") ^ 2
  other <- outer(survey, survey, "!=")
  diag(difference) <- NA
  partners_other <- rowSums(other)
  data.frame(
    partners = length(y) - 1,
    pair_any = rowMeans(difference, na.rm = TRUE),
    partners_other_survey = partners_other,
    pair_other_survey = ifelse(partners_other > 0,
                               rowSums(difference * other, na.rm = TRUE) /
                                 partners_other,
                               NA_real_)
  )
}


# clusters ----------------------------------------------------------------------

# Connected components of a bipartite graph: items (here pixel-years) linked
# through shared labels (pixels and surveys). Union-find on the labels, with
# path halving; returns a component index per item.
link_components <- function(item, label) {
  label <- as.integer(factor(label))
  parent <- seq_len(max(label))
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  # each item's labels are joined to the item's first label
  for (rows in split(seq_along(item), item)) {
    root <- find(label[rows[1]])
    for (j in rows[-1]) {
      other <- find(label[j])
      if (other != root) parent[other] <- root
    }
  }
  roots <- vapply(seq_along(parent), find, integer(1))
  first_label <- label[match(unique(item), item)]
  setNames(as.integer(factor(roots[first_label])), unique(item))
}

# Clusters for the survey-structured bootstrap (Table 2): pixel-years that share
# a pixel or a survey are in one cluster. Pixel-years at one pixel share the
# map's error there, and the two assays of a pair share their surveys' errors
# with every other assay of those surveys. `assays` holds both assays of every
# pair (an assay with a different-survey partner is that partner's partner), so
# its surveys are all the pairs' surveys. Returns a cluster per row of `assays`.
survey_pixel_clusters <- function(assays) {
  links <- bind_rows(
    assays %>% distinct(pixel_year, label = paste("pixel", cell)),
    assays %>% distinct(pixel_year, label = paste("survey", survey))
  )
  components <- link_components(links$pixel_year, links$label)
  unname(components[as.character(assays$pixel_year)])
}


# n* --------------------------------------------------------------------------

# R, the map's error variance in population mortality per unit of one
# bioassay's noise variance W:
#
#   W = D_b / 2,   R = (D_m - W) / W,   n* = 1 / R
#
# with D_b the mean squared difference between two assays at a pixel-year and
# D_m the map's mean squared error against the held-out assay. The held-out
# assay's own noise W is in both, so it cancels without an estimate of rho.
# n pooled bioassays have noise variance W / n, so n* is the number whose mean
# predicts a held-out assay as accurately as the map.
#
# The model-floor variant treats the fitted u + p as error every assay at the
# pixel-year shares, of variance U, so that pooled bioassays approach the
# population mortality plus u + p rather than the population mortality:
#
#   R_model = (D_m - W - 2 U) / W
#
# The interval comes from bootstrapping R, whose denominator is always
# positive, and inverting: n* in [1 / R_upper, 1 / R_lower], with an infinite
# upper limit when R_lower <= 0. Bootstrapping n* directly fails in small
# cells, where the map's net error can be <= 0 in a draw.
excess_ratio <- function(pair, map_error, model_floor = 0) {
  w <- mean(pair) / 2
  (mean(map_error) - w - 2 * mean(model_floor)) / w
}

# n* from R, infinite when the map's net error is <= 0: no number of bioassays
# matches it
invert_ratio <- function(ratio) {
  ifelse(ratio > 0, 1 / ratio, Inf)
}

# bootstrap replicates of `statistic` (a function of a subset of `data`
# returning a named numeric vector), resampling whole clusters with
# replacement. As pixel_bootstrap() in R/validation_scoring.R, for any
# cluster column, and always a matrix with one column per statistic (t() of
# replicate() loses the column name of a single statistic).
cluster_bootstrap <- function(data, cluster, statistic, n_bootstrap = 2000) {
  ids <- unique(cluster)
  rows <- split(seq_len(nrow(data)), cluster)
  do.call(rbind, replicate(n_bootstrap, {
    picked <- sample(ids, length(ids), replace = TRUE)
    statistic(data[unlist(rows[as.character(picked)], use.names = FALSE), ])
  }, simplify = FALSE))
}

# n* point estimate and interval from a ratio and its bootstrap replicates
n_star_interval <- function(ratio, replicates) {
  ratio_lower <- unname(quantile(replicates, 0.025, na.rm = TRUE))
  ratio_upper <- unname(quantile(replicates, 0.975, na.rm = TRUE))
  c(n_star = invert_ratio(ratio),
    n_star_lower = invert_ratio(ratio_upper),
    n_star_upper = invert_ratio(ratio_lower))
}

