# Bioassay or map: which better estimates mortality at a pixel-year with data?
# (#30)
#
#   Rscript R/bioassay_vs_map.R
#
# Run after R/validation_metrics.R, which writes the per-record scores read
# here (outputs/cv_scores.csv). Nothing is refitted or rescored.
#
# The quantity is n*: the number of bioassays from the same pixel-year whose
# mean predicts a held-out bioassay as accurately as the map. On held-out
# pixel-years (cell x year x insecticide) with at least two assays,
#
#   D_b = mean squared difference between an assay and its partners at the
#         pixel-year, averaged per held-out assay
#   D_m = mean squared difference between the map and the held-out assay
#   W   = D_b / 2, one assay's noise variance
#   n*  = W / (D_m - W)
#
# The held-out assay's noise is in both D_b and D_m, so no rho enters (see
# excess_ratio() in R/bioassay_vs_map_functions.R). The map is the map's point
# prediction (`map`, #32): for the two-stage model the posterior mean of
# plogis(m + omega + xi), without u and p; for the dynamical model the mean of
# its draws.
#
# Duplicate records (#31) are removed from the held-out records before pairing:
# copies of one assay held by two databases pair at zero difference.
#
# Writes, rows by model x CV experiment x insecticide class (plus all classes):
#
#   outputs/cv_bioassay_vs_map_independent.csv (Table 1)
#       n* from every held-out assay with a partner, with a bootstrap over
#       pixels: held-out pixel-years concentrate in few pixels, and those at one
#       pixel share the map's error there. Also n* with the model floor
#       (supplementary: u + p as error shared by every assay at the pixel-year,
#       n* = W / (D_m - W - 2 U)), and n* from rho-hat on all held-out assays,
#       F / (MSE - F) with F the data floor, for comparison. The latter is what
#       #36's caption quotes, and Table 1 is its empirical check.
#   outputs/cv_bioassay_vs_map_surveys.csv (Table 2, a robustness check)
#       n* from different-survey partners only, with a bootstrap over clusters
#       of pixel-years linked by a shared pixel or survey. A partner from the
#       held-out assay's own survey shares its survey error, which favours the
#       bioassay. Different-survey pairs are also more often both at 0% or 100%,
#       so a difference between the tables is not the survey effect alone.
#
# Columns `map_exact` (share of the model's records scored against the map
# itself) and `indicative` mark rows whose prediction is not yet the map: the
# September two-stage folds store no map draws, so there `map` falls back to
# the CV mean, which adds fresh u and p. Rows with few pixels (Table 1) or few
# clusters (Table 2) are flagged in `note`; cluster-bootstrap intervals are
# unreliable with few clusters (Cameron & Miller 2015, J Human Resources
# 50:317-372).

source("R/validation_scoring.R")
source("R/two_stage_correction.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/bioassay_vs_map_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(terra)
})

set.seed(2026 - 10 - 5)
n_bootstrap <- 2000
n_posterior <- 4000

models <- c("two_stage", "dynamical")

# flag cells below these sizes: Table 1 by pixels, Table 2 by clusters
min_pixels <- 10
min_clusters <- 20


# per-record scores ------------------------------------------------------------

scores <- read_csv("outputs/cv_scores.csv", show_col_types = FALSE,
                   col_select = any_of(c(
                     "model", "experiment", "fold", "insecticide_type",
                     "insecticide_class", "country_name", "year_start", "cell",
                     "died", "mosquito_number", "observed", "predicted",
                     "rho_external", "map", "map_exact", "floor",
                     "model_floor"))) %>%
  filter(model %in% models)

# the map, the data floor and the two-stage model floor per record (#32)
stopifnot(all(c("map", "map_exact", "floor", "model_floor") %in% names(scores)))

scores <- scores %>%
  mutate(fold = as.character(fold)) %>%
  group_by(model, experiment, fold) %>%
  mutate(record = row_number()) %>%
  ungroup()

# every model must describe the same records in the same order, or the map
# columns would be swapped between assays
reference <- scores %>% filter(model == models[1])
for (m in models[-1]) {
  other <- scores %>% filter(model == m)
  stopifnot(identical(paste(reference$experiment, reference$fold,
                            record_key(reference)),
                      paste(other$experiment, other$fold, record_key(other))))
}


# held-out records, with surveys, less #31's duplicates -------------------------

data <- subset_modelled_bioassays(
  readRDS("data/clean/all_gambiae_complex_data.RDS"),
  rast("data/clean/raster_mask.tif"))

records <- reference %>%
  select(experiment, fold, record, insecticide_type, insecticide_class,
         country_name, year_start, cell, died, mosquito_number, observed,
         floor) %>%
  group_by(experiment, fold) %>%
  group_modify(~ attach_bioassay_records(.x, data)) %>%
  ungroup()

unmatched <- sum(is.na(records$data_row))
cat(sprintf("held-out records: %i, of which %i not found in the bioassay data\n",
            nrow(records), unmatched))
stopifnot(unmatched / nrow(records) < 0.01)
records <- records %>% filter(!is.na(data_row))

records$survey <- survey_id(records$citation, records$country_name,
                            records$year_start,
                            project_km(records$longitude, records$latitude))

duplicates <- records %>%
  group_by(experiment, fold) %>%
  group_modify(~ tibble(data_row = cross_database_duplicates(.x))) %>%
  ungroup() %>%
  mutate(duplicate = TRUE)
cat("held-out records dropped as cross-database duplicates (#31):\n")
print(count(duplicates, experiment, fold))
records <- records %>%
  left_join(duplicates, by = c("experiment", "fold", "data_row")) %>%
  filter(is.na(duplicate)) %>%
  select(-duplicate)


# replicate partners ------------------------------------------------------------

records <- records %>%
  group_by(experiment, cell, year_start, insecticide_type) %>%
  mutate(pixel_year = cur_group_id(), n_pixel_year = n()) %>%
  ungroup()

partners <- records %>%
  filter(n_pixel_year >= 2) %>%
  group_by(experiment, pixel_year) %>%
  group_modify(~ bind_cols(select(.x, fold, record),
                           replicate_partners(.x$observed, .x$survey))) %>%
  ungroup()

records <- records %>%
  left_join(partners, by = c("experiment", "pixel_year", "fold", "record"))

# one row per model and held-out record
assays <- records %>%
  inner_join(scores %>% select(model, experiment, fold, record, map, map_exact,
                               model_floor),
             by = c("experiment", "fold", "record"),
             relationship = "one-to-many") %>%
  mutate(map_error = (observed - map) ^ 2) %>%
  # only what the tables use, which keeps the bootstraps fast
  select(model, experiment, insecticide_class, cell, pixel_year, survey,
         observed, floor, map_error, map_exact, model_floor, pair_any,
         pair_other_survey, partners_other_survey)

# the rows of every table: each class, and all classes
by_class <- function(data) {
  bind_rows(data, mutate(data, insecticide_class = "all"))
}
class_order <- c(sort(unique(records$insecticide_class)), "all")


# Table 1: all partners, pixel bootstrap ----------------------------------------

table_independent <- assays %>%
  by_class() %>%
  group_by(model, experiment, insecticide_class) %>%
  group_modify(function(data, key) {

    replicated <- data %>% filter(!is.na(pair_any))
    if (nrow(replicated) == 0) return(tibble())

    ratio <- function(d) {
      c(data = excess_ratio(d$pair_any, d$map_error),
        model = excess_ratio(d$pair_any, d$map_error,
                             ifelse(is.na(d$model_floor), 0, d$model_floor)))
    }
    point <- ratio(replicated)
    drawn <- cluster_bootstrap(replicated, replicated$cell, ratio, n_bootstrap)

    # n* from rho-hat on all held-out assays: F / (MSE - F), the level at which
    # #36's bar meets the score of n pooled bioassays
    floor_ratio <- function(d) {
      floor <- mean(d$floor, na.rm = TRUE)
      c(rho = (mean(d$map_error) - floor) / floor)
    }
    point_rho <- floor_ratio(data)
    drawn_rho <- cluster_bootstrap(data, data$cell, floor_ratio, n_bootstrap)

    has_model_floor <- any(!is.na(replicated$model_floor))
    w <- mean(replicated$pair_any) / 2

    bind_cols(
      tibble(
        map_exact = mean(data$map_exact),
        pixel_years = n_distinct(replicated$pixel_year),
        pixels = n_distinct(replicated$cell),
        assays = nrow(replicated),
        one_bioassay = w,
        map_error = mean(replicated$map_error) - w),
      as_tibble(as.list(n_star_interval(point[["data"]], drawn[, "data"]))),
      tibble(model_floor = if (has_model_floor)
               mean(replicated$model_floor) else NA_real_),
      as_tibble(as.list(setNames(
        if (has_model_floor)
          n_star_interval(point[["model"]], drawn[, "model"]) else rep(NA, 3),
        c("n_star_model_floor", "n_star_model_floor_lower",
          "n_star_model_floor_upper")))),
      tibble(assays_all = nrow(data),
             pixels_all = n_distinct(data$cell),
             data_floor_all = mean(data$floor, na.rm = TRUE),
             map_error_all = mean(data$map_error) -
               mean(data$floor, na.rm = TRUE)),
      as_tibble(as.list(setNames(
        n_star_interval(point_rho[["rho"]], drawn_rho[, "rho"]),
        c("n_star_rho", "n_star_rho_lower", "n_star_rho_upper"))))
    )
  }) %>%
  ungroup() %>%
  mutate(indicative = map_exact < 1,
         note = ifelse(pixels < min_pixels,
                       sprintf("fewer than %i pixels", min_pixels), ""),
         insecticide_class = factor(insecticide_class, class_order)) %>%
  arrange(model, experiment, insecticide_class)

write.csv(table_independent, "outputs/cv_bioassay_vs_map_independent.csv",
          row.names = FALSE)


# Table 2: different-survey partners, cluster bootstrap --------------------------

table_surveys <- assays %>%
  filter(!is.na(pair_other_survey)) %>%
  by_class() %>%
  group_by(model, experiment, insecticide_class) %>%
  group_modify(function(data, key) {

    cluster <- survey_pixel_clusters(data)
    ratio <- function(d) {
      c(data = excess_ratio(d$pair_other_survey, d$map_error))
    }
    point <- ratio(data)
    drawn <- cluster_bootstrap(data, cluster, ratio, n_bootstrap)
    w <- mean(data$pair_other_survey) / 2

    bind_cols(
      tibble(
        map_exact = mean(data$map_exact),
        pairs = sum(data$partners_other_survey) / 2,
        surveys = n_distinct(data$survey),
        clusters = n_distinct(cluster),
        pixel_years = n_distinct(data$pixel_year),
        pixels = n_distinct(data$cell),
        assays = nrow(data),
        one_bioassay = w,
        map_error = mean(data$map_error) - w),
      as_tibble(as.list(n_star_interval(point[["data"]], drawn[, "data"])))
    )
  }) %>%
  ungroup() %>%
  mutate(indicative = map_exact < 1,
         note = ifelse(clusters < min_clusters,
                       sprintf("fewer than %i clusters: interval not interpreted",
                               min_clusters), ""),
         insecticide_class = factor(insecticide_class, class_order)) %>%
  arrange(model, experiment, insecticide_class)

write.csv(table_surveys, "outputs/cv_bioassay_vs_map_surveys.csv",
          row.names = FALSE)


# print ---------------------------------------------------------------------------

interval <- function(estimate, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          estimate, lower, upper)
}

cat("\nTable 1, all partners (pixel bootstrap):\n")
print(as.data.frame(table_independent %>%
  transmute(model, experiment, class = insecticide_class, pixel_years, pixels,
            one_bioassay = round(one_bioassay, 3),
            map_error = round(map_error, 3),
            n_star = interval(n_star, n_star_lower, n_star_upper),
            n_star_model_floor = interval(n_star_model_floor,
                                          n_star_model_floor_lower,
                                          n_star_model_floor_upper),
            n_star_rho = interval(n_star_rho, n_star_rho_lower,
                                  n_star_rho_upper),
            indicative, note)), row.names = FALSE)

cat("\nTable 2, different-survey partners (cluster bootstrap):\n")
print(as.data.frame(table_surveys %>%
  transmute(model, experiment, class = insecticide_class, pairs, surveys,
            clusters, one_bioassay = round(one_bioassay, 3),
            map_error = round(map_error, 3),
            n_star = interval(n_star, n_star_lower, n_star_upper),
            indicative, note)), row.names = FALSE)
