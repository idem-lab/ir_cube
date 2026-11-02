# The panels of the main cross-validation figure (#36, A-E;
# fig_variance_explained.R): % of variance explained in single held-out
# bioassays, bioassay-based estimates against the models.
#
#   A  in sample      the interpolation test assays, predicted by the full-data
#                     fit (outputs/full_fit_maps_at_assays.csv), which saw them
#   B  interpolation  the same assays, held out (spatial_interpolation)
#   C  extrapolation  the held-out block assays (spatial_blocks, folds pooled)
#   D  forecast       the held-out assays of 2014-2018, fitted to the data
#                     before 2014 (temporal_forecasting_2014)
#   E  forecast       the same, 2018-2022 from the data before 2018
#
# A and B share their assays, so the step from A to B changes only whether the
# model saw the local data.
#
# Bars, all scored as in R/validation_skill.R:
#   local      A only: expected score of an independent second bioassay at
#              the same pixel-year, from rho-hat: 100 - 2 x noise share; and
#              local_2, local_3, of the mean of 2 or 3 (the pooled figure)
#   nearest    the nearest_neighbour null: the single nearest training assay
#              of the same insecticide from the two most recent years
#              (R/null_models.R)
#   dynamical, two-stage
#              the map's point prediction (`map`, #32)
#
# Every bar carries n*_CV (#40), the number of bioassays at the pixel-year it
# is worth (n_star_cv() in R/validation_skill.R). Also the same bars per
# insecticide for B and C, for the per-insecticide figure and the
# supplementary n*_CV table, and the median over held-out assays of the
# distance to the bioassay the nearest bar used, and the years from it
# (nearest_km and nearest_year in outputs/cv_scores.csv), for the panel
# titles.
#
# Reads outputs/cv_scores.csv (R/validation_metrics.R) and
# outputs/full_fit_maps_at_assays.csv (R/full_fit_maps_at_assays.R).
# Writes outputs/cv_level_skill.csv,
# outputs/cv_level_skill_by_insecticide.csv and outputs/cv_n_star.csv (the
# supplementary table: n*_CV [95% interval] of the nearest bioassay and the
# two models, per insecticide and over all of them, for B and C).

source("R/validation_skill.R")

set.seed(2026 - 10 - 5)

scores <- require_map_columns(read.csv("outputs/cv_scores.csv",
                                   colClasses = c(fold = "character")))

experiments <- c(interpolation = "spatial_interpolation",
                 extrapolation = "spatial_blocks",
                 forecast_2014 = "temporal_forecasting_2014",
                 forecast_2018 = "temporal_forecasting_2018")
models <- c(nearest = "nearest_neighbour", dynamical = "dynamical",
            two_stage = "two_stage")

# One row per held-out assay, one prediction column per bar. The models'
# records are matched on position, checked on record_key(), since identical
# records (#31) make the key not unique
held_out <- function(experiment) {
  data <- scores %>% filter(.data$experiment == !!experiment)
  reference <- data %>% filter(model == "dynamical")
  out <- reference %>%
    transmute(cell, year, insecticide_type, insecticide_class, fold,
              died, mosquito_number, observed, floor)
  for (bar in names(models)) {
    rows <- data %>% filter(model == models[[bar]])
    stopifnot(identical(record_key(rows), record_key(reference)),
              identical(rows$fold, reference$fold))
    out[[bar]] <- rows$map
    if (bar == "nearest") {
      out$nearest_km <- rows$nearest_km
      out$nearest_years <- rows$year - rows$nearest_year
    }
  }
  stopifnot(!anyNA(out[c(names(models), "nearest_km", "nearest_years")]))
  out
}

interpolation <- held_out(experiments[["interpolation"]])
extrapolation <- held_out(experiments[["extrapolation"]])
forecast_2014 <- held_out(experiments[["forecast_2014"]])
forecast_2018 <- held_out(experiments[["forecast_2018"]])

# panel A: the full-data fit's maps at the same assays, as two more columns,
# so that A and B are scored on the same bootstrap resamples and share their
# ceiling and local bar
full_fit_file <- "outputs/full_fit_maps_at_assays.csv"
full_fit <- interpolation %>%
  select(cell, year, insecticide_type) %>%
  left_join(read.csv(full_fit_file),
            by = c("cell", "year", "insecticide_type"),
            relationship = "many-to-one")
if (anyNA(full_fit$dynamical) || anyNA(full_fit$two_stage)) {
  stop(sum(is.na(full_fit$dynamical)), " interpolation test assays have ",
       "no full-fit map in ", full_fit_file,
       "; it must come from the same data as outputs/cv_scores.csv")
}
interpolation$full_dynamical <- full_fit$dynamical
interpolation$full_two_stage <- full_fit$two_stage

rho_draws <- rho_replicates()

# the numbers of pooled bioassays scored in sample
pools <- 1:3

panel_bars <- function(data, bars, pools = 1) {
  bootstrap_summary(data, level_statistic(bars, pools), rho_draws,
                    n_star = TRUE) %>%
    mutate(assays = nrow(data), pixels = n_distinct(data$cell),
           variance = mean((data$observed - mean(data$observed)) ^ 2))
}

interpolation_bars <- panel_bars(interpolation,
                                 c(names(models), "full_dynamical",
                                   "full_two_stage"),
                                 pools)
level <- bind_rows(
  interpolation_bars %>%
    filter(bar %in% c("ceiling", local_names(pools), "full_dynamical",
                      "full_two_stage")) %>%
    mutate(panel = "A", source = "full-data fit (outputs/two_stage/maps)",
           bar = sub("^full_", "", bar)),
  interpolation_bars %>%
    filter(!grepl("^full_", bar)) %>%
    mutate(panel = "B", source = "spatial_interpolation"),
  panel_bars(extrapolation, names(models)) %>%
    mutate(panel = "C", source = "spatial_blocks 1, 2"),
  panel_bars(forecast_2014, names(models)) %>%
    mutate(panel = "D", source = experiments[["forecast_2014"]]),
  panel_bars(forecast_2018, names(models)) %>%
    mutate(panel = "E", source = experiments[["forecast_2018"]])
) %>%
  # the local bioassays only in sample, where they exist
  filter(!grepl("^local", bar) | panel == "A") %>%
  relocate(panel, source, .before = 1)


# distances to the training data ----------------------------------------------------

# the medians over held-out assays of the distance and years to the bioassay
# the nearest bar used
nearest_medians <- function(data) {
  tibble(median_km = median(data$nearest_km),
         median_years = median(data$nearest_years))
}
distances <- bind_rows(
  B = nearest_medians(interpolation),
  C = nearest_medians(extrapolation),
  D = nearest_medians(forecast_2014),
  E = nearest_medians(forecast_2018),
  .id = "panel")

level <- level %>% left_join(distances, by = "panel")


write.csv(level, "outputs/cv_level_skill.csv", row.names = FALSE)


# per insecticide, B and C ----------------------------------------------------------

# Each insecticide bootstrapped on its own pixels, after the panels, so that
# their resamples are unchanged. Every insecticide with held-out assays at two
# or more pixels and some spread in them is scored. The figure draws only
# those whose comparison can be read (`shown`), on two counts that are not
# substitutes:
#
#   ceiling >= min_ceiling
#       more than half the observed spread in held-out mortality has to be real
#       variation in resistance rather than assay noise. Where it is not, the
#       ratio is two small numbers divided by each other: Bendiocarb in the
#       interpolation fold has 154 assays in 67 pixels and a predictable
#       standard deviation of about 6 points of mortality, so a model would
#       have to land inside 6 points to score above zero. A sample-size rule
#       cannot see this.
#
#   interval width <= max_ci_width
#       the estimate has to be resolvable, for the nearest bioassay and the
#       dynamical model, or the comparison the panel exists to make cannot be
#       read. This is where the pixel count enters, but only through its effect
#       on precision.
min_ceiling <- 50
max_ci_width <- 100
width_bars <- c("nearest", "dynamical")

by_insecticide <- bind_rows(lapply(
  list(B = interpolation, C = extrapolation), function(data) {
    bind_rows(lapply(split(data, data$insecticide_type), function(type) {
      if (n_distinct(type$cell) < 2 ||
          mean((type$observed - mean(type$observed)) ^ 2) == 0) {
        return(NULL)
      }
      panel_bars(type, c("nearest", "dynamical", "two_stage")) %>%
        mutate(insecticide_type = type$insecticide_type[1], .before = 1)
    }))
  }), .id = "panel") %>%
  group_by(panel, insecticide_type) %>%
  mutate(
    ceiling = estimate[bar == "ceiling"],
    ci_width = max((upper - lower)[bar %in% width_bars]),
    drop_reason = case_when(
      ceiling < min_ceiling & ci_width > max_ci_width ~ "no signal, imprecise",
      ceiling < min_ceiling                           ~ "no signal",
      ci_width > max_ci_width                         ~ "imprecise",
      TRUE                                            ~ NA_character_),
    shown = is.na(drop_reason)) %>%
  ungroup()

write.csv(by_insecticide, "outputs/cv_level_skill_by_insecticide.csv",
          row.names = FALSE)


# the supplementary n*_CV table -----------------------------------------------------

n_star_bars <- c(nearest = "nearest bioassay", dynamical = "dynamical model",
                 two_stage = "two-stage model")
n_star_table <- bind_rows(
  by_insecticide,
  level %>%
    filter(panel %in% c("B", "C")) %>%
    mutate(insecticide_type = "all", drop_reason = NA_character_)) %>%
  filter(bar %in% names(n_star_bars)) %>%
  transmute(experiment = c(B = "interpolation", C = "extrapolation")[panel],
            insecticide_type, assays, pixels,
            bar = n_star_bars[bar],
            n_star = sprintf("%.2f [%.2f, %.2f]", n_star, n_star_lower,
                             n_star_upper),
            note = coalesce(drop_reason, "")) %>%
  pivot_wider(names_from = bar, values_from = n_star) %>%
  relocate(note, .after = last_col()) %>%
  mutate(insecticide_type = factor(insecticide_type,
                                   c(sort(unique(by_insecticide$insecticide_type)),
                                     "all"))) %>%
  arrange(experiment != "interpolation", insecticide_type)

write.csv(n_star_table, "outputs/cv_n_star.csv", row.names = FALSE)

cat("\n% of variance explained in single held-out bioassays:\n")
print(as.data.frame(level %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(panel, assays, pixels, median_km, bar, value) %>%
  pivot_wider(names_from = bar, values_from = value)), row.names = FALSE)

cat("\nn*_CV, bioassays at the pixel-year [95% interval]:\n")
print(as.data.frame(n_star_table), row.names = FALSE)

cells <- by_insecticide %>% distinct(panel, insecticide_type, assays, pixels,
                                     ceiling, ci_width, shown, drop_reason)
cat(sprintf("\nby insecticide: %i panel-insecticides, %i shown\n",
            nrow(cells), sum(cells$shown)))
print(as.data.frame(cells %>%
  mutate(across(c(ceiling, ci_width), ~ round(.x, 1)),
         shown = ifelse(shown, "yes", drop_reason)) %>%
  select(-drop_reason)), row.names = FALSE)
