# The level panels of the main cross-validation figure (#36, A-C;
# fig_variance_explained.R): % of variance explained in single held-out
# bioassays, bioassay-based estimates against the models.
#
#   A  in sample      the interpolation test assays, predicted by the full-data
#                     fit (outputs/full_fit_maps_at_assays.csv), which saw them
#   B  interpolation  the same assays, held out (spatial_interpolation)
#   C  extrapolation  the held-out block assays (spatial_blocks, folds pooled)
#
# A and B share their assays, so the step from A to B changes only whether the
# model saw the local data.
#
# Bars, all scored as in R/validation_skill.R:
#   local      expected score of an independent second bioassay at the same
#              pixel-year, from rho-hat: 100 - 2 x noise share
#   best K     the nearest_neighbour_oracle null: the pooled K nearest
#              training assays of the same insecticide from the two most
#              recent years, at the K minimising its own MSE on the held-out
#              assays (R/null_models.R)
#   nearest    the nearest_neighbour null: the single nearest such assay
#   dynamical, two-stage
#              the map's point prediction (`map`, #32)
#
# Every bar carries n*_CV (#40), the number of bioassays at the pixel-year it
# is worth (n_star_cv() in R/validation_skill.R). Also the same bars per
# insecticide for B and C, for the per-insecticide figure and the
# supplementary n*_CV table, and the median over held-out assays of the
# distance to the nearest training assay of the same insecticide
# (distance_same, R/validation_geometry.R), for the panel titles.
#
# Reads outputs/cv_scores.csv (R/validation_metrics.R),
# outputs/full_fit_maps_at_assays.csv (R/full_fit_maps_at_assays.R) and
# outputs/cv_geometry.csv (R/validation_geometry.R).
# Writes outputs/cv_level_skill.csv,
# outputs/cv_level_skill_by_insecticide.csv and outputs/cv_n_star.csv (the
# supplementary table: n*_CV [95% interval] of the nearest bioassay and the
# two models, per insecticide and over all of them, for B and C).

source("R/validation_skill.R")

set.seed(2026 - 10 - 5)

scores <- require_map_columns(read.csv("outputs/cv_scores.csv",
                                   colClasses = c(fold = "character")))

experiments <- c(interpolation = "spatial_interpolation",
                 extrapolation = "spatial_blocks")
models <- c(best_k = "nearest_neighbour_oracle", nearest = "nearest_neighbour",
            dynamical = "dynamical", two_stage = "two_stage")

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
  }
  stopifnot(!anyNA(out[names(models)]))
  out
}

interpolation <- held_out(experiments[["interpolation"]])
extrapolation <- held_out(experiments[["extrapolation"]])

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

panel_bars <- function(data, bars) {
  bootstrap_summary(data, level_statistic(bars), rho_draws, n_star = TRUE) %>%
    mutate(assays = nrow(data), pixels = n_distinct(data$cell),
           variance = mean((data$observed - mean(data$observed)) ^ 2))
}

interpolation_bars <- panel_bars(interpolation,
                                 c(names(models), "full_dynamical",
                                   "full_two_stage"))
level <- bind_rows(
  interpolation_bars %>%
    filter(bar %in% c("ceiling", "local", "full_dynamical",
                      "full_two_stage")) %>%
    mutate(panel = "A", source = "full-data fit (outputs/two_stage/maps)",
           bar = sub("^full_", "", bar)),
  interpolation_bars %>%
    filter(!grepl("^full_", bar)) %>%
    mutate(panel = "B", source = "spatial_interpolation"),
  panel_bars(extrapolation, names(models)) %>%
    mutate(panel = "C", source = "spatial_blocks 1, 2")
) %>%
  relocate(panel, source, .before = 1)


# distances to the training data ----------------------------------------------------

# per held-out assay, distance_same at its pixel-year-insecticide in its fold
geometry <- read.csv("outputs/cv_geometry.csv",
                     colClasses = c(fold = "character"))
median_km <- function(data, experiment) {
  km <- data %>%
    select(fold, cell, year, insecticide_type) %>%
    left_join(geometry %>%
                filter(.data$experiment == !!experiment) %>%
                select(fold, cell, year, insecticide_type,
                       distance_same),
              by = c("fold", "cell", "year", "insecticide_type"),
              relationship = "many-to-one")
  stopifnot(!anyNA(km$distance_same))
  median(km$distance_same)
}
distances <- tibble(
  panel = c("B", "C"),
  median_km = c(median_km(interpolation, experiments[["interpolation"]]),
                median_km(extrapolation, experiments[["extrapolation"]])))

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
