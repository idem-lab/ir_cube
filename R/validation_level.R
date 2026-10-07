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
# Also the median over held-out assays of the distance to the nearest training
# assay of the same insecticide (as R/validation_geometry.R's distance_same,
# but per assay rather than per pixel-year-insecticide), for the panel titles.
#
# Reads outputs/cv_scores.csv (R/validation_metrics.R),
# and outputs/full_fit_maps_at_assays.csv (R/full_fit_maps_at_assays.R; panel A
# is a placeholder without it).
# Writes outputs/cv_level_skill.csv.

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
    transmute(cell, year_start, insecticide_type, insecticide_class, fold,
              died, mosquito_number, observed, floor)
  for (bar in names(models)) {
    rows <- data %>% filter(model == models[[bar]])
    stopifnot(identical(record_key(rows), record_key(reference)),
              identical(rows$fold, reference$fold))
    out[[bar]] <- rows$map
    if (bar == "two_stage") out$two_stage_map_exact <- rows$map_exact
  }
  stopifnot(!anyNA(out[names(models)]))
  out
}

interpolation <- held_out(experiments[["interpolation"]])
extrapolation <- held_out(experiments[["extrapolation"]])

# panel A: the full-data fit's maps at the same assays, as two more columns,
# so that A and B are scored on the same bootstrap resamples and share their
# ceiling and local bar. Without them A's model bars are placeholders
full_fit_file <- "outputs/full_fit_maps_at_assays.csv"
if (file.exists(full_fit_file)) {
  full_fit <- interpolation %>%
    select(cell, year_start, insecticide_type) %>%
    left_join(read.csv(full_fit_file),
              by = c("cell", "year_start", "insecticide_type"),
              relationship = "many-to-one")
  if (anyNA(full_fit$dynamical) || anyNA(full_fit$two_stage)) {
    stop(sum(is.na(full_fit$dynamical)), " interpolation test assays have ",
         "no full-fit map in ", full_fit_file,
         "; it must come from the same data as outputs/cv_scores.csv")
  }
  interpolation$full_dynamical <- full_fit$dynamical
  interpolation$full_two_stage <- full_fit$two_stage
} else {
  warning("no ", full_fit_file, "; panel A's model bars are placeholders")
  interpolation$full_dynamical <- NA_real_
  interpolation$full_two_stage <- NA_real_
}

rho_draws <- rho_replicates()

panel_bars <- function(data, bars) {
  bootstrap_summary(data, level_statistic(bars), rho_draws) %>%
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
  mutate(
    # a bar with no estimate is a placeholder, and drawn as one
    source = ifelse(is.na(estimate), "placeholder", source),
    map_exact = case_when(
      bar != "two_stage" ~ NA,
      panel == "B" ~ all(interpolation$two_stage_map_exact),
      panel == "C" ~ all(extrapolation$two_stage_map_exact),
      TRUE ~ TRUE)) %>%
  relocate(panel, source, .before = 1)


# distances to the training data ----------------------------------------------------

# the fold definitions (validation_blocks.R sources validation_folds.R)
suppressMessages({
  sink("/dev/null")
  source("R/validation_blocks.R")
  sink()
})

# per held-out assay, the great-circle distance between pixel centres to the
# nearest training assay of the same insecticide, any year. The distance
# depends only on the pixel and insecticide, so it is computed once for each
nearest_training_km <- function(training, test) {
  keys <- test %>% distinct(cell, insecticide_type)
  keys$km <- NA_real_
  for (type in unique(keys$insecticide_type)) {
    index <- which(keys$insecticide_type == type)
    training_cells <- unique(training$cell[training$insecticide_type == type])
    distance <- fields::rdist.earth(terra::xyFromCell(mask,
                                                      keys$cell[index]),
                                    terra::xyFromCell(mask,
                                                      training_cells),
                                    miles = FALSE)
    keys$km[index] <- apply(distance, 1, min)
  }
  keys$km[match(paste(test$cell, test$insecticide_type),
                paste(keys$cell, keys$insecticide_type))]
}

distances <- bind_rows(
  tibble(panel = "B",
         km = nearest_training_km(spatial_interpolation$training,
                                  spatial_interpolation$test)),
  bind_rows(lapply(spatial_blocks, function(block) {
    tibble(panel = "C", km = nearest_training_km(block$training, block$test))
  }))
) %>%
  group_by(panel) %>%
  summarise(median_km = median(km), .groups = "drop")
# the held-out sets must be the scored ones
stopifnot(nrow(spatial_interpolation$test) == nrow(interpolation),
          sum(vapply(spatial_blocks, function(b) nrow(b$test),
                     numeric(1))) == nrow(extrapolation))

level <- level %>% left_join(distances, by = "panel")


write.csv(level, "outputs/cv_level_skill.csv", row.names = FALSE)

cat("\n% of variance explained in single held-out bioassays:\n")
print(as.data.frame(level %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(panel, assays, pixels, median_km, bar, value) %>%
  pivot_wider(names_from = bar, values_from = value)), row.names = FALSE)
