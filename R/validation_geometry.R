# How predictive skill depends on how far the held-out record sits from the
# training data.
#
# Each cross-validation experiment reduces to a single number, but the folds
# differ enormously in difficulty: a held-out record two cells from a training
# record is a different problem from one 800 km away. Reporting skill against
# distance and against local data volume turns each fold set into a continuum of
# difficulty rather than one summary, makes the arbitrary geometry of the folds
# matter less, and answers directly the question a user of the map has — at what
# separation, and at what data density, should the mechanistic model be trusted
# over local interpolation (#12 review).
#
# Three axes, for every held-out record:
#
#   distance to the nearest training record of the same insecticide, in km
#   number of training records of the same insecticide within `radius_km`
#   years since the most recent training observation at that same pixel
#
# Distances are between pixel centroids, because the pixel is the unit the model
# predicts at. Sourcing this file needs only the fold definitions and the saved
# scores; nothing is refitted.

source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/validation_blocks.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

radius_km <- 100

# breaks chosen so that no bin is so thin that its excess MSE is noise: the
# distance breaks straddle the ~13 km minimum separation the interpolation fold
# guarantees, and the ~5 km grid resolution

cell_coordinates <- function(cells) {
  xy <- terra::xyFromCell(mask, cells)
  colnames(xy) <- c("longitude", "latitude")
  xy
}

# the three axes for one fold. Computed on distinct (cell, year, insecticide)
# keys rather than per assay, since that is all they depend on
fold_geometry <- function(training, test, experiment, fold) {

  keys <- test %>%
    distinct(cell, year, insecticide_type)

  distance <- fields::rdist.earth(cell_coordinates(keys$cell),
                                  cell_coordinates(training$cell),
                                  miles = FALSE)

  same_insecticide <- outer(keys$insecticide_type,
                            training$insecticide_type,
                            FUN = "==")
  distance_same <- distance
  distance_same[!same_insecticide] <- Inf

  # the most recent training year at the same pixel. Negative where the nearest
  # training observation at that pixel is from a later year than the held-out
  # record, which should not happen in a forecasting experiment and is how the
  # training set leak showed itself
  last_year <- tapply(training$year, training$cell, max)

  # computed outside the mutate: a column created earlier in a mutate() call
  # shadows a variable of the same name later in it, which silently turned the
  # distance matrix into the vector of minima
  nearest_any <- apply(distance, 1, min)
  nearest_same <- apply(distance_same, 1, min)
  within_radius <- rowSums(distance_same <= radius_km)

  keys %>%
    mutate(experiment = experiment,
           fold = fold,
           distance_any = nearest_any,
           distance_same = nearest_same,
           n_within_radius = within_radius,
           years_since_cell = year -
             as.numeric(last_year[as.character(cell)]),
           .before = everything())

}

geometry <- bind_rows(
  bind_rows(
    lapply(seq_along(spatial_blocks), function(index) {
      fold_geometry(spatial_blocks[[index]]$training,
                    spatial_blocks[[index]]$test,
                    "spatial_blocks",
                    as.character(index))
    })
  ),
  fold_geometry(spatial_interpolation$training,
                spatial_interpolation$test,
                "spatial_interpolation",
                "all"),
  bind_rows(
    lapply(seq_along(temporal_forecasting_folds), function(index) {
      fold <- temporal_forecasting_folds[[index]]
      fold_geometry(fold$training,
                    fold$test,
                    paste0("temporal_forecasting_", fold$cut_year),
                    names(temporal_forecasting_folds)[index])
    })
  )
)

write.csv(geometry, "outputs/cv_geometry.csv", row.names = FALSE)

cat("\ndistance from each held-out pixel to the nearest training record",
    "of the same insecticide (km):\n")
print(geometry %>%
        group_by(experiment) %>%
        summarise(keys = n(),
                  min = min(distance_same),
                  q25 = quantile(distance_same, 0.25),
                  median = median(distance_same),
                  q75 = quantile(distance_same, 0.75),
                  max = max(distance_same),
                  .groups = "drop") %>%
        mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
        as.data.frame())

cat("\nyears since the last training observation at the same pixel",
    "(negative means a later year, which is a training set leak):\n")
print(geometry %>%
        group_by(experiment) %>%
        summarise(keys = n(),
                  pixel_never_sampled = sum(is.na(years_since_cell)),
                  later_year_in_training = sum(years_since_cell < 0,
                                               na.rm = TRUE),
                  median = median(years_since_cell, na.rm = TRUE),
                  max = suppressWarnings(max(years_since_cell, na.rm = TRUE)),
                  .groups = "drop") %>%
        as.data.frame())
