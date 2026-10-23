# Fit the final two-stage model (#21, R/two_stage_correction.R) inside one
# cross-validation fold of #12, per insecticide type, and save its held-out
# predictive draws in the #12 fold format, so they are scored by the same code.
#
#   Rscript R/run_two_stage_folds.R <experiment> <fold> <type index>
#   Rscript R/run_two_stage_folds.R <experiment> <fold> assemble
#
# e.g. Rscript R/run_two_stage_folds.R spatial_interpolation all 3
#
# Both read the dynamical fold's logit draws at its training and held-out
# assays (outputs/cv_draws/assays/dynamical__<experiment>__<fold>.rds, saved
# by R/run_one_fold.R; fold_assay_draws() in R/fit_validation_fold.R), not the
# fold itself. The per-type fits run in parallel, one process each:
#
#   <k>        for insecticide type k, fit the final model to the training
#              assays with the posterior mean logit as m_ref, and draw at the
#              type's test assays with the paired dynamical draws, so every
#              draw sits on its own dynamical draw (the cut posterior):
#              m + omega + xi plus fresh draws of the noise terms u (per
#              held-out pixel-year), p (per held-out pixel) and sy (per
#              held-out study-year), never their fitted values; and the
#              same draws without the noise, the map plogis(m + omega + xi),
#              whose posterior mean is what the maps show and what the point
#              metrics score (#32); and, for the forecasting folds, the
#              map's posterior mean at the type's before-window assays
#              (training assays in the window before the cut), for the
#              change score (#36). All three from one set of field draws.
#              Saved to outputs/two_stage/parts/;
#   assemble   combine the types into
#              outputs/cv_draws/two_stage__<experiment>__<fold>.rds and the
#              fold's rows of outputs/two_stage/fit_summary.csv. Every type
#              needs its part: a fit that fails or does not converge stops
#              its type's step, and the fold is not assembled without it.
#
# The saved object has the #12 fields (model, experiment, fold, test_df,
# p_draws) plus map_draws (draws x held-out assays of the map), rho_type, the
# per-type rho of the fit (a data frame of insecticide_type and rho, from the
# table rho_lookup() scores against) and fit_summary (which holds each type's
# tau, sigma_p and sigma_sy, the model floor's inputs). The forecasting folds
# also carry map_before, the map's posterior mean at each assay of the
# dynamical fold's before_df, in its order.
#
# Run with OpenBLAS (reference BLAS is ~10x slower), e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=3 nice -n 10 Rscript ...

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3)
experiment_name <- arguments[1]
fold_name <- arguments[2]
step <- arguments[3]

source("R/two_stage_helpers.R")
output_dir <- "outputs/two_stage"
parts_dir <- file.path(output_dir, "parts")
dir.create(parts_dir, showWarnings = FALSE, recursive = TRUE)
fold_key <- sprintf("%s__%s", experiment_name, fold_name)
part_file <- function(k) file.path(parts_dir, sprintf("%s__%i.rds", fold_key, k))

suppressMessages(library(dplyr))
assays <- read_assay_draws(file.path("outputs/cv_draws",
                                     sprintf("dynamical__%s.rds", fold_key)))
types <- assays$types
# the forecasting folds' before window: training assays, in before_df's order
stopifnot((experiment_name == "temporal_forecasting") ==
            any(assays$training$before))
# the per-type replicate-based rho every model is scored against
rho_type <- setNames(rho_for_record(tibble(insecticide_type = types),
                                    rho_lookup()), types)


# assemble the types -----------------------------------------------------------

if (step == "assemble") {
  test_df <- assays$test_df
  # draws x held-out assays, filled type by type
  p_out <- matrix(NA_real_, nrow(assays$logit_test), ncol(assays$logit_test))
  map_out <- p_out
  # the before window's types, in before_df's order
  before <- assays$training$before
  before_type <- assays$training$type_id[before]
  map_before <- if (!any(before)) NULL else rep(NA_real_, sum(before))
  experiment_label <- assays$experiment
  rm(assays)
  summaries <- list()
  for (k in seq_along(types)) {
    part <- readRDS(part_file(k))
    test_rows <- which(test_df$type_id == k)
    p_out[, test_rows] <- part$p_draws
    map_out[, test_rows] <- part$map_draws
    map_before[before_type == k] <- part$map_before
    summaries[[k]] <- part$summary
  }
  stopifnot(!anyNA(p_out), !anyNA(map_out), !anyNA(map_before))
  summary_table <- bind_rows(summaries) %>%
    mutate(experiment = experiment_name, fold = fold_name, .before = 1)
  draws_file <- file.path("outputs/cv_draws",
                          sprintf("two_stage__%s.rds", fold_key))
  saveRDS(list(model = "two_stage",
               experiment = experiment_label,
               fold = fold_name,
               test_df = test_df,
               p_draws = p_out,
               map_draws = map_out,
               map_before = map_before,
               rho_type = tibble(insecticide_type = types, rho = rho_type),
               fit_summary = summary_table),
          draws_file)
  report("saved %s", draws_file)

  # one row per experiment x fold x type; a rerun replaces its fold's rows
  summary_file <- file.path(output_dir, "fit_summary.csv")
  if (file.exists(summary_file) && file.size(summary_file) > 0) {
    previous <- read.csv(summary_file, colClasses = c(fold = "character")) %>%
      filter(!(experiment == experiment_name & fold == fold_name))
    summary_table <- bind_rows(previous, summary_table)
  }
  write.csv(summary_table, summary_file, row.names = FALSE)
  report("wrote %s", summary_file)
  quit(save = "no")
}


# Pixel-centre coordinates for rows of `cell` (and `year`), as the published
# map uses (R/two_stage_maps.R)
at_pixel_centres <- function(rows, file = "data/clean/raster_mask.tif") {
  xy <- terra::xyFromCell(terra::rast(file), rows$cell)
  tibble(lon = xy[, 1], lat = xy[, 2], year = rows$year, cell = rows$cell)
}


# one insecticide type ---------------------------------------------------------

source("R/two_stage_correction.R")
k <- as.integer(step)
stopifnot(!is.na(k), k >= 1, k <= length(types))
type <- types[k]
train_rows <- which(assays$training$type_id == k)
test_rows <- which(assays$test_df$type_id == k)
stopifnot(length(train_rows) > 0, length(test_rows) > 0)
m_draws_train <- assays$logit_train[, train_rows, drop = FALSE]
m_draws_test <- assays$logit_test[, test_rows, drop = FALSE]
# the type's before-window assays (forecasting folds only), training assays in
# before_df's order
before_rows <- which(assays$training$type_id == k & assays$training$before)
m_draws_before <- if (length(before_rows) == 0) NULL else
  assays$logit_train[, before_rows, drop = FALSE]
before_k <- if (length(before_rows) == 0) NULL else
  at_pixel_centres(assays$training[before_rows, ])
train_k <- mutate(select(assays$training[train_rows, ], -before),
                  m = colMeans(m_draws_train), rho = rho_type[[type]])
test_k <- with(assays$test_df[test_rows, ],
               tibble(lon = longitude, lat = latitude, year = year,
                      cell = cell, study = study))
# the map is evaluated where the published map is: at the pixel centre, not
# the assay's own coordinates (R/two_stage_maps.R)
test_centres_k <- at_pixel_centres(select(test_k, year, cell))
t0 <- assays$t0
rm(assays)
invisible(gc())

# m_ref is the posterior mean logit at the training assays; T defaults to the
# type's last training year: xi is represented up to it and forecast to later
# test years. The meshes cover the training sites and the prediction mask. The
# seed is the fold's and the type's
set.seed(string_seed(paste(experiment_name, fold_name, type, sep = "__")))
time_fit <- system.time({
  meshes <- build_correction_meshes(coords_km(train_k),
                                    prediction_mask_coords())
  fit <- fit_correction(train_k, t0 = t0, meshes = meshes)
})[["elapsed"]]

# One set of joint draws of the correction fields, mixed over the CCD points
# and each on its own dynamical draw (the cut posterior), projected to: the
# held-out assays' coordinates with fresh noise (the predictive draws); their
# pixel centres without noise (the map); and, for the forecasting folds, the
# before-window assays without noise (the map's posterior mean there, for the
# change score, #36: training assays, so omega + xi are fitted there)
map_before <- NULL
time_predict <- system.time({
  fields <- correction_node_draws(
    fit, unique(c(test_k$year, before_k$year)), nrow(m_draws_train),
    m_draws_train = m_draws_train)
  p_draws <- plogis(t(t(m_draws_test) +
                        project_correction(fit, fields, test_k, noise = TRUE)))
  map_draws <- plogis(t(t(m_draws_test) +
                          project_correction(fit, fields, test_centres_k)))
  if (!is.null(m_draws_before)) {
    map_before <- rowMeans(plogis(t(m_draws_before) +
                                    project_correction(fit, fields, before_k)))
  }
  rm(fields)
  stopifnot(!anyNA(p_draws), !anyNA(map_draws), !anyNA(map_before))
})[["elapsed"]]

summary <- bind_cols(
  tibble(insecticide_type = type, rho = rho_type[[type]],
         n_test = length(test_rows),
         share_test_pixel_in_train = mean(test_k$cell %in% train_k$cell),
         max_test_horizon = max(c(0, test_k$year - fit$T))),
  fit_summary(fit),
  tibble(time_fit_s = time_fit, time_predict_s = time_predict,
         peak_memory_gb = peak_memory_gb())
)
saveRDS(list(summary = summary, p_draws = p_draws, map_draws = map_draws,
             map_before = map_before),
        part_file(k))
report(paste("%s %-18s n=%5i test=%4i range_omega=%.0f sigma_omega=%.2f",
             "range_eta=%.0f sigma_eta=%.3f phi=%.2f tau=%.3f sigma_p=%.2f",
             "sigma_sy=%.2f; fit %.0f s, predict %.0f s, peak %.1f GB"),
       fold_key, type, nrow(train_k), length(test_rows),
       fit$hyper$range_omega, fit$hyper$sigma_omega, fit$hyper$range_eta,
       fit$hyper$sigma_eta, fit$hyper$phi, fit$hyper$tau, fit$hyper$sigma_p,
       fit$hyper$sigma_sy, time_fit, time_predict, peak_memory_gb())
