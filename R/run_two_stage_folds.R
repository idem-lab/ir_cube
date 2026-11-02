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
#              held-out pixel-year) and p (per held-out pixel), never their
#              fitted values; and the same draws without u and p, the map
#              plogis(m + omega + xi), whose posterior mean is what the maps
#              show and what the point metrics score (#32); and, for the
#              forecasting folds, the map's posterior mean at the type's
#              before-window assays (training assays in the window before the
#              cut), for the change score (#36). Saved to
#              outputs/two_stage/parts/;
#   assemble   combine the types into
#              outputs/cv_draws/two_stage__<experiment>__<fold>.rds and the
#              fold's rows of outputs/two_stage/fit_summary.csv. A type with
#              no fit keeps the dynamical draws, which are also its map, and
#              is listed in fallback_types.
#
# The saved object has the #12 fields (model, experiment, fold, test_df,
# p_draws) plus map_draws (draws x held-out assays of the map), rho_type, the
# per-type rho of the fit (a data frame of insecticide_type and rho, from the
# table rho_lookup() scores against), fallback_types and fit_summary (which
# holds each type's tau and sigma_p, the model floor's inputs). The
# forecasting folds also carry map_before, the map's posterior mean at each
# assay of the dynamical fold's before_df, in its order (the dynamical
# posterior mean for fallback types). Parts saved before map_before was added
# are an error: rerun them.
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
  p_out <- plogis(assays$logit_test)
  # the dynamical model has no u or p, so its draws are also its map
  map_out <- p_out
  # the before window's types, in before_df's order
  before <- assays$training$before
  before_type <- assays$training$type_id[before]
  map_before <- if (!any(before)) NULL else
    colMeans(plogis(assays$logit_train[, before, drop = FALSE]))
  experiment_label <- assays$experiment
  rm(assays)
  summaries <- list()
  for (k in seq_along(types)) {
    part <- if (file.exists(part_file(k))) readRDS(part_file(k)) else
      list(summary = tibble(insecticide_type = types[k],
                            error = "no part file"))
    test_rows <- which(test_df$type_id == k)
    if (!is.null(part$p_draws)) {
      stopifnot(!is.null(part$map_draws))
      p_out[, test_rows] <- part$p_draws
      map_out[, test_rows] <- part$map_draws
      before_rows <- which(before_type == k)
      if (length(before_rows) > 0) {
        stopifnot(!is.null(part$map_before))
        map_before[before_rows] <- part$map_before
      }
    }
    summaries[[k]] <- mutate(part$summary, fallback_dynamical =
                               length(test_rows) > 0 && is.null(part$p_draws))
  }
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
               fallback_types = summary_table$insecticide_type[
                 summary_table$fallback_dynamical],
               fit_summary = summary_table),
          draws_file)
  report("saved %s (fallback: %s)", draws_file,
         paste(summary_table$insecticide_type[summary_table$fallback_dynamical],
               collapse = ", "))

  # one row per experiment x fold x type; a rerun replaces its fold's rows
  summary_file <- file.path(output_dir, "fit_summary.csv")
  if (file.exists(summary_file)) {
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
                      cell = cell))
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
error <- NA_character_
fit <- NULL
time_fit <- system.time(
  if (length(train_rows) == 0) {
    error <- "no training assays"
  } else {
    meshes <- build_correction_meshes(coords_km(train_k),
                                      prediction_mask_coords())
    fit <- tryCatch(fit_correction(train_k, t0 = t0, meshes = meshes),
                    error = function(e) {
                      error <<- conditionMessage(e)
                      NULL
                    })
  }
)[["elapsed"]]

p_draws <- NULL
map_draws <- NULL
time_predict <- system.time(
  if (!is.null(fit) && fit$opt$convergence == 0 && length(test_rows) > 0) {
    # held-out predictive draws at the assays' coordinates, then the map at
    # their pixel centres in a second call, so the predictive draws are
    # unchanged by it
    lambda <- predict_correction(fit, test_k, m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = nrow(m_draws_train))
    stopifnot(!anyNA(lambda))
    p_draws <- plogis(lambda)
    lambda <- predict_correction(fit, test_centres_k,
                                 m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = nrow(m_draws_train), noise = FALSE)
    stopifnot(!anyNA(lambda))
    map_draws <- plogis(lambda)
    rm(lambda)
  }
)[["elapsed"]]

# The map's posterior mean at the before-window assays, for the change score
# (#36): training assays, so omega + xi are fitted there. Drawn after the
# held-out draws, so those are unchanged
map_before <- NULL
if (!is.null(p_draws) && !is.null(m_draws_before)) {
  map_before <- colMeans(plogis(predict_correction(
    fit, before_k, m_draws_train = m_draws_train,
    m_draws_new = m_draws_before, n_draws = nrow(m_draws_train),
    noise = FALSE)))
  stopifnot(!anyNA(map_before))
}

summary <- bind_cols(
  tibble(insecticide_type = type, rho = rho_type[[type]],
         n_test = length(test_rows),
         share_test_pixel_in_train = mean(test_k$cell %in% train_k$cell),
         max_test_horizon = if (is.null(fit)) NA else
           max(c(0, test_k$year - fit$T))),
  if (is.null(fit)) NULL else fit_summary(fit),
  tibble(error = error, time_fit_s = time_fit, time_predict_s = time_predict,
         peak_memory_gb = peak_memory_gb())
)
saveRDS(list(summary = summary, p_draws = p_draws, map_draws = map_draws,
             map_before = map_before),
        part_file(k))
report("%s %-18s n=%5i test=%4i %s; fit %.0f s, predict %.0f s, peak %.1f GB",
       fold_key, type, nrow(train_k), length(test_rows),
       if (is.null(fit)) paste("FAILED:", error) else
         sprintf(paste("conv=%i range_omega=%.0f sigma_omega=%.2f",
                       "range_eta=%.0f sigma_eta=%.3f phi=%.2f tau=%.3f",
                       "sigma_p=%.2f refit=%s"),
                 fit$opt$convergence, fit$hyper$range_omega,
                 fit$hyper$sigma_omega, fit$hyper$range_eta,
                 fit$hyper$sigma_eta, fit$hyper$phi, fit$hyper$tau,
                 fit$hyper$sigma_p, fit$stage_b$refit),
       time_fit, time_predict, peak_memory_gb())
