# Fit the final two-stage model (#21, R/two_stage_correction.R) inside one
# cross-validation fold of #12, per insecticide type, and save its held-out
# predictive draws in the #12 fold format, so they are scored by the same code.
#
#   Rscript R/run_two_stage_folds.R <experiment> <fold> stage_one
#   Rscript R/run_two_stage_folds.R <experiment> <fold> <type index>
#   Rscript R/run_two_stage_folds.R <experiment> <fold> assemble
#
# e.g. Rscript R/run_two_stage_folds.R spatial_interpolation all stage_one
#      Rscript R/run_two_stage_folds.R spatial_interpolation all 3
#
# Three steps, so that the memory-heavy one runs alone and the per-type fits
# can run in parallel (one process each):
#
#   stage_one  rebuild the fold's training and test sets exactly as
#              run_one_fold.R did and check the test set against the saved
#              dynamical fold's; recompute the dynamical model's 2000 paired
#              posterior logit draws at the training and test assays
#              (R/dynamical_predictions.R), checking the test draws against
#              the saved p_draws, and, for the forecasting folds, at the
#              before-window assays (training assays in the window before the
#              cut, the fold's before_df), checked against p_draws_before;
#              save them to
#              outputs/two_stage/stage_one__<experiment>__<fold>.rds. Loading a
#              saved fold takes 4-8 GB;
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
#              before-window assays, for the change score (#36). Saved to
#              outputs/two_stage/parts/;
#   assemble   combine the types into
#              outputs/cv_draws/two_stage__<experiment>__<fold>.rds and the
#              fold's rows of outputs/two_stage/fit_summary.csv. A type with
#              no fit keeps the dynamical draws, which are also its map, and
#              is listed in fallback_types.
#
# The saved object has the #12 fields (model, experiment, fold, test_df,
# p_draws) plus map_draws (draws x held-out assays of the map, NA in the
# columns of a type whose part predates #32 and so has none), rho_type, the
# per-type rho of the fit (a data frame of insecticide_type and rho, from the
# table rho_lookup() scores against), fallback_types and fit_summary (which
# holds each type's tau and sigma_p, the model floor's inputs). The
# forecasting folds also carry before_df (the dynamical fold's, row for row)
# and map_before, the map's posterior mean at each of its assays (NA for
# types whose parts have none; the dynamical posterior mean for fallback
# types); both are absent where the stage-one cache predates them.
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
cache_file <- file.path(output_dir, sprintf("stage_one__%s.rds", fold_key))
part_file <- function(k) file.path(parts_dir, sprintf("%s__%i.rds", fold_key, k))


# stage one: the fold and its paired dynamical draws ---------------------------

if (step == "stage_one") {

  wait_for_memory()
  fold <- readRDS(file.path("outputs/cv_draws",
                            sprintf("dynamical__%s.rds", fold_key)))
  # the design matrix is rebuilt with the fold's own model options
  # (validation_covariates.R), so that it matches the draws
  model_options <- fold$options

  suppressMessages({
    sink("/dev/null")
    source("R/validation_folds.R")
    source("R/validation_covariates.R")
    if (experiment_name == "spatial_blocks") source("R/validation_blocks.R")
    sink()
  })
  source("R/dynamical_predictions.R")

  fold_sets <- switch(
    experiment_name,
    spatial_interpolation = spatial_interpolation,
    spatial_blocks = spatial_blocks[[as.integer(fold_name)]],
    temporal_forecasting = temporal_forecasting_folds[[fold_name]]
  )
  training <- fold_sets$training
  test <- fold_sets$test
  stopifnot(!is.null(training), !is.null(test))
  # the before window, for the change score (#36): forecasting folds only
  before_df <- fold$before_df
  stopifnot(is.null(before_df) == is.null(fold$p_draws_before))
  if (!is.null(before_df)) {
    before <- fold_sets$before
    stopifnot(nrow(before) == nrow(before_df),
              identical(as.integer(before$cell_id),
                        as.integer(before_df$cell_id)),
              identical(as.integer(before$type_id),
                        as.integer(before_df$type_id)),
              identical(as.integer(before$year_id),
                        as.integer(before_df$year_id)),
              identical(as.numeric(before$died), as.numeric(before_df$died)))
  }
  report("%s: %i training, %i held-out assays", fold_key, nrow(training),
         nrow(test))

  # the rebuilt held-out set must be the fold's, row for row, or the draws
  # cannot be paired
  test_df <- fold$test_df
  stopifnot(
    nrow(test) == nrow(test_df),
    identical(as.integer(test$cell_id), as.integer(test_df$cell_id)),
    identical(as.integer(test$type_id), as.integer(test_df$type_id)),
    identical(as.integer(test$year_id), as.integer(test_df$year_id)),
    identical(as.numeric(test$died), as.numeric(test_df$died)),
    identical(as.numeric(test$mosquito_number),
              as.numeric(test_df$mosquito_number))
  )

  # the stored test draws, thinned by the rule the pairing uses (thin_draws())
  draw_index <- paired_draw_index(fold)
  p_saved <- thin_draws(fold$p_draws)
  stopifnot(nrow(p_saved) == length(draw_index))
  p_before_saved <- if (is.null(before_df)) NULL else
    thin_draws(fold$p_draws_before)
  parameters <- dynamical_parameter_draws(fold, classes_index, types,
                                          draw_index)
  experiment_label <- fold$experiment
  rm(fold)
  invisible(gc())

  logit_test <- dynamical_logit(parameters, test_df, df, x_cell_years,
                                cell_years_index)
  logit_train <- dynamical_logit(parameters, training, df, x_cell_years,
                                 cell_years_index)
  logit_before <- if (is.null(before_df)) NULL else
    dynamical_logit(parameters, before_df, df, x_cell_years, cell_years_index)
  rm(parameters)

  # The recomputed test draws must reproduce the saved ones, compared on the
  # logit scale, where the model is additive. The saved draws are
  # probabilities, which round to exactly 0 or 1 beyond |logit| ~ 37, so both
  # sides are truncated at +-27.6 (p = 1e-12) for the comparison. 1e-6
  # because the recursion over long forecast windows accumulates rounding
  # (4.9e-8 on temporal_forecasting 2014)
  truncate <- function(x) pmin(pmax(x, -27.6), 27.6)
  max_logit_diff <- max(abs(truncate(logit_test) - truncate(qlogis(p_saved))))
  report("recomputed vs saved test draws: max |logit diff| = %.3g",
         max_logit_diff)
  stopifnot(max_logit_diff < 1e-6)
  if (!is.null(before_df)) {
    max_before_diff <- max(abs(truncate(logit_before) -
                                 truncate(qlogis(p_before_saved))))
    report("recomputed vs saved before-window draws: max |logit diff| = %.3g",
           max_before_diff)
    stopifnot(max_before_diff < 1e-6)
  }

  training <- transmute(training, type_id, lon = longitude, lat = latitude,
                        year = year_start, cell, died, mosquito_number)
  saveRDS(list(experiment_label = experiment_label, types = types,
               t0 = baseline_year, training = training, test_df = test_df,
               p_saved = p_saved, logit_train = logit_train,
               logit_test = logit_test, before_df = before_df,
               p_before_saved = p_before_saved, logit_before = logit_before),
          cache_file)
  report("saved %s; peak memory %.1f GB", cache_file, peak_memory_gb())
  quit(save = "no")
}

suppressMessages(library(dplyr))
cache <- readRDS(cache_file)
types <- cache$types
# the per-type replicate-based rho every model is scored against
rho_type <- setNames(rho_for_record(tibble(insecticide_type = types),
                                    rho_lookup()), types)


# assemble the types -----------------------------------------------------------

if (step == "assemble") {
  test_df <- cache$test_df
  p_out <- cache$p_saved
  # the dynamical model has no u or p, so its draws are also its map
  map_out <- p_out
  before_df <- cache$before_df
  map_before <- if (is.null(before_df)) NULL else
    colMeans(cache$p_before_saved)
  summaries <- list()
  for (k in seq_along(types)) {
    part <- if (file.exists(part_file(k))) readRDS(part_file(k)) else
      list(summary = tibble(insecticide_type = types[k],
                            error = "no part file"))
    test_rows <- which(test_df$type_id == k)
    if (!is.null(part$p_draws)) {
      p_out[, test_rows] <- part$p_draws
      map_out[, test_rows] <- if (is.null(part$map_draws)) NA_real_ else
        part$map_draws
      if (!is.null(before_df)) {
        before_rows <- which(before_df$type_id == k)
        map_before[before_rows] <- if (is.null(part$map_before)) NA_real_ else
          part$map_before
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
               experiment = cache$experiment_label,
               fold = fold_name,
               test_df = test_df,
               p_draws = p_out,
               map_draws = map_out,
               before_df = before_df,
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
train_rows <- which(cache$training$type_id == k)
test_rows <- which(cache$test_df$type_id == k)
m_draws_train <- cache$logit_train[, train_rows, drop = FALSE]
m_draws_test <- cache$logit_test[, test_rows, drop = FALSE]
# the type's before-window assays (forecasting folds), where the cache has
# them
before_rows <- if (is.null(cache$logit_before)) integer(0) else
  which(cache$before_df$type_id == k)
m_draws_before <- if (length(before_rows) == 0) NULL else
  cache$logit_before[, before_rows, drop = FALSE]
before_k <- if (length(before_rows) == 0) NULL else
  at_pixel_centres(with(cache$before_df[before_rows, ],
                        tibble(year = year_start, cell = cell)))
train_k <- mutate(cache$training[train_rows, ], m = colMeans(m_draws_train),
                  rho = rho_type[[type]])
test_k <- with(cache$test_df[test_rows, ],
               tibble(lon = longitude, lat = latitude, year = year_start,
                      cell = cell))
# the map is evaluated where the published map is: at the pixel centre, not
# the assay's own coordinates (R/two_stage_maps.R)
test_centres_k <- at_pixel_centres(select(test_k, year, cell))
t0 <- cache$t0
rm(cache)
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
    # unchanged by it; its noise draws are discarded
    lambda <- predict_correction(fit, test_k, m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = nrow(m_draws_train))
    stopifnot(!anyNA(lambda))
    p_draws <- plogis(lambda)
    lambda <- predict_correction(fit, test_centres_k,
                                 m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = nrow(m_draws_train), map = TRUE)
    stopifnot(!anyNA(lambda$map))
    map_draws <- plogis(lambda$map)
    rm(lambda)
  }
)[["elapsed"]]

# The map's posterior mean at the before-window assays, for the change score
# (#36): training assays, so omega + xi are fitted there. Drawn after the
# held-out draws, so those are unchanged; the noise draws it also makes are
# discarded
map_before <- NULL
if (!is.null(p_draws) && !is.null(m_draws_before)) {
  map_before <- colMeans(plogis(predict_correction(
    fit, before_k, m_draws_train = m_draws_train,
    m_draws_new = m_draws_before, n_draws = nrow(m_draws_train),
    map = TRUE)$map))
  stopifnot(!anyNA(map_before))
}

# the pull of the predictive mean towards 0.5 that averaging fresh u and p
# after the inverse logit causes: predictive mean less map, per held-out assay
pull <- if (is.null(p_draws)) NA else colMeans(p_draws) - colMeans(map_draws)

summary <- bind_cols(
  tibble(insecticide_type = type, rho = rho_type[[type]],
         n_test = length(test_rows),
         share_test_pixel_in_train = mean(test_k$cell %in% train_k$cell),
         max_test_horizon = if (is.null(fit)) NA else
           max(c(0, test_k$year - fit$T))),
  if (is.null(fit)) NULL else fit_summary(fit),
  tibble(mean_pull = mean(pull), max_abs_pull = max(abs(pull)),
         error = error, time_fit_s = time_fit, time_predict_s = time_predict,
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
