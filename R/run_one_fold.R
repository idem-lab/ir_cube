# Fit the dynamical model to a single cross-validation fold and save the draws.
#
#   Rscript R/run_one_fold.R <experiment> <fold> [n_chains|default] [threads]
#
# e.g. Rscript R/run_one_fold.R spatial_blocks 1 4 4
#      Rscript R/run_one_fold.R temporal_forecasting 2014 4 4
#
# One fold per process, so that greta's warmup and sampling progress goes to
# this fold's own log rather than being buffered out of sight. Dispatched by
# run_validation_folds.R, but can be run directly to redo a single fold.

arguments <- commandArgs(trailingOnly = TRUE)
experiment_name <- arguments[1]
fold_name <- arguments[2]
threads <- if (length(arguments) >= 4) as.integer(arguments[4]) else 4L
# overrides of the sampler settings (dynamical_mcmc_settings(), in
# R/dynamical_model.R): the number of chains, and for smoke-testing the path
# without a real fit, warmup and samples
setting_overrides <- list()
if (length(arguments) >= 3 && arguments[3] != "default") {
  setting_overrides$n_chains <- as.integer(arguments[3])
}
if (length(arguments) >= 5) {
  setting_overrides$warmup <- as.integer(arguments[5])
}
if (length(arguments) >= 6) {
  setting_overrides$n_samples <- as.integer(arguments[6])
}

# greta first: TensorFlow's thread count is fixed once it starts, and python
# has to start before terra and sf are attached (R/greta_setup.R)
source("R/greta_setup.R")
start_greta(threads = threads)

source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/validation_covariates.R")
source("R/fit_validation_fold.R")
settings <- do.call(dynamical_mcmc_settings, setting_overrides)

# find the requested fold
before <- NULL
experiment_label <- experiment_name
stopifnot(experiment_name %in% validation_experiments)
if (experiment_name == "spatial_interpolation") {
  training <- spatial_interpolation$training
  test <- spatial_interpolation$test
} else if (experiment_name == "spatial_blocks") {
  source("R/validation_blocks.R")
  index <- as.integer(fold_name)
  stopifnot(!is.na(index), index >= 1, index <= length(spatial_blocks))
  training <- spatial_blocks[[index]]$training
  test <- spatial_blocks[[index]]$test
} else if (experiment_name == "temporal_forecasting") {
  # `fold` is the cut year: training is everything strictly before it, the
  # holdout is the window starting at it. The origins are defined in
  # validation_folds.R.
  stopifnot(fold_name %in% names(temporal_forecasting_folds))
  fold <- temporal_forecasting_folds[[fold_name]]
  training <- fold$training
  test <- fold$test
  # predictions in the before window are needed for the change-based score, and
  # can only be requested at fitting time
  before <- fold$before
  # Each origin is scored as its own experiment: pooling a 2014 forecast with
  # a 2018 one would average over holdout windows whose true rates of decline
  # differ by a factor of two. The file name keeps the plain experiment name,
  # with the origin in the fold field.
  experiment_label <- paste0("temporal_forecasting_", fold$cut_year)
} else {
  stop("unknown experiment: ", experiment_name)
}

# CV_DRAWS_DIR lets a smoke test write somewhere harmless, so a five-sample
# wiring check cannot be mistaken for a fold and skipped later
draws_dir <- Sys.getenv("CV_DRAWS_DIR", "outputs/cv_draws")
destination <- file.path(draws_dir,
                         sprintf("dynamical__%s__%s.rds",
                                 experiment_name, fold_name))

# resume by file: a fold already on disk is not refitted, so the dispatcher can
# be restarted without checking anything itself
if (file.exists(destination)) {
  cat(sprintf("%s | %s / %s already fitted, skipping\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              experiment_name, fold_name))
  cat("FOLD COMPLETE\n")
  quit(save = "no", status = 0)
}

cat(sprintf("%s | %s / %s | %i training, %i held out | %i chains, %i threads\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            experiment_name, fold_name, nrow(training), nrow(test),
            settings$n_chains, threads))
flush(stdout())

elapsed <- system.time(
  fit <- fit_fold(
    train_df = training,
    test_df = test,
    before_df = before,
    x_cell_years = x_cell_years,
    cell_years_index = cell_years_index,
    df = df,
    classes_index = classes_index,
    types = types,
    # dynamical_model_options(), set in validation_covariates.R
    options = model_options,
    x_cells_init = x_cells_init,
    settings = settings
  )
)

cat(sprintf("%s | fit complete in %.1f hours\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            elapsed[["elapsed"]] / 3600))

dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

# Everything fit_fold() returns, under the three labels that identify the fold.
# That includes the draws object and the greta arrays the predictions came from,
# so that calculate() can be used on a reloaded fold to predict a quantity that
# was not asked for at fitting time, though the scoring reads only the matrices
# and they are the bulk of each file's size.
saveRDS(
  c(list(model = "dynamical",
         experiment = experiment_label,
         fold = fold_name),
    fit),
  destination
)

cat(sprintf("%s | saved. p ESS median %.0f, rho ESS median %.0f, worst Rhat %.3f\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            median(fit$ess_p, na.rm = TRUE),
            median(fit$ess_rho, na.rm = TRUE),
            max(fit$convergence[, 1], na.rm = TRUE)))
cat("FOLD COMPLETE\n")
