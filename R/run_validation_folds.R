# Fit the dynamical model to each cross-validation training fold and save
# posterior predictive draws for the held-out data.
#
# This replaces the earlier point-estimate extraction
# (#10). The MCMC is unchanged; what is saved is the posterior draws of the
# predicted population fraction and the observation overdispersion, rather than
# their means, so that held-out data can be scored against the full posterior
# predictive distribution. Scoring and plotting are separate, in
# validation_metrics.R and fig_predictive_validation.R, so metrics can be
# revised without refitting.
#
# The null models are cheap and are run here too, so that every candidate is
# stored in the same format and scored by the same code.

# This process fits the null models and dispatches the dynamical folds as
# separate processes; it never calls fit_fold() itself. So it needs neither
# greta's python session nor the covariate extraction, both of which used to be
# paid for here for nothing. run_one_fold.R does that setup, in the order it
# has to be done in (#12 review).
source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/null_models.R")

draws_dir <- "outputs/cv_draws"
dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

n_draws <- 1000

# store one fold of one model in a consistent format: the draws, and the
# held-out records they correspond to
save_fold <- function(fit, model, experiment, fold, label = experiment) {
  object <- list(
    model = model,
    # the experiment this fold is pooled under when scored, which for the
    # forecasting folds is the origin rather than the bare experiment name; the
    # file keeps the plain name so all the forecasting folds sit together
    experiment = label,
    fold = fold,
    p_draws = fit$p_draws,
    # no fitted overdispersion: every model is scored at the external
    # replicate-based estimate, and a null that fits its own only buys coverage
    # by being vague (#12 review)
    test_df = fit$test_df
  )
  file <- file.path(draws_dir,
                    sprintf("%s__%s__%s.rds", model, experiment, fold))
  saveRDS(object, file)
  cat(sprintf("saved %s (%i held-out assays)\n", file, nrow(fit$test_df)))
  invisible(file)
}

# The nearest neighbour null is no longer tuned. It used to read a neighbour
# count per experiment from outputs/optimal_nn.csv, chosen by grid search on an
# internal test set of 100 records sampled at random from the training data -
# which sat a median 0 km from their nearest usable neighbour, against 36 to 194
# km for the records actually held out, so it chose too few neighbours and
# handicapped the baseline. The lookup also had no row for the block folds,
# where it returned numeric(0) and reduced the null to Beta(0.5, 0.5) without
# erroring. Both problems go away with the table (#12).
#
# In its place, two specifications, neither of them a chosen value: the practice
# baseline at one neighbour and the most recent two available years, and the
# oracle bound at whichever neighbour count minimises the null's own error on
# the held-out records. See nn_null_draws() and nn_oracle_draws().

# The folds to fit: the sub-national block folds (the test of spatial skill),
# the interpolation fold, and the two five-year forecast origins. The six
# leave-one-country-out folds and the 2020 three-year forecasting fold are
# defunct and are not refitted. Nulls are rebuilt for every fold either way,
# since they cost minutes.
source("R/validation_blocks.R")

forecast_folds <- temporal_forecasting_folds[as.character(forecast_cuts)]

folds <- c(
  lapply(
    seq_along(spatial_blocks),
    function(i) list(experiment = "spatial_blocks",
                     fold = as.character(i),
                     label = "spatial_blocks",
                     training = spatial_blocks[[i]]$training,
                     test = spatial_blocks[[i]]$test)
  ),
  list(
    list(experiment = "spatial_interpolation",
         fold = "all",
         label = "spatial_interpolation",
         training = spatial_interpolation$training,
         test = spatial_interpolation$test)
  ),
  # The rolling forecast origins. The nulls need no special handling for the
  # horizon any more: their year window is anchored at prediction time, so a
  # forecasting fold reads the last two training years for every held-out
  # record, however far ahead it sits.
  lapply(
    names(forecast_folds),
    function(name) {
      fold <- forecast_folds[[name]]
      list(experiment = "temporal_forecasting",
           label = paste0("temporal_forecasting_", fold$cut_year),
           fold = name,
           training = fold$training,
           test = fold$test)
    }
  )
)


# null models --------------------------------------------------------------

for (fold in folds) {

  null_file <- function(model) {
    file.path(draws_dir,
              sprintf("%s__%s__%s.rds", model, fold$experiment, fold$fold))
  }

  if (!file.exists(null_file("intercept"))) {
    save_fold(
      intercept_null_draws(fold$training, fold$test, n_draws = n_draws),
      model = "intercept",
      experiment = fold$experiment,
      fold = fold$fold,
      label = fold$label
    )
  }

  if (!file.exists(null_file("nearest_neighbour"))) {
    save_fold(
      nn_null_draws(fold$training,
                    fold$test,
                    n_neighbours = 1,
                    n_years_prior = 1,
                    n_draws = n_draws),
      model = "nearest_neighbour",
      experiment = fold$experiment,
      fold = fold$fold,
      label = fold$label
    )
  }

  if (file.exists(null_file("nearest_neighbour_oracle"))) next

  oracle <- nn_oracle_draws(fold$training,
                            fold$test,
                            n_years_prior = 1,
                            n_draws = n_draws)
  cat(sprintf("  oracle neighbour count for %s / %s: %i\n",
              fold$experiment, fold$fold, oracle$n_neighbours))
  save_fold(
    oracle,
    model = "nearest_neighbour_oracle",
    experiment = fold$experiment,
    fold = fold$fold,
    label = fold$label
  )

}


# dynamical model ----------------------------------------------------------

# Each fold is fitted by its own process, via run_one_fold.R, so that greta's
# warmup and sampling progress goes to that fold's log and can be watched while
# it runs. Two folds at a time, four chains each.
#
# Why this split: TensorFlow vectorises chains into a single op rather than
# running one per core, and that op scales poorly beyond about four threads, so
# a fold confined to four threads loses little. Four chains rather than two
# because greta pools information across chains when adapting during warmup, and
# two chains adapt poorly — the Kenya fold reached Rhat 7.6 that way.
#
# Folds whose draws are already on disk are skipped, so the run resumes.
#
# Two folds at a time, four threads each. With the closed-form recursion (#25)
# and the sampler settings of doc/cv_run_plan.md section 3, a fold at four
# threads takes about 4 h (2014 forecasting fold, 3.0 s per iteration) to 6 h
# (interpolation fold), against 62 h with the greta.dynamics loop. Each uses
# about 3.5 GB; check free memory before raising n_concurrent.
# The number of chains and the rest of the sampler settings are
# dynamical_mcmc_settings() (R/dynamical_model.R), which run_one_fold.R uses
# when passed "default".
n_concurrent <- 2
threads_per_fold <- 4

log_dir <- "outputs/cv_logs"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

pending <- Filter(
  function(fold) {
    !file.exists(file.path(
      draws_dir,
      sprintf("dynamical__%s__%s.rds", fold$experiment, fold$fold)))
  },
  folds
)

cat(sprintf("\n%i folds to fit, %i at a time, %i threads each\n",
            length(pending), n_concurrent, threads_per_fold))
cat(sprintf("progress logs: %s/\n\n", log_dir))

running <- list()

launch <- function(fold) {
  log_file <- file.path(log_dir,
                        sprintf("%s__%s.log", fold$experiment, fold$fold))
  cat(sprintf("%s | launching %s / %s -> %s\n",
              format(Sys.time(), "%H:%M:%S"),
              fold$experiment, fold$fold, log_file))
  process <- processx::process$new(
    "Rscript",
    c("R/run_one_fold.R", fold$experiment, fold$fold,
      "default", as.character(threads_per_fold)),
    stdout = log_file,
    stderr = "2>&1"
  )
  list(fold = fold, process = process, log = log_file)
}

queue <- pending

while (length(queue) > 0 || length(running) > 0) {

  # start jobs while there is room
  while (length(running) < n_concurrent && length(queue) > 0) {
    running <- c(running, list(launch(queue[[1]])))
    queue <- queue[-1]
  }

  Sys.sleep(60)

  # reap anything that has finished
  still_running <- list()
  for (job in running) {
    if (job$process$is_alive()) {
      still_running <- c(still_running, list(job))
    } else {
      status <- job$process$get_exit_status()
      cat(sprintf("%s | finished %s / %s (exit %s)\n",
                  format(Sys.time(), "%H:%M:%S"),
                  job$fold$experiment, job$fold$fold, status))
      if (!identical(status, 0L)) {
        cat(sprintf("  NON-ZERO EXIT: see %s\n", job$log))
      }
    }
  }
  running <- still_running

}

cat("\nall folds finished\n")
