# Fit the dynamical model to one cross-validation training fold and return
# posterior predictive draws for the held-out data; fold_assay_draws() adds
# the logit draws at the training and held-out assays, for the two-stage model.
#
# The model is build_dynamical_model() in R/dynamical_model.R, which fit_model.R
# also uses: the likelihood here is restricted to the training fold. Rather than
# collapsing the posterior to a mean predicted fraction and a mean
# overdispersion, this returns the draws themselves, so that the held-out data
# can be scored against the full posterior predictive distribution.
#
# The predicted fraction and the overdispersion are drawn in a single
# calculate() call, so that each pair comes from the same posterior sample.
#
# Because the greta arrays are local to this function they go out of scope when
# it returns, so the model does not need to be purged by hand between folds.

# dynamical_logit() and paired_draw_index(), also for fold_assay_draws(), and
# check_chains() (R/chain_checks.R)
source("R/dynamical_predictions.R")

fit_fold <- function(train_df,
                     test_df,
                     before_df = NULL,
                     x_cell_years,
                     cell_years_index,
                     df,
                     classes_index,
                     types,
                     options = dynamical_model_options(),
                     x_cells_init = NULL,
                     settings = dynamical_mcmc_settings(),
                     stored_draws = 2000) {

  # the model, with the likelihood over the training fold (R/dynamical_model.R)
  built <- build_dynamical_model(train_df = train_df,
                                 df = df,
                                 x_cell_years = x_cell_years,
                                 cell_years_index = cell_years_index,
                                 classes_index = classes_index,
                                 types = types,
                                 options = options,
                                 x_cells_init = x_cells_init)

  # start from the initial values in git (dynamical_inits_file), with the
  # sampler settings (R/dynamical_model.R)
  inits_one <- dynamical_inits(readRDS(dynamical_inits_file), built)
  draws <- run_dynamical_mcmc(built$model, built$variables, inits_one,
                              settings)

  # a fixed-length run, not one topped up to an ESS target: folds compared
  # with each other must share their sampling settings (#12 review)
  report <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...), "\n")
    flush(stdout())
  }

  sampled <- settings$n_samples
  ess <- coda::effectiveSize(draws)
  report("sampled %d per chain | raw parameter ESS min %.0f median %.0f",
         sampled, min(ess, na.rm = TRUE), median(ess, na.rm = TRUE))

  # whether the chains agree (R/chain_checks.R)
  chain_check <- check_chains(draws, chain_log_density(built$model, draws),
                              built$variable_labels)
  for (line in chain_check_lines(chain_check)) report("%s", line)

  # Predictions at the held-out data, by calculate() on the draws without
  # nsim, which keeps the MCMC order (nsim resamples), so their effective
  # sample size can be measured.
  population_mortality_vec_test <- built$mortality(test_df)
  # the overdispersion of each insecticide type
  rho_types <- built$terms$rho_types

  # Optionally, predictions at a second set of cell-years: the forecasting
  # experiment is scored on the change in mortality between the window before
  # the cut and the holdout window (#12 review 5.1).
  report("computing predictions at %d held-out assays%s", nrow(test_df),
         if (is.null(before_df)) "" else
           sprintf(" and %d before-window records", nrow(before_df)))
  prediction_draws <- calculate(
    population_mortality_vec_test = population_mortality_vec_test,
    rho_types = rho_types,
    values = draws
  )

  # effective sample size of the quantities the validation metrics actually
  # consume, rather than of the raw model parameters, over the chains that
  # R/run_one_fold.R keeps (drop_fold_chains() drops the stuck ones; the
  # predictions are only stored thinned, so it cannot recompute this)
  moving <- setdiff(seq_along(draws), stuck_chains(draws))
  if (length(moving) == 0) moving <- seq_along(draws)
  ess_prediction <- coda::effectiveSize(prediction_draws[moving])
  ess_p <- ess_prediction[grep("population_mortality_vec_test\\[",
                               names(ess_prediction))]
  ess_rho <- ess_prediction[grep("rho_types", names(ess_prediction))]

  report("prediction ESS: p median %.0f min %.0f | rho median %.0f min %.0f",
         median(ess_p, na.rm = TRUE), min(ess_p, na.rm = TRUE),
         median(ess_rho, na.rm = TRUE), min(ess_rho, na.rm = TRUE))

  # flatten the mcmc.list to a draws x quantity matrix, preserving order
  prediction_matrix <- as.matrix(prediction_draws)
  p_columns <- grep("population_mortality_vec_test\\[",
                    colnames(prediction_matrix))
  rho_columns <- grep("rho_types", colnames(prediction_matrix))
  p_draws <- prediction_matrix[, p_columns, drop = FALSE]
  rho_type_draws <- prediction_matrix[, rho_columns, drop = FALSE]
  rm(prediction_matrix, prediction_draws)
  invisible(gc())

  # The before-window predictions come from a second calculate() on the same
  # draws, to save memory. Without nsim calculate() is deterministic in the
  # draws, so draw i here is paired with draw i above.
  p_draws_before <- NULL
  if (!is.null(before_df) && nrow(before_df) > 0) {
    population_mortality_vec_before <- built$mortality(before_df)
    before_draws <- calculate(
      population_mortality_vec_before = population_mortality_vec_before,
      values = draws
    )
    p_draws_before <- as.matrix(before_draws)
    rm(before_draws)
    invisible(gc())
  }

  # Thin the stored prediction draws, on one set of indices so that p, rho and
  # the before-window p stay paired. The scoring thins to 2,000 draws anyway;
  # ESS is measured above, on the unthinned draws.
  keep_draws <- if (nrow(p_draws) > stored_draws) {
    round(seq(1, nrow(p_draws), length.out = stored_draws))
  } else {
    seq_len(nrow(p_draws))
  }
  p_draws <- p_draws[keep_draws, , drop = FALSE]
  rho_type_draws <- rho_type_draws[keep_draws, , drop = FALSE]
  if (!is.null(p_draws_before)) {
    p_draws_before <- p_draws_before[keep_draws, , drop = FALSE]
  }

  convergence <- coda::gelman.diag(draws,
                                   multivariate = FALSE,
                                   autoburnin = FALSE)$psrf
  report("Rhat worst %.3f, %d of %d parameters above 1.01",
         max(convergence[, 1], na.rm = TRUE),
         sum(convergence[, 1] > 1.01, na.rm = TRUE),
         nrow(convergence))

  list(# The draws object, and the greta arrays the predictions were computed
       # from, so calculate(values = draws) can predict more from a reloaded
       # fold. Sampling cannot be continued after the session ends
       # (extra_samples() fails on a reloaded draws object).
       draws = draws,
       prediction_arrays = list(p = population_mortality_vec_test,
                                rho = rho_types),
       p_draws = p_draws,
       # the overdispersion is shared by every assay of an insecticide type, so
       # it is stored by type with the index needed to expand it, rather than
       # as one column per held-out assay
       rho_type_draws = rho_type_draws,
       type_id = test_df$type_id,
       options = built$options,
       x_cells_init = x_cells_init,
       test_df = test_df,
       p_draws_before = p_draws_before,
       before_df = before_df,
       convergence = convergence,
       chain_check = chain_check,
       ess = ess,
       ess_p = ess_p,
       ess_rho = ess_rho,
       n_sampled = sampled,
       n_chains = settings$n_chains,
       settings = settings)

}

# The dynamical model's logit draws at a fold's training and held-out assays,
# which the two-stage model (R/run_two_stage_folds.R) pairs its draws with,
# from fit_fold()'s result `fit` and its training set `train_df`. At the draws
# p_draws was thinned to, by dynamical_logit() (R/dynamical_predictions.R):
# on the logit scale the predictions do not round to 0 or 1, as p_draws does
# beyond |logit| ~ 37. The before window (forecasting folds) is the training
# assays in `before_years`, flagged `before`, in the order of fit$before_df.
# Returns
#   training     the training assays: type_id, lon, lat, year, cell, died,
#                mosquito_number, before
#   test_df      fit$test_df
#   logit_train  draws x training assays
#   logit_test   draws x held-out assays
#   max_difference  max |plogis(logit) - p| against p_draws and
#                p_draws_before; stops above 1e-6. On the probability scale,
#                since qlogis(p) loses precision as p nears 0 or 1
fold_assay_draws <- function(fit, train_df, before_years = NULL, df,
                             x_cell_years, cell_years_index, classes_index,
                             types) {

  # the draws p_draws holds, less any chains in fit$stuck_chains
  index <- paired_draw_index(fit)
  stopifnot(length(index) == nrow(fit$p_draws))
  parameters <- dynamical_parameter_draws(fit, classes_index, types, index)
  logit_train <- dynamical_logit(parameters, train_df, df, x_cell_years,
                                 cell_years_index)
  logit_test <- dynamical_logit(parameters, fit$test_df, df, x_cell_years,
                                cell_years_index)
  rm(parameters)

  training <- transmute(train_df, type_id, lon = longitude, lat = latitude,
                        year, cell, died, mosquito_number, study,
                        before = year %in% before_years)
  before_rows <- train_df[training$before, ]
  stopifnot(nrow(before_rows) == NROW(fit$before_df),
            identical(paste(before_rows$cell_id, before_rows$type_id,
                            before_rows$year_id),
                      paste(fit$before_df$cell_id, fit$before_df$type_id,
                            fit$before_df$year_id)))

  # the recursion over long forecast windows accumulates rounding (4.9e-8 on
  # the logit scale on temporal_forecasting 2014)
  difference <- max(
    abs(plogis(logit_test) - fit$p_draws),
    if (!is.null(fit$p_draws_before)) {
      abs(plogis(logit_train[, training$before, drop = FALSE]) -
            fit$p_draws_before)
    })
  cat(sprintf(paste("%s | assay draws: %i x (%i training + %i held out),",
                    "max |difference| from p_draws %.3g\n"),
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(index),
              ncol(logit_train), ncol(logit_test), difference))
  stopifnot(difference < 1e-6)

  list(training = training,
       test_df = fit$test_df,
       logit_train = logit_train,
       logit_test = logit_test,
       max_difference = difference)
}

