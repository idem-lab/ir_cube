# Fit the dynamical model to one cross-validation training fold and return
# posterior predictive draws for the held-out data.
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
                     stored_draws = 2000,
                     inits_file = dynamical_inits_file) {

  # the model, with the likelihood over the training fold (R/dynamical_model.R)
  built <- build_dynamical_model(train_df = train_df,
                                 df = df,
                                 x_cell_years = x_cell_years,
                                 cell_years_index = cell_years_index,
                                 classes_index = classes_index,
                                 types = types,
                                 options = options,
                                 x_cells_init = x_cells_init)

  # use cached posterior means as inits, and the sampler settings
  # (R/dynamical_model.R)
  inits_one <- dynamical_inits(readRDS(inits_file), built$variables,
                               levels = built$lookups$levels,
                               columns = colnames(x_cell_years))
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
  # consume, rather than of the raw model parameters
  ess_prediction <- coda::effectiveSize(prediction_draws)
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
       ess = ess,
       ess_p = ess_p,
       ess_rho = ess_rho,
       n_sampled = sampled,
       n_chains = settings$n_chains,
       settings = settings)

}

