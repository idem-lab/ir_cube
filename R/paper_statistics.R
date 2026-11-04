# Numbers for the paper from the final run's outputs (#69), printed to the
# log; nothing is written.
#
#   Rscript R/paper_statistics.R
#
# Run after R/fit_model.R, R/fig_illustrate_bioassay_variability.R,
# R/two_stage_maps.R and R/validation_metrics.R. Prints:
#   1. the full fit's sampler settings, and its convergence over every element
#      of its draws (the chains kept, as `draws`): rank-normalised Rhat, and
#      bulk and tail ESS (posterior);
#   2. the external overdispersion rho per insecticide, with the replicate-pair
#      check of it (R/fig_illustrate_bioassay_variability.R);
#   3. the second-stage hyperparameters per insecticide: the mode, and the
#      posterior mean and 95% interval from the CCD points (below);
#   4. cross-validation coverage of the central 50% and 95% predictive
#      intervals and mean PIT, per model and test (outputs/cv_summary.csv).
# The mechanistic model's parameters are summarised by R/summarise_posterior.R,
# in-sample calibration by R/in_sample_calibration.R, and the data by
# R/summarise_data.R.

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(posterior)
})
options(width = 200)

# 1. sampler and convergence -----------------------------------------------------

fit <- new.env()
load("temporary/fitted_model.RData", envir = fit)
cat("1. full fit: sampler settings (dynamical_mcmc_settings())\n")
str(fit$settings)
draws_array <- as_draws_array(coda::as.mcmc.list(
  lapply(fit$draws, function(x) coda::mcmc(as.matrix(x)))))
convergence <- summarise_draws(draws_array, rhat, ess_bulk, ess_tail)
cat(sprintf(paste("chains kept %i of %i; %i elements; log density Rhat %.3f;",
                  "max Rhat %.3f (%s); min bulk ESS %.0f (%s);",
                  "min tail ESS %.0f (%s)\n"),
            length(fit$draws), nrow(fit$chain_check$chains),
            nrow(convergence), fit$chain_check$lp_rhat,
            max(convergence$rhat), convergence$variable[which.max(convergence$rhat)],
            min(convergence$ess_bulk),
            convergence$variable[which.min(convergence$ess_bulk)],
            min(convergence$ess_tail),
            convergence$variable[which.min(convergence$ess_tail)]))
rm(fit, draws_array)
invisible(gc())


# 2. rho -----------------------------------------------------------------------------

cat("\n2. external overdispersion rho per insecticide (posterior mean and 95%",
    "interval), and its replicate-pair check (ratio of the pair-based noise",
    "to rho's, with its 95% interval)\n")
rho <- read.csv("outputs/bioassay_rho_hierarchical.csv")
pairs <- read.csv("outputs/bioassay_rho_pair_check.csv")
print(rho %>%
        select(insecticide_type, rho, rho_lower, rho_upper, n_assays,
               worst_rhat) %>%
        left_join(select(pairs, insecticide_type, pair_assays = assays,
                         ratio, ratio_lower = lower, ratio_upper = upper),
                  by = "insecticide_type") %>%
        mutate(across(where(is.double), ~ signif(.x, 3))))


# 3. second-stage hyperparameters ------------------------------------------------

# The CCD points (ccd_points(), R/two_stage_correction.R) and their weights
# reproduce the posterior's mean and covariance of theta (the optimiser's
# scale) if it is Gaussian. The 95% interval is that of the Gaussian with the
# weighted points' mean and SD, mapped to the hyperparameter's scale; the
# posterior mean is the weighted mean of the hyperparameter at the points (as
# the _ccd columns of hyperparameters.csv)
hyper_names <- c("sigma_omega", "range_omega", "sigma_eta", "range_eta",
                 "phi", "tau", "sigma_p", "sigma_sy")
hyper_transform <- list(exp, exp, exp, exp, plogis, exp, exp, exp)
types <- read.csv("outputs/two_stage/maps/hyperparameters.csv")$insecticide_type
hyper <- bind_rows(lapply(types, function(type) {
  fit <- readRDS(file.path("outputs/two_stage/maps", type, "fit.rds"))
  theta <- fit$ccd$theta
  weight <- fit$ccd$weight
  stopifnot(ncol(theta) == length(hyper_names))
  mean_theta <- colSums(theta * weight)
  sd_theta <- sqrt(colSums(sweep(theta, 2, mean_theta) ^ 2 * weight))
  tibble(insecticide_type = type,
         hyperparameter = hyper_names,
         mode = unlist(fit$hyper[hyper_names]),
         mean = colSums(as.matrix(fit$ccd$hyper[hyper_names]) * weight),
         lower = mapply(function(f, m, s) f(m - qnorm(0.975) * s),
                        hyper_transform, mean_theta, sd_theta),
         upper = mapply(function(f, m, s) f(m + qnorm(0.975) * s),
                        hyper_transform, mean_theta, sd_theta))
}))
cat("\n3. second-stage hyperparameters (ranges in km; phi the AR(1)",
    "coefficient of xi's increments)\n")
print(hyper %>% mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf)


# 4. cross-validation calibration ------------------------------------------------

cat("\n4. cross-validation: coverage of the central 50% and 95% predictive",
    "intervals, and mean PIT, per model and test\n")
print(read.csv("outputs/cv_summary.csv") %>%
        select(experiment, model, n, coverage_50, coverage_95, mean_pit) %>%
        arrange(experiment, model) %>%
        mutate(across(where(is.double), ~ signif(.x, 3))))
