# Functions for assessing out-of-sample posterior predictive ability of the
# bioassay models.
#
# The quantity these all describe is the posterior predictive distribution for a
# held out bioassay i: a finite mixture, over the S posterior draws from a
# training fold fit, of beta-binomial distributions
#
#   f_i(y) = (1 / S) sum_s dbetabinom(y; size = n_i, p = p_i^(s), rho = rho_i^(s))
#
# where n_i is the number of mosquitoes tested, p_i^(s) is the predicted
# population susceptible fraction under draw s, and rho_i^(s) the corresponding
# draw of the observation overdispersion for that insecticide class.
#
# The beta-binomial functions here are thin wrappers around extraDistr, which
# the rest of the validation code already uses. What they add is the p / rho
# parameterisation used by betabinomial_p_rho() in R/functions.R, rather than
# the shape parameters extraDistr takes; see R/check_validation_functions.R for
# checks that the reparameterisation is right.


# beta-binomial primitives ------------------------------------------------

# convert a mean proportion p and overdispersion rho (on the unit interval, 0
# being no overdispersion) to the beta distribution shape parameters, solving
#   p = a / (a + b), rho = 1 / (a + b + 1)
bb_shape <- function(p, rho, epsilon = 1e-10) {
  p <- pmin(pmax(p, epsilon), 1 - epsilon)
  rho <- pmin(pmax(rho, epsilon), 1 - epsilon)
  total <- 1 / rho - 1
  list(a = p * total,
       b = (1 - p) * total)
}

# beta-binomial probability mass function, in the p / rho parameterisation
dbetabinom <- function(y, size, p, rho, log = FALSE) {
  shape <- bb_shape(p, rho)
  extraDistr::dbbinom(y, size, alpha = shape$a, beta = shape$b, log = log)
}

# beta-binomial cumulative distribution function
pbetabinom <- function(q, size, p, rho) {
  shape <- bb_shape(p, rho)
  extraDistr::pbbinom(q, size, alpha = shape$a, beta = shape$b)
}

# draw beta-binomial variates
rbetabinom <- function(n, size, p, rho) {
  shape <- bb_shape(p, rho)
  extraDistr::rbbinom(n, size, alpha = shape$a, beta = shape$b)
}


# posterior predictive summaries -------------------------------------------

# Core single-pass summary of the posterior predictive distribution at a set of
# held-out observations.
#
#   died, mosquito_number: vectors of length n_obs, the observed data
#   p_draws:   n_draws x n_obs matrix of posterior draws of the population
#              susceptible fraction (e.g. draws of population_mortality_vec_test)
#   rho_draws: n_draws x n_obs matrix of the matching draws of the observation
#              overdispersion, or a vector of length n_draws if a single
#              overdispersion applies to all observations
#
# Returns a data frame with one row per observation, holding everything the
# metrics below need: the log predictive density at the observation, the
# mixture cdf just below the observation and the mixture mass at it (the two
# ingredients of a randomised quantile residual), and the predictive mean.
ppd_summary <- function(died, mosquito_number, p_draws, rho_draws,
                        chunk_size = 1000) {

  n_obs <- length(died)
  stopifnot(
    length(mosquito_number) == n_obs,
    ncol(p_draws) == n_obs
  )

  if (is.null(dim(rho_draws))) {
    rho_draws <- matrix(rho_draws,
                        nrow = nrow(p_draws),
                        ncol = n_obs)
  }
  stopifnot(all(dim(rho_draws) == dim(p_draws)))

  n_draws <- nrow(p_draws)
  cdf_below <- numeric(n_obs)
  pmf_at <- numeric(n_obs)

  # Both ingredients are means over draws of a single extraDistr evaluation:
  # the mixture mass at the observation is dbbinom at y, and the mixture cdf
  # just below it is pbbinom at y - 1, which is 0 at y = 0 as it should be.
  # There is no need to evaluate the pmf over the whole support below y.
  #
  # Observations are taken in chunks purely to bound peak memory at
  # n_draws x chunk_size rather than n_draws x n_obs; the largest held-out
  # fold has ~10,000 records against 2,000 draws.
  for (index in split(seq_len(n_obs), ceiling(seq_len(n_obs) / chunk_size))) {

    shape <- bb_shape(as.vector(p_draws[, index, drop = FALSE]),
                      as.vector(rho_draws[, index, drop = FALSE]))
    y <- rep(died[index], each = n_draws)
    size <- rep(mosquito_number[index], each = n_draws)

    pmf_at[index] <- colMeans(
      matrix(extraDistr::dbbinom(y, size, alpha = shape$a, beta = shape$b),
             nrow = n_draws)
    )
    cdf_below[index] <- colMeans(
      matrix(extraDistr::pbbinom(y - 1, size, alpha = shape$a, beta = shape$b),
             nrow = n_draws)
    )

  }

  data.frame(
    died = died,
    mosquito_number = mosquito_number,
    observed = died / mosquito_number,
    predicted = colMeans(p_draws),
    log_score = log(pmf_at),
    cdf_below = cdf_below,
    pmf_at = pmf_at
  )

}

# simulate replicate observations from the posterior predictive distribution,
# returning an n_draws x n_obs matrix of counts
ppd_simulate <- function(mosquito_number, p_draws, rho_draws) {

  if (is.null(dim(rho_draws))) {
    rho_draws <- matrix(rho_draws,
                        nrow = nrow(p_draws),
                        ncol = ncol(p_draws))
  }

  size <- rep(mosquito_number, each = nrow(p_draws))
  sims <- rbetabinom(length(p_draws),
                     size = size,
                     p = as.vector(p_draws),
                     rho = as.vector(rho_draws))
  matrix(sims, nrow = nrow(p_draws), ncol = ncol(p_draws))

}


# randomised quantile residuals --------------------------------------------

# Randomised probability integral transform values for discrete observations
# (Dunn & Smyth 1996):
#   u_i = F_i(y_i - 1) + v_i f_i(y_i),   v_i ~ U(0, 1)
# which is standard uniform if the predictive distribution is correct. Returns
# an n_obs x n_rep matrix; averaging summaries over the replicates keeps
# conclusions from depending on a single draw of v.
ppd_pit <- function(summary, n_rep = 100) {
  n_obs <- nrow(summary)
  v <- matrix(runif(n_obs * n_rep), nrow = n_obs, ncol = n_rep)
  summary$cdf_below + v * summary$pmf_at
}

# residual z scores, for plotting against covariates and space
pit_to_z <- function(pit) {
  qnorm(pmin(pmax(pit, 1e-10), 1 - 1e-10))
}


# calibration --------------------------------------------------------------

# Empirical coverage of central predictive intervals, at each nominal level.
# Computed from the randomised PIT values, for which an observation falls in
# the central interval of level `level` exactly when the PIT lies within
# ((1 - level) / 2, (1 + level) / 2). Doing it this way avoids the conservatism
# that discreteness introduces into quantile-based intervals, so a calibrated
# model sits on the diagonal.
coverage_curve <- function(pit, levels = seq(0.1, 0.95, by = 0.05)) {
  pit <- as.matrix(pit)
  covered <- vapply(
    levels,
    function(level) {
      lower <- (1 - level) / 2
      upper <- (1 + level) / 2
      mean(pit > lower & pit < upper)
    },
    numeric(1)
  )
  data.frame(nominal = levels,
             empirical = covered)
}

# Cramer-von Mises criterion for deviation from a standard uniform. Expectation
# 0 for a calibrated model; corresponds to the PS2 statistic of Taggart (2022),
# which decomposes into over/under-prediction and over/under-dispersion
cvm_stat <- function(u) {
  u <- sort(u)
  n <- length(u)
  1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n)) ^ 2)
}

# apply a statistic of uniformity across PIT randomisation replicates and
# average, to remove dependence on any single draw of v
pit_statistic <- function(pit, statistic = cvm_stat) {
  pit <- as.matrix(pit)
  mean(apply(pit, 2, statistic))
}

# Null distribution of a uniformity statistic at a given sample size. This
# assumes the PIT values are independent, which holds for the simulated data in
# check_validation_functions.R but not for held-out records scored against a
# shared posterior: there the draws of p are common to every record, so the PIT
# values are dependent and this band is too narrow. It is therefore used only by
# the check suite, and the reported uniformity statistics are read as an
# ordering between models rather than as a test (#12 review)
pit_null_band <- function(n_obs,
                          statistic = cvm_stat,
                          n_sim = 1000,
                          probs = c(0.025, 0.5, 0.975)) {
  sims <- replicate(n_sim, statistic(runif(n_obs)))
  quantile(sims, probs)
}


# scores -------------------------------------------------------------------

# CRPS of the posterior predictive distribution at each observation, on the
# mortality proportion scale. scoringRules takes one row per observation, so
# the draws are transposed
ppd_crps <- function(died, mosquito_number, sims) {
  proportion_sims <- sweep(sims, 2, mosquito_number, FUN = "/")
  scoringRules::crps_sample(y = died / mosquito_number,
                            dat = t(proportion_sims))
}


# the noise floor ----------------------------------------------------------

# The irreducible component of mean squared error on the mortality proportion
# scale. Decomposing the error of a prediction m of an assay with true
# population fraction p:
#
#   E[(y / n - m) ^ 2] = p (1 - p) (1 + (n - 1) rho) / n  +  (p - m) ^ 2
#
# the first term is assay noise that no model can remove. It needs p(1 - p)
# rather than p, and that is recoverable from the data without knowing p:
#
#   E[yhat (1 - yhat)] = p (1 - p) [1 - (1 + (n - 1) rho) / n]
#
# so the observed proportion gives an unbiased estimate of p(1 - p) after
# dividing by the known factor. Averaged over many assays the noise in that
# estimate washes out.
noise_floor_mse <- function(died, mosquito_number, rho) {
  mean(noise_floor_record(died, mosquito_number, rho), na.rm = TRUE)
}

# The per-assay terms noise_floor_mse() averages: each is unbiased for that
# assay's sampling variance but noisy (exactly 0 at 0% or 100%), so only means
# over many assays are meaningful. NA for an assay that carries no information
# about p(1 - p)
noise_floor_record <- function(died, mosquito_number, rho) {
  yhat <- died / mosquito_number
  inflation <- (1 + (mosquito_number - 1) * rho) / mosquito_number
  # An assay of a single mosquito has inflation exactly 1, so the estimator
  # below divides by zero: a single observation carries no information about
  # p(1 - p). Such assays are dropped from the floor rather than allowed to
  # make it undefined; there are a handful in the held-out data.
  usable <- mosquito_number > 1 & is.finite(inflation) & inflation < 1
  out <- rep(NA_real_, length(yhat))
  pq <- pmax(yhat[usable] * (1 - yhat[usable]) / (1 - inflation[usable]), 0)
  out[usable] <- pq * inflation[usable]
  out
}


# the model floor -----------------------------------------------------------

# The two-stage model (#21) treats u (per pixel-year), p (per pixel) and sy
# (per study-year, #55) as noise shared by the assays of a site or a survey:
# an assay samples the population mortality q = plogis(logit(m) + e),
# e ~ N(0, sd^2), around the map value m, with sd^2 = tau^2 + sigma_p^2 +
# sigma_sy^2 from the fit. The model floor U (#32) is the squared error this
# adds to the map's, as the model attributes it, for a level score:
#
#   U(m, sd) = E[(q - m) ^ 2]
#
# The expectation is the plain mean over 1000 evenly spaced normal quantiles,
# within 4e-5 of 60-node Gauss-Hermite quadrature for map values 0.02 to 0.999
# and SDs up to 1.5. Elementwise over map; m of exactly 0 or 1 gives U = 0.
model_floor_mse <- function(map, sd) {
  stopifnot(length(sd) %in% c(1, length(map)))
  e <- qnorm(ppoints(1000))
  sd <- rep_len(sd, length(map))
  deviation <- plogis(outer(qlogis(map), rep(1, length(e))) +
                        outer(sd, e)) - map
  rowMeans(deviation ^ 2)
}

# root mean squared error on the proportion scale. (The definition previously
# used in the validation scripts squared the mean error rather than the errors,
# which measures bias rather than error; see issue #11)
rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted) ^ 2))
}


# reliability --------------------------------------------------------------

# bin observations by predicted mortality and compare the mean prediction with
# the mean observation in each bin. Averaging within a bin estimates the
# population-level quantity with much less noise than any single assay, so this
# is the direct check on whether predictions are right on average
reliability_bins <- function(predicted, observed, n_bins = 10) {
  breaks <- quantile(predicted,
                     probs = seq(0, 1, length.out = n_bins + 1),
                     na.rm = TRUE)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  bin <- cut(predicted, breaks = breaks, labels = FALSE)
  out <- lapply(
    sort(unique(bin)),
    function(b) {
      keep <- bin == b
      data.frame(bin = b,
                 n = sum(keep),
                 predicted = mean(predicted[keep]),
                 observed = mean(observed[keep]))
    }
  )
  do.call(rbind, out)
}

# The scatter a perfect model would still show in a reliability bin, obtained by
# simulation. An analytic envelope - the standard error of the mean of k assays
# of size n at fraction p - was also computed here once, and is gone: it assumes
# the assays in a bin are independent and uses a single assay size for all of
# them; more importantly it conditions on the model's predicted fraction being
# the truth, so it cannot say how much of a reliability gap posterior
# uncertainty in that fraction would produce on its own.
#
# This instead generates replicate held-out datasets from the model's own
# posterior predictive distribution and re-runs the identical binning on each.
# Under those replicates the model is the truth by construction, so the spread
# of the gaps is exactly the gap a correct model would show, with assay noise,
# overdispersion, posterior uncertainty and the finite bin size all included.
# Anything outside it is the model's own error.
#
# This matters for a reason found the hard way: binning instead on the observed
# mortality induces regression to the mean, and a perfectly calibrated model
# then shows an apparent bias of +0.15 in the lowest decile and -0.07 in the
# highest, with none of it real. Reliability must be conditioned on the
# prediction, and even then needs this envelope to be read.
reliability_ppc <- function(predicted, mosquito_number, p_draws, rho,
                            n_bins = 10, n_rep = 200,
                            probs = c(0.025, 0.5, 0.975)) {

  gaps <- function(observed) {
    bins <- reliability_bins(predicted, observed, n_bins = n_bins)
    bins$observed - bins$predicted
  }

  replicates <- t(replicate(n_rep, {
    draw <- sample(nrow(p_draws), 1)
    simulated <- rbetabinom(length(mosquito_number), mosquito_number,
                            p_draws[draw, ], rho)
    gaps(simulated / mosquito_number)
  }))

  quantiles <- t(apply(replicates, 2, quantile, probs = probs))
  colnames(quantiles) <- paste0("ppc_", c("lower", "median", "upper"))
  as.data.frame(quantiles)

}


# aggregation --------------------------------------------------------------

# Posterior predictive distribution for a pooled group of assays. For each
# draw, the group total is the sum of independent beta-binomials with that
# draw's per-assay fractions, so heterogeneity in the true fraction within the
# group is carried by the model's own predictions rather than assumed away.
# This makes the comparison valid at any level of aggregation.
#
# `group` is a vector of group labels, one per observation. Returns a data
# frame with one row per group holding the observed pooled mortality and
# summaries of the pooled predictive distribution.
ppd_aggregate <- function(died, mosquito_number, group, sims,
                          probs = c(0.025, 0.5, 0.975)) {

  groups <- unique(group)
  out <- lapply(
    groups,
    function(g) {
      keep <- group == g
      total_tested <- sum(mosquito_number[keep])
      # pooled mortality under each posterior predictive draw
      pooled_sims <- rowSums(sims[, keep, drop = FALSE]) / total_tested
      pooled_observed <- sum(died[keep]) / total_tested
      quantiles <- unname(quantile(pooled_sims, probs))
      data.frame(group = g,
                 n_assays = sum(keep),
                 n_tested = total_tested,
                 observed = pooled_observed,
                 predicted = mean(pooled_sims),
                 lower = quantiles[1],
                 median = quantiles[2],
                 upper = quantiles[3],
                 # rank of the observation within the predictive draws, a
                 # PIT-like calibration check at the group level
                 pit = mean(pooled_sims < pooled_observed) +
                   runif(1) * mean(pooled_sims == pooled_observed))
    }
  )
  do.call(rbind, out)

}

# bioassay overdispersion per insecticide type ----------------------------

# rho sets the noise floor on every score, so every script has to use the same
# one: the hierarchical fit over replicate bioassay groups, types nested in
# class, fitted by MCMC in fig_illustrate_bioassay_variability.R. There is no
# fallback - scoring against a different overdispersion than the rest of the
# pipeline is worse than not scoring at all.
rho_lookup <- function(file = "outputs/bioassay_rho_hierarchical.csv") {
  if (!file.exists(file)) {
    stop("no per-type overdispersion at ", file,
         "; run R/fig_illustrate_bioassay_variability.R first")
  }
  hierarchical <- utils::read.csv(file)
  if (!all(hierarchical$worst_rhat < 1.05)) {
    stop("the per-type overdispersion fit has not converged (worst Rhat ",
         round(max(hierarchical$worst_rhat), 3), ")")
  }
  list(source = "per insecticide type (hierarchical, MCMC)",
       key = "insecticide_type",
       table = data.frame(key = hierarchical$insecticide_type,
                          rho = hierarchical$rho,
                          rho_lower = hierarchical$rho_lower,
                          rho_upper = hierarchical$rho_upper))
}

# rho per record. `data` must carry insecticide_type.
rho_for_record <- function(data, lookup) {
  # without this the missing column gives match(NULL, ...) -> integer(0), a
  # zero-length result, and the check below passes on nothing at all
  stopifnot(lookup$key %in% names(data))
  index <- match(data[[lookup$key]], lookup$table$key)
  out <- lookup$table$rho[index]
  stopifnot(!any(is.na(out)))
  out
}

# a seed from a string: a polynomial string hash, exact in double precision
string_seed <- function(key) {
  Reduce(function(h, x) (h * 31 + x) %% 2147483647, utf8ToInt(key), 0)
}
