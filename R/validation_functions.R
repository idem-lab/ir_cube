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

# probability mass of a beta-binomial over the values in `y`, for each of the
# draws implied by vectors of shape parameters `a` and `b`. Returns a matrix
# with one row per draw and one column per element of `y`
bb_pmf_matrix <- function(y, size, a, b) {
  n_draws <- length(a)
  n_y <- length(y)
  density <- extraDistr::dbbinom(
    x = rep(y, each = n_draws),
    size = size,
    alpha = rep(a, times = n_y),
    beta = rep(b, times = n_y)
  )
  matrix(density, nrow = n_draws, ncol = n_y)
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
ppd_summary <- function(died, mosquito_number, p_draws, rho_draws) {

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

  log_score <- numeric(n_obs)
  cdf_below <- numeric(n_obs)
  pmf_at <- numeric(n_obs)

  # loop over observations, vectorising over draws and over the support below
  # each observation. This keeps memory to n_draws x (died + 1) at a time
  for (i in seq_len(n_obs)) {

    shape <- bb_shape(p_draws[, i], rho_draws[, i])
    support <- 0:died[i]
    # mixture pmf over 0:died[i], averaging the per-draw pmfs
    mixture <- colMeans(
      bb_pmf_matrix(support, mosquito_number[i], shape$a, shape$b)
    )

    pmf_at[i] <- mixture[length(mixture)]
    cdf_below[i] <- sum(mixture) - pmf_at[i]
    log_score[i] <- log(pmf_at[i])

  }

  data.frame(
    died = died,
    mosquito_number = mosquito_number,
    observed = died / mosquito_number,
    predicted = colMeans(p_draws),
    log_score = log_score,
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

# Kolmogorov-Smirnov statistic for deviation from a standard uniform
ks_stat <- function(u) {
  n <- length(u)
  u <- sort(u) - (0:(n - 1)) / n
  max(c(u, 1 / n - u))
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

# null distribution of a uniformity statistic at a given sample size, for
# drawing reference bands
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
  yhat <- died / mosquito_number
  inflation <- (1 + (mosquito_number - 1) * rho) / mosquito_number
  pq <- pmax(yhat * (1 - yhat) / (1 - inflation), 0)
  mean(pq * inflation)
}

# expected CRPS of an oracle that knows the true population fraction, which is
# half the mean absolute difference between two draws from the assay
# distribution
noise_floor_crps <- function(mosquito_number, p, rho, n_sim = 1000) {
  n <- length(mosquito_number)
  p <- rep_len(p, n)
  rho <- rep_len(rho, n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    x <- rbetabinom(n_sim, mosquito_number[i], p[i], rho[i]) / mosquito_number[i]
    x_prime <- rbetabinom(n_sim, mosquito_number[i], p[i], rho[i]) / mosquito_number[i]
    out[i] <- mean(abs(x - x_prime)) / 2
  }
  mean(out)
}

# proportion of the explainable error that a model removes: 0 for the null
# model, 1 at the noise floor
mse_skill <- function(mse_model, mse_null, mse_floor) {
  (mse_null - mse_model) / (mse_null - mse_floor)
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

# the scatter a perfect model would still show in a reliability bin: the
# standard error of the mean of k assays of size n at fraction p, with
# overdispersion rho. Drawn as an envelope around the 1:1 line, this separates
# measurement noise from model error, and narrows visibly as k grows
reliability_envelope <- function(p, mosquito_number, k, rho) {
  variance <- p * (1 - p) * (1 + (mosquito_number - 1) * rho) / mosquito_number
  sqrt(variance / k)
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
