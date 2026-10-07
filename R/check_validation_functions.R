# Checks on R/validation_functions.R, run against simulated data where the
# truth is known. These verify that the beta-binomial primitives are correct,
# and that the calibration metrics behave as intended: flat for a model whose
# predictive distribution is the data-generating one, and deviating in the
# expected direction for models that are biased, overconfident or
# underconfident.

source("R/validation_functions.R")

set.seed(2026 - 8 - 31)

report <- function(label, pass, detail = "") {
  cat(sprintf("%-52s %s %s\n", label, ifelse(pass, "PASS", "FAIL"), detail))
  invisible(pass)
}


# beta-binomial primitives ------------------------------------------------

size <- 40
p <- 0.7
rho <- 0.12
support <- 0:size
pmf <- dbetabinom(support, size, p, rho)

report("pmf sums to one", isTRUE(all.equal(sum(pmf), 1)))
report("pmf mean matches size * p",
       isTRUE(all.equal(sum(support * pmf), size * p)))

# variance of a beta-binomial is n p (1 - p) [1 + (n - 1) rho]
variance_expected <- size * p * (1 - p) * (1 + (size - 1) * rho)
variance_observed <- sum((support - size * p) ^ 2 * pmf)
report("pmf variance matches n p q [1 + (n-1) rho]",
       isTRUE(all.equal(variance_observed, variance_expected)))

# as rho tends to zero the beta-binomial tends to the binomial
report("binomial limit as rho -> 0",
       isTRUE(all.equal(dbetabinom(support, size, p, 1e-9),
                        dbinom(support, size, p),
                        tolerance = 1e-5)))

report("cdf reaches one", isTRUE(all.equal(pbetabinom(size, size, p, rho), 1)))
report("cdf matches cumulative pmf",
       isTRUE(all.equal(pbetabinom(0:size, size, p, rho), cumsum(pmf))))

# sampled moments match
draws <- rbetabinom(2e5, size, p, rho)
report("rbetabinom mean", abs(mean(draws) - size * p) < 0.05,
       sprintf("(%.3f vs %.3f)", mean(draws), size * p))
report("rbetabinom variance",
       abs(var(draws) / variance_expected - 1) < 0.05,
       sprintf("(%.1f vs %.1f)", var(draws), variance_expected))

# ppd_crps against a brute force evaluation of the CRPS double sum
sims_test <- matrix(rbinom(300, 100, 0.4), ncol = 1)
observed_test <- 45
brute <- mean(abs(sims_test / 100 - observed_test / 100)) -
  0.5 * mean(abs(outer(sims_test / 100, sims_test / 100, "-")))
report("ppd_crps matches brute force double sum",
       isTRUE(all.equal(as.numeric(ppd_crps(observed_test, 100, sims_test)),
                        brute, tolerance = 1e-6)))


# behaviour of the calibration metrics ------------------------------------

n_obs <- 500
n_draws <- 300
rho_true <- 0.12
mosquito_number <- rep(100, n_obs)

# true population fractions, spanning the range seen in the data
p_true <- rbeta(n_obs, 6, 2)
died <- rbetabinom(n_obs, mosquito_number, p_true, rho_true)

# four candidate models, expressed as posterior draws. The first has the
# data-generating distribution as its predictive distribution
make_draws <- function(p, rho) {
  list(p = matrix(p, nrow = n_draws, ncol = n_obs, byrow = TRUE),
       rho = matrix(rho, nrow = n_draws, ncol = n_obs))
}

candidates <- list(
  calibrated = make_draws(p_true, rho_true),
  biased = make_draws(pmin(p_true + 0.08, 0.999), rho_true),
  overconfident = make_draws(p_true, rho_true / 4),
  underconfident = make_draws(p_true, rho_true * 4)
)

results <- lapply(
  candidates,
  function(candidate) {
    summary <- ppd_summary(died, mosquito_number, candidate$p, candidate$rho)
    pit <- ppd_pit(summary, n_rep = 20)
    coverage <- coverage_curve(pit, levels = c(0.5, 0.95))
    data.frame(
      mean_pit = mean(pit),
      coverage_50 = coverage$empirical[1],
      coverage_95 = coverage$empirical[2],
      cvm = pit_statistic(pit),
      elpd = mean(summary$log_score),
      mse = mean((summary$observed - summary$predicted) ^ 2)
    )
  }
)
results <- do.call(rbind, results)

cat("\n")
print(round(results, 4))
cat("\n")

null_band <- pit_null_band(n_obs, n_sim = 500)

report("calibrated: mean PIT near 0.5",
       abs(results["calibrated", "mean_pit"] - 0.5) < 0.03)
report("calibrated: 95% coverage near nominal",
       abs(results["calibrated", "coverage_95"] - 0.95) < 0.03)
report("calibrated: 50% coverage near nominal",
       abs(results["calibrated", "coverage_50"] - 0.5) < 0.06)
report("calibrated: CvM within null band",
       results["calibrated", "cvm"] < null_band[3],
       sprintf("(%.3f vs upper %.3f)", results["calibrated", "cvm"], null_band[3]))
report("biased: mean PIT shifted",
       abs(results["biased", "mean_pit"] - 0.5) > 0.05)
report("biased: CvM outside null band",
       results["biased", "cvm"] > null_band[3])
report("overconfident: coverage below nominal",
       results["overconfident", "coverage_95"] < 0.9)
report("underconfident: coverage above nominal",
       results["underconfident", "coverage_95"] > 0.98)
report("calibrated model has the best elpd",
       which.max(results$elpd) == 1)

# dispersion errors are invisible to a point summary: the biased model is the
# only one a mean squared error comparison can distinguish, which is the reason
# for scoring the distribution rather than the point
report("MSE cannot separate the dispersion errors",
       isTRUE(all.equal(results["overconfident", "mse"],
                        results["underconfident", "mse"],
                        tolerance = 1e-8)))


# the randomised PIT where the predictive distribution has a large atom --------

# The distortion a discrete predictive distribution can introduce is driven by
# the size of the atom at the observation, so it is invisible in the simulation
# above, where the fractions span the range and no single count is very likely.
# It bites in the real held-out data, where about 30% of assays record 100%
# mortality, so the randomised PIT is checked here on data generated with
# fractions close to one. (The pipeline stores one randomisation replicate; the
# mean over replicates is the mid-P value, which is not a PIT, and the figures
# once read it as one - #12 review.)
n_atom <- 4000
p_atom <- rbeta(n_atom, 20, 1.2)
size_atom <- rep(100, n_atom)
died_atom <- rbetabinom(n_atom, size_atom, p_atom, rho_true)

pit_atom <- ppd_pit(
  ppd_summary(died_atom, size_atom,
              matrix(p_atom, nrow = 200, ncol = n_atom, byrow = TRUE),
              matrix(rho_true, nrow = 200, ncol = n_atom)),
  n_rep = 100
)
randomised <- pit_atom[, 1]

cat(sprintf("\n  %.0f%% of these assays are at 100%% mortality\n",
            100 * mean(died_atom == size_atom)))

report("randomised PIT gives nominal coverage",
       abs(mean(randomised > 0.025 & randomised < 0.975) - 0.95) < 0.02,
       sprintf("(%.3f)", mean(randomised > 0.025 & randomised < 0.975)))
report("randomised PIT is uniform under a large atom",
       abs(var(randomised) - 1 / 12) < 0.005,
       sprintf("(variance %.4f, uniform is %.4f)", var(randomised), 1 / 12))


# the noise floor ----------------------------------------------------------

# with predictions equal to the truth, the mean squared error is entirely
# irreducible, so the floor estimator should recover it
mse_at_truth <- mean((died / mosquito_number - p_true) ^ 2)
floor_estimate <- noise_floor_mse(died, mosquito_number, rho_true)
report("noise floor recovers irreducible MSE",
       abs(floor_estimate / mse_at_truth - 1) < 0.1,
       sprintf("(%.5f vs %.5f)", floor_estimate, mse_at_truth))

# the pooled noise floor ---------------------------------------------------

# the single-assay estimator is the one-assay case of the pooled one
one_assay <- noise_floor_var_pooled(died[1], mosquito_number[1], rho_true)
report("pooled floor reduces to the single-assay floor",
       abs(one_assay - noise_floor_mse(died[1], mosquito_number[1], rho_true)) <
         1e-12)

# and over many replicate groups it should recover the variance of the pooled
# proportion around the fraction that generated them
group_size <- 8
n_group <- 3000
p_group <- rbeta(n_group, 6, 2)
group_died <- matrix(NA_real_, n_group, group_size)
for (j in seq_len(group_size)) {
  group_died[, j] <- rbetabinom(n_group, rep(100, n_group), p_group, rho_true)
}
pooled_observed <- rowSums(group_died) / (group_size * 100)
pooled_floor <- mean(vapply(
  seq_len(n_group),
  function(i) noise_floor_var_pooled(group_died[i, ], rep(100, group_size),
                                     rho_true),
  numeric(1)
))
report("pooled floor recovers the variance of a pooled proportion",
       abs(pooled_floor / mean((pooled_observed - p_group) ^ 2) - 1) < 0.1,
       sprintf("(%.5f vs %.5f)", pooled_floor,
               mean((pooled_observed - p_group) ^ 2)))

report("pooling lowers the floor",
       pooled_floor < mean(vapply(
         seq_len(n_group),
         function(i) noise_floor_var_pooled(group_died[i, 1], 100, rho_true),
         numeric(1))),
       sprintf("(%.5f pooled over %i assays)", pooled_floor, group_size))


# reliability and aggregation ----------------------------------------------

sims <- ppd_simulate(mosquito_number, candidates$calibrated$p,
                     candidates$calibrated$rho)

reliability <- reliability_bins(colMeans(candidates$calibrated$p),
                                died / mosquito_number)
report("calibrated model is reliable across bins",
       max(abs(reliability$predicted - reliability$observed)) < 0.05,
       sprintf("(max deviation %.3f)",
               max(abs(reliability$predicted - reliability$observed))))

# pooled groups of assays: the group-level PIT should also be uniform
group <- rep(seq_len(n_obs / 5), each = 5)
aggregated <- ppd_aggregate(died, mosquito_number, group, sims)
report("aggregated predictive is calibrated",
       abs(mean(aggregated$pit) - 0.5) < 0.06,
       sprintf("(mean group PIT %.3f)", mean(aggregated$pit)))

# and pooling should shrink the scatter around the truth
scatter_single <- sd(died / mosquito_number - p_true)
scatter_pooled <- sd(aggregated$observed - aggregated$predicted)
report("pooling reduces scatter",
       scatter_pooled < scatter_single,
       sprintf("(%.4f vs %.4f, ratio %.2f)",
               scatter_pooled, scatter_single, scatter_pooled / scatter_single))


# the model floor (#32) ----------------------------------------------------

# the per-record floor terms average to the floor
report("per-record floor terms average to the floor",
       isTRUE(all.equal(mean(noise_floor_record(died, mosquito_number,
                                                rho_true), na.rm = TRUE),
                        noise_floor_mse(died, mosquito_number, rho_true))))

# U and the pull against Monte Carlo, at map values across the range and the
# site noise SDs the September two-stage fits found (combined 0.2 to 0.85)
map_values <- c(0.02, 0.2, 0.5, 0.8, 0.98)
noise_sds <- c(0.2, 0.5, 0.85)
grid <- expand.grid(map = map_values, sd = noise_sds)
e <- rnorm(1e6)
monte_carlo <- t(mapply(function(m, s) {
  deviation <- plogis(qlogis(m) + s * e) - m
  c(pull = mean(deviation), mse = mean(deviation ^ 2))
}, grid$map, grid$sd))
quadrature <- site_noise_moments(grid$map, grid$sd)
report("model floor matches Monte Carlo",
       max(abs(quadrature$mse - monte_carlo[, "mse"])) < 1e-4,
       sprintf("(max |diff| %.1e)",
               max(abs(quadrature$mse - monte_carlo[, "mse"]))))
report("pull towards 0.5 matches Monte Carlo",
       max(abs(quadrature$pull - monte_carlo[, "pull"])) < 1e-4,
       sprintf("(max |diff| %.1e; %.4f at 80%%, sd 0.5)",
               max(abs(quadrature$pull - monte_carlo[, "pull"])),
               quadrature$pull[grid$map == 0.8 & grid$sd == 0.5]))
report("pull is towards 0.5, and zero at 0.5",
       all(sign(quadrature$pull[grid$map != 0.5]) ==
             sign(0.5 - grid$map[grid$map != 0.5])) &&
         max(abs(quadrature$pull[grid$map == 0.5])) < 1e-12)
report("model floor combines tau and sigma_p in quadrature",
       isTRUE(all.equal(model_floor_mse(0.8, 0.3, 0.4),
                        site_noise_moments(0.8, 0.5)$mse)))

# the change floor: independent u at the two pixel-years
change <- (plogis(qlogis(0.8) + 0.5 * e) - 0.8) -
  (plogis(qlogis(0.4) + 0.5 * rev(e)) - 0.4)
report("change model floor matches Monte Carlo",
       abs(model_floor_change_mse(0.4, 0.8, 0.5) - mean(change ^ 2)) < 1e-4,
       sprintf("(%.5f vs %.5f)", model_floor_change_mse(0.4, 0.8, 0.5),
               mean(change ^ 2)))
