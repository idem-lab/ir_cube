# evaluate validation metrics against bias and error for different true
# population-level susceptibility fractions
rbetabinom <- function(n, probability, mosquito_number = 100, rho = 0.14) {

  # reparameterise from prediction and overdispersion ot the beta parameters
  a <- probability * (1 / rho - 1)
  b <- a * (1 - probability) / probability
  
  extraDistr::rbbinom(n = n,
                      size = mosquito_number,
                      alpha = a,
                      beta = b)
  
}

betabinom_dev <- function(died, mosquito_number, predicted, rho = 0.14) {
  
  # reparameterise from prediction and overdispersion ot the beta parameters
  a <- predicted * (1 / rho - 1)
  b <- a * (1 - predicted) / predicted
  
  log_probs <- extraDistr::dbbinom(x = died,
                                   size = mosquito_number,
                                   alpha = a,
                                   beta = b,
                                   log = TRUE)
  -2 * sum(log_probs)
}

# binomial deviance
binom_dev <- function(died, mosquito_number, predicted) {
  
  log_probs <- dbinom(x = died,
                      size = mosquito_number,
                      prob = predicted,
                      log = TRUE)
  -2 * sum(log_probs)
}

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted) ^ 2))
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}

crps <- function(died, mosquito_number, predicted, rho = 0.14, n_sims = 1000) {
  
  ratio_observed <- died / mosquito_number
  n <- length(ratio_observed)
  
  if (length(mosquito_number) == 1) {
    mosquito_number <- rep(mosquito_number, n)
  }
  
  if (length(predicted) == 1) {
    predicted <- rep(predicted, n)
  }
  
  if (length(rho) == 1) {
    rho <- rep(rho, n)
  }
  
  # reparameterise from prediction and overdispersion to the beta parameters
  a <- predicted * (1 / rho - 1)
  b <- a * (1 - predicted) / predicted
  
  # run sims for each observation
  sim_bbinom <- function(size, alpha, beta) {
    extraDistr::rbbinom(
      n = n_sims,
      size = size,
      alpha = alpha,
      beta = beta
    )    
  }
  # simulate ratios  
  t_died_sims <- mapply(sim_bbinom,
                      size = mosquito_number,
                      alpha = a,
                      beta = b)
  died_sims <- t(t_died_sims)
  ratio_sims <- sweep(died_sims,
                      2,
                      mosquito_number,
                      FUN = "/")
  
  # reshape into matrix (rows are sims for one observation)
  scoringRules::crps_sample(y = ratio_observed,
                            dat = ratio_sims)
}

mean_crps <- function(died, mosquito_number, predicted, rho = 0.14, n_sims = 1000) {
  mean(
    crps(died = died,
         mosquito_number = mosquito_number,
         predicted = predicted,
         rho = rho,
         n_sims = n_sims)
  )
}

# CDF of the betabinomial distribution
pbetabinom <- function(q, size, prob, rho = 0.14) {
  
  # reparameterise from prediction and overdispersion ot the beta parameters
  a <- prob * (1 / rho - 1)
  b <- a * (1 - prob) / prob
  
  # evaluate the CDF  
  extraDistr::pbbinom(q = q,
                      size = size,
                      alpha = a,
                      beta = b)  
}

# Kolmogorov-Smornov statistic for a standard uniform distribution
ks_stat <- function(u) {
  n <- length(u)
  u <- sort(u) - (0:(n - 1))/n
  max(c(u, 1/n - u))
}

# Kolmogorov-Smirnov test statistic of deviation from the betabinomial distribution
ks_pred <- function(died, mosquito_number, probability, rho = 0.14) {
  # compute the CDF of the distribution at these data, given the prediction and
  # the data-assumed rho
  u <- pbetabinom(q = died,
                  size = mosquito_number,
                  prob = probability,
                  rho = rho)
  # compute and return the Kolmogorov-Smirnov statistic for deviation of this
  # from a standard uniform
  ks_stat(u)
}

# Cramer von Mises test criterion
cvm_stat <- function(u) {
  u <- sort(u)
  n <- length(u)
  1 / (12 * n) + sum((u - (2 * seq_len(n) - 1) / (2 * n) ) ^ 2)
}

# Cramer von Mises test criterion of deviation from the betabinomial
# distribution
cvm_pred <- function(died, mosquito_number, probability, rho = 0.14) {
  # compute the CDF of the distribution at these data, given the prediction and
  # the data-assumed rho
  u <- pbetabinom(q = died,
                  size = mosquito_number,
                  prob = probability,
                  rho = rho)
  # compute and return the Cramer von Mises criterion for deviation of this from
  # a standard uniform
  cvm_stat(u)
}

library(extraDistr)
library(tidyverse)



# Compute metrics against each:
# betabinomial deviance, RMSE, bias
epsilon <- 0.001

error_amount <- 0.1
bias <- error_amount
# this is the sd of a half normal that has expectation error_amount
noise_sd <- error_amount / sqrt(2 / pi)

# set a vector of true prevalences
data <- expand_grid(
  true_fraction = seq(error_amount, 1 - error_amount, length.out = 100),
  observation = 1:100,
  scenario = c("truth", "noisy", "biased")
) %>%
  # simulate betabinomial data from each
  mutate(
    mosquito_number = 100,
    died = rbetabinom(n = n(),
                      probability = true_fraction,
                      mosquito_number = mosquito_number)
  ) %>%
  mutate(
    # add error to some predictions
    error = case_when(
      scenario == "truth" ~ 0,
      scenario == "biased" ~ bias,
      scenario == "noisy" ~ rnorm(n(), 0, noise_sd)
    ),
    # clamp predicted probabilities to (0, 1)
    prediction = pmax(epsilon,
                      pmin(1 - epsilon,
                           true_fraction + error)),
  )

# Compute metrics against each:
# betabinomial deviance, RMSE, bias
# and plot against true prevalence

metrics <- data %>%
  mutate(
    observed_fraction = died / mosquito_number
  ) %>%
  group_by(
    true_fraction, scenario
  ) %>%
  summarise(
    rmse = rmse(observed_fraction, prediction),
    crps = mean_crps(died = died,
                     mosquito_number = mosquito_number,
                     predicted = prediction),
    deviance = betabinom_dev(died = died,
                             mosquito_number = mosquito_number,
                             predicted = prediction),
    deviance0 = binom_dev(died = died,
                          mosquito_number = mosquito_number,
                          predicted = prediction),
    ks = ks_pred(died = died,
                 mosquito_number = mosquito_number,
                 probability = prediction),
    cvm = cvm_pred(died = died,
                 mosquito_number = mosquito_number,
                 probability = prediction),
    bias = mean(prediction - observed_fraction),
    .groups = "drop"
  )

# plot these against true prevalence
metrics %>%
  pivot_longer(
    cols = all_of(c("rmse", "deviance", "deviance0", "bias", "crps", "ks", "cvm")),
    names_to = "metric",
    values_to = "value"
  ) %>%
  ggplot(
    aes(x = true_fraction,
        y = value)
  ) +
  geom_line() +
  facet_grid(metric ~ scenario,
             # ncol = 5,
             scales = "free_y") +
  coord_cartesian(
    xlim = c(0, 1)
  ) +
  theme_minimal()

# try computing ratios of scores against that of the true prevalence
relative_metrics <- metrics %>%
  pivot_longer(
    cols = all_of(c("rmse", "deviance", "deviance0", "bias", "crps", "ks", "cvm")),
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = value
  ) %>%
  mutate(
    biased = biased / mean(truth),
    noisy = noisy / mean(truth),
    truth = truth / mean(truth)
  ) %>%
  pivot_longer(
    cols = all_of(c("biased", "noisy", "truth")),
    names_to = "scenario",
    values_to = "value"
  )

relative_metrics %>%
  ggplot(
    aes(x = true_fraction,
        y = value,
        colour = scenario)
  ) +
  geom_line() +
  facet_wrap(~metric,
             scales = "free_y") +
  theme_minimal()

relative_metrics %>%
  group_by(metric, scenario) %>%
  summarise(
    mean = mean(value),
    sd = sd(value)
  )



# try computing ratios of scores against that of the true prevalence
relative_metrics_pointwise <- metrics %>%
  pivot_longer(
    cols = all_of(c("rmse", "deviance", "deviance0", "bias", "crps", "ks", "cvm")),
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = value
  ) %>%
  mutate(
    biased = biased / truth,
    noisy = noisy / truth,
    truth = truth / truth
  ) %>%
  pivot_longer(
    cols = all_of(c("biased", "noisy", "truth")),
    names_to = "scenario",
    values_to = "value"
  )

relative_metrics_pointwise %>%
  ggplot(
    aes(x = true_fraction,
        y = value,
        colour = scenario)
  ) +
  geom_line() +
  facet_wrap(~metric,
             scales = "free_y") +
  theme_minimal()

# bias clearly highlights where predictions are biased when shown as a direct
# measure

# relative deviance clearly highlights the effects of noisiness, and this is
# larger when the noise is large relative to the fraction

# deviance under a binomial distribution works just as well as betabinomial here
# - so just use this instead?

# rmse is more affected by bias as a direct measure, and relative difference are
# not clear

# cvm, ks, and crps are all better than rmse, but not better than deviance.

# Note the PS2 metric in Taggart 2022
# (http://www.bom.gov.au/research/publications/researchreports/BRR-064.pdf)
# corresponds to cvm, but with decomposition into over/underprediction and
# incorrect dispersion (over/underconfidence). As seen here, a well-calibrated
# model has expected value 0



# Behaviour of the posterior predictive metrics ----------------------------

# The comparison above varies the location of a prediction, but not its
# dispersion, so it cannot show what the newer metrics are for (#10): a model
# can be wrong by predicting the right average with the wrong spread, and a
# point summary such as mean squared error is blind to that. This adds
# over- and under-dispersed scenarios, and marks the noise floor that no model
# can beat.

source("R/validation_functions.R")

set.seed(2026 - 8 - 31)

n_per_fraction <- 200
rho_true <- 0.155
mosquito_number <- 100

scenarios <- c("correct", "biased", "overdispersed", "underdispersed")

behaviour <- lapply(
  seq(0.1, 0.95, length.out = 25),
  function(true_fraction) {

    died <- rbetabinom(n_per_fraction, mosquito_number, true_fraction, rho_true)

    lapply(
      scenarios,
      function(scenario) {

        predicted <- switch(
          scenario,
          correct = true_fraction,
          biased = min(true_fraction + 0.1, 0.99),
          overdispersed = true_fraction,
          underdispersed = true_fraction
        )
        rho <- switch(
          scenario,
          correct = rho_true,
          biased = rho_true,
          overdispersed = rho_true * 3,
          underdispersed = rho_true / 3
        )

        # a predictive distribution with no parameter uncertainty, so the
        # comparison isolates the effect of getting location or spread wrong
        p_draws <- matrix(predicted, nrow = 1, ncol = n_per_fraction)
        rho_draws <- matrix(rho, nrow = 1, ncol = n_per_fraction)

        summary <- ppd_summary(died, rep(mosquito_number, n_per_fraction),
                               p_draws, rho_draws)
        pit <- ppd_pit(summary, n_rep = 20)
        sims <- ppd_simulate(rep(mosquito_number, n_per_fraction),
                             matrix(predicted, nrow = 500, ncol = n_per_fraction),
                             matrix(rho, nrow = 500, ncol = n_per_fraction))

        data.frame(
          true_fraction = true_fraction,
          scenario = scenario,
          coverage_95 = coverage_curve(pit, levels = 0.95)$empirical,
          mean_pit = mean(pit),
          cvm = pit_statistic(pit, cvm_stat),
          crps = mean(ppd_crps(died, rep(mosquito_number, n_per_fraction), sims)),
          elpd = mean(summary$log_score),
          mse = mean((summary$observed - summary$predicted) ^ 2)
        )

      }
    ) %>%
      bind_rows()

  }
) %>%
  bind_rows()

# the floor: the error that remains when the prediction is exactly right
floor_by_fraction <- behaviour %>%
  filter(scenario == "correct") %>%
  select(true_fraction, floor_mse = mse, floor_crps = crps, floor_elpd = elpd)

behaviour_plot <- behaviour %>%
  left_join(floor_by_fraction, by = "true_fraction") %>%
  pivot_longer(
    cols = c("coverage_95", "mean_pit", "cvm", "crps", "elpd", "mse"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      coverage_95 = "95% interval coverage",
      mean_pit = "mean PIT",
      cvm = "Cramer-von Mises",
      crps = "CRPS",
      elpd = "log predictive density",
      mse = "mean squared error"
    ),
    metric = factor(metric,
                    levels = c("95% interval coverage", "mean PIT",
                               "Cramer-von Mises", "CRPS",
                               "log predictive density",
                               "mean squared error"))
  ) %>%
  ggplot(
    aes(x = true_fraction,
        y = value,
        colour = scenario)
  ) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(
    values = c(correct = "#2166AC",
               biased = "#B2182B",
               overdispersed = "#E08214",
               underdispersed = "#4D9221"),
    name = ""
  ) +
  labs(
    x = "true population susceptible fraction",
    y = "",
    title = "What each validation metric can and cannot detect",
    subtitle = paste("mean squared error is identical for the correct,",
                     "over- and under-dispersed predictions:",
                     "\nonly the distributional metrics separate them")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/validation_metric_behaviour.png",
       behaviour_plot,
       bg = "white",
       width = 10,
       height = 7)

# confirm the point the figure is making
behaviour %>%
  group_by(scenario) %>%
  summarise(
    mse = mean(mse),
    coverage_95 = mean(coverage_95),
    cvm = mean(cvm),
    .groups = "drop"
  ) %>%
  print()
