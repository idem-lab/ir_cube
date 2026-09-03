# Score the saved cross-validation posterior predictive draws.
#
# Reads everything written by run_validation_folds.R and writes tidy tables of
# per-record scores and per-experiment summaries, which the figures then read.
# Keeping scoring separate from fitting means metrics can be revised without
# refitting anything (#10).
#
# Three questions are asked of each model, one measure each:
#
#   is the predictive distribution the right width and shape?
#       coverage of central predictive intervals, summarised by the
#       Cramer-von Mises statistic on the randomised PIT values
#   is the whole distribution close to the data?
#       CRPS, in mortality units, against the noise floor
#   is it right on average, and conditional on the prediction?
#       reliability bins, and the mean PIT
#
# Scores are also computed on pooled groups of assays. A single bioassay is a
# noisy measurement of the population fraction the model is predicting, so
# aggregation is what brings the comparison to bear on that quantity: the
# reference in every case is the model's own aggregated predictive
# distribution, so heterogeneity in the true fraction within a group is carried
# by the model's predictions rather than assumed away.

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

draws_dir <- "outputs/cv_draws"
n_pit_reps <- 100
coverage_levels <- seq(0.1, 0.95, by = 0.05)

set.seed(2026 - 8 - 31)

# externally estimated overdispersion, from replicate bioassays in the same
# pixel, year and insecticide (estimate_bioassay_rho.R). This sets the noise
# floor, and is independent of any of the models being scored
rho_external <- read.csv("outputs/bioassay_rho.csv")
rho_for_class <- function(insecticide_class) {
  index <- match(insecticide_class, rho_external$insecticide_class)
  # fall back to the pooled estimate for any class not fitted separately
  pooled <- rho_external$rho[rho_external$insecticide_class == "all"]
  ifelse(is.na(index), pooled, rho_external$rho[index])
}

files <- list.files(draws_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("no draws found in ", draws_dir, "; run R/run_validation_folds.R first")
}

cat(sprintf("scoring %i saved folds\n", length(files)))


# per record ---------------------------------------------------------------

score_fold <- function(file) {

  fold <- readRDS(file)
  test <- fold$test_df

  summary <- ppd_summary(test$died,
                         test$mosquito_number,
                         fold$p_draws,
                         fold$rho_draws)

  pit <- ppd_pit(summary, n_rep = n_pit_reps)
  sims <- ppd_simulate(test$mosquito_number, fold$p_draws, fold$rho_draws)

  scores <- summary %>%
    mutate(
      model = fold$model,
      experiment = fold$experiment,
      fold = fold$fold,
      insecticide_type = test$insecticide_type,
      insecticide_class = test$insecticide_class,
      country_name = test$country_name,
      year_start = test$year_start,
      cell = test$cell,
      # averaged over randomisation replicates, for plotting and for the mean;
      # the uniformity statistics use the replicates individually
      pit = rowMeans(pit),
      crps = ppd_crps(test$died, test$mosquito_number, sims),
      rho_external = rho_for_class(test$insecticide_class),
      .before = everything()
    )

  list(scores = scores,
       pit = pit,
       sims = sims,
       fold = fold)

}

scored <- lapply(files, score_fold)
names(scored) <- basename(files)

all_scores <- bind_rows(lapply(scored, `[[`, "scores"))

write.csv(all_scores, "outputs/cv_scores.csv", row.names = FALSE)


# per experiment -----------------------------------------------------------

# summarise one model in one experiment, pooling its folds
summarise_experiment <- function(scores, pit_list) {

  pit <- do.call(rbind, pit_list)
  n_obs <- nrow(scores)

  coverage <- coverage_curve(pit, levels = c(0.5, 0.95))
  floor_mse <- noise_floor_mse(scores$died,
                               scores$mosquito_number,
                               scores$rho_external)

  data.frame(
    n = n_obs,
    mean_pit = mean(pit),
    coverage_50 = coverage$empirical[1],
    coverage_95 = coverage$empirical[2],
    cvm = pit_statistic(pit, cvm_stat),
    cvm_null_upper = pit_null_band(n_obs, cvm_stat, n_sim = 200)[3],
    ks = pit_statistic(pit, ks_stat),
    crps = mean(scores$crps),
    elpd = mean(scores$log_score),
    bias = mean(scores$predicted - scores$observed),
    mse = mean((scores$observed - scores$predicted) ^ 2),
    mse_floor = floor_mse
  )

}

keys <- bind_rows(lapply(scored, function(x) {
  data.frame(model = x$fold$model, experiment = x$fold$experiment)
}))

summaries <- lapply(
  split(seq_along(scored), paste(keys$model, keys$experiment)),
  function(index) {
    summarise_experiment(
      bind_rows(lapply(scored[index], `[[`, "scores")),
      lapply(scored[index], `[[`, "pit")
    ) %>%
      mutate(model = keys$model[index[1]],
             experiment = keys$experiment[index[1]],
             .before = everything())
  }
)
summaries <- bind_rows(summaries)

# skill against the nearest neighbour null, anchored at the noise floor: 0 is
# the null model, 1 is as good as bioassay noise allows
summaries <- summaries %>%
  group_by(experiment) %>%
  mutate(
    mse_null = mse[model == "nearest_neighbour"],
    skill = mse_skill(mse, mse_null, mse_floor)
  ) %>%
  ungroup()

write.csv(summaries, "outputs/cv_summary.csv", row.names = FALSE)

cat("\nsummary by experiment and model:\n")
print(summaries %>%
        select(experiment, model, n, coverage_95, mean_pit, crps, skill, cvm) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# coverage curves ----------------------------------------------------------

coverage_curves <- lapply(
  split(seq_along(scored), paste(keys$model, keys$experiment)),
  function(index) {
    pit <- do.call(rbind, lapply(scored[index], `[[`, "pit"))
    coverage_curve(pit, levels = coverage_levels) %>%
      mutate(model = keys$model[index[1]],
             experiment = keys$experiment[index[1]],
             .before = everything())
  }
)
coverage_curves <- bind_rows(coverage_curves)
write.csv(coverage_curves, "outputs/cv_coverage.csv", row.names = FALSE)


# reliability --------------------------------------------------------------

reliability <- all_scores %>%
  group_by(model, experiment) %>%
  group_modify(~ {
    bins <- reliability_bins(.x$predicted, .x$observed, n_bins = 10)
    # the scatter a perfect model would still show in each bin, from the
    # external overdispersion, the assay sizes and the number of assays pooled
    bins$envelope <- reliability_envelope(
      p = bins$predicted,
      mosquito_number = median(.x$mosquito_number),
      k = bins$n,
      rho = mean(.x$rho_external)
    )
    bins
  }) %>%
  ungroup()

write.csv(reliability, "outputs/cv_reliability.csv", row.names = FALSE)


# aggregated scores --------------------------------------------------------

# Two ladders. The first pools assays sharing a pixel, year and insecticide:
# the model asserts a single fraction there, so these are replicates of one
# population quantity, and this is the strictest test of the model at its own
# unit of inference. The second pools by country, year and insecticide: many
# more assays per group, so assay noise falls further, at the cost of testing
# an aggregate rather than any single pixel.
aggregate_fold <- function(entry, grouping) {

  test <- entry$fold$test_df
  group <- switch(
    grouping,
    pixel_year = paste(test$cell, test$year_start, test$insecticide_type),
    country_year = paste(test$country_name, test$year_start,
                         test$insecticide_type)
  )

  ppd_aggregate(test$died, test$mosquito_number, group, entry$sims) %>%
    mutate(model = entry$fold$model,
           experiment = entry$fold$experiment,
           grouping = grouping,
           .before = everything())

}

aggregated <- bind_rows(
  unname(
    c(lapply(scored, aggregate_fold, grouping = "pixel_year"),
      lapply(scored, aggregate_fold, grouping = "country_year"))
  )
)

write.csv(aggregated, "outputs/cv_aggregate.csv", row.names = FALSE)

aggregate_summary <- aggregated %>%
  group_by(grouping, experiment, model) %>%
  summarise(
    groups = n(),
    mean_assays = mean(n_assays),
    coverage_95 = mean(observed >= lower & observed <= upper),
    mean_pit = mean(pit),
    bias = mean(predicted - observed),
    rmse = rmse(observed, predicted),
    .groups = "drop"
  )

cat("\naggregated calibration:\n")
print(aggregate_summary %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())

write.csv(aggregate_summary, "outputs/cv_aggregate_summary.csv",
          row.names = FALSE)


# the model's overdispersion against the external estimate ------------------

# a fitted overdispersion larger than the replicate-based estimate would mean
# the model is absorbing process misfit into the observation process, which
# would also show as over-coverage
rho_comparison <- bind_rows(lapply(scored, function(entry) {
  test <- entry$fold$test_df
  data.frame(
    model = entry$fold$model,
    experiment = entry$fold$experiment,
    fold = entry$fold$fold,
    insecticide_class = test$insecticide_class,
    rho_fitted = colMeans(entry$fold$rho_draws)
  )
})) %>%
  group_by(model, experiment, insecticide_class) %>%
  summarise(rho_fitted = mean(rho_fitted), .groups = "drop") %>%
  mutate(rho_external = rho_for_class(insecticide_class))

write.csv(rho_comparison, "outputs/cv_rho_comparison.csv", row.names = FALSE)

cat("\nfitted against externally estimated overdispersion:\n")
print(rho_comparison %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# WHO threshold check ------------------------------------------------------

# The same comparison in the terms the results are used in: WHO classifies a
# population as showing confirmed resistance below 90% mortality, and possible
# resistance between 90% and 98%. This asks whether the model predicts the right
# proportion of held-out bioassays in each category, which needs no
# distributional vocabulary to read.
who_thresholds <- bind_rows(unname(lapply(scored, function(entry) {

  test <- entry$fold$test_df
  observed_proportion <- test$died / test$mosquito_number
  simulated_proportion <- sweep(entry$sims, 2, test$mosquito_number, FUN = "/")

  categorise <- function(x) {
    c(confirmed = mean(x < 0.9),
      possible = mean(x >= 0.9 & x < 0.98),
      susceptible = mean(x >= 0.98))
  }

  observed <- categorise(observed_proportion)
  # the same summary under each posterior predictive replicate dataset
  simulated <- t(apply(simulated_proportion, 1, categorise))

  data.frame(
    model = entry$fold$model,
    experiment = entry$fold$experiment,
    fold = entry$fold$fold,
    category = names(observed),
    observed = as.numeric(observed),
    predicted = colMeans(simulated),
    lower = apply(simulated, 2, quantile, 0.025),
    upper = apply(simulated, 2, quantile, 0.975)
  )

})))

write.csv(who_thresholds, "outputs/cv_who_thresholds.csv", row.names = FALSE)

cat("\nWHO resistance categories, observed against predicted:\n")
print(who_thresholds %>%
        group_by(model, experiment, category) %>%
        summarise(across(c(observed, predicted, lower, upper), mean),
                  .groups = "drop") %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())
