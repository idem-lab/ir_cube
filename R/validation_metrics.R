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
#       coverage of central predictive intervals, with the Cramer-von Mises
#       statistic on the randomised PIT values as a single scalar summary
#   is the whole distribution close to the data?
#       CRPS, in mortality units, and mean squared error against the noise floor
#   is it right on average, and conditional on the prediction?
#       reliability bins, and the mean PIT
#
# Coverage, mean PIT and the Cramer-von Mises statistic are all functionals of
# one PIT distribution, so only these are kept: the Kolmogorov-Smirnov
# statistic added a fourth view of the same object, and a null band for the
# Cramer-von Mises statistic assumed independent PIT values, which does not hold
# for records scored against a shared posterior (#12 review).
#
# Point metrics (MSE, excess, bias, reliability bins, pooled means) score the
# map, the point prediction a reader would take from the published maps:
# for the two-stage model the posterior mean without its site noise u and p,
# for every other model its predictive mean (#32; score_fold()). The
# distributional metrics score each model's own predictive distribution.
#
# Every model is scored at the same overdispersion, the external
# replicate-based estimate, rather than at its own fitted value. Letting each
# model choose made coverage a comparison of dispersion rather than of
# prediction: a model can look calibrated by being vague, and the intercept null
# reached nominal coverage that way. The dispersion each model's own residuals
# imply is reported separately, in cv_rho_comparison.csv.
#
# Scores are also computed on pooled groups of assays. A single bioassay is a
# noisy measurement of the population fraction the model is predicting, so
# aggregation is what brings the comparison to bear on that quantity: the
# reference in every case is the model's own aggregated predictive
# distribution, so heterogeneity in the true fraction within a group is carried
# by the model's predictions rather than assumed away. Only the country-year
# rung is reported; pooling within pixel, year and insecticide averaged 1.3
# assays per group, so it was indistinguishable from the unpooled scores.

source("R/validation_scoring.R")

set.seed(2026 - 8 - 31)
coverage_levels <- seq(0.1, 0.95, by = 0.05)

# externally estimated overdispersion, from replicate bioassays in the same
# pixel, year and insecticide. This sets the noise floor, and is independent of
# any of the models being scored. See rho_lookup() in validation_functions.R.
rho_source <- rho_lookup()
cat("overdispersion:", rho_source$source, "\n")

# every fold on disk, the two-stage model's among them (two_stage__*.rds)
files <- draws_files()
cat(sprintf("scoring %i saved folds\n", length(files)))


# per record ---------------------------------------------------------------

# scored one at a time, with an explicit collection between folds: reading a
# model fold means holding its 1.8 GB saved object briefly
scored <- lapply(seq_along(files), function(i) {
  cat(sprintf("%s | scoring %s\n", format(Sys.time(), "%H:%M:%S"),
              basename(files[i])))
  flush(stdout())
  on.exit(gc(verbose = FALSE))
  with_fold_stream(files[i], "score",
                   score_fold(files[i], rho_source))
})
names(scored) <- basename(files)

all_scores <- bind_rows(lapply(scored, `[[`, "scores"))

write.csv(all_scores, "outputs/cv_scores.csv", row.names = FALSE)


# per experiment -----------------------------------------------------------

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

# Variance in the population fraction that each model explains, anchored on the
# intercept null and on the noise floor: 0 is the no-information baseline, 1 is
# as good as bioassay noise allows. The anchor was previously the nearest
# neighbour null, which pinned an informative baseline at zero by construction
# and hid that it is itself worse than a global per-insecticide mean under
# spatial extrapolation (#12 review).
#
# `excess` is reported alongside every ratio: it is mean squared error above the
# noise floor in absolute mortality-squared units, so the conclusion does not
# rest entirely on the floor. `rms_p` is its square root, an error in the
# population fraction itself.
#
# The two-stage model also attributes part of its excess to its site noise u
# and p: `model_floor`, the model floor U (#32), labelled as attributed by the
# model; `map_error` is the rest, the map's own error. NA for the other models
summaries <- summaries %>%
  group_by(experiment) %>%
  mutate(
    excess = mse - mse_floor,
    rms_p = sqrt(pmax(excess, 0)),
    map_error = excess - model_floor
  ) %>%
  ungroup()

write.csv(summaries, "outputs/cv_summary.csv", row.names = FALSE)

cat("\nsummary by experiment and model:\n")
print(summaries %>%
        select(experiment, model, n, coverage_95, coverage_95_map, mean_pit,
               crps, bias, bias_predictive, mse, mse_floor, excess,
               model_floor, rms_p, cvm) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# coverage curves ----------------------------------------------------------

coverage_curves <- lapply(
  split(seq_along(scored), paste(keys$model, keys$experiment)),
  function(index) {
    # the exact expectation over the PIT randomisation, as in the summaries
    scores <- bind_rows(lapply(scored[index], `[[`, "scores"))
    data.frame(nominal = coverage_levels,
               empirical = vapply(coverage_levels, function(level) {
                 mean(expected_coverage(scores$cdf_below, scores$pmf_at,
                                        level))
               }, numeric(1))) %>%
      mutate(model = keys$model[index[1]],
             experiment = keys$experiment[index[1]],
             .before = everything())
  }
)
coverage_curves <- bind_rows(coverage_curves)
write.csv(coverage_curves, "outputs/cv_coverage.csv", row.names = FALSE)


# reliability --------------------------------------------------------------

# Binned on the map (score_fold()), never on the observation. Binning on the observed
# mortality induces regression to the mean and makes a calibrated model look
# badly biased at both extremes; conditioning on the prediction is the question
# actually of interest — when the model says 60%, is the average outcome 60%.
#
# The envelope comes from the model's own posterior predictive distribution, so
# a gap outside it is the model's error rather than the diagnostic's; for the
# two-stage model that includes its fresh u and p, so the envelope is centred
# off the map by the pull towards 0.5. That
# requires the posterior draws, so it is computed per fold and then pooled by
# experiment, weighting each fold by its held-out records.
reliability <- all_scores %>%
  group_by(model, experiment) %>%
  group_modify(~ reliability_bins(.x$map, .x$observed, n_bins = 10)) %>%
  ungroup()

reliability_checks <- bind_rows(lapply(scored, function(entry) {
  bins <- reliability_bins(entry$scores$map,
                           entry$scores$observed,
                           n_bins = 10)
  with_fold_stream(entry$fold$file, "reliability", bind_cols(
    data.frame(model = entry$fold$model,
               experiment = entry$fold$experiment,
               fold = entry$fold$fold),
    bins,
    reliability_ppc(predicted = entry$scores$map,
                    mosquito_number = entry$scores$mosquito_number,
                    p_draws = entry$p_draws,
                    rho = entry$rho_scoring,
                    n_bins = 10,
                    n_rep = 200)
  ))
})) %>%
  mutate(gap = observed - predicted,
         beyond_ppc = gap < ppc_lower | gap > ppc_upper)

# and the same envelope carried onto the pooled bins the figure draws, each
# fold weighted by its held-out records. An analytic envelope also used to sit
# on this table; it conditioned on the prediction being the truth, so it could
# not express posterior uncertainty, and reliability_envelope() is gone with it.
reliability <- reliability %>%
  left_join(
    reliability_checks %>%
      group_by(model, experiment, bin) %>%
      summarise(ppc_lower = weighted.mean(ppc_lower, n),
                ppc_upper = weighted.mean(ppc_upper, n),
                .groups = "drop"),
    by = c("model", "experiment", "bin")
  )
stopifnot(!anyNA(reliability$ppc_lower))

write.csv(reliability, "outputs/cv_reliability.csv", row.names = FALSE)
write.csv(reliability_checks, "outputs/cv_reliability_ppc.csv",
          row.names = FALSE)

cat("\nreliability against the model's own posterior predictive envelope",
    "(lowest and highest predicted decile):\n")
print(reliability_checks %>%
        filter(bin %in% c(1, 10)) %>%
        select(experiment, fold, model, bin, n, predicted, observed, gap,
               ppc_lower, ppc_upper, beyond_ppc) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# aggregated scores --------------------------------------------------------

# pooled by country, year and insecticide (aggregate_fold())
aggregated <- bind_rows(
  unname(lapply(scored, function(entry) {
    with_fold_stream(entry$fold$file, "aggregate",
                     aggregate_fold(entry, grouping = "country_year"))
  }))
)

write.csv(aggregated, "outputs/cv_aggregate.csv", row.names = FALSE)

aggregate_summary <- aggregated %>%
  group_by(grouping, experiment, model) %>%
  summarise(
    groups = n(),
    mean_assays = mean(n_assays),
    coverage_95 = mean(observed >= lower & observed <= upper),
    mean_pit = mean(pit),
    # the pooled map; the predictive mean's bias alongside
    bias = mean(map - observed),
    bias_predictive = mean(predicted - observed),
    rmse = rmse(observed, map),
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
# sampling diagnostics per fold, so that the convergence caveat travels with
# the results. Null model folds are analytic and have none
sampling_diagnostics <- bind_rows(lapply(scored, function(entry) {
  fold <- entry$fold
  if (is.null(fold$ess_p)) return(NULL)
  data.frame(
    model = fold$model,
    experiment = fold$experiment,
    fold = fold$fold,
    n_chains = fold$n_chains,
    n_sampled = fold$n_sampled,
    ess_p_median = median(fold$ess_p, na.rm = TRUE),
    ess_p_min = min(fold$ess_p, na.rm = TRUE),
    ess_rho_median = median(fold$ess_rho, na.rm = TRUE),
    rhat_worst = max(fold$convergence[, 1], na.rm = TRUE),
    rhat_above_1.01 = sum(fold$convergence[, 1] > 1.01, na.rm = TRUE)
  )
}))

if (nrow(sampling_diagnostics) > 0) {
  write.csv(sampling_diagnostics, "outputs/cv_sampling_diagnostics.csv",
            row.names = FALSE)
  cat("\nsampling diagnostics:\n")
  print(sampling_diagnostics %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          as.data.frame())
}


rho_comparison <- bind_rows(lapply(scored, function(entry) {
  entry$scores %>%
    select(model, experiment, fold, insecticide_type, insecticide_class,
           rho_fitted)
})) %>%
  # by type as well as class: the models fit one overdispersion per class, but
  # the external estimate is now per type, so the comparison is made at the
  # finer of the two
  group_by(model, experiment, insecticide_class, insecticide_type) %>%
  summarise(rho_fitted = mean(rho_fitted), .groups = "drop")
rho_comparison$rho_external <- rho_for_record(rho_comparison, rho_source)

write.csv(rho_comparison, "outputs/cv_rho_comparison.csv", row.names = FALSE)

cat("\nfitted against externally estimated overdispersion:\n")
print(rho_comparison %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# the data floor, checked directly ------------------------------------------

# The floor every score uses is noise_floor_mse() at the external rho. On the
# held-out pixel-years with replicate assays it can be checked against the
# replicates' own differences, which need no rho (floor_pair_check(), #32).
# A ratio near 1 supports the floor; the replicated subset need not be
# representative of all held-out assays (share_replicated, rho_floor_all)
floor_check <- floor_pair_check(all_scores %>% filter(model == "dynamical"))
write.csv(floor_check, "outputs/cv_floor_pair_check.csv", row.names = FALSE)

cat("\ndata floor: replicate pairs against rho, same held-out assays:\n")
print(floor_check %>%
        mutate(across(where(is.double), ~ round(.x, 4))) %>%
        as.data.frame())


# per fold and per year ----------------------------------------------------

# The pooled numbers hide which folds carry the result, and master reported a
# per-country and a per-lead-year breakdown that the first version of this
# pipeline dropped. `excess` is again mean squared error above the noise floor.
# A per-group share of the intercept null's excess used to sit here too; it was
# the intercept-referenced, floor-corrected variance explained under another
# name, and variance_explained.R carries the one definition of that (#12
# review). by_group() is in validation_scoring.R.

by_fold <- by_group(all_scores, fold)
write.csv(by_fold, "outputs/cv_by_fold.csv", row.names = FALSE)

cat("\nby fold:\n")
print(by_fold %>%
        select(experiment, fold, model, n, bias, coverage_95, excess,
               model_floor) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())

by_year <- by_group(
  all_scores %>% filter(startsWith(experiment, "temporal_forecasting")),
  year_start
)
write.csv(by_year, "outputs/cv_by_year.csv", row.names = FALSE)

if (nrow(by_year) > 0) {
  cat("\nforecasting, by lead year:\n")
  print(by_year %>%
          select(experiment, year_start, model, n, mean_observed,
                 mean_predicted, bias, mean_pit, excess) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          as.data.frame())
}


# uncertainty --------------------------------------------------------------

# The paired pixel-cluster bootstrap on excess MSE that used to sit here is
# gone with the intercept-null skill metric it supported. variance_explained.R
# carries the one bootstrap and the one definition of variance explained; two
# of each, differing in reference and in floor treatment, was the confusion the
# review objected to. It also pivoted models wide on value columns, which would
# have collapsed genuinely duplicate assays had any two matched exactly.
