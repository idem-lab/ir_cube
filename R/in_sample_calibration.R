# In-sample calibration of both models (#69, #70): coverage of the central 50%
# and 95% predictive intervals and the mean PIT, at every modelled bioassay,
# per insecticide and pooled.
#
#   Rscript R/in_sample_calibration.R
#
# Run after R/two_stage_maps.R. The predictive draws are assay_draws()'s
# (R/two_stage_predictions.R): a new assay at each modelled bioassay's
# pixel-year, from the mechanistic model alone (ilogit(m)) and from the full
# model (with fresh pixel-year, pixel and study-year noise). Both are scored as
# the cross-validation scores held-out bioassays (score_fold(),
# summarise_experiment()): beta-binomial at the external per-type
# overdispersion, coverage as its expectation over the PIT randomisation, and
# the mean PIT over 100 randomisations.
#
# Prints the table, and writes the coverage curve at the cross-validation's
# levels to outputs/in_sample_coverage.csv (as outputs/cv_coverage.csv), for
# R/fig_predictive_validation.R.

source("R/packages.R")
source("R/two_stage_predictions.R")
source("R/validation_scoring.R")

fit_env <- new.env()
load("temporary/fitted_model.RData", envir = fit_env)
df <- fit_env$df
types <- fit_env$types
rm(fit_env)
invisible(gc())

# the coverage curve's levels, as R/validation_metrics.R
coverage_levels <- seq(0.1, 0.95, by = 0.05)

sims <- assay_draws(df, types)
rho <- rho_for_record(df, rho_lookup())

models <- c(dynamical = "dynamical", two_stage = "two_stage")
set.seed(string_seed("in-sample calibration"))
scored <- lapply(models, function(model) {
  summary <- ppd_summary(df$died, df$mosquito_number, sims[[model]],
                         matrix(rho, nrow(sims[[model]]), nrow(df),
                                byrow = TRUE))
  list(summary = summary, pit = ppd_pit(summary, n_rep = n_pit_reps))
})

calibration <- function(model, rows) {
  summary <- scored[[model]]$summary[rows, ]
  tibble(n = length(rows),
         coverage_50 = mean(expected_coverage(summary$cdf_below,
                                              summary$pmf_at, 0.5)),
         coverage_95 = mean(expected_coverage(summary$cdf_below,
                                              summary$pmf_at, 0.95)),
         mean_pit = mean(scored[[model]]$pit[rows, ]))
}

groups <- c(list(All = seq_len(nrow(df))),
            split(seq_len(nrow(df)), df$insecticide_type))
table <- bind_rows(lapply(names(models), function(model) {
  bind_rows(lapply(names(groups), function(group) {
    calibration(model, groups[[group]]) %>%
      mutate(model = model, insecticide_type = group, .before = everything())
  }))
}))
cat("\nIn-sample calibration at the modelled bioassays (coverage: share of",
    "bioassays in the central interval; mean PIT: 0.5 if unbiased)\n")
print(table, n = Inf, width = Inf)

coverage_curves <- bind_rows(lapply(names(models), function(model) {
  summary <- scored[[model]]$summary
  tibble(model = model, experiment = "in_sample",
         nominal = coverage_levels,
         empirical = vapply(coverage_levels, function(level) {
           mean(expected_coverage(summary$cdf_below, summary$pmf_at, level))
         }, numeric(1)))
}))
write.csv(coverage_curves, "outputs/in_sample_coverage.csv", row.names = FALSE)
