# Scoring the saved cross-validation folds: the functions behind
# R/validation_metrics.R, R/variance_explained.R and R/two_stage_metrics.R.
# Functions and settings only; sourcing this file scores nothing. Sources
# R/validation_functions.R, which holds the primitives (the beta-binomial
# predictive distribution, PIT, CRPS, the noise floor and rho_lookup()).

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

draws_dir <- "outputs/cv_draws"
n_pit_reps <- 100

# only the three experiments that survive; the leave-one-country-out folds
# confound spatial prediction with the country initial condition, and the
# three-year 2020 forecast fold leaked its cut year into the holdout. Draws for
# both are parked in outputs/cv_draws_defunct.
scored_experiments <- c("spatial_interpolation", "spatial_blocks",
                        "temporal_forecasting")

# Folds fitted by MCMC hold up to 20,000 draws (4 chains x 5,000). The scoring
# functions build an n_draws x (died + 1) matrix per observation, so that is
# twenty times the work of the 1,000 draws the null models carry, for no gain:
# at roughly 50 draws per effective sample, every tenth draw retains almost all
# the information. Thinning is applied to the model folds only; the null models
# are analytic and already independent.
max_draws <- 2000

thin_draws <- function(x, maximum = max_draws) {
  if (nrow(x) <= maximum) return(x)
  keep <- round(seq(1, nrow(x), length.out = maximum))
  x[keep, , drop = FALSE]
}


# the saved folds ------------------------------------------------------------

# Every fold on disk, <model>__<experiment>__<fold>.rds, named by model. Stops
# on anything outside the scored experiments rather than silently scoring it
draws_files <- function(dir = draws_dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("no draws found in ", dir, "; run R/run_validation_folds.R first")
  }
  parts <- strsplit(sub("\\.rds$", "", basename(files)), "__")
  stopifnot(all(lengths(parts) == 3),
            all(vapply(parts, `[`, "", 2) %in% scored_experiments))
  setNames(files, vapply(parts, `[`, "", 1))
}


# random number streams ---------------------------------------------------------

# The models of #12 draw their PIT randomisations and simulated replicates from
# the one stream seeded at the top of validation_metrics.R, in file order, as
# they did when their published tables were made. Any other model draws from a
# stream of its own per fold and per use, seeded from the fold's file name and
# the use, and leaves the main stream where it was: adding a model then changes
# no number for the others, and its folds do not share uniforms.
main_stream_models <- c("dynamical", "intercept", "nearest_neighbour",
                        "nearest_neighbour_oracle")

with_fold_stream <- function(file, use, expr) {
  if (sub("__.*$", "", basename(file)) %in% main_stream_models) return(expr)
  key <- paste(basename(file), use)
  saved <- .Random.seed
  on.exit(assign(".Random.seed", saved, envir = globalenv()))
  set.seed(string_seed(key))
  expr
}

# a seed from a string: a polynomial string hash, exact in double precision
string_seed <- function(key) {
  Reduce(function(h, x) (h * 31 + x) %% 2147483647, utf8ToInt(key), 0)
}


# per record ----------------------------------------------------------------------

# Score one saved fold, every record at the external, replicate-based
# overdispersion of rho_source (rho_lookup()). Returns only what the summaries
# need: a saved model fold carries the greta draws and arrays, most of a
# 1.8 GB file, and keeping all of them held 26 GB.
score_fold <- function(file, rho_source, n_rep = n_pit_reps) {

  fold <- readRDS(file)
  # the model is reported under its saved name, which is the file's prefix
  stopifnot(identical(fold$model, sub("__.*$", "", basename(file))))
  test <- fold$test_df

  # predictions come from the saved object. For model folds these were produced
  # by greta's calculate(values = draws) in MCMC order (or by the two-stage
  # model, R/run_two_stage_folds.R); for the null models they are analytic.
  # Either way they are not recomputed here
  p_draws <- thin_draws(fold$p_draws)

  # the overdispersion the fit implies, kept for the diagnostic table: a
  # posterior for the dynamical model. The nulls no longer fit their own - a
  # null earns its place on point prediction, and letting it choose its own
  # dispersion made coverage a comparison of vagueness (#12 review) - and the
  # two-stage model fixes it at the external value, so they read as missing
  rho_fitted <- if (!is.null(fold$rho_type_draws)) {
    # saved compactly as draws by insecticide type, expanded here
    thin_draws(fold$rho_type_draws)[, fold$type_id, drop = FALSE]
  } else {
    matrix(NA_real_, nrow = nrow(p_draws), ncol = nrow(test))
  }

  # every model is scored at the external, replicate-based estimate
  rho_scoring <- rho_for_record(test, rho_source)
  rho_draws <- matrix(rho_scoring, nrow = nrow(p_draws), ncol = nrow(test),
                      byrow = TRUE)

  summary <- ppd_summary(test$died, test$mosquito_number, p_draws, rho_draws)
  pit <- ppd_pit(summary, n_rep = n_rep)
  sims <- ppd_simulate(test$mosquito_number, p_draws, rho_draws)

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
      # One randomisation replicate, not the mean of them. Averaging over
      # replicates converges to the mid-P value `cdf_below + 0.5 * pmf_at`,
      # which is not uniform under calibration for discrete data: with ~30% of
      # held-out assays at 100% mortality a perfectly calibrated model reads
      # 0.969 at nominal 0.95. The uniformity statistics use the full matrix
      # and were never affected; this column feeds the figures (#12 review)
      pit = pit[, 1],
      crps = ppd_crps(test$died, test$mosquito_number, sims),
      rho_external = rho_scoring,
      rho_fitted = colMeans(rho_fitted),
      .before = everything()
    )

  list(scores = scores,
       pit = pit,
       sims = sims,
       p_draws = p_draws,
       rho_scoring = rho_scoring,
       fold = list(file = file,
                   model = fold$model,
                   experiment = fold$experiment,
                   fold = fold$fold,
                   test_df = test,
                   convergence = fold$convergence,
                   ess_p = fold$ess_p,
                   ess_rho = fold$ess_rho,
                   n_sampled = fold$n_sampled,
                   n_chains = fold$n_chains))
}

# The expected coverage of each record's central interval at `level`, over the
# PIT randomisation u = F(y - 1) + v f(y), v ~ U(0, 1): the share of
# [cdf_below, cdf_below + pmf_at] inside ((1 - level) / 2, (1 + level) / 2).
# The mean over records is the limit of coverage_curve() as the replicates
# grow, without the Monte Carlo noise, so differences between models can be
# bootstrapped record by record
expected_coverage <- function(cdf_below, pmf_at, level) {
  lower <- (1 - level) / 2
  upper <- (1 + level) / 2
  overlap <- pmax(0, pmin(upper, cdf_below + pmf_at) - pmax(lower, cdf_below))
  ifelse(pmf_at > 0, overlap / pmf_at,
         as.numeric(cdf_below > lower & cdf_below < upper))
}


# summaries -------------------------------------------------------------------------

# one model in one experiment, pooling its folds. Coverage is the exact
# expectation over the PIT randomisation (expected_coverage()), as in every
# other table
summarise_experiment <- function(scores, pit_list) {
  pit <- do.call(rbind, pit_list)
  data.frame(
    n = nrow(scores),
    mean_pit = mean(pit),
    coverage_50 = mean(expected_coverage(scores$cdf_below, scores$pmf_at,
                                         0.5)),
    coverage_95 = mean(expected_coverage(scores$cdf_below, scores$pmf_at,
                                         0.95)),
    cvm = pit_statistic(pit, cvm_stat),
    crps = mean(scores$crps),
    elpd = mean(scores$log_score),
    bias = mean(scores$predicted - scores$observed),
    mse = mean((scores$observed - scores$predicted) ^ 2),
    mse_floor = noise_floor_mse(scores$died, scores$mosquito_number,
                                scores$rho_external)
  )
}

# Assays pooled by country, year and insecticide: 17-18 assays per group, so
# assay noise falls roughly seventeen-fold and the comparison is nearly purely
# about the population fraction, at the cost of testing an aggregate rather
# than any single pixel. Pooling within pixel, year and insecticide was also
# reported, and dropped: it averaged 1.3 assays per group, so it was the
# unpooled comparison under another name (#12 review).
aggregate_fold <- function(entry, grouping) {
  test <- entry$fold$test_df
  group <- switch(
    grouping,
    country_year = paste(test$country_name, test$year_start,
                         test$insecticide_type)
  )
  ppd_aggregate(test$died, test$mosquito_number, group, entry$sims) %>%
    mutate(model = entry$fold$model,
           experiment = entry$fold$experiment,
           grouping = grouping,
           .before = everything())
}

# per-record scores summarised by the grouping columns in `...`; `excess` is
# mean squared error above the noise floor
by_group <- function(scores, ...) {
  scores %>%
    group_by(experiment, ..., model) %>%
    summarise(
      n = n(),
      mean_observed = mean(observed),
      mean_predicted = mean(predicted),
      bias = mean(predicted - observed),
      mean_pit = mean(pit),
      coverage_95 = mean(expected_coverage(cdf_below, pmf_at, 0.95)),
      crps = mean(crps),
      mse = mean((observed - predicted) ^ 2),
      mse_floor = noise_floor_mse(died, mosquito_number, rho_external),
      .groups = "drop"
    ) %>%
    mutate(excess = mse - mse_floor)
}


# pairing records across models ---------------------------------------------------

# the key identifying a held-out assay, so that predictions from different
# models are matched on the record rather than on row order
record_key <- function(x) {
  columns <- c("cell", "year_start", "insecticide_type", "died",
               "mosquito_number")
  # paste() drops a NULL argument silently, so a missing column would make two
  # keys agree on fewer fields, or on none at all
  stopifnot(all(columns %in% names(x)))
  do.call(paste, x[columns])
}

# Bootstrap a statistic of a data frame by resampling its pixels (`cell`), all
# of a pixel's rows together: bioassays cluster hard by pixel, and the assay
# count badly overstates the information a fold holds. With one column per
# model in `data`, every model is scored on the same resample, so differences
# formed within a replicate are paired. Returns n_bootstrap x
# length(statistic(data))
pixel_bootstrap <- function(data, statistic, n_bootstrap = 2000) {
  cells <- unique(data$cell)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell)
  t(replicate(n_bootstrap, {
    picked <- sample(cells, length(cells), replace = TRUE)
    statistic(data[unlist(rows_by_cell[as.character(picked)],
                          use.names = FALSE), ])
  }))
}
