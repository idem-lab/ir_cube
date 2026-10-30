# Skill scores for the main cross-validation figure (#36): bioassay-based
# estimates against the models, for the level of mortality and for its change.
# Functions only; sourced by R/validation_level.R and R/validation_change.R.
# Sources R/validation_scoring.R.
#
#   level   100 (1 - MSE / Var(y)) against single held-out bioassays y: the %
#           of variance explained
#   change  100 (1 - sum w (P - A)^2 / sum w A^2), the skill against no change
#           (the persistence index; Kitanidis & Bras 1980), where A is the
#           change between a pair of single bioassays at one pixel and P its
#           prediction, both divided by the years between the two bioassays.
#           0 is as good as predicting no change; negative is worse
#
# Each bar also carries quantities from the data floor F (the expected
# squared error of a perfect prediction of population mortality, from the
# assays' beta-binomial noise at rho-hat; noise_floor_record()):
#
#   ceiling  the score of a perfect prediction: 100 (1 - F / Var(y)) for level,
#            100 (1 - sum w F_pair / sum w A^2) for change, with F_pair the
#            pair's two floors over its gap squared. The dotted line
#   local    the expected score of an independent second bioassay at the same
#            pixel-year (level: 100 (1 - 2 F / Var(y)), the target's noise and
#            its own), or an independent second pair at the same pixel (change).
#            With `pools`, also local_<n>: the mean of n such bioassays or
#            pairs, whose noise is 1/n of one's (level: 100 (1 - (1 + 1/n) F /
#            Var(y)), the ceiling less (100 - ceiling) / n)
#   n*_CV    for level, the number of bioassays the bar is worth (n_star_cv())
#
# Intervals resample pixels (cluster_bootstrap(), in R/validation_scoring.R).
# Replicate i also takes posterior draw i of rho by type, so the intervals on
# the ceiling and the local bars carry the uncertainty in rho as well as in the
# sample.

source("R/validation_scoring.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

n_bootstrap <- 2000


# The per-record scores carry the map and the data floor (#32); stop on
# scores written before that
require_map_columns <- function(scores) {
  stopifnot(all(c("map", "floor") %in% names(scores)))
  scores
}


# The joint posterior of rho by type (fig_illustrate_bioassay_variability.R),
# one row per bootstrap replicate
rho_replicates <- function(n = n_bootstrap,
                           file = "outputs/bioassay_rho_type_draws.rds") {
  if (!file.exists(file)) {
    stop("no rho draws at ", file,
         "; run R/fig_illustrate_bioassay_variability.R first")
  }
  draws <- readRDS(file)
  draws[sample(nrow(draws), n, replace = nrow(draws) < n), , drop = FALSE]
}

# The statistic at the point estimate and over pixel-bootstrap replicates, the
# replicate's rho draw passed as `rho` (a named vector by type; NULL at the
# point estimate). Returns estimate, lower and upper per element, and with
# n_star = TRUE (level statistics) n*_CV and its interval (n_star_cv())
bootstrap_summary <- function(data, statistic, rho_draws, n_star = FALSE) {
  point <- statistic(data, NULL)
  replicates <- cluster_bootstrap(data, statistic, rho_draws = rho_draws)
  out <- tibble(bar = names(point), estimate = unname(point),
                lower = apply(replicates, 2, quantile, 0.025, na.rm = TRUE),
                upper = apply(replicates, 2, quantile, 0.975, na.rm = TRUE))
  if (n_star) out <- bind_cols(out, n_star_cv(point, replicates))
  out
}

# the data floor per record, at rho-hat or at a draw of rho by type
record_floor <- function(data, rho = NULL) {
  if (is.null(rho)) return(data$floor)
  noise_floor_record(data$died, data$mosquito_number,
                     unname(rho[data$insecticide_type]))
}


# level ----------------------------------------------------------------------------

# the local bars' names: local for one bioassay or pair, local_<n> for the
# mean of n
local_names <- function(pools) {
  ifelse(pools == 1, "local", paste0("local_", pools))
}

# Scores of every column of `data` named in `predictions` (bar label =
# column), plus the ceiling and the local bars for each number of bioassays in
# `pools`. The floor's mean is over the assays that carry one
# (noise_floor_mse())
level_statistic <- function(predictions, pools = 1) {
  function(data, rho) {
    variance <- mean((data$observed - mean(data$observed)) ^ 2)
    floor <- mean(record_floor(data, rho), na.rm = TRUE)
    scores <- vapply(predictions, function(column) {
      100 * (1 - mean((data$observed - data[[column]]) ^ 2) / variance)
    }, numeric(1))
    c(ceiling = 100 * (1 - floor / variance),
      setNames(100 * (1 - (1 + 1 / pools) * floor / variance),
               local_names(pools)),
      scores)
  }
}

# n*_CV (#40): the number of bioassays at a held-out bioassay's pixel-year
# whose mean would predict it as accurately as a bar does, F / (D_m - F), with
# F the data floor and D_m the bar's mean squared error. From the level scores
# that is (100 - ceiling) / (ceiling - score). The out-of-sample check on the
# maps' n* (n_star() in R/two_stage_map_functions.R), with the same bioassay
# noise (rho-hat), but in bioassays of the held-out ones' sizes rather than of
# 100 mosquitoes.
#
# The interval inverts that of R = (ceiling - score) / (100 - ceiling), whose
# denominator is always positive: n* in [1 / R_upper, 1 / R_lower], infinite
# where R <= 0 (no number of bioassays matches the bar). Inverting a
# replicate's n* fails where its R crosses 0. NA for the ceiling, the perfect
# prediction; the local bar is 1 by construction
n_star_cv <- function(point, replicates) {
  excess <- function(x) (x[, "ceiling"] - x) / (100 - x[, "ceiling"])
  invert <- function(r) ifelse(r > 0, 1 / r, Inf)
  r <- excess(replicates)
  out <- tibble(
    n_star = invert(excess(t(point))[1, ]),
    n_star_lower = invert(apply(r, 2, quantile, 0.975, na.rm = TRUE)),
    n_star_upper = invert(apply(r, 2, quantile, 0.025, na.rm = TRUE)))
  out[names(point) == "ceiling", ] <- NA
  out
}


# change ---------------------------------------------------------------------------

# A pair's data floor: the two assays' floors over the gap squared
pair_floor <- function(pairs, rho = NULL) {
  if (is.null(rho)) return(pairs$floor)
  floor_b <- noise_floor_record(pairs$died_b, pairs$mosquito_number_b,
                                unname(rho[pairs$insecticide_type]))
  floor_a <- noise_floor_record(pairs$died_a, pairs$mosquito_number_a,
                                unname(rho[pairs$insecticide_type]))
  (floor_b + floor_a) / pairs$gap ^ 2
}

# The expected squared error of the mean of n independent pairs at the same
# pixel, as a prediction of the target pair: their noise, each pair's taken to
# be the target's two assays' noise over a gap drawn from the targets' gaps
# (the mean of 1 / gap^2 under weights w), over n, plus the target's own
# noise. Rewriting the target's floor F_pair = (F_b + F_a) / g^2 that way
# gives F_pair (1 + g^2 E[1 / gap^2] / n)
local_pair_error <- function(pairs, floor, w, n = 1) {
  inverse_square_gap <- sum(w / pairs$gap ^ 2) / sum(w)
  floor * (1 + pairs$gap ^ 2 * inverse_square_gap / n)
}

# Scores of every rate column of `pairs` named in `predictions` (a predicted
# change over the pair's gap, per year), plus the ceiling and the local pairs
# for each number of pairs in `pools`, each pair weighted by its column
# `weight` (giving each pixel-insecticide equal weight)
change_statistic <- function(predictions, weight = "w", pools = 1) {
  function(pairs, rho) {
    w <- pairs[[weight]]
    floor <- pair_floor(pairs, rho)
    total <- sum(w * pairs$rate ^ 2)
    scores <- vapply(predictions, function(column) {
      100 * (1 - sum(w * (pairs[[column]] - pairs$rate) ^ 2) / total)
    }, numeric(1))
    local <- vapply(pools, function(n) {
      100 * (1 - sum(w * local_pair_error(pairs, floor, w, n)) / total)
    }, numeric(1))
    c(ceiling = 100 * (1 - sum(w * floor) / total),
      setNames(local, local_names(pools)),
      scores)
  }
}
