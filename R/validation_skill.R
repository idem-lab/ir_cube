# Skill scores for the main cross-validation figure (#36): bioassay-based
# estimates against the models, for the level of mortality and for its change.
# Functions only; sourced by R/validation_level.R and R/validation_change.R.
# Sources R/validation_scoring.R.
#
#   level   100 (1 - MSE / Var(y)) against single held-out bioassays y: the %
#           of variance explained, as in R/variance_explained.R
#   change  100 (1 - sum w (P - A)^2 / sum w A^2), the skill against no change
#           (the persistence index; Kitanidis & Bras 1980), where A is the
#           change between a pair of single bioassays at one pixel and P its
#           prediction, both divided by the years between the two bioassays.
#           0 is as good as predicting no change; negative is worse
#
# Each bar also carries two quantities from the data floor F (the expected
# squared error of a perfect prediction of population mortality, from the
# assays' beta-binomial noise at rho-hat; noise_floor_record()):
#
#   ceiling  the score of a perfect prediction: 100 (1 - F / Var(y)) for level,
#            100 (1 - sum w F_pair / sum w A^2) for change, with F_pair the
#            pair's two floors over its gap squared. The dotted line
#   local    the expected score of an independent second bioassay at the same
#            pixel-year (level: 100 (1 - 2 F / Var(y)), the target's noise and
#            its own), or an independent second pair at the same pixel (change)
#
# Intervals resample pixels (pixel_bootstrap(), as R/variance_explained.R
# does). Replicate i also takes posterior draw i of rho by type, so the
# intervals on the ceiling and the local bars carry the uncertainty in rho as
# well as in the sample.

source("R/validation_scoring.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

n_bootstrap <- 2000


# The per-record scores carry the map and the data floor (#32); stop on
# scores written before that
require_map_columns <- function(scores) {
  stopifnot(all(c("map", "map_exact", "floor") %in% names(scores)))
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
# point estimate). Returns estimate, lower and upper per element
bootstrap_summary <- function(data, statistic, rho_draws) {
  point <- statistic(data, NULL)
  cells <- unique(data$cell)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell)
  replicates <- t(vapply(seq_len(nrow(rho_draws)), function(i) {
    picked <- sample(cells, length(cells), replace = TRUE)
    statistic(data[unlist(rows_by_cell[as.character(picked)],
                          use.names = FALSE), ], rho_draws[i, ])
  }, numeric(length(point))))
  tibble(bar = names(point), estimate = unname(point),
         lower = apply(replicates, 2, quantile, 0.025, na.rm = TRUE),
         upper = apply(replicates, 2, quantile, 0.975, na.rm = TRUE))
}

# the data floor per record, at rho-hat or at a draw of rho by type
record_floor <- function(data, rho = NULL) {
  if (is.null(rho)) return(data$floor)
  noise_floor_record(data$died, data$mosquito_number,
                     unname(rho[data$insecticide_type]))
}


# level ----------------------------------------------------------------------------

# Scores of every column of `data` named in `predictions` (bar label =
# column), plus the ceiling and the local bar. The floor's mean is over the
# assays that carry one (noise_floor_mse())
level_statistic <- function(predictions) {
  function(data, rho) {
    variance <- mean((data$observed - mean(data$observed)) ^ 2)
    floor <- mean(record_floor(data, rho), na.rm = TRUE)
    scores <- vapply(predictions, function(column) {
      100 * (1 - mean((data$observed - data[[column]]) ^ 2) / variance)
    }, numeric(1))
    c(ceiling = 100 * (1 - floor / variance),
      local = 100 * (1 - 2 * floor / variance),
      scores)
  }
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

# The expected squared error of an independent second pair at the same pixel,
# as a prediction of the target pair: the second pair's noise, taken to be
# the target's two assays' noise over a gap drawn from the targets' gaps (the
# weighted mean of 1 / gap^2), plus the target's own noise. Rewriting the
# target's floor F_pair = (F_b + F_a) / g^2 that way gives
# F_pair (1 + g^2 E[1 / gap^2])
local_pair_error <- function(pairs, floor) {
  inverse_square_gap <- sum(pairs$w / pairs$gap ^ 2) / sum(pairs$w)
  floor * (1 + pairs$gap ^ 2 * inverse_square_gap)
}

# Scores of every rate column of `pairs` named in `predictions` (a predicted
# change over the pair's gap, per year), plus the ceiling and the local pair.
# Weights w give each pixel equal weight. A column that is all NA scores NA
change_statistic <- function(predictions) {
  function(pairs, rho) {
    floor <- pair_floor(pairs, rho)
    total <- sum(pairs$w * pairs$rate ^ 2)
    scores <- vapply(predictions, function(column) {
      100 * (1 - sum(pairs$w * (pairs[[column]] - pairs$rate) ^ 2) / total)
    }, numeric(1))
    c(ceiling = 100 * (1 - sum(pairs$w * floor) / total),
      local = 100 * (1 - sum(pairs$w * local_pair_error(pairs, floor)) /
                       total),
      scores)
  }
}
