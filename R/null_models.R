# Null model definitions used in the cross-validation experiments: the
# insecticide-type intercept model and the nearest neighbour heuristic.
#
# Separated from the scoring so that both the point-prediction
# scoring and the posterior predictive scoring (#10) use one definition of each
# null model. The lower half of this file adds predictive distributions for the
# nulls, because a proper scoring rule needs a distribution rather than a point
# from every candidate. Every model is scored at the same external,
# replicate-based overdispersion (see validation_metrics.R).
#
# The nulls do not fit an overdispersion of their own. They earn their place on
# point prediction - mean squared error and variance explained - and that needs
# only the predicted fraction, while the dynamical model's calibration is judged
# against held-out data directly. Fitting one rho per null made coverage a
# comparison of vagueness rather than of prediction: the intercept null reached
# nominal coverage by inflating its rho to 0.45 against an external estimate of
# 0.12-0.22 (#12 review).

# approximate non-zero and non-one proportions by applying the empirical logit
# transform to the data and then the inverse logit to provide a proportion
emplog_prop <- function(died, mosquito_number) {
  emplog <- log((died + 0.5) / (mosquito_number - died + 0.5))
  plogis(emplog)
}

# the nearest neighbour null, and predictive distributions for the nulls -----

# For each held-out record, average over the k nearest training records of the
# same insecticide from the most recent two years available.

# Pooled counts over the nearest neighbours of each held-out record - died and
# tested summed over them, rather than a single proportion, so that uncertainty
# in the predicted fraction can be represented by a beta posterior on those
# counts.
#
# Several neighbour counts at once. Neither the distance matrix nor the
# per-record ordering of those distances depends on k, so both are computed once
# and every k is read off the same cumulative sums. nn_oracle_draws() used to
# call this once per candidate k and pay for the full distance matrix and a sort
# per record fifteen times over (#12 review).
#
# `k_values` is a vector; the returned matrices have one row per held-out record
# and one column per element of it, in that order.
nn_counts_grid <- function(latitude,
                           longitude,
                           year,
                           insecticide_type,
                           training_data,
                           k_values,
                           n_years_prior = 1) {

  # A missing row in outputs/optimal_nn.csv returns numeric(0) from the lookup
  # in run_validation_folds.R, and every step below then degrades silently:
  # sort(x)[integer(0)] is numeric(0), `<= numeric(0)` is logical(0),
  # which(logical(0)) is integer(0), and sum() of no elements is 0. The result
  # is a pooled count of 0 died out of 0 tested for every held-out record, whose
  # Jeffreys posterior is Beta(0.5, 0.5) — so the "nearest neighbour null"
  # becomes a random number generator with no data in it, and nothing errors.
  # This is exactly what happened to the two spatial block folds, whose
  # experiment name had no row in that file: their reported score of -0.86 was
  # the score of the prior, not of a nearest neighbour model.
  stopifnot(
    is.numeric(k_values),
    length(k_values) >= 1,
    all(is.finite(k_values)),
    all(k_values >= 1)
  )

  training_coords <- as.matrix(training_data[, c("longitude", "latitude")])
  test_coords <- cbind(longitude, latitude)

  dists <- fields::rdist.earth(test_coords, training_coords,
                               miles = FALSE)

  # The year window is anchored at prediction time, not at the record's own
  # year: the most recent `n_years_prior` + 1 years of data that exist when the
  # prediction is made. For the spatial experiments training spans every year,
  # so the anchor is the record's own year and this is the familiar
  # {year, year - 1}. For a forecast the training data stop at the cut, so the
  # anchor is the last training year and every held-out record draws on the same
  # window — which is the situation a person forecasting from that cut is
  # actually in.
  #
  # Anchoring on the record's own year instead made the null weaker the further
  # ahead it had to predict, for an indexing reason rather than an information
  # one: a record one year past the cut saw the whole window, one five years past
  # saw a single year of it, or none at all.
  anchor <- pmin(year, max(training_data$year))
  year_diff <- -1 * seq(0, n_years_prior)

  n_test <- nrow(test_coords)
  total_died <- matrix(NA_real_, nrow = n_test, ncol = length(k_values))
  total_tested <- matrix(NA_real_, nrow = n_test, ncol = length(k_values))

  # which training records a record may draw on depends only on its anchor year
  # and its insecticide, so the masks are built once per distinct pair rather
  # than once per record
  mask_key <- paste(anchor, insecticide_type)
  unique_key <- unique(mask_key)
  key_index <- match(mask_key, unique_key)
  valid_masks <- lapply(match(unique_key, mask_key), function(j) {
    training_data$year %in% (anchor[j] + year_diff) &
      training_data$insecticide_type == insecticide_type[j]
  })

  training_died <- training_data$died
  training_tested <- training_data$mosquito_number

  for (i in seq_len(n_test)) {

    distance_vec <- dists[i, ]
    valid <- valid_masks[[key_index[i]]]
    masked_distance_vec <- ifelse(valid, distance_vec, Inf)

    # If no training record matches this record's insecticide and year window,
    # every masked distance is Inf, the threshold is Inf, and `<= threshold`
    # then selects the entire training set — every insecticide, every year. That
    # is not a nearest neighbour prediction, it is a global mean wearing one,
    # and it is silent. It bites when the holdout runs further past the cut than
    # `n_years_prior` reaches back: with a five-year forecast window and
    # `n_years_prior = 3`, the fourth and fifth lead years have no valid
    # training year at all. Fail here instead, so the caller has to set
    # `n_years_prior` to at least the window length.
    if (!any(valid)) {
      valid_years <- anchor[i] + year_diff
      stop("no training records in ", paste(range(valid_years), collapse = "-"),
           " for ", insecticide_type[i], ": the nearest neighbour null has ",
           "nothing to predict from")
    }

    # the k-th smallest masked distance is the threshold, and every record at or
    # inside it is pooled - so ties at the threshold all count, and a k larger
    # than the number of valid records gives an infinite threshold and pools the
    # whole training set, exactly as the per-k version did. Sorting once and
    # taking cumulative sums along that order makes each k an index lookup
    order_index <- order(masked_distance_vec)
    sorted_distance <- masked_distance_vec[order_index]
    cumulative_died <- cumsum(training_died[order_index])
    cumulative_tested <- cumsum(training_tested[order_index])

    # findInterval() on an ascending vector returns the number of elements at or
    # below the threshold, which is the size of the pooled set
    n_pooled <- findInterval(sorted_distance[k_values], sorted_distance)

    total_died[i, ] <- cumulative_died[n_pooled]
    total_tested[i, ] <- cumulative_tested[n_pooled]

  }

  list(total_died = total_died,
       total_tested = total_tested)

}

# the same for a single neighbour count, as a data frame of one column each
predict_null_fixed_nn_counts <- function(latitude,
                                         longitude,
                                         year,
                                         insecticide_type,
                                         training_data,
                                         n_nearest_neighbours,
                                         n_years_prior = 1) {

  stopifnot(length(n_nearest_neighbours) == 1)

  counts <- nn_counts_grid(
    latitude = latitude,
    longitude = longitude,
    year = year,
    insecticide_type = insecticide_type,
    training_data = training_data,
    k_values = n_nearest_neighbours,
    n_years_prior = n_years_prior
  )

  data.frame(total_died = counts$total_died[, 1],
             total_tested = counts$total_tested[, 1])

}

# posterior draws of the predicted fraction from pooled counts, under a
# Jeffreys beta prior
pooled_count_draws <- function(total_died, total_tested, n_draws) {
  n_obs <- length(total_died)
  draws <- rbeta(n_draws * n_obs,
                 rep(total_died + 0.5, each = n_draws),
                 rep(total_tested - total_died + 0.5, each = n_draws))
  matrix(draws, nrow = n_draws, ncol = n_obs)
}

# Predictive distribution of the insecticide-type intercept null: the fraction
# for each type has a beta posterior from the pooled training counts for that
# type
intercept_null_draws <- function(training_data, test_data, n_draws = 1000) {

  pooled <- aggregate(
    cbind(died, mosquito_number) ~ insecticide_type,
    data = training_data,
    FUN = sum
  )

  index <- match(test_data$insecticide_type, pooled$insecticide_type)
  p_draws <- pooled_count_draws(pooled$died[index],
                                pooled$mosquito_number[index],
                                n_draws)

  list(p_draws = p_draws,
       test_df = test_data)

}

# Predictive distribution of the nearest neighbour null.
#
# This is not a competitor model to be optimised, it is a stand-in for what a
# person would do without a model: read a site's value off the nearby recent
# surveys. Tuning it separately for each experiment would build a series of new
# models to compare against the one model under analysis, so it is not tuned at
# all. Two specifications are reported instead, and neither involves a chosen
# value:
#
#   the practice baseline - `n_neighbours = 1`, `n_years_prior = 1`: the single
#     nearest record in the most recent two years available at prediction time.
#     Short recency is justified a priori, since anyone running this
#     surveillance knows resistance is moving fast, and the k-curves agree.
#
#   the oracle bound - nn_oracle_draws() below: the same method at whichever
#     number of neighbours minimises its own error on the held-out records. It
#     answers a different question, how well any such rule could have done in
#     this test, and a bound is properly per-test.
#
# This replaces a tuned `n_neighbours` read from outputs/optimal_nn.csv. That
# was chosen on an internal test set of 100 records sampled at random from the
# training data, which sat a median of 0 km from their nearest usable neighbour
# - 26% of them at the same site - while the held-out records sit at 36 km
# (interpolation), 139 km (blocks) and 194 km (countries). Optimal k rises with
# the distance that has to be reached, so tuning at 0 km chose it too small and
# handicapped the baseline. The lookup also had no row for the block folds,
# which silently reduced the null there to its prior (#12).
nn_null_draws <- function(training_data, test_data, n_neighbours = 1,
                          n_years_prior = 1, n_draws = 1000) {

  counts <- predict_null_fixed_nn_counts(
    latitude = test_data$latitude,
    longitude = test_data$longitude,
    year = test_data$year,
    insecticide_type = test_data$insecticide_type,
    training_data = training_data,
    n_nearest_neighbours = n_neighbours,
    n_years_prior = n_years_prior
  )

  p_draws <- pooled_count_draws(counts$total_died,
                                counts$total_tested,
                                n_draws)

  list(p_draws = p_draws,
       test_df = test_data)

}


# The oracle bound: the nearest neighbour null at whichever number of
# neighbours minimises its own mean squared error on the held-out records.
#
# This is deliberately given hindsight the dynamical model is not given, so that
# no choice of neighbour count can be said to have handicapped the baseline. The
# selection is one scalar over thousands of records, so the optimism it buys is
# small, and it runs in the conservative direction for any claim that the
# dynamical model beats the null.
#
# Minimised on mean squared error because that is what the comparison is
# reported on - excess MSE and variance explained - so the bound is a bound on
# the statistic actually quoted. The coverage and CRPS of the returned fold are
# therefore at the MSE-optimal k, not at their own optima.
nn_oracle_draws <- function(training_data, test_data,
                            n_years_prior = 1, n_draws = 1000,
                            k_grid = c(1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 120,
                                       200, 300, 500, 800)) {

  observed <- test_data$died / test_data$mosquito_number

  # one pass over the grid: the distances and their ordering are shared by every
  # candidate k, and the winning k's counts are already in hand, so no further
  # neighbour search is needed to draw from it
  counts <- nn_counts_grid(
    latitude = test_data$latitude,
    longitude = test_data$longitude,
    year = test_data$year,
    insecticide_type = test_data$insecticide_type,
    training_data = training_data,
    k_values = k_grid,
    n_years_prior = n_years_prior
  )

  predicted <- emplog_prop(counts$total_died, counts$total_tested)
  # mean() per column rather than colMeans(): the two accumulate in different
  # orders and so can differ in the last bits, and these numbers are compared
  # against the per-k version this replaced
  mse <- apply((observed - predicted) ^ 2, 2, mean)

  best_index <- which.min(mse)
  best <- k_grid[best_index]
  if (best == max(k_grid)) {
    warning("the oracle neighbour count is at the top of the grid (", best,
            "); widen k_grid so the minimum is interior")
  }

  p_draws <- pooled_count_draws(counts$total_died[, best_index],
                                counts$total_tested[, best_index],
                                n_draws)

  list(p_draws = p_draws,
       test_df = test_data,
       n_neighbours = best,
       k_grid = k_grid,
       k_mse = mse)

}
