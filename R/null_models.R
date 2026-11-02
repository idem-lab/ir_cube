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

# the nearest neighbour null, and predictive distributions for the nulls -----

# The year the nearest neighbour null's window ends at, for a record of `year`.
# The window is anchored at prediction time, not at the record's own year: the
# most recent `n_years_prior` + 1 years of data that exist when the prediction
# is made. For the spatial experiments training spans every year, so the anchor
# is the record's own year and the window is the familiar {year, year - 1}.
# For a forecast the training data stop at the cut, so the anchor is the last
# training year and every held-out record draws on the same window — which is
# the situation a person forecasting from that cut is actually in.
#
# Anchoring on the record's own year instead made the null weaker the further
# ahead it had to predict, for an indexing reason rather than an information
# one: a record one year past the cut saw the whole window, one five years past
# saw a single year of it, or none at all.
nn_anchor_year <- function(year, training_years) {
  pmin(year, max(training_years))
}

# The nearest bioassay: for each held-out record, the single training record of
# the same insecticide in its year window (nn_anchor_year(), `n_years_prior` + 1
# years) at the smallest distance between assay coordinates. Ties - records at
# the same site - go to the most recent year, and any left (one site and year)
# to a pick drawn from the stream seeded by `seed`, so a fold always picks the
# same record. Returns, per held-out record, the chosen training row and its
# distance in km.
nearest_bioassay <- function(training_data, test_data, n_years_prior = 1,
                             seed) {

  dists <- fields::rdist.earth(
    as.matrix(test_data[, c("longitude", "latitude")]),
    as.matrix(training_data[, c("longitude", "latitude")]),
    miles = FALSE)

  # which training records a record may draw on depends only on its anchor year
  # and its insecticide, so the candidates are found once per distinct pair
  anchor <- nn_anchor_year(test_data$year, training_data$year)
  key <- paste(anchor, test_data$insecticide_type)
  candidates <- lapply(split(seq_along(key), key), function(rows) {
    j <- rows[1]
    which(training_data$year >= anchor[j] - n_years_prior &
            training_data$year <= anchor[j] &
            training_data$insecticide_type == test_data$insecticide_type[j])
  })

  saved <- if (exists(".Random.seed", globalenv())) .Random.seed
  on.exit(if (is.null(saved)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", saved, envir = globalenv()))
  set.seed(seed)

  chosen <- integer(nrow(test_data))
  distance <- numeric(nrow(test_data))
  for (i in seq_len(nrow(test_data))) {
    valid <- candidates[[key[i]]]
    # with no record of the insecticide in the window the null has nothing to
    # predict from: stop, rather than fall back to a wider set
    if (length(valid) == 0) {
      stop("no training records in ", anchor[i] - n_years_prior, "-",
           anchor[i], " for ", test_data$insecticide_type[i],
           ": the nearest bioassay null has nothing to predict from")
    }
    d <- dists[i, valid]
    nearest <- valid[d == min(d)]
    nearest <- nearest[training_data$year[nearest] ==
                         max(training_data$year[nearest])]
    if (length(nearest) > 1) {
      nearest <- nearest[sample.int(length(nearest), 1)]
    }
    chosen[i] <- nearest
    distance[i] <- min(d)
  }

  list(row = chosen, distance_km = distance)
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
# person would do without a model: read a site's value off the nearest recent
# bioassay, the single nearest record in the most recent two years available at
# prediction time (nearest_bioassay()). Tuning it would build a new model to
# compare against the one under analysis, so it is not tuned at all. Its
# fraction has a beta posterior from that bioassay's counts. Also returns the
# chosen bioassay's distance (km) and year per held-out record, which the CV
# figure's panel titles report.
nn_null_draws <- function(training_data, test_data, n_years_prior = 1,
                          n_draws = 1000, seed) {

  nearest <- nearest_bioassay(training_data, test_data,
                              n_years_prior = n_years_prior, seed = seed)

  p_draws <- pooled_count_draws(training_data$died[nearest$row],
                                training_data$mosquito_number[nearest$row],
                                n_draws)

  list(p_draws = p_draws,
       test_df = test_data,
       nearest_km = nearest$distance_km,
       nearest_year = training_data$year[nearest$row])

}
