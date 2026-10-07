# Out-of-sample variance explained, and the share of observed variance that
# bioassay sampling makes unexplainable, for every fitted fold and model.
#
# Two quantities, deliberately kept apart because they carry different kinds of
# uncertainty and the first must not depend on the second:
#
#   explained_i = 100 (1 - MSE_model / Var(y))
#       the fraction of observed variance in held-out mortality the model
#       accounts for. No noise floor enters it, so it does not inherit the
#       uncertainty in rho. Interval by resampling pixels, since bioassays
#       cluster hard by pixel and the assay count badly overstates the
#       information a fold holds.
#
#   noise = 100 floor / Var(y)
#       the share of that variance attributable to beta-binomial sampling in
#       the assay itself, which no model can explain. Interval from the
#       posterior of rho, estimated per insecticide type from replicate
#       bioassays (fig_illustrate_bioassay_variability.R).
#
#   ceiling = 100 - noise
#       the most any model could explain.
#
#   model floor = 100 U / Var(y)
#       for the two-stage model only, the part of its remainder below the
#       ceiling that it attributes to its site noise u and p (#32;
#       model_floor_mse() in validation_functions.R). The rest is the map's
#       own error.
#
# MSE is the map's (#32): every model's point prediction as the maps show it,
# from outputs/cv_scores.csv, so R/validation_metrics.R runs first.
#
# The floor is the sampling variance of the observed proportion, p(1-p) k_i,
# with k_i = (1 + (m_i - 1) rho_t) / m_i the beta-binomial design effect over
# the assay size, and it is noise_floor_mse() from validation_functions.R - the
# same estimator every other script uses. See the note above that function for
# why: it is exactly unbiased for a subset's mean sampling variance and assumes
# nothing about how p is distributed, which matters because p is bimodal for
# DDT and Alpha-cypermethrin.
#
# A constant rho overstates assay noise as mortality approaches 100%, so for the
# near-saturated insecticides the floor can exceed the observed variance however
# p is estimated - Fenitrothion in the interpolation fold does. Those cells are
# excluded from the per-insecticide figure by the rule below rather than papered
# over.
#
# Writes outputs/cv_variance_explained.csv (pooled per experiment),
# outputs/cv_variance_explained_by_fold.csv,
# outputs/cv_variance_explained_by_horizon.csv (forecasting, origins pooled)
# and outputs/cv_variance_explained_by_insecticide.csv.

source("R/validation_scoring.R")

set.seed(2026 - 9 - 24)
n_bootstrap <- 2000
n_posterior <- 4000

# the three out-of-sample experiments worth reporting. Leave-one-country-out is
# excluded: it measures the difficulty of an entirely unsampled country, which
# the deployed model never faces, and its held-out bias tracks the country's own
# fitted effect at r = -0.94 (#12 review).
experiments <- list(
  list(label = "spatial interpolation", experiment = "spatial_interpolation",
       folds = "all"),
  list(label = "spatial extrapolation", experiment = "spatial_blocks",
       folds = c("1", "2")),
  # the two five-year rolling origins, pooled into one bar as the two spatial
  # blocks are. The 2020 three-year fold they replace both leaked training
  # records into the holdout and landed on the one pause in twenty years of
  # decline, which was also the thinnest window.
  list(label = "temporal change", experiment = "temporal_forecasting",
       folds = c("2014", "2018"))
)

# each contrast is the first model less the second, in percentage points of the
# observed variance
contrasts <- list(
  c("two_stage", "dynamical"),
  c("dynamical", "nearest_neighbour"),
  c("dynamical", "nearest_neighbour_oracle"),
  c("dynamical", "intercept"),
  c("nearest_neighbour", "intercept"))

models <- c(dynamical = "dynamical model",
            two_stage = "two-stage model",
            nearest_neighbour = "nearest recent survey",
            nearest_neighbour_oracle = "nearest surveys, best k",
            intercept = "insecticide mean")

# the same per-type overdispersion the rest of the pipeline scores against.
# rho_table, not rho_lookup: the latter is the function that supplies it.
rho_spec <- rho_lookup()
rho_key <- rho_spec$key
rho_table <- rho_spec$table

# the joint posterior of rho by type, saved by the fit that produced it
rho_draws_file <- "outputs/bioassay_rho_type_draws.rds"
if (!file.exists(rho_draws_file)) {
  stop("no rho draws at ", rho_draws_file,
       "; run R/fig_illustrate_bioassay_variability.R first")
}
rho_draws <- readRDS(rho_draws_file)
stopifnot(all(rho_table$key %in% colnames(rho_draws)))
cat("overdispersion used for the noise share:", rho_spec$source, "\n")


# assemble the held-out records, with one prediction column per model ---------

# The predictions are each model's map (#32): the point prediction
# validation_metrics.R writes to outputs/cv_scores.csv as `map`, which for the
# two-stage model leaves out its site noise u and p and for every other model
# is its predictive mean. Read from there rather than from the saved folds,
# which take 4-8 GB each to load, so run validation_metrics.R first. It is the
# mean of the thinned draws every other score uses (thin_draws()), where the
# folds' means used all of them: on the September folds that moved the
# dynamical model's pooled spatial scores by at most 0.002 points.
#
# Predictions are matched on the record (record_key(), in
# validation_scoring.R). Every model must describe the reference's records (the
# dynamical model's), in the same order. They do, being built from one fold
# definition, but nothing else enforces it and a silent misalignment would swap
# predictions between assays (#12 review)
scores <- read.csv("outputs/cv_scores.csv", colClasses = c(fold = "character"))
stopifnot(all(c("map", "map_exact", "model_floor") %in% names(scores)))

read_fold <- function(model, experiment, fold, reference = NULL) {
  # the forecasting folds are scored under one experiment per origin
  label <- if (experiment == "temporal_forecasting") {
    paste(experiment, fold, sep = "_")
  } else {
    experiment
  }
  x <- scores[scores$model == model & scores$experiment == label &
                scores$fold == fold, ]
  stopifnot(nrow(x) > 0,
            is.null(reference) ||
              identical(record_key(x), record_key(reference)))
  list(test_df = x, predicted = x$map, model_floor = x$model_floor,
       map_exact = x$map_exact)
}

records <- bind_rows(lapply(experiments, function(spec) {
  bind_rows(lapply(spec$folds, function(fold) {

    dynamical <- read_fold("dynamical", spec$experiment, fold)
    reference <- dynamical$test_df
    two_stage <- read_fold("two_stage", spec$experiment, fold, reference)
    predictions <- lapply(setNames(names(models), names(models)), function(m) {
      if (m == "dynamical") return(dynamical$predicted)
      if (m == "two_stage") return(two_stage$predicted)
      read_fold(m, spec$experiment, fold, reference)$predicted
    })

    reference %>%
      transmute(experiment = spec$label, fold = fold,
                cell, insecticide_type, insecticide_class, year_start,
                died, mosquito_number,
                observed = died / mosquito_number) %>%
      bind_cols(as_tibble(setNames(predictions,
                                   paste0("p_", names(predictions))))) %>%
      # the two-stage model's floor U and whether its map is exact
      mutate(model_floor = two_stage$model_floor,
             map_exact = two_stage$map_exact)
  }))
}))

# every model must have a prediction for every record, or a bar would be drawn
# from a different denominator than its neighbours
prediction_columns <- paste0("p_", names(models))
stopifnot(all(prediction_columns %in% names(records)),
          !any(is.na(records[, prediction_columns])))

records <- records %>%
  left_join(rho_table %>% select(key, rho_type = rho),
            by = setNames("key", rho_key))
stopifnot(!any(is.na(records$rho_type)))

cat(sprintf("%i held-out records across %i experiments\n",
            nrow(records), n_distinct(records$experiment)))


# the two quantities ---------------------------------------------------------

# The noise floor is noise_floor_mse() from validation_functions.R, the same
# estimator the rest of the pipeline uses. It is y(1-y) k / (1-k), and since
# E[y(1-y)] = p(1-p)(1-k) the division makes it exactly unbiased for the mean
# sampling variance of a subset, with no assumption about how p is distributed.
#
# It replaces an empirical Bayes floor fitted here: a Beta(a_t, b_t) prior per
# insecticide type with each assay treated as 1/k independent draws. That was
# introduced because the plug-in is exactly zero for an assay reading 0% or
# 100%, which a large share of records for the near-saturated insecticides do.
# But the zero is a per-record artefact and only the subset mean is ever used,
# where the de-attenuation on the unsaturated records compensates for it
# exactly - which is what makes the plug-in unbiased. Simulation put the Beta
# version 1.4 to 9.6% high depending on rho, and it is the more fragile of the
# two where p is bimodal, as it is for DDT and Alpha-cypermethrin (#12 review).

# with, last, the two-stage model's floor U as a share of the same variance
# (#32), so that it is resampled with the scores
explained_for <- function(data) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  c(vapply(prediction_columns, function(column) {
    100 * (1 - mean((data$observed - data[[column]]) ^ 2) / variance)
  }, numeric(1)),
  model_floor = 100 * mean(data$model_floor) / variance)
}

# Posterior draws of the noise share, propagating rho only: this carries the
# uncertainty in the overdispersion but treats Var(y) as known, which is where
# essentially all of the floor's uncertainty sits.
#
# The draws come from the hierarchical fit itself rather than being rebuilt as
# independent logit-normals from each type's interval. The hierarchy correlates
# rho between types, so independent draws understate the uncertainty in a share
# averaged over types (#12 review).
noise_posterior <- function(data) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  index <- match(data$insecticide_type, colnames(rho_draws))
  stopifnot(!any(is.na(index)))
  vapply(sample(nrow(rho_draws), n_posterior, replace = nrow(rho_draws) < n_posterior),
         function(i) {
           100 * noise_floor_mse(data$died, data$mosquito_number,
                                 rho_draws[i, index]) / variance
         }, numeric(1))
}

summarise_subset <- function(data, label, stratum = NA_character_) {

  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  point <- explained_for(data)
  # resampling pixels (pixel_bootstrap(), in validation_scoring.R)
  replicates <- pixel_bootstrap(data, explained_for, n_bootstrap)
  noise <- noise_posterior(data)
  noise_point <- 100 * noise_floor_mse(data$died, data$mosquito_number,
                                       data$rho_type) / variance

  bind_rows(
    bind_rows(lapply(prediction_columns, function(column) {
      data.frame(
        quantity = models[sub("^p_", "", column)],
        kind = "model",
        estimate = point[[column]],
        lower = quantile(replicates[, column], 0.025, na.rm = TRUE),
        upper = quantile(replicates[, column], 0.975, na.rm = TRUE))
    })),
    # Paired contrasts, formed within each bootstrap replicate. Every model is
    # scored on the same resample, so the per-model intervals are strongly
    # correlated and whether two of them overlap says nothing about the
    # difference; this is the only thing that can settle a comparison. The
    # denominator cancels within a replicate, so the contrast is just
    # 100 (MSE_reference - MSE_model) / Var(y).
    bind_rows(lapply(contrasts, function(pair) {
      a <- paste0("p_", pair[1])
      b <- paste0("p_", pair[2])
      drawn <- replicates[, a] - replicates[, b]
      data.frame(quantity = paste(models[[pair[1]]], "-", models[[pair[2]]]),
                 kind = "contrast",
                 estimate = point[[a]] - point[[b]],
                 lower = quantile(drawn, 0.025, na.rm = TRUE),
                 upper = quantile(drawn, 0.975, na.rm = TRUE))
    })),
    data.frame(quantity = "bioassay variability", kind = "noise",
               estimate = noise_point,
               lower = quantile(noise, 0.025), upper = quantile(noise, 0.975)),
    # The two-stage model's floor (#32): the share of observed variance it
    # attributes to its site noise u and p, part of its remainder above the
    # noise share, as the model attributes it rather than as measured. In the
    # tables only, not the bar figure (#36). Its interval is the pixel
    # bootstrap's, with tau and sigma_p held at their fitted values. The share
    # of records whose map is exact (map_exact) says whether U was taken at
    # the map or at the predictive mean
    data.frame(quantity = "two-stage model floor (u + p)",
               kind = "model floor",
               estimate = point[["model_floor"]],
               lower = quantile(replicates[, "model_floor"], 0.025),
               upper = quantile(replicates[, "model_floor"], 0.975),
               share_map_exact = mean(data$map_exact))
  ) %>%
    mutate(experiment = label, stratum = stratum,
           assays = nrow(data), pixels = n_distinct(data$cell),
           variance = variance, .before = everything())
}

pooled <- bind_rows(lapply(split(records, records$experiment), function(data) {
  summarise_subset(data, data$experiment[1])
}))
row.names(pooled) <- NULL
write.csv(pooled, "outputs/cv_variance_explained.csv", row.names = FALSE)

# and per fold, for experiments built from more than one, so that a pooled bar
# can be checked against the folds it pools
by_fold <- bind_rows(lapply(
  split(records, paste(records$experiment, records$fold)), function(data) {
    if (n_distinct(records$fold[records$experiment == data$experiment[1]]) < 2) {
      return(NULL)
    }
    summarise_subset(data, data$experiment[1], data$fold[1])
  }))
row.names(by_fold) <- NULL
write.csv(by_fold, "outputs/cv_variance_explained_by_fold.csv", row.names = FALSE)

cat("\nper experiment, % of observed variance in held-out mortality:\n")
print(as.data.frame(pooled %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(experiment, assays, pixels, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)


cat("\npaired contrasts, within bootstrap replicates:\n")
print(as.data.frame(pooled %>%
  filter(kind == "contrast") %>%
  mutate(value = sprintf("%6.1f [%6.1f, %6.1f]", estimate, lower, upper)) %>%
  select(experiment, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)

cat("\nper fold within experiment:\n")
print(as.data.frame(by_fold %>%
  filter(kind != "contrast") %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(experiment, fold = stratum, assays, pixels, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)

cat("\npaired contrasts per fold:\n")
print(as.data.frame(by_fold %>%
  filter(kind == "contrast") %>%
  mutate(value = sprintf("%6.1f [%6.1f, %6.1f]", estimate, lower, upper)) %>%
  select(experiment, fold = stratum, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)


# and by insecticide ---------------------------------------------------------

# Every insecticide with held-out records is scored, and nothing is dropped
# from the table. What the figure shows is decided afterwards, by whether the
# cell can carry a reading at all, on two counts that are not substitutes:
#
#   ceiling >= min_ceiling
#       more than half the observed spread in held-out mortality has to be real
#       variation in resistance rather than assay noise. Where it is not, the
#       ratio is two small numbers divided by each other: Bendiocarb in the
#       interpolation fold has 154 assays in 67 pixels and a predictable
#       standard deviation of about 6 points of mortality, so a model would
#       have to land inside 6 points to score above zero. A sample-size rule
#       cannot see this, which is why the 20-assay threshold this replaces let
#       that cell through at -158%.
#
#   interval width <= max_ci_width
#       the estimate has to be resolvable, for both plotted models, or the
#       comparison the panel exists to make cannot be read. This is where the
#       pixel count enters, but only through its effect on precision:
#       Lambda-cyhalothrin in the forecast fold has 25 points of predictable
#       variation and 18 pixels to estimate it from.
min_ceiling <- 50
max_ci_width <- 100
shown_models <- c("dynamical model", "nearest recent survey")

by_insecticide <- bind_rows(lapply(
  split(records, paste(records$experiment, records$insecticide_type)),
  function(data) {
    # a single pixel leaves no between-cluster variation to bootstrap, and a
    # constant holdout leaves no denominator
    if (n_distinct(data$cell) < 2) return(NULL)
    if (mean((data$observed - mean(data$observed)) ^ 2) == 0) return(NULL)
    summarise_subset(data, data$experiment[1], data$insecticide_type[1])
  }))

by_insecticide <- by_insecticide %>%
  group_by(experiment, stratum) %>%
  mutate(
    ceiling = 100 - estimate[kind == "noise"],
    ci_width = max((upper - lower)[quantity %in% shown_models]),
    drop_reason = case_when(
      ceiling < min_ceiling & ci_width > max_ci_width ~ "no signal, imprecise",
      ceiling < min_ceiling                           ~ "no signal",
      ci_width > max_ci_width                         ~ "imprecise",
      TRUE                                            ~ NA_character_),
    shown = is.na(drop_reason)) %>%
  ungroup()
row.names(by_insecticide) <- NULL
write.csv(by_insecticide,
          "outputs/cv_variance_explained_by_insecticide.csv", row.names = FALSE)

cells <- by_insecticide %>% distinct(experiment, stratum, assays, pixels,
                                     variance, ceiling, ci_width, shown,
                                     drop_reason)
cat(sprintf("\nby insecticide: %i experiment-insecticide combinations, %i shown\n",
            nrow(cells), sum(cells$shown)))
cat(sprintf("thresholds: ceiling >= %i%%, 95%% interval width <= %i points\n",
            min_ceiling, max_ci_width))
print(as.data.frame(cells %>%
  transmute(experiment, insecticide = stratum, assays, pixels,
            sd_pp = round(100 * sqrt(variance), 1),
            ceiling = round(ceiling, 1),
            predictable_sd_pp = round(100 * sqrt(variance) *
                                        sqrt(pmax(ceiling, 0) / 100), 1),
            ci_width = round(ci_width, 1),
            shown = ifelse(shown, "yes", drop_reason)) %>%
  arrange(desc(shown == "yes"), experiment, insecticide)), row.names = FALSE)


# and by forecast horizon --------------------------------------------------

# The two forecasting origins pooled, by years ahead of the last training year
# (training is year_start < the fold's cut year, so the cut year is one year
# ahead). A pixel's records from both origins move together in the bootstrap.
# Last in the script, so the bootstraps above draw the random numbers they did
# before it was added
forecasts <- records %>%
  filter(experiment == "temporal change") %>%
  mutate(horizon = year_start - as.integer(fold) + 1)
by_horizon <- bind_rows(lapply(split(forecasts, forecasts$horizon),
                               function(data) {
  summarise_subset(data, data$experiment[1], as.character(data$horizon[1]))
}))
row.names(by_horizon) <- NULL
write.csv(by_horizon, "outputs/cv_variance_explained_by_horizon.csv",
          row.names = FALSE)
