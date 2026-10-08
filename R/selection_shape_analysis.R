# Paired-difference analysis of the shape of the covariate effects on selection
# (#23).
#
# The model's selection recursion is additive on the logit scale:
#   logit q_t2 - logit q_t1 = -sum_{r = t1 + 1}^{t2} log w_r,
#   w_r = 1 + x_r' b,  b >= 0,
# where q is the fraction susceptible (= expected bioassay mortality). So for
# repeat bioassays at the same pixel and insecticide, the per-year change in
# logit resistance is the mean of log w_r over the interval between them, with
# no dependence on the initial state. This fits that relationship directly to
# empirical logits, with the model's covariates, and asks whether adding
# saturating hinge bases min(x, k) (non-negative coefficients, so the effect
# stays non-decreasing and concave) improves cross-validated fit over the
# linear terms the model uses now.
#
# Outputs:
#   outputs/selection_shape_cv.csv      cross-validated comparison of variants
#   outputs/selection_shape_coefs.csv   coefficients of the variants, all data
#   outputs/selection_shape_net_slope.csv  univariate net-use slope check
#   outputs/selection_shape_net_profile.csv    scores with the net-use
#     coefficient fixed at fractions of the saved model's value
#   outputs/selection_shape_net_bootstrap.csv  block bootstrap of the net-use
#     coefficient, on all pairs and subsets
#   outputs/selection_shape_net_confounding.csv  how much net-use variation
#     the pairs carry, and its collinearity with the other columns

library(tidyverse)
library(terra)
source("R/functions.R")

# build df and the cell-year covariates exactly as fit_model.R does, by
# evaluating its code up to the design matrix (skipping its sourcing of
# R/packages.R and R/functions.R)
fit_model_exprs <- parse("R/fit_model.R")
fit_model_text <- vapply(fit_model_exprs,
                         function(e) paste(deparse(e), collapse = " "),
                         "")
last_expr <- which(startsWith(fit_model_text, "x_cell_years <-"))
stopifnot(length(last_expr) == 1)
for (i in seq_len(last_expr)) {
  if (grepl("^source\\(\"R/(packages|functions).R\"\\)", fit_model_text[i])) {
    next
  }
  eval(fit_model_exprs[[i]], envir = globalenv())
}

# the trended designs (selection_design()): population by the encounter
# transform of density (the default), d / (d + d_half), raw or log, times
# g_dom, and the crops times g_ag, both linear, 0 in 1995 and 1 in 2025; d_half
# 50, the value of this analysis
x_trend <- map(c("encounter", "saturating", "raw", "log"), function(pop) {
  selection_design_matrix(unique_cells, baseline_year, final_data_year,
                          selection_design(pop = pop,
                                           pop_d_half = 50))$x_cell_years
})
x_trend <- do.call(cbind, c(x_trend[1],
                            map(x_trend[-1], ~ .x[, 3, drop = FALSE])))
trend_pop_names <- c("pop_enc:g_dom", "pop_sat:g_dom", "pop:g_dom",
                     "log_pop:g_dom")
trend_crop_names <- grep(":g_ag$", colnames(x_trend), value = TRUE)
stopifnot(length(trend_crop_names) == 10,
          all(trend_pop_names %in% colnames(x_trend)))
x_trend <- x_trend[, c(trend_pop_names, trend_crop_names)]

# the design of the saved fit, without the trend
selection <- selection_design_matrix(unique_cells, baseline_year,
                                     final_data_year,
                                     selection_design_untrended())
stopifnot(identical(selection$cell_years_index, cell_years_index))
x_cell_years <- selection$x_cell_years
rm(selection)

# covariates for all cell-years with data, in the model's column order
covs <- bind_cols(cell_years_index, as_tibble(x_cell_years),
                  as_tibble(x_trend))
crop_names <- setdiff(colnames(x_cell_years), c("nets", "irs", "pop"))

# log population density, min-max scaled over the whole cube (2000-2030) as
# prep_rasters.R would with its log transform switched on. This is the column
# proposed to replace pop. The 2031-2049 projections that prep_rasters.R also
# scales over are not in data/clean; they move the scaling slightly, not the
# ranking.
pop_raw <- c(rast("data/clean/pop_cube.tif"),
             rast("data/clean/pop_cube_future.tif"))
log_pop_range <- range(log(global(pop_raw, "range", na.rm = TRUE)))
pop_raw <- pre_pad_cube(pop_raw[[paste0("pop_", 2000:final_data_year)]],
                        baseline_year)
pop_raw_extract <- terra::extract(pop_raw, unique_cells) %>%
  mutate(cell_id = seq_along(unique_cells)) %>%
  pivot_longer(-cell_id,
               names_prefix = "pop_",
               names_to = "year",
               values_to = "pop_raw") %>%
  mutate(year_id = as.numeric(year) - baseline_year + 1) %>%
  select(cell_id, year_id, pop_raw)
covs <- covs %>%
  left_join(pop_raw_extract, by = c("cell_id", "year_id")) %>%
  mutate(
    log_pop = (log(pop_raw) - log_pop_range[1]) / diff(log_pop_range),
    # a constant, for selection not tied to any covariate
    one = 1
  )
stopifnot(!anyNA(covs$log_pop))

# check the raw population reproduces the model's scaled column (a linear map,
# up to the per-year relevelling in prep_rasters.R)
stopifnot(cor(covs$pop_raw, covs$pop) > 0.999)


# data: one record per pixel, insecticide type and year

rho_type <- read_csv("outputs/bioassay_rho_hierarchical.csv",
                     show_col_types = FALSE) %>%
  select(insecticide_type, rho)

# empirical logit of mortality with a continuity correction, and its
# approximate variance: binomial variance on the logit scale, inflated for
# beta-binomial overdispersion within each assay, and combining the assays
# pooled into the record
records <- df %>%
  left_join(rho_type, by = "insecticide_type") %>%
  # one location and country per pixel
  group_by(cell_id) %>%
  mutate(latitude = first(latitude),
         longitude = first(longitude),
         country_name = first(country_name)) %>%
  group_by(cell_id, cell, latitude, longitude, country_name,
           insecticide_class, insecticide_type, type_id, year_id) %>%
  summarise(
    died = sum(died),
    n = sum(mosquito_number),
    # sum of n_i (1 + (n_i - 1) rho) over pooled assays
    n_eff_factor = sum(mosquito_number * (1 + (mosquito_number - 1) * rho)),
    .groups = "drop"
  ) %>%
  mutate(
    p_hat = (died + 0.5) / (n + 1),
    logit_mort = qlogis(p_hat),
    # var(logit p_hat) ~= var(p_hat) / (p (1 - p))^2, where
    # var(p_hat) = p (1 - p) sum n_i (1 + (n_i - 1) rho) / n^2
    var_logit = n_eff_factor / (n ^ 2 * p_hat * (1 - p_hat))
  )
stopifnot(!anyDuplicated(records[, c("cell_id", "type_id", "year_id")]))

# consecutive pairs of records at the same pixel and insecticide type. Using
# only consecutive pairs keeps each record in at most two differences.
pairs <- records %>%
  arrange(cell_id, type_id, year_id) %>%
  group_by(cell_id, type_id) %>%
  mutate(
    year_1 = lag(year_id),
    logit_mort_1 = lag(logit_mort),
    var_logit_1 = lag(var_logit),
    p_hat_1 = lag(p_hat)
  ) %>%
  ungroup() %>%
  filter(!is.na(year_1)) %>%
  transmute(
    cell_id, latitude, longitude, country_name,
    insecticide_class, insecticide_type, type_id,
    year_1,
    year_2 = year_id,
    dt = year_2 - year_1,
    p_hat_1,
    p_hat_2 = p_hat,
    # per-year increase in logit resistance = decrease in logit mortality
    y = (logit_mort_1 - logit_mort) / dt,
    var_y = (var_logit_1 + var_logit) / dt ^ 2
  ) %>%
  mutate(pair_id = row_number())

# cell-years in each interval: the fitness applied between the two records is
# that of years year_1 + 1, ..., year_2
pair_years <- pairs %>%
  select(pair_id, cell_id, year_1, year_2) %>%
  mutate(year_id = map2(year_1 + 1, year_2, seq)) %>%
  unnest(year_id) %>%
  left_join(covs, by = c("cell_id", "year_id")) %>%
  select(-year_1, -year_2)
stopifnot(!anyNA(pair_years))

# interval means of the covariates, for the univariate checks
pair_means <- pair_years %>%
  group_by(pair_id) %>%
  summarise(across(c(nets, irs, pop, log_pop, all_of(crop_names)), mean))
pairs <- pairs %>%
  left_join(pair_means, by = "pair_id")

pairs %>%
  count(insecticide_class, insecticide_type) %>%
  print()


# net-use slope checks. The figure quoted in #23 (-0.01 [-0.37, +0.34], model
# +0.61) used the 2014 forecasting fold's held-out pyrethroid assays paired with
# training assays at the same pixel, and that fold's posterior; it needs the
# fold's draws, so is not recomputed here. These are the same regression on all
# the data: first pyrethroid pairs 3-8 years apart (all pairs, not only
# consecutive ones), unweighted regression of the per-year change in logit
# resistance on net use averaged over the years y1, ..., y2, with a bootstrap
# over pixels

net_slope_pairs <- records %>%
  filter(insecticide_class == "Pyrethroids") %>%
  select(cell_id, type_id, year_id, logit_mort, var_logit) %>%
  inner_join(., ., by = c("cell_id", "type_id"),
             suffix = c("_1", "_2"),
             relationship = "many-to-many") %>%
  mutate(dt = year_id_2 - year_id_1) %>%
  filter(dt >= 3, dt <= 8) %>%
  mutate(
    y = (logit_mort_1 - logit_mort_2) / dt,
    var_y = (var_logit_1 + var_logit_2) / dt ^ 2,
    pair_id = row_number()
  )
net_slope_pairs <- net_slope_pairs %>%
  select(pair_id, cell_id, year_id_1, year_id_2) %>%
  mutate(year_id = map2(year_id_1, year_id_2, seq)) %>%
  unnest(year_id) %>%
  left_join(covs, by = c("cell_id", "year_id")) %>%
  group_by(pair_id) %>%
  summarise(nets = mean(nets)) %>%
  right_join(net_slope_pairs, by = "pair_id")

# the same regression on the consecutive pairs used below, weighted by the
# inverse approximate variance, with net use averaged over the interval the
# model applies (years y1 + 1, ..., y2)
net_slope_consecutive <- pairs %>%
  filter(insecticide_class == "Pyrethroids")

cluster_boot_slope <- function(data, weighted, n_boot = 1000, seed = 11) {
  fit_slope <- function(d) {
    w <- if (weighted) 1 / d$var_y else NULL
    unname(coef(lm(y ~ nets, data = d, weights = w))[2])
  }
  set.seed(seed)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell_id)
  boots <- replicate(n_boot, {
    cells_boot <- sample(names(rows_by_cell), replace = TRUE)
    fit_slope(data[unlist(rows_by_cell[cells_boot]), ])
  })
  tibble(
    slope = fit_slope(data),
    lower = quantile(boots, 0.025),
    upper = quantile(boots, 0.975),
    n_pairs = nrow(data),
    n_cells = n_distinct(data$cell_id)
  )
}

net_slope <- bind_rows(
  cluster_boot_slope(net_slope_pairs, weighted = FALSE) %>%
    mutate(pairs = "all pairs 3-8 years apart, unweighted"),
  cluster_boot_slope(net_slope_pairs, weighted = TRUE) %>%
    mutate(pairs = "all pairs 3-8 years apart, inverse-variance weighted"),
  cluster_boot_slope(net_slope_consecutive, weighted = FALSE) %>%
    mutate(pairs = "consecutive pairs, unweighted"),
  cluster_boot_slope(net_slope_consecutive, weighted = TRUE) %>%
    mutate(pairs = "consecutive pairs, inverse-variance weighted")
) %>%
  relocate(pairs)
print(net_slope)


# model-implied per-year change for the same consecutive pairs, under the
# saved full-data fit (posterior mean over draws of the implied rate)

fitted <- new.env()
load("temporary/fitted_model.RData", envir = fitted)
stopifnot(identical(fitted$types, types),
          identical(fitted$classes_index, classes_index),
          identical(colnames(fitted$x_cell_years), colnames(x_cell_years)))
# the selection effects of 500 draws, from the fit's variables in whichever
# parameterisation it sampled (dynamical_parameter_draws())
source("R/dynamical_predictions.R")
fold <- list(draws = fitted$draws, options = fitted$model_options)
n_total <- nrow(as.matrix(fitted$draws))
parameters <- dynamical_parameter_draws(
  fold, classes_index, types,
  draw_index = round(seq(1, n_total, length.out = 500)))
n_draws <- parameters$n_draws
n_covs <- ncol(x_cell_years)
n_types <- length(types)

x_pair_years <- as.matrix(pair_years[, colnames(x_cell_years)])
pair_type <- pairs$type_id[pair_years$pair_id]
model_rate <- matrix(NA, n_draws, nrow(pairs))
model_effect <- matrix(0, n_covs, n_types,
                       dimnames = list(colnames(x_cell_years), types))
for (s in seq_len(n_draws)) {
  effect_type <- matrix(parameters$effect_type[s, , ], n_covs, n_types)
  log_w <- log1p(x_pair_years %*% effect_type)
  log_w <- log_w[cbind(seq_along(pair_type), pair_type)]
  model_rate[s, ] <- rowsum(log_w, pair_years$pair_id)[, 1] / pairs$dt
  model_effect <- model_effect + effect_type / n_draws
}
pairs$model_rate <- colMeans(model_rate)
rm(fitted, parameters, model_rate)

net_slope <- bind_rows(
  net_slope,
  cluster_boot_slope(pairs %>%
                       filter(insecticide_class == "Pyrethroids") %>%
                       mutate(y = model_rate),
                     weighted = TRUE) %>%
    mutate(pairs = "consecutive pairs, weighted, full-data model's implied rate")
)
print(net_slope)
write_csv(net_slope, "outputs/selection_shape_net_slope.csv")


# fit w = 1 + f(x)' b per insecticide class to the paired differences

# knots for the hinge bases min(x, k): quartiles of the covariate over the
# pair-years (for IRS, of its non-zero values; two thirds of pair-years have
# none)
knot_probs <- c(0.25, 0.5, 0.75)
knots <- list(
  nets = unname(quantile(pair_years$nets, knot_probs)),
  irs = unname(quantile(pair_years$irs[pair_years$irs > 0], knot_probs)),
  pop = unname(quantile(pair_years$pop, knot_probs))
)
print(knots)

# variants: linear terms, hinge covariates and knot indices, and whether a
# constant (non-positive) change in logit resistance per year is allowed, as
# the reversion in #24 would add. log_pop spans 0.48-0.97 at the data, so it
# acts mostly as a constant; the _const variants add an explicit constant
# column to separate the two.
base_pop <- c("nets", "irs", "pop", crop_names)
base_log_pop <- c("nets", "irs", "log_pop", crop_names)
base_const <- c("one", "nets", "irs", "log_pop", crop_names)
variants <- list(
  linear_pop = list(linear = base_pop),
  linear_pop_const = list(linear = c("one", base_pop)),
  const_no_pop = list(linear = c("one", "nets", "irs", crop_names)),
  pop_k123 = list(linear = base_pop, hinge = list(pop = 1:3)),
  pop_k123_nets_k123 = list(linear = base_pop,
                            hinge = list(pop = 1:3, nets = 1:3)),
  linear = list(linear = base_log_pop),
  linear_const = list(linear = base_const),
  nets_k1 = list(hinge = list(nets = 1)),
  nets_k2 = list(hinge = list(nets = 2)),
  nets_k3 = list(hinge = list(nets = 3)),
  nets_k123 = list(hinge = list(nets = 1:3)),
  irs_k1 = list(hinge = list(irs = 1)),
  irs_k2 = list(hinge = list(irs = 2)),
  irs_k3 = list(hinge = list(irs = 3)),
  irs_k123 = list(hinge = list(irs = 1:3)),
  both_k123 = list(hinge = list(nets = 1:3, irs = 1:3)),
  const_nets_k1 = list(linear = base_const, hinge = list(nets = 1)),
  const_nets_k123 = list(linear = base_const, hinge = list(nets = 1:3)),
  const_irs_k123 = list(linear = base_const, hinge = list(irs = 1:3)),
  linear_reversion = list(reversion = TRUE),
  const_reversion = list(linear = base_const, reversion = TRUE),
  nets_k123_reversion = list(hinge = list(nets = 1:3), reversion = TRUE),
  # the trended designs: population and crops only as products with the
  # linear trend
  trend_pop_enc = list(linear = c("nets", "irs", "pop_enc:g_dom",
                                  trend_crop_names)),
  trend_pop_sat = list(linear = c("nets", "irs", "pop_sat:g_dom",
                                  trend_crop_names)),
  trend_pop_raw = list(linear = c("nets", "irs", "pop:g_dom",
                                  trend_crop_names)),
  trend_pop_log = list(linear = c("nets", "irs", "log_pop:g_dom",
                                  trend_crop_names))
)
variants <- map(variants, function(v) {
  list(linear = v$linear %||% base_log_pop,
       hinge = v$hinge %||% list(),
       reversion = v$reversion %||% FALSE)
})

# basis matrix over the pair-year rows
make_basis <- function(variant) {
  x <- as.matrix(pair_years[, variant$linear])
  for (cov in names(variant$hinge)) {
    for (k in knots[[cov]][variant$hinge[[cov]]]) {
      x <- cbind(x, pmin(pair_years[[cov]], k))
      colnames(x)[ncol(x)] <- sprintf("min(%s, %.3f)", cov, k)
    }
  }
  x
}

# weighted least squares for the interval mean of log w, with b >= 0 and the
# reversion constant <= 0, by L-BFGS-B with an analytic gradient. `rows`
# indexes pairs.
fit_pairs <- function(x, rows, weights, reversion, fixed = NULL) {
  keep <- pair_years$pair_id %in% rows
  x <- x[keep, , drop = FALSE]
  group <- match(pair_years$pair_id[keep], rows)
  dt <- pairs$dt[rows]
  y <- pairs$y[rows]
  w <- weights[rows]
  n_b <- ncol(x)
  predict_rate <- function(par) {
    rate <- rowsum(log1p(x %*% par[seq_len(n_b)]), group)[, 1] / dt
    if (reversion) rate <- rate + par[n_b + 1]
    rate
  }
  objective <- function(par) sum(w * (y - predict_rate(par)) ^ 2)
  gradient <- function(par) {
    resid <- y - predict_rate(par)
    d_rate <- rowsum(x / c(1 + x %*% par[seq_len(n_b)]), group) / dt
    if (reversion) d_rate <- cbind(d_rate, 1)
    -2 * colSums(w * resid * d_rate)
  }
  n_par <- n_b + reversion
  lower <- c(rep(0, n_b), if (reversion) -Inf)
  upper <- c(rep(Inf, n_b), if (reversion) 0)
  # coefficients fixed at given values, as named elements of `fixed`
  fixed_idx <- match(names(fixed), colnames(x))
  lower[fixed_idx] <- upper[fixed_idx] <- fixed
  starts <- list(rep(0.05, n_par), rep(0.5, n_par))
  fits <- map(starts, function(start) {
    if (reversion) start[n_par] <- 0
    start[fixed_idx] <- fixed
    optim(start, objective, gradient,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = list(maxit = 5000, factr = 1e5))
  })
  best <- fits[[which.min(map_dbl(fits, "value"))]]
  par <- best$par
  names(par) <- c(colnames(x), if (reversion) "reversion")
  list(par = par, convergence = best$convergence)
}

# predict for any pairs from fitted parameters
predict_pairs <- function(x, par, rows, reversion) {
  keep <- pair_years$pair_id %in% rows
  group <- match(pair_years$pair_id[keep], rows)
  n_b <- ncol(x)
  rate <- rowsum(log1p(x[keep, , drop = FALSE] %*% par[seq_len(n_b)]),
                 group)[, 1] / pairs$dt[rows]
  if (reversion) rate <- rate + par[n_b + 1]
  rate
}

class_rows <- split(pairs$pair_id, pairs$insecticide_class)

# variance of each pair: the approximate sampling variance scaled by phi, plus
# an extra variance sigma2 per record for between-sample variation that the
# beta-binomial rho does not cover. phi < 1 is expected: the logit variance
# approximation overstates the scatter of records at or near 100% mortality.
# Both are estimated per class by maximum likelihood under the linear variant,
# alternating with the coefficients, then held fixed for every variant.
pair_var <- function(phi, sigma2) {
  phi * pairs$var_y + 2 * sigma2 / pairs$dt ^ 2
}
x_linear <- make_basis(variants$linear)
dispersion_class <- map_dfr(class_rows, function(rows) {
  par <- c(0, 0.1)
  for (iter in 1:5) {
    v <- pair_var(par[1], par[2])
    fit <- fit_pairs(x_linear, rows, 1 / v, reversion = FALSE)
    resid <- pairs$y[rows] - predict_pairs(x_linear, fit$par, rows, FALSE)
    par <- optim(par, function(p) {
      v <- pair_var(p[1], p[2])[rows]
      -sum(dnorm(resid, 0, sqrt(v), log = TRUE))
    }, method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(100, 20))$par
  }
  tibble(phi = par[1], sigma2 = par[2])
}, .id = "insecticide_class")
print(dispersion_class)
pair_dispersion <- dispersion_class[match(pairs$insecticide_class,
                                          dispersion_class$insecticide_class), ]
pairs$var_total <- pair_var(pair_dispersion$phi, pair_dispersion$sigma2)
pair_weights <- 1 / pairs$var_total

# full-data fits of every variant, per class
bases <- map(variants, make_basis)
full_fits <- imap(variants, function(variant, name) {
  map(class_rows, function(rows) {
    fit_pairs(bases[[name]], rows, pair_weights, variant$reversion)
  })
})
stopifnot(all(unlist(map_depth(full_fits, 2, "convergence")) == 0))

# mean squared standardised residual under the linear variant (close to 1
# if the variance model is calibrated)
imap_dbl(class_rows, function(rows, class) {
  pred <- predict_pairs(bases$linear, full_fits$linear[[class]]$par, rows,
                        FALSE)
  mean((pairs$y[rows] - pred) ^ 2 / pairs$var_total[rows])
}) %>%
  print()

coefs <- imap_dfr(full_fits, function(fits, name) {
  imap_dfr(fits, function(fit, class) {
    tibble(variant = name, insecticide_class = class,
           term = names(fit$par), estimate = unname(fit$par))
  })
})
write_csv(coefs, "outputs/selection_shape_coefs.csv")


# cross-validation by blocks of pixels: 1-degree squares, all pairs in a
# square (every insecticide) held out together, 10 folds, 5 repeats

pairs$block <- paste(floor(pairs$latitude), floor(pairs$longitude))
blocks <- unique(pairs$block)
n_folds <- 10
n_repeats <- 5

log_density <- function(rows, pred) {
  dnorm(pairs$y[rows], pred, sqrt(pairs$var_total[rows]), log = TRUE)
}

set.seed(2026)
pair_folds <- map(seq_len(n_repeats), function(rep) {
  block_fold <- sample(rep_len(seq_len(n_folds), length(blocks)))
  block_fold[match(pairs$block, blocks)]
})
cv_scores <- map_dfr(seq_len(n_repeats), function(rep) {
  pair_fold <- pair_folds[[rep]]
  imap_dfr(variants, function(variant, name) {
    map_dfr(seq_len(n_folds), function(fold) {
      imap_dfr(class_rows, function(rows, class) {
        train <- rows[pair_fold[rows] != fold]
        test <- rows[pair_fold[rows] == fold]
        fit <- fit_pairs(bases[[name]], train, pair_weights, variant$reversion)
        pred <- predict_pairs(bases[[name]], fit$par, test, variant$reversion)
        tibble(repeat_id = rep, variant = name, insecticide_class = class,
               pair_id = test, pred = pred,
               log_density = log_density(test, pred))
      })
    })
  })
})

# in-sample log density of the saved full-data model's implied rates, as a
# reference (not cross-validated, so favoured)
model_reference <- pairs %>%
  mutate(log_density = log_density(pair_id, model_rate)) %>%
  group_by(insecticide_class) %>%
  summarise(log_density = sum(log_density))
print(model_reference)

# summarise: held-out log density summed over pairs (mean over repeats),
# difference from the linear variant, and its standard error from the spread
# of per-block differences
cv_pair <- cv_scores %>%
  group_by(variant, insecticide_class, pair_id) %>%
  summarise(log_density = mean(log_density), .groups = "drop") %>%
  left_join(pairs %>% select(pair_id, block), by = "pair_id")

summarise_cv <- function(data) {
  base <- data %>%
    filter(variant == "linear") %>%
    select(pair_id, base = log_density)
  data %>%
    left_join(base, by = "pair_id") %>%
    group_by(variant, block) %>%
    summarise(diff = sum(log_density - base),
              log_density = sum(log_density),
              .groups = "drop_last") %>%
    summarise(log_density = sum(log_density),
              diff_vs_linear = sum(diff),
              se = sqrt(n()) * sd(diff),
              .groups = "drop")
}
cv_table <- bind_rows(
  cv_pair %>%
    summarise_cv() %>%
    mutate(insecticide_class = "all"),
  cv_pair %>%
    split(.$insecticide_class) %>%
    imap_dfr(~ summarise_cv(.x) %>% mutate(insecticide_class = .y))
) %>%
  mutate(variant = factor(variant, names(variants))) %>%
  arrange(insecticide_class, variant) %>%
  left_join(
    pairs %>%
      count(insecticide_class, name = "n_pairs") %>%
      bind_rows(tibble(insecticide_class = "all", n_pairs = nrow(pairs))),
    by = "insecticide_class"
  )
print(cv_table, n = Inf)
write_csv(cv_table, "outputs/selection_shape_cv.csv")


# robustness of the fitted net-use coefficient

# the saved model's net-use effect per class: posterior mean effect per type,
# averaged over the types in the class weighted by their numbers of pairs
model_net_class <- pairs %>%
  count(insecticide_class, insecticide_type) %>%
  mutate(effect = model_effect["nets", insecticide_type]) %>%
  group_by(insecticide_class) %>%
  summarise(model_net = weighted.mean(effect, n))
model_net <- setNames(model_net_class$model_net,
                      model_net_class$insecticide_class)
print(model_net)

# profile: in-sample and held-out log density with the net-use coefficient
# fixed at fractions of the model's value and the others refitted, for the
# current columns (raw pop) and for log pop
profile_bases <- c("linear_pop", "linear")
net_fractions <- c(0, 0.25, 0.5, 1)
net_profile_pairs <- expand_grid(variant = profile_bases,
                                 fraction = net_fractions) %>%
  pmap_dfr(function(variant, fraction) {
    imap_dfr(class_rows, function(rows, class) {
      fixed <- c(nets = fraction * model_net[[class]])
      x <- bases[[variant]]
      fit <- fit_pairs(x, rows, pair_weights, FALSE, fixed)
      in_sample <- log_density(rows, predict_pairs(x, fit$par, rows, FALSE))
      held_out <- map(pair_folds, function(pair_fold) {
        score <- numeric(length(rows))
        for (fold in seq_len(n_folds)) {
          test_idx <- which(pair_fold[rows] == fold)
          train <- rows[-test_idx]
          test <- rows[test_idx]
          fit <- fit_pairs(x, train, pair_weights, FALSE, fixed)
          score[test_idx] <- log_density(test,
                                         predict_pairs(x, fit$par, test, FALSE))
        }
        score
      })
      tibble(variant = variant, fraction = fraction,
             insecticide_class = class, net_coef = fixed[["nets"]],
             pair_id = rows, in_sample = in_sample,
             held_out = reduce(held_out, `+`) / length(held_out))
    })
  })

# differences from a zero net-use coefficient, with standard errors from the
# spread of per-block differences
summarise_profile <- function(data) {
  data %>%
    left_join(pairs %>% select(pair_id, block), by = "pair_id") %>%
    group_by(variant, pair_id) %>%
    mutate(in_sample_vs_0 = in_sample - in_sample[fraction == 0],
           held_out_vs_0 = held_out - held_out[fraction == 0]) %>%
    group_by(variant, fraction, block) %>%
    summarise(across(c(in_sample, held_out, in_sample_vs_0, held_out_vs_0),
                     sum),
              .groups = "drop_last") %>%
    summarise(se_held_out_vs_0 = sqrt(n()) * sd(held_out_vs_0),
              in_sample = sum(in_sample),
              held_out = sum(held_out),
              in_sample_vs_0 = sum(in_sample_vs_0),
              held_out_vs_0 = sum(held_out_vs_0),
              .groups = "drop")
}
net_profile <- bind_rows(
  net_profile_pairs %>%
    split(.$insecticide_class) %>%
    imap_dfr(~ summarise_profile(.x) %>%
               mutate(insecticide_class = .y,
                      net_coef = fraction * model_net[[.y]])),
  summarise_profile(net_profile_pairs) %>%
    mutate(insecticide_class = "all")
) %>%
  relocate(variant, insecticide_class, fraction, net_coef) %>%
  relocate(se_held_out_vs_0, .after = held_out_vs_0)
print(net_profile, n = Inf, width = Inf)
write_csv(net_profile, "outputs/selection_shape_net_profile.csv")

# block bootstrap of the net-use coefficient: resample 1-degree blocks with
# replacement, weighting each pair by its block's multiplicity
bootstrap_net <- function(rows, variant, n_boot = 200, seed = 1) {
  x <- bases[[variant]]
  estimate <- fit_pairs(x, rows, pair_weights, FALSE)$par[["nets"]]
  rows_by_block <- split(rows, pairs$block[rows])
  set.seed(seed)
  boots <- parallel::mclapply(seq_len(n_boot), function(i) {
    sampled <- sample(names(rows_by_block), replace = TRUE)
    mult <- table(unlist(rows_by_block[sampled]))
    boot_rows <- as.integer(names(mult))
    w <- pair_weights
    w[boot_rows] <- w[boot_rows] * as.vector(mult)
    fit_pairs(x, boot_rows, w, FALSE)$par[["nets"]]
  }, mc.cores = 6)
  boots <- unlist(boots)
  tibble(variant = variant,
         n_pairs = length(rows),
         n_blocks = length(rows_by_block),
         estimate = estimate,
         lower = quantile(boots, 0.025),
         upper = quantile(boots, 0.975),
         share_zero = mean(boots < 1e-8),
         model_net = NA_real_)
}

# per-pair net use at the two records, for the subset spanning an increase
pairs <- pairs %>%
  select(-any_of(c("nets_1", "nets_2"))) %>%
  left_join(covs %>% select(cell_id, year_1 = year_id, nets_1 = nets),
            by = c("cell_id", "year_1")) %>%
  left_join(covs %>% select(cell_id, year_2 = year_id, nets_2 = nets),
            by = c("cell_id", "year_2"))

pyrethroid <- pairs$insecticide_class == "Pyrethroids"
top_countries <- pairs %>%
  filter(insecticide_class == "Pyrethroids") %>%
  count(country_name, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(country_name)
subsets <- c(
  map(class_rows, identity) %>%
    set_names(paste(names(class_rows), "all pairs", sep = ": ")),
  list(
    "Pyrethroids: first record 2010 or later" =
      pairs$pair_id[pyrethroid & pairs$year_1 + baseline_year - 1 >= 2010],
    "Pyrethroids: net use rose by > 0.2 over the interval" =
      pairs$pair_id[pyrethroid & pairs$nets_2 - pairs$nets_1 > 0.2]
  ),
  map(top_countries, function(country) {
    pairs$pair_id[pyrethroid & pairs$country_name == country]
  }) %>%
    set_names(paste("Pyrethroids:", top_countries))
)
net_bootstrap <- expand_grid(subset = names(subsets),
                             variant = profile_bases) %>%
  pmap_dfr(function(subset, variant) {
    bootstrap_net(subsets[[subset]], variant) %>%
      mutate(subset = subset, .before = everything())
  }) %>%
  mutate(model_net = unname(.env$model_net[str_remove(subset, ":.*$")]))
print(net_bootstrap, n = Inf, width = Inf)
write_csv(net_bootstrap, "outputs/selection_shape_net_bootstrap.csv")


# how much net-use variation the pairs carry, and how collinear it is with
# the constant-like columns: interval-mean net use across pairs, its
# correlation with log pop and raw pop, and the share of its variance that a
# regression on a constant, log pop, IRS and the crops explains
net_confounding <- pairs %>%
  mutate(group = insecticide_class) %>%
  bind_rows(pairs %>% mutate(group = "all")) %>%
  group_by(group) %>%
  group_modify(function(d, key) {
    others <- lm(reformulate(c("log_pop", "irs", sprintf("`%s`", crop_names)),
                             "nets"),
                 data = d)
    tibble(
      n_pairs = nrow(d),
      nets_q05 = quantile(d$nets, 0.05),
      nets_q25 = quantile(d$nets, 0.25),
      nets_median = median(d$nets),
      nets_q75 = quantile(d$nets, 0.75),
      nets_q95 = quantile(d$nets, 0.95),
      nets_sd = sd(d$nets),
      share_below_0.18 = mean(d$nets < knots$nets[1]),
      log_pop_mean = mean(d$log_pop),
      log_pop_sd = sd(d$log_pop),
      cor_nets_log_pop = cor(d$nets, d$log_pop),
      cor_nets_pop = cor(d$nets, d$pop),
      r2_nets_on_others = summary(others)$r.squared
    )
  }) %>%
  ungroup()
print(net_confounding, width = Inf)
write_csv(net_confounding, "outputs/selection_shape_net_confounding.csv")
