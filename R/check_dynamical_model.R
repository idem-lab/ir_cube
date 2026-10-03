# Check that the greta model (build_dynamical_model(), with the closed-form op)
# and the plain-R predictions (R/dynamical_predictions.R) agree, at a random
# free state, for a given set of model options:
#   - predicted mortality at every training assay, and rho per type
#   - the log likelihood: the model's log density over all the data, less that
#     of the same model with one assay in the likelihood (which shares every
#     prior and Jacobian term), against the plain-R betabinomial log likelihood
#     of the other assays
# and print the log density itself, for regression checks between versions.
#
#   Rscript R/check_dynamical_model.R '<options>' [seed] [sd]
# (the free state is N(0, sd^2), sd 0.5 by default; a smaller sd avoids states
# where p rounds to 1 at assays with survivors, and the log density is NaN)
# e.g.
#   Rscript R/check_dynamical_model.R 'dynamical_model_options(reversion = FALSE)'
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
options_text <- if (length(arguments) >= 1) arguments[1] else
  "dynamical_model_options()"
seed <- if (length(arguments) >= 2) as.integer(arguments[2]) else 1L
free_sd <- if (length(arguments) >= 3) as.numeric(arguments[3]) else 0.5

source("R/greta_setup.R")
start_greta(threads = 4)
source("R/dynamical_model.R")
model_options <- eval(parse(text = options_text))
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/validation_functions.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

cat("options:", options_text, "\n")

build <- function(train_df) {
  build_dynamical_model(train_df = train_df,
                        df = df,
                        x_cell_years = x_cell_years,
                        cell_years_index = cell_years_index,
                        classes_index = classes_index,
                        types = types,
                        options = model_options,
                        x_cells_init = x_cells_init)
}

built <- build(df)
built_one <- build(df[1, ])
# the options as built, with the centre of the initial-state covariates
model_options <- built$options

log_density <- function(model, free) {
  f <- model$dag$generate_log_prob_function(which = "adjusted")
  free_tf <- tensorflow::tf$constant(matrix(free, nrow = 1),
                                     dtype = tensorflow::tf$float64)
  as.numeric(f(free_tf))
}

n_free <- length(unlist(built$model$dag$example_parameters(free = TRUE)))
stopifnot(n_free == length(unlist(
  built_one$model$dag$example_parameters(free = TRUE))))
set.seed(seed)
free <- rnorm(n_free, 0, free_sd)

# the free-state elements of a variable
free_columns <- function(model, name) {
  columns <- free_state_columns(model)
  columns[[attr(columns, "targets")[[name]]]]
}

# A random reversion rate is on the scale of its free state, around 1 per year,
# which drives p to 1 at most assays by the end of the series; so it is set to
# 0.05 per year (its free state is the log of the rate) and checked
if (identical(model_options$reversion, "estimated")) {
  columns <- free_columns(built$model, "reversion_rate")
  free[columns] <- log(0.05)
  trace <- built$model$dag$trace_values(matrix(free, nrow = 1))
  stopifnot(isTRUE(all.equal(
    unname(trace[1, grep("^reversion_rate", colnames(trace))]),
    rep(0.05, length(columns)))))
}

ld_all <- log_density(built$model, free)
ld_one <- log_density(built_one$model, free)

# the variables' values at this free state, as a one-draw fold for the plain-R
# path
trace <- built$model$dag$trace_values(matrix(free, nrow = 1))
stopifnot(identical(unname(trace),
                    unname(built_one$model$dag$trace_values(
                      matrix(free, nrow = 1)))))
fold <- list(draws = coda::mcmc.list(coda::mcmc(trace)),
             options = model_options, x_cells_init = x_cells_init)
parameters <- dynamical_parameter_draws(fold, classes_index, types, df,
                                        draw_index = 1)

# greta, at these values: calculate() on a one-draw greta_mcmc_list, as on a
# fitted model's draws
values <- greta:::as_greta_mcmc_list(
  coda::mcmc.list(coda::mcmc(trace)),
  list(raw_draws = coda::mcmc.list(coda::mcmc(matrix(free, nrow = 1))),
       model = built$model))
greta_values <- calculate(p = built$population_mortality_vec,
                          rho = built$terms$rho_types,
                          values = values)
greta_values <- as.matrix(greta_values)
p_greta <- greta_values[1, grep("^p\\[", colnames(greta_values))]
rho_greta <- greta_values[1, grep("^rho\\[", colnames(greta_values))]

# plain R
p_r <- plogis(c(dynamical_logit(parameters, df, df, x_cell_years,
                                cell_years_index)))
rho_r <- c(parameters$rho_types)

# the betabinomial log likelihood, parameterised as in betabinomial_p_rho() and
# not clamped (dbetabinom() clamps p away from 0 and 1)
a <- p_r * (1 / rho_r[df$type_id] - 1)
b <- a * (1 - p_r) / p_r
loglik_r <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                beta = b, log = TRUE)

cat(sprintf("free parameters %d, log density %.10g\n", n_free, ld_all))
cat(sprintf("p: max abs diff %.3g, max logit diff %.3g (range %.3g-%.3g)\n",
            max(abs(p_greta - p_r)),
            max(abs(qlogis(pmin(pmax(p_greta, 1e-12), 1 - 1e-12)) -
                      qlogis(pmin(pmax(p_r, 1e-12), 1 - 1e-12)))),
            min(p_r), max(p_r)))
cat(sprintf("rho: max abs diff %.3g\n", max(abs(rho_greta - rho_r))))
cat(sprintf("log likelihood of assays 2..n: greta %.10g, plain R %.10g, diff %.3g\n",
            ld_all - ld_one, sum(loglik_r[-1]),
            ld_all - ld_one - sum(loglik_r[-1])))

# the map path (R/two_stage_map_functions.R), at the data cells in 2000, 2012
# and 2024, against dynamical_logit()
map_years <- c(2000, 2012, 2024)
map_rows <- df %>%
  distinct(cell, cell_id) %>%
  mutate(country_name = countries[
    built$lookups$cell_country_lookup[cell_id]])
covariates <- map_covariates(map_rows$cell, baseline_year, max(map_years),
                             model_options$selection_columns)

# the map path's covariates are x_cell_years at the data cells
n_fit_years <- max(cell_years_index$year_id)
x_map <- matrix(aperm(map_x(covariates, seq_len(nrow(map_rows)), n_fit_years),
                      c(2, 1, 3)),
                ncol = ncol(x_cell_years))
x_fit <- x_cell_years[match(paste(rep(map_rows$cell_id, each = n_fit_years),
                                  rep(seq_len(n_fit_years), nrow(map_rows))),
                            paste(cell_years_index$cell_id,
                                  cell_years_index$year_id)), ]
cat(sprintf("map covariates vs x_cell_years, %d cell-years x %d columns: max abs diff %.3g\n",
            nrow(x_fit), ncol(x_fit), max(abs(x_map - x_fit))))
stopifnot(identical(dim(x_map), dim(x_fit)), max(abs(x_map - x_fit)) == 0)

# the trends' product columns are 0 in the baseline year, unless a trend is
# given as a region x year matrix (e.g. g(1995) = 0.37)
trend_columns <- grep(":g_(dom|ag)$", colnames(x_cell_years))
if (length(trend_columns) > 0) {
  x_baseline <- x_cell_years[cell_years_index$year_id == 1, trend_columns]
  cat(sprintf("%d trend product columns in %d: %d non-zero values\n",
              length(trend_columns), baseline_year, sum(x_baseline != 0)))
  if (any(x_baseline != 0)) {
    design <- model_options$selection_columns
    if (!is.matrix(design$trend_pop) && !is.matrix(design$trend_crops)) {
      stop("trend product columns are not 0 in the baseline year")
    }
    warning("trend product columns are not 0 in the baseline year")
  }
}

logit_init_all <- map_logit_init(parameters, countries, regions, df)
cell_country_index <- match(map_rows$country_name,
                            dimnames(logit_init_all)[[2]])
x_years <- map_x(covariates, seq_len(nrow(map_rows)),
                 max(map_years) - baseline_year + 1)
clamp <- function(l) pmin(pmax(l, qlogis(1e-12)), qlogis(1 - 1e-12))
map_difference <- 0
for (k in seq_along(types)) {
  dyn <- dynamical_logit_cells(parameters, k,
                               matrix(logit_init_all[, cell_country_index, k],
                                      1),
                               x_years, map_years - baseline_year + 1,
                               x_init = covariates$init)
  for (y in map_years) {
    rows <- tibble(cell_id = map_rows$cell_id, type_id = k,
                   year_id = y - baseline_year + 1)
    l_rows <- c(dynamical_logit(parameters, rows, df, x_cell_years,
                                cell_years_index))
    l_map <- c(dyn[[as.character(y - baseline_year + 1)]])
    map_difference <- max(map_difference, abs(clamp(l_map) - clamp(l_rows)))
  }
}
cat(sprintf("map path vs plain R, logit, %d cells x %d types x %d years: max abs diff %.3g\n",
            nrow(map_rows), length(types), length(map_years), map_difference))
