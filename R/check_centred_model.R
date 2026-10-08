# Check the centred parameterisation (centred_options(), R/dynamical_model.R)
# against the non-centred model. Centring changes the coordinates HMC samples
# in, not the model, so at corresponding points (centre_variables(),
# noncentre_variables()):
#   - the log densities differ by the log Jacobian of the change of
#     variables: for each centred effect b = m + s z, sampled as b rather
#     than z, the centred log density is lower by log s;
#   - the predicted mortality at every assay and rho per type are the same;
#   - mapping a point to the other coordinates and back returns it.
# At three points: a random free state of each model, and the cached initial
# values (dynamical_inits()). Run R/check_dynamical_model.R with the centred
# options too, for the greta and plain-R predictions.
#
#   IR_CUBE_MODEL_OPTIONS='<options>' Rscript R/check_centred_model.R \
#     [seed] [sd] [threads]
# The centring is that of IR_CUBE_MODEL_OPTIONS, or if it has none,
# centred_options_data_informed(); the non-centred model has the same options
# without it. The free state is N(0, sd^2), sd 0.5 by default.
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
seed <- if (length(arguments) >= 1) as.integer(arguments[1]) else 1L
free_sd <- if (length(arguments) >= 2) as.numeric(arguments[2]) else 0.5
threads <- if (length(arguments) >= 3) as.integer(arguments[3]) else 4L

source("R/greta_setup.R")
start_greta(threads = threads)
source("R/dynamical_model.R")
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/dynamical_predictions.R")

centred_model_options <- model_options
if (identical(model_options$centred, centred_options())) {
  centred_model_options$centred <- centred_options_data_informed()
}
noncentred_model_options <- model_options
noncentred_model_options$centred <- centred_options()
cat("centred:", deparse(centred_model_options$centred), "\n")

build <- function(options) {
  build_dynamical_model(train_df = df,
                        df = df,
                        x_cell_years = x_cell_years,
                        cell_years_index = cell_years_index,
                        classes_index = classes_index,
                        types = types,
                        options = options,
                        x_cells_init = x_cells_init)
}
noncentred <- build(noncentred_model_options)
centred <- build(centred_model_options)
# the options as built, with the centre of the initial-state covariates
centred_model_options <- centred$options

log_density <- function(model, free) {
  f <- model$dag$generate_log_prob_function(which = "adjusted")
  free_tf <- tensorflow::tf$constant(matrix(free, nrow = 1),
                                     dtype = tensorflow::tf$float64)
  as.numeric(f(free_tf))
}

# the values of a model's variables at free state `free`, as a named list of
# arrays
variable_values <- function(model, free) {
  trace <- model$dag$trace_values(matrix(free, nrow = 1))
  names <- unique(sub("\\[.*$", "", colnames(trace)))
  v <- lapply(setNames(nm = names), extract_parameter, draws_matrix = trace)
  variables_at_draw(v, 1)
}

# the free state of a model at the values `values` of its variables, by
# greta's own transforms; greta's free state is row-major within a variable
free_state <- function(model, values) {
  columns <- free_state_columns(model)
  targets <- attr(columns, "targets")
  free <- rep(NA_real_, sum(lengths(columns)))
  for (name in names(targets)) {
    node <- greta:::get_node(model$target_greta_arrays[[name]])
    value <- array(values[[name]], node$dim)
    free[columns[[targets[[name]]]]] <-
      greta:::flatten_rowwise(greta:::to_free(node, value))
  }
  stopifnot(!anyNA(free))
  free
}

# the log Jacobian of the move from the non-centred to the centred
# coordinates: log s for each centred effect, s its prior sd
log_jacobian <- function(values) {
  rows <- which(centred_selection_rows(centred_model_options))
  n_classes <- ncol(values$beta_class_raw)
  n_types <- ncol(values$beta_type_raw)
  out <- n_classes * sum(log(values$sigma_overall[rows])) +
    n_types * sum(log(values$sigma_class[rows]))
  if (centred_rho_type(centred_model_options)) {
    out <- out + n_types * log(c(values$rho_sigma_type))
  }
  out
}

# predicted mortality at every assay and rho per type, by calculate()
predictions <- function(built, free) {
  trace <- built$model$dag$trace_values(matrix(free, nrow = 1))
  values <- greta:::as_greta_mcmc_list(
    coda::mcmc.list(coda::mcmc(trace)),
    list(raw_draws = coda::mcmc.list(coda::mcmc(matrix(free, nrow = 1))),
         model = built$model))
  out <- calculate(p = built$population_mortality_vec,
                   rho = built$terms$rho_types, values = values)
  c(as.matrix(out))
}

# compare the models at the non-centred free state `free_noncentred`
compare <- function(label, free_noncentred) {
  v_noncentred <- variable_values(noncentred$model, free_noncentred)
  v_centred <- centre_variables(v_noncentred, classes_index,
                                centred_model_options)
  free_centred <- free_state(centred$model, v_centred)
  back <- noncentre_variables(variable_values(centred$model, free_centred),
                              classes_index, centred_model_options)
  free_back <- free_state(noncentred$model, back)
  ld_noncentred <- log_density(noncentred$model, free_noncentred)
  ld_centred <- log_density(centred$model, free_centred)
  jacobian <- log_jacobian(v_noncentred)
  ld_diff <- ld_centred - (ld_noncentred - jacobian)
  round_trip <- max(abs(free_back - free_noncentred))
  prediction_diff <- max(abs(predictions(centred, free_centred) -
                               predictions(noncentred, free_noncentred)))
  cat(sprintf(paste0("%s: log density non-centred %.10g, centred %.10g, ",
                     "log Jacobian %.6g, difference %.3g; round trip %.3g; ",
                     "predictions %.3g\n"),
              label, ld_noncentred, ld_centred, jacobian, ld_diff, round_trip,
              prediction_diff))
  stopifnot(is.finite(ld_noncentred), abs(ld_diff) < 1e-8,
            round_trip < 1e-10, prediction_diff < 1e-10)
}

n_free <- length(unlist(noncentred$model$dag$example_parameters(free = TRUE)))
stopifnot(n_free ==
            length(unlist(centred$model$dag$example_parameters(free = TRUE))))
cat("free parameters", n_free, "\n")

# A random reversion rate is on the scale of its free state, around 1 per year,
# which drives p to 1 at most assays by the end of the series; so it is set to
# 0.05 per year, as in R/check_dynamical_model.R
set_reversion <- function(model, free) {
  if (identical(centred_model_options$reversion, "estimated")) {
    columns <- free_state_columns(model)
    name <- attr(columns, "targets")[["reversion_rate"]]
    free[columns[[name]]] <- log(0.05)
  }
  free
}

# 1. a random free state of the non-centred model
set.seed(seed)
free <- set_reversion(noncentred$model, rnorm(n_free, 0, free_sd))
compare("random non-centred state", free)

# 2. a random free state of the centred model, moved to the non-centred one
free_centred <- set_reversion(centred$model, rnorm(n_free, 0, free_sd))
v_centred <- variable_values(centred$model, free_centred)
free <- free_state(noncentred$model,
                   noncentre_variables(v_centred, classes_index,
                                       centred_model_options))
compare("random centred state", free)

# 3. the cached initial values, given to each model as dynamical_inits() gives
# them
inits <- function(built) {
  unclass(dynamical_chain_inits(dynamical_inits_files()[1], built$variables,
                                levels = built$lookups$levels,
                                columns = colnames(x_cell_years),
                                n_chains = 1, options = built$options,
                                classes_index = classes_index)[[1]])
}
inits_noncentred <- inits(noncentred)
inits_centred <- inits(centred)
free <- free_state(noncentred$model, inits_noncentred)
stopifnot(max(abs(free_state(centred$model, inits_centred) -
                    free_state(centred$model,
                               centre_variables(variable_values(
                                 noncentred$model, free),
                                 classes_index, centred_model_options)))) <
            1e-10)
compare("cached initial values", free)

# the draws of the centred model, moved to the non-centred variables as
# fit_model.R caches them (noncentred_draws()), against the same draws of the
# non-centred model
draws_centred <- lapply(setNames(nm = names(centred$variables)), function(name) {
  array(inits_centred[[name]], c(1, dim(centred$variables[[name]])))
})
draws_moved <- noncentred_draws(draws_centred, classes_index,
                                centred_model_options)
moved_diff <- max(vapply(names(inits_noncentred), function(name) {
  max(abs(c(draws_moved[[name]]) - c(inits_noncentred[[name]])))
}, numeric(1)))
cat(sprintf("noncentred_draws() of the centred inits: max abs diff %.3g\n",
            moved_diff))
stopifnot(setequal(names(draws_moved), names(inits_noncentred)),
          moved_diff < 1e-12)
cat("all checks passed\n")
