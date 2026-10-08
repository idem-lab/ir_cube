# Recompute the dynamical model's predictions, as logit draws, at arbitrary
# cell-years from its saved posterior draws, in plain R.
#
# A saved draws object cannot be resumed to ask greta for more predictions, so
# the two-stage model (paired draw for draw with a fold's stored test
# predictions), R/predict.R and the figure scripts recompute them from the
# sampled parameters, with the model's own transforms (dynamical_terms() in
# R/dynamical_model.R) and its options (fold$options). The recursion is
# solved in closed form on the logit scale, as in closed_form_states():
#
#   logit q_t = logit q_0 - sum_{s <= t} log w_s - t kappa
#
# the state for year index t having had the fitness of years 1..t applied.
# The logit draws are of predicted mortality: q_t, or with a mortality floor
# f, f + (1 - f) q_t.

source("R/dynamical_model.R")
# thin_draws() and max_draws
source("R/validation_scoring.R")

# Draw indices into as.matrix(fold$draws) that the stored predictions
# correspond to: newer folds thin p_draws to `maximum` at fitting time with
# this rule, older folds are thinned by the same rule at scoring time. The
# rule is thin_draws() (R/validation_scoring.R). Chains recorded as stuck
# (fold$stuck_chains, set by R/drop_stuck_chains.R, which drops their rows
# from the stored predictions) are left out.
paired_draw_index <- function(fold, maximum = max_draws) {
  n_total <- sum(vapply(fold$draws, nrow, integer(1)))
  index <- thin_draws(matrix(seq_len(n_total)), maximum)[, 1]
  if (length(fold$stuck_chains) > 0) {
    index <- index[!draw_chain(fold$draws)[index] %in% fold$stuck_chains]
  }
  index
}

# The chain of each row of as.matrix(draws).
draw_chain <- function(draws) {
  rep(seq_along(draws), vapply(draws, nrow, integer(1)))
}

# The chains of `draws` that stopped moving: those whose share of distinct
# draws is below `min_distinct`. With one step size for all chains
# (windowed_hmc()), a chain can sit where every trajectory at that step size
# fails, and then repeats one draw: in the October 2026 run, one chain of
# spatial blocks fold 1 kept a single draw for all 3,000 samples.
stuck_chains <- function(draws, min_distinct = 0.5) {
  which(vapply(draws, function(chain) {
    chain <- as.matrix(chain)
    nrow(unique(chain)) / nrow(chain) < min_distinct
  }, logical(1)))
}

# `draws` (a greta_mcmc_list) with only the chains `keep`, in the raw
# free-state draws that calculate() reads as well.
drop_chains <- function(draws, keep) {
  out <- draws[keep]
  model_info <- attr(draws, "model_info")
  model_info$raw_draws <- model_info$raw_draws[keep]
  attr(out, "model_info") <- model_info
  class(out) <- class(draws)
  out
}

# One named parameter of a draws x parameter matrix, as a draws x
# dim(parameter) array. greta names elements by their R (column-major) index,
# e.g. "beta_type_raw[13,9]", so the index is parsed from the names
extract_parameter <- function(draws_matrix, name) {
  # a scalar has one column, named without an index
  if (name %in% colnames(draws_matrix)) {
    return(draws_matrix[, name, drop = FALSE])
  }
  columns <- grep(paste0("^", name, "\\["), colnames(draws_matrix))
  stopifnot(length(columns) > 0)
  index <- colnames(draws_matrix)[columns] %>%
    str_remove(paste0("^", name, "\\[")) %>%
    str_remove("\\]$") %>%
    str_split(",", simplify = TRUE)
  index <- matrix(as.integer(index), nrow = length(columns))
  dims <- apply(index, 2, max)
  # vectors are stored as column vectors, [i,1]
  if (length(dims) == 2 && dims[2] == 1) dims <- dims[1]
  linear <- if (length(dims) == 1) index[, 1] else
    index[, 1] + (index[, 2] - 1) * dims[1]
  out <- array(NA_real_, dim = c(nrow(draws_matrix), prod(dims)))
  out[, linear] <- draws_matrix[, columns]
  dim(out) <- c(nrow(draws_matrix), dims)
  out
}

# The columns of each variable's free state in the raw draws of greta model
# `model`, which concatenate them in the dag's variable order, named by the
# variables' TensorFlow names; attribute "targets" maps the model's target
# names to those names.
free_state_columns <- function(model) {
  dag <- model$dag
  sizes <- vapply(dag$example_parameters(free = TRUE), length, integer(1))
  columns <- Map(function(end, size) (end - size + 1):end, cumsum(sizes),
                 sizes)
  node_names <- vapply(dag$node_list, function(node) node$unique_name,
                       character(1))
  targets <- vapply(model$target_greta_arrays,
                    function(x) greta:::get_node(x)$unique_name,
                    character(1))
  attr(columns, "targets") <- setNames(
    dag$get_tf_names()[match(targets, node_names)], names(targets))
  columns
}

# dynamical_terms() for every draw of the variables `v` (a named list of
# draws x dim arrays), stacked as draws x dim(term) arrays, for the terms
# `terms`.
dynamical_terms_draws <- function(v, classes_index, types, terms, options) {
  n_draws <- nrow(v[[1]])
  out <- NULL
  for (i in seq_len(n_draws)) {
    v_i <- variables_at_draw(v, i)
    terms_i <- dynamical_terms(v_i, classes_index = classes_index,
                               types = types, options = options)[terms]
    if (is.null(out)) {
      out <- lapply(terms_i, function(x) {
        array(NA_real_, c(n_draws, dim(as.matrix(x))))
      })
    }
    for (name in terms) {
      out[[name]][i, , ] <- terms_i[[name]]
    }
  }
  out
}

# The parameters the predictions need, for the draws in `draw_index`:
#   effect_type          draws x n_covs x n_types, exp(beta_type)
#   logit_init_relative  draws x n_countries x n_types, the logit relative
#                        initial state (above init_frac_min) of each fitted
#                        country, to which a cell's initial-state covariate
#                        effects are added (dynamical_logit_cells())
#   init_coef            draws x n_init_covs x n_types, the coefficients of
#                        the initial-state covariates (NULL for none)
#   rho_types            draws x n_types, the observation overdispersion
#   mortality_floor      draws (NULL for none)
#   kappa_type           draws x n_types, the reversion kappa (<= 0; NULL for
#                        none, see reversion_kappa())
#   init_min             n_types, init_frac_min
#   x_cells_init         the fit's initial-state covariates, one row per
#                        cell_id (NULL for none)
#   variables            every variable, as draws x dim arrays, from which
#                        map_logit_init() draws countries without data
# and the options, types and classes_index.
dynamical_parameter_draws <- function(fold,
                                      classes_index,
                                      types,
                                      draw_index = paired_draw_index(fold),
                                      options = fold$options) {

  # as.matrix() on an mcmc.list stacks chains in order, which is how fit_fold()
  # flattened the calculate() output that p_draws came from
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  n_draws <- nrow(draws_matrix)
  n_types <- length(types)

  variable_names <- unique(sub("\\[.*$", "", colnames(draws_matrix)))
  v <- lapply(setNames(nm = variable_names), extract_parameter,
              draws_matrix = draws_matrix)
  stopifnot(!is.null(options),
            identical(dim(v$logit_init_mean), c(n_draws, n_types)))

  reversion <- !isFALSE(options$reversion)
  terms <- dynamical_terms_draws(
    v, classes_index, types,
    terms = c("beta_type", "logit_init_relative", "rho_types",
              if (reversion) "kappa_type"),
    options = options)

  list(effect_type = exp(terms$beta_type),
       logit_init_relative = terms$logit_init_relative,
       init_coef = if (!is.null(options$init_covariates)) {
         array(v$init_coef,
               c(n_draws, length(options$init_covariates), n_types),
               dimnames = list(NULL, options$init_covariates, types))
       },
       rho_types = matrix(terms$rho_types, n_draws),
       mortality_floor = if (isTRUE(options$mortality_floor)) {
         c(v$mortality_floor)
       },
       kappa_type = if (reversion) matrix(terms$kappa_type, n_draws),
       init_min = init_frac_constants(types)$min,
       x_cells_init = select_init_covariates(fold$x_cells_init, options),
       variables = v,
       options = options,
       types = types,
       classes_index = classes_index,
       n_draws = n_draws)
}

# `parameters` (dynamical_parameter_draws()) for its draws `draws` only
subset_draws <- function(parameters, draws) {
  rows <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) return(x[draws])
    do.call(`[`, c(list(x, draws), rep(list(TRUE), length(dim(x)) - 1),
                   drop = FALSE))
  }
  for (name in c("effect_type", "logit_init_relative", "init_coef",
                 "rho_types", "mortality_floor", "kappa_type")) {
    parameters[name] <- list(rows(parameters[[name]]))
  }
  parameters$variables <- lapply(parameters$variables, rows)
  parameters$n_draws <- length(draws)
  parameters
}

# Logit q_0 at cells, draws x cells, for type k, from the logit relative
# initial state at each cell's country (`logit_init`, draws x cells) and the
# cells' initial-state covariates x_init (cells x covariates, named columns;
# for fits with them), as logit_init_relative_rows() and closed_form_states()
# in the model
cell_logit_init <- function(parameters, k, logit_init, x_init = NULL) {
  init_coef <- parameters$init_coef
  if (!is.null(init_coef)) {
    stopifnot(!is.null(x_init), nrow(x_init) == ncol(logit_init))
    logit_init <- logit_init +
      matrix(init_coef[, , k], nrow = parameters$n_draws) %*%
      t(x_init[, dimnames(init_coef)[[2]], drop = FALSE])
  }
  floored_logit(logit_init, parameters$init_min[k])
}

# The recursion for cells of insecticide type k:
#   parameters  dynamical_parameter_draws() (or subset_draws())
#   logit_init  draws x cells, the logit relative initial state at each cell's
#               country (parameters$logit_init_relative, or map_logit_init())
#   x           cells x years x n_covs covariates, year index 1 the baseline
#               year, covariates in the column order of x_cell_years
#   years_keep  year indices to return
#   x_init      cells x initial-state covariates (named columns), for fits
#               with them
# Returns a list named by years_keep of draws x cells logit mortality.
dynamical_logit_cells <- function(parameters, k, logit_init, x, years_keep,
                                  x_init = NULL) {
  n_cells <- dim(x)[1]
  n_draws <- parameters$n_draws
  stopifnot(ncol(logit_init) == n_cells, nrow(logit_init) == n_draws,
            max(years_keep) <= dim(x)[2])

  logit_init <- cell_logit_init(parameters, k, logit_init, x_init)

  effect <- matrix(parameters$effect_type[, , k], nrow = n_draws)
  kappa <- parameters$kappa_type[, k]
  cumulative <- 0
  out <- list()
  for (t in seq_len(max(years_keep))) {
    # log1p because the selection term can be tiny
    cumulative <- cumulative +
      log1p(effect %*% t(matrix(x[, t, ], nrow = n_cells)))
    # reversion: - t kappa in year t (reversion_kappa()), one per draw
    if (!is.null(kappa)) cumulative <- cumulative + kappa
    if (t %in% years_keep) {
      out[[as.character(t)]] <- floored_logit(
        logit_init - cumulative, parameters$mortality_floor)
    }
  }
  out
}

# Draws of the dynamical model's logit predicted mortality at (cell, type,
# year) rows (cell_id, type_id, year_id, as in `df`), with the covariates of
# x_cell_years (one row per (cell_id, year_id) of cell_years_index). The
# initial condition of a cell is that of the country of its first record in
# the full `df`, as in the model (dynamical_lookups()), whatever the rows'
# country_id. Returns a draws x nrow(rows) matrix, paired with
# thin_draws(fold$p_draws) when `parameters` are at paired_draw_index().
dynamical_logit <- function(parameters, rows, df, x_cell_years,
                            cell_years_index, max_block = 2.5e7) {

  n_draws <- parameters$n_draws
  n_times <- max(cell_years_index$year_id)
  if (any(rows$year_id > n_times | rows$year_id < 1)) {
    stop("rows fall outside the covariate years 1..", n_times)
  }
  cell_country <- dynamical_lookups(df)$cell_country_lookup
  stopifnot(!anyNA(cell_country[rows$cell_id]))

  # row of x_cell_years for each (cell, year)
  x_row <- matrix(NA_integer_, max(cell_years_index$cell_id), n_times)
  x_row[cbind(cell_years_index$cell_id, cell_years_index$year_id)] <-
    seq_len(nrow(cell_years_index))

  # assays sharing a (cell, type, year) share a prediction, computed once
  keys <- paste(rows$cell_id, rows$type_id, rows$year_id)
  unique_rows <- rows[!duplicated(keys), c("cell_id", "type_id", "year_id")]
  result <- matrix(NA_real_, n_draws, nrow(unique_rows))
  chunk_size <- max(1, floor(max_block / (n_times * n_draws)))

  for (k in sort(unique(unique_rows$type_id))) {
    cells_k <- sort(unique(unique_rows$cell_id[unique_rows$type_id == k]))
    for (cells in split(cells_k, ceiling(seq_along(cells_k) / chunk_size))) {
      target <- which(unique_rows$type_id == k &
                        unique_rows$cell_id %in% cells)
      years_keep <- sort(unique(unique_rows$year_id[target]))
      x_index <- x_row[cells, seq_len(max(years_keep)), drop = FALSE]
      stopifnot(!anyNA(x_index))
      x <- array(x_cell_years[as.vector(x_index), , drop = FALSE],
                 c(length(cells), max(years_keep), ncol(x_cell_years)))
      logit <- dynamical_logit_cells(
        parameters, k,
        matrix(parameters$logit_init_relative[, cell_country[cells], k],
               nrow = n_draws),
        x, years_keep,
        x_init = parameters$x_cells_init[cells, , drop = FALSE])
      for (t in years_keep) {
        at_t <- target[unique_rows$year_id[target] == t]
        result[, at_t] <- logit[[as.character(t)]][
          , match(unique_rows$cell_id[at_t], cells), drop = FALSE]
      }
    }
  }

  result[, match(keys, keys[!duplicated(keys)]), drop = FALSE]
}
