# Fit the dynamical model to one cross-validation training fold and return
# posterior predictive draws for the held-out data.
#
# The three cross-validation experiments in dynamic_predictive_validation.R
# previously each carried their own copy of the model definition, so a change to
# the model had to be made in three places. This is that definition, once (#10).
#
# The model itself is unchanged from fit_model.R: the likelihood is simply
# restricted to the training fold. What differs from the previous validation
# code is what comes back. Rather than collapsing the posterior to a mean
# predicted fraction and a mean overdispersion, this returns the draws
# themselves, so that the held-out data can be scored against the full posterior
# predictive distribution.
#
# The predicted fraction and the overdispersion are drawn in a single
# calculate() call, so that each pair comes from the same posterior sample.
# Drawing them separately, as the earlier code did, breaks that pairing; it did
# not matter when only their means were used, but a predictive distribution
# needs them coupled.
#
# Because the greta arrays are local to this function they go out of scope when
# it returns, so the model does not need to be purged by hand between folds.

fit_fold <- function(train_df,
                     test_df,
                     x_cell_years,
                     df,
                     classes_index,
                     types,
                     n_covs,
                     n_times,
                     n_unique_cells,
                     n_classes,
                     n_types,
                     n_regions,
                     n_countries,
                     n_chains = 4,
                     warmup = 2000,
                     n_samples = 500,
                     n_sim = 1000,
                     Lmax = 30,
                     target_ess = 1000,
                     max_samples = 5000,
                     batch_samples = 500,
                     threads = 0,
                     inits_file = "temporary/inits.RDS") {

  # TensorFlow parallelises one op across cores rather than running a chain per
  # core, and scales poorly beyond about four threads, so confining each fold to
  # a couple of threads costs it little and leaves room to run several folds at
  # once. greta exposes no interface for this; it has to be set through
  # reticulate, before any ops are created
  if (threads > 0) {
    tensorflow_module <- reticulate::import("tensorflow")
    tensorflow_module$config$threading$set_intra_op_parallelism_threads(
      as.integer(threads))
    tensorflow_module$config$threading$set_inter_op_parallelism_threads(
      as.integer(threads))
  }


  # doubly hierarchical version

  # between-classes
  beta_overall <- normal(0, 1, dim = n_covs)
  sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
  sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))

  # between-types
  beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
  beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")

  # betas for types
  beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
  beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index] + beta_type_sigma

  # multiply through to get relative fitness of resistance for each insecticide
  # type at each cell

  effect_type <- exp(beta_type)


  selection_cell_years <- x_cell_years %*% effect_type
  fitness_cell_years <- 1 + selection_cell_years

  # reformat this in to a 3D array with dimensions:
  #   n_times x n_unique_cells x n_types x 1
  # to solve dynamics with time-varying fitness (time must be first, then other
  # two must match state variable, which has a trailing dimension of size 1)
  fitness_array <- fitness_cell_years
  dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)

  # prior and minimum values for the initial fractions susceptible
  init_frac_prior <- ifelse(types == "DDT", 0.9, 0.95)
  init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)

  # mean logit proportion of the distance from the minimum to 1, for each type
  init_frac_relative_prior <- (init_frac_prior - init_frac_min) / (1 - init_frac_min)
  logit_init_mean <- normal(qlogis(init_frac_relative_prior), 1, dim = n_types)

  # variance at each level, for each insecticide
  init_region_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
  init_country_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)

  # unscaled deviation parameters (hierarchical decentring)
  init_region_raw <- normal(0, 1, dim = c(n_regions, dim = n_types))
  init_country_raw <- normal(0, 1, dim = c(n_countries, dim = n_types))

  # combine to get the regional and country-level deviation from the prior logit
  # mean
  init_region_effect <- sweep(init_region_raw, 2, init_region_sd, FUN = "*")
  init_country_effect <- sweep(init_country_raw, 2, init_country_sd, FUN = "*")

  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)

  # combine these together into the country-level logit-mean initial value
  # logit_init_mean <- qlogis(init_frac_mean)
  init_country_overall_effect <- init_country_effect + init_region_effect[country_region_index, ]
  logit_init_country <- sweep(init_country_overall_effect,
                              2,
                              logit_init_mean,
                              FUN = "+")

  # convert from relative (0-1) to the constrained scale (above init_frac_min)
  init_country_relative <- ilogit(logit_init_country)
  init_range <- 1 - init_frac_min
  init_country_magnitude <- sweep(init_country_relative, 2, init_range, FUN = "*")
  init_country <- sweep(init_country_magnitude, 2, init_frac_min, FUN = "+")

  cell_country_lookup <- df %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cell_id, country_id) %>%
    arrange(cell_id) %>%
    pull(country_id)

  # expand out to all observations
  init_array <- init_country[cell_country_lookup, ]

  # add a trailing dimension to match greta.dynamics interface
  dim(init_array) <- c(dim(init_array), 1)

  # iterate through time to get fraction without resistant allele for all years
  # at all cells with data
  dynamic_cells <- iterate_dynamic_function(
    transition_function = haploid_next,
    initial_state = init_array,
    niter = n_times,
    w = fitness_array,
    parameter_is_time_varying = c("w"),
    tol = 0)

  # pull out the values at the cells, insecticides and years corresponding to
  # the training data
  index <- cbind(train_df$cell_id, train_df$type_id, train_df$year_id)
  fraction_susceptible_vec <- dynamic_cells$all_states[index]

  # resistant phenotype expression rate = expected mortality in bioassays
  population_mortality_vec <- fraction_susceptible_vec

  # define observation model
  rho_classes <- normal(0, 0.025,
                        truncation = c(0, 1),
                        dim = n_classes)

  distribution(train_df$died) <- betabinomial_p_rho(
    N = train_df$mosquito_number,
    p = population_mortality_vec,
    rho = rho_classes[train_df$class_id])

  m <- model(
    # initial fractions susceptible
    init_region_sd,
    init_country_sd,
    init_region_raw,
    init_country_raw,
    # hierarchical regression coefficients
    beta_overall,
    beta_class_raw,
    beta_type_raw,
    sigma_overall,
    sigma_class,
    # dispersion parameters
    rho_classes
  )

  # use cached posterior means as inits
  inits_one <- readRDS(inits_file)
  inits <- replicate(n_chains,
                     inits_one,
                     simplify = FALSE)

  Lmin <- round(Lmax / 2)

  draws <- mcmc(m,
                chains = n_chains,
                initial_values = inits,
                warmup = warmup,
                sampler = hmc(Lmin = Lmin, Lmax = Lmax),
                n_samples = n_samples)

  # Extend sampling until the effective sample size is adequate, rather than
  # guessing a sample count up front. How many draws are needed per effective
  # sample cannot be known before the sampler has adapted, so take an initial
  # batch, measure, and top up. Each batch reports its effective sample size, so
  # progress is visible in the fold's log rather than only on completion.
  report <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...), "\n")
    flush(stdout())
  }

  sampled <- n_samples
  ess <- coda::effectiveSize(draws)
  report("sampled %d per chain | ESS min %.0f median %.0f | target %d",
         sampled, min(ess, na.rm = TRUE), median(ess, na.rm = TRUE), target_ess)

  while (min(ess, na.rm = TRUE) < target_ess && sampled < max_samples) {
    draws <- extra_samples(draws,
                           n_samples = batch_samples,
                           verbose = TRUE)
    sampled <- sampled + batch_samples
    ess <- coda::effectiveSize(draws)
    report("sampled %d per chain | ESS min %.0f median %.0f | target %d",
           sampled, min(ess, na.rm = TRUE), median(ess, na.rm = TRUE), target_ess)
  }

  # Predictions at the held-out data.
  #
  # These come from greta's calculate() applied to the draws object, which is
  # the supported way to predict from a fitted greta model. Note there is no
  # nsim argument: with nsim, calculate() returns an independent resample of the
  # posterior, which is a valid posterior sample but destroys the MCMC ordering,
  # so effective sample size cannot be recovered from it. Without nsim it
  # returns the draws in order, as an mcmc.list, and the predicted fractions can
  # be diagnosed like any other monitored quantity.
  index_test <- cbind(test_df$cell_id, test_df$type_id, test_df$year_id)
  population_mortality_vec_test <- dynamic_cells$all_states[index_test]

  report("computing predictions at %d held-out assays", nrow(test_df))
  prediction_draws <- calculate(population_mortality_vec_test,
                                rho_classes,
                                values = draws)

  # effective sample size of the quantities the validation metrics actually
  # consume, rather than of the raw model parameters
  ess_prediction <- coda::effectiveSize(prediction_draws)
  ess_p <- ess_prediction[grep("population_mortality_vec_test",
                               names(ess_prediction))]
  ess_rho <- ess_prediction[grep("rho_classes", names(ess_prediction))]

  report("prediction ESS: p median %.0f min %.0f | rho median %.0f min %.0f",
         median(ess_p, na.rm = TRUE), min(ess_p, na.rm = TRUE),
         median(ess_rho, na.rm = TRUE), min(ess_rho, na.rm = TRUE))

  # flatten the mcmc.list to a draws x quantity matrix, preserving order
  prediction_matrix <- as.matrix(prediction_draws)
  p_columns <- grep("population_mortality_vec_test", colnames(prediction_matrix))
  rho_columns <- grep("rho_classes", colnames(prediction_matrix))
  p_draws <- prediction_matrix[, p_columns, drop = FALSE]
  rho_class_draws <- prediction_matrix[, rho_columns, drop = FALSE]
  # expand the class-level overdispersion out to one column per held-out assay
  rho_draws <- rho_class_draws[, test_df$class_id, drop = FALSE]

  convergence <- coda::gelman.diag(draws,
                                   multivariate = FALSE,
                                   autoburnin = FALSE)$psrf
  report("Rhat worst %.3f, %d of %d parameters above 1.01",
         max(convergence[, 1], na.rm = TRUE),
         sum(convergence[, 1] > 1.01, na.rm = TRUE),
         nrow(convergence))

  list(# the draws object itself: greta's calculate() needs this to predict, so
       # it is what must be kept
       draws = draws,
       p_draws = p_draws,
       rho_draws = rho_draws,
       test_df = test_df,
       n_train = nrow(train_df),
       convergence = convergence,
       ess = ess,
       ess_p = ess_p,
       ess_rho = ess_rho,
       n_sampled = sampled,
       n_chains = n_chains,
       draws_per_ess = (sampled * n_chains) / median(ess, na.rm = TRUE))

}


# dynamic iteration function used by the model above
haploid_next <- function(state, iter, w) {
  # fraction susceptible
  q <- state
  # fraction resistant
  p <- 1 - q
  # number resistant grows relative to resistant by ratio w, re-normalise to get
  # fraction remaining susceptible
  q / (q + p * w)
}
