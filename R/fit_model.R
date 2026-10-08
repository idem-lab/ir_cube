# fit model

#   Rscript R/fit_model.R [threads]

# load packages and functions
# greta first, so python starts before terra and sf are attached; TensorFlow's
# thread count has to be set before then too
arguments <- commandArgs(trailingOnly = TRUE)
threads <- if (length(arguments) >= 1) as.integer(arguments[1]) else NULL
source("R/greta_setup.R")
start_greta(threads = threads)
source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/dynamical_model.R")
source("R/chain_floor_modes.R")

# set the start of the timeseries considered in modelling (the start of
# non-negligible levels of resistance) - assume it's before the mass-rollout of
# nets
baseline_year <- 1995

# set the final year of data (insufficient and spatially biased data for 2025)
final_data_year <- 2024

# the modelled subset: the insecticide types with at least 1000 unique
# places/times, and alpha-cypermethrin (914 unique) because of its use in LLINs,
# each at its modal concentration, in the modelled years and inside the mask,
# with indices to the classes, types, regions, countries and cells, and the
# class of each type (R/bioassay_subset.R)
list2env(modelled_bioassays(baseline_year, final_data_year), environment())
years <- baseline_year - 1 + sort(unique(df$year_id))

# the model terms (R/dynamical_model.R), or an R expression for them in
# IR_CUBE_MODEL_OPTIONS (docker/run_pod_job.sh)
model_options <- eval(str2lang(Sys.getenv("IR_CUBE_MODEL_OPTIONS",
                                          "dynamical_model_options()")))

# create design matrix at all unique cells and for all years, as the model
# options' selection design asks (R/model_covariates.R)
selection <- selection_design_matrix(unique_cells, baseline_year,
                                     final_data_year,
                                     model_options$selection_columns)
cell_years_index <- selection$cell_years_index
x_cell_years <- selection$x_cell_years
rm(selection)

# the initial-state covariates (#19) at each cell, one row per cell_id
x_cells_init <- init_covariate_matrix(unique_cells,
                                      model_options$selection_columns)

# dimensions of things in the fitting stage
n_covs <- ncol(x_cell_years)
n_obs <- nrow(df)
n_unique_cells <- length(unique_cells)
n_times <- max(df$year_start) - min(df$year_start) + 1
n_classes <- length(classes)
n_types <- length(types)
n_regions <- length(regions)
n_countries <- length(countries)

# build the model, with the likelihood over all the data (R/dynamical_model.R)
built <- build_dynamical_model(train_df = df,
                               df = df,
                               x_cell_years = x_cell_years,
                               cell_years_index = cell_years_index,
                               classes_index = classes_index,
                               types = types,
                               options = model_options,
                               x_cells_init = x_cells_init)
m <- built$model
# the options as built, with the centre of the initial-state covariates
model_options <- built$options

# the variables and derived quantities, as named objects in the saved image,
# which the figure and prediction scripts read
list2env(built$variables, globalenv())
list2env(built$terms, globalenv())
effect_type <- exp(beta_type)
# the states at every cell, type and year, as cells x types x years (created
# after model(), so they are not computed while sampling)
dynamic_cells <- list(all_states = built$all_states())
population_mortality_vec <- built$population_mortality_vec
country_region_index <- built$lookups$country_region_index
cell_country_lookup <- built$lookups$cell_country_lookup

# the sampler settings (R/dynamical_model.R), as for the folds: more chains
# cost about proportionally more per iteration or worse, and did not give more
# effective samples per core-second (doc/cv_run_plan.md, section 3)
settings <- eval(str2lang(Sys.getenv("IR_CUBE_MCMC_SETTINGS",
                                     "dynamical_mcmc_settings()")))

# used cached posterior means as inits: those of dynamical_inits_file, or
# with IR_CUBE_INITS, of each of its files for a share of the chains
inits <- dynamical_chain_inits(dynamical_inits_files(), built$variables,
                               levels = built$lookups$levels,
                               columns = colnames(x_cell_years),
                               n_chains = settings$n_chains,
                               options = model_options,
                               classes_index = classes_index)

# the time is kept in the saved image, to compare runs on effective samples
# per hour (R/centred_pilot.R)
sampling_time <- system.time(
  draws <- run_dynamical_mcmc(m, built$variables, inits, settings)
)
sampling_time

# check convergence
rhats <- coda::gelman.diag(draws,
                           autoburnin = FALSE,
                           multivariate = FALSE)
summary(rhats$psrf)

# the mortality-floor mode and log posterior of each chain (#37)
chain_modes <- chain_floor_modes_safely(m, draws)

# save fitted model to use for plotting and predictions
save.image(file = "temporary/fitted_model.RData")


# save posterior means as initial values for a future model run

# these have to match the arguments to model(), and need to be greta variable
# nodes (not operation nodes)
posts <- do.call(calculate,
                 c(built$variables, list(values = draws, nsim = 100)))
# those of the non-centred model, if some levels were centred: the cache is
# always in the non-centred variables, which dynamical_inits() moves to any
# model's
posts <- noncentred_draws(posts, classes_index, model_options)

post_means <- lapply(posts,
                     function(x) {
                       apply(x, 2:3, mean)
                     })
attr(post_means, "columns") <- colnames(x_cell_years)
attr(post_means, "levels") <- built$lookups$levels
saveRDS(post_means, dynamical_inits_file)


