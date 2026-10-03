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

# set the start of the timeseries considered in modelling (the start of
# non-negligible levels of resistance) - assume it's before the mass-rollout of
# nets
baseline_year <- 1995

# set the final year of data (insufficient and spatially biased data for 2025)
final_data_year <- 2024

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")

# the modelled subset: the insecticide types with at least 1000 unique
# places/times, and alpha-cypermethrin (914 unique) because of its use in LLINs,
# each at its modal concentration, in the modelled years and inside the mask
# (R/bioassay_subset.R)
df <- subset_modelled_bioassays(ir_africa,
                                mask,
                                baseline_year = baseline_year,
                                final_data_year = final_data_year)

# indices to the classes, types, regions, countries and cells, and the class
# of each type (R/bioassay_subset.R)
list2env(index_bioassays(df), environment())
years <- baseline_year - 1 + sort(unique(df$year_id))

# the model terms (R/dynamical_model.R)
model_options <- dynamical_model_options()

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
settings <- dynamical_mcmc_settings()

# used cached posterior means as inits
inits_one <- dynamical_inits(readRDS(dynamical_inits_file), built$variables,
                             levels = built$lookups$levels,
                             columns = colnames(x_cell_years))

system.time(
  draws <- run_dynamical_mcmc(m, built$variables, inits_one, settings)
)

# check convergence
rhats <- coda::gelman.diag(draws,
                           autoburnin = FALSE,
                           multivariate = FALSE)
summary(rhats$psrf)

# save fitted model to use for plotting and predictions
save.image(file = "temporary/fitted_model.RData")


# save posterior means as initial values for a future model run

# these have to match the arguments to model(), and need to be greta variable
# nodes (not operation nodes)
posts <- do.call(calculate,
                 c(built$variables, list(values = draws, nsim = 100)))

post_means <- lapply(posts,
                     function(x) {
                       apply(x, 2:3, mean)
                     })
attr(post_means, "columns") <- colnames(x_cell_years)
attr(post_means, "levels") <- built$lookups$levels
saveRDS(post_means, dynamical_inits_file)


