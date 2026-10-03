# illustrate the predictive distribution validation metric

# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# get posterior predictive simulations of observations
died_sim <- betabinomial_p_rho(N = df$mosquito_number,
                               p = population_mortality_vec,
                               rho = rho_types[df$type_id])
mortality_sim <- died_sim / df$mosquito_number

# summarise fit to data
died <- calculate(died_sim, values = draws, nsim = 1e3)[[1]][, , 1]

# create a dharma object to compute randomised quantile residuals and
# corresponding residual z scores
dharma <- DHARMa::createDHARMa(
  simulatedResponse = t(died),
  observedResponse = df$died,
  integerResponse = TRUE
)

# Display the GoF statistic (Anderson Darling) for each of these histograms,
# versus a uniform distribution

