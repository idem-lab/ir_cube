# plot estimated baseline susceptibility

# load packages and functions
# greta first, so python starts before terra and sf are attached (a fit that
# did not pass logit_init_mean to model() needs its saved dag)
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

# posterior draws of the parameters
fit_env <- new.env()
load("temporary/fitted_model.RData", envir = fit_env)
types <- fit_env$types
fold <- list(draws = fit_env$draws, options = fit_env$model_options)
parameters <- dynamical_parameter_draws(fold, fit_env$classes_index, types,
                                        fit_env$df)

country_borders <- readRDS("data/clean/country_borders.RDS")

# logit initial fraction susceptible for every country, with countries and
# regions without data drawn from the hierarchical prior, as in R/predict.R
# (same draws and seed, so the same initial states as the maps). With
# initial-state covariates, this is the initial state at the covariates' mean
# (they are standardised)
set.seed(1)
logit_init_all <- map_logit_init(parameters, fit_env$countries,
                                 fit_env$regions, fit_env$df)
rm(fit_env)
logit_init_all <- sweep(logit_init_all, 3, parameters$init_min,
                        FUN = floored_logit)
stopifnot(all(country_borders$country_name %in% dimnames(logit_init_all)[[2]]))

init_all <- plogis(logit_init_all)
country_init_sims <- expand_grid(
  country_name = dimnames(init_all)[[2]],
  insecticide = types
) %>%
  mutate(
    susc_mean = as.vector(t(apply(init_all, 2:3, mean))),
    susc_lower = as.vector(t(apply(init_all, 2:3, quantile, 0.025))),
    susc_upper = as.vector(t(apply(init_all, 2:3, quantile, 0.975)))
  )

init_polys <- country_borders %>%
  left_join(
    country_init_sims,
    by = join_by(country_name)
  )

ggplot() +
  geom_sf(
    aes(
      fill = susc_mean
    ),
    data = init_polys
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    high = "palegreen",
    name = "Initial\nsusceptibility",
    limits = c(0.5, 1),
    na.value = "transparent") +
  facet_wrap(~insecticide, ncol = 3) +
  theme_ir_maps()

ggsave("figures/initial_susceptibility_map.png",
       bg = "white",
       width = 8,
       height = 8,
       dpi = 300)
