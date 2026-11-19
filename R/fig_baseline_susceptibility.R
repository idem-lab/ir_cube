# map the estimated baseline (1995) susceptibility per cell, for both models:
#   figures/initial_susceptibility_map.png            the dynamical model: the
#     country level plus the effects of the initial-state covariates
#   figures/initial_susceptibility_map_two_stage.png  the two-stage model,
#     ilogit(m + omega) (xi is 0 in 1995), from R/two_stage_maps.R
# The two are not quite like for like: the first is q0, before any selection,
# while m in the second is the 1995 state, after that year's selection and
# reversion (dynamical_logit_cells()). They differ by omega plus that one
# year, which is small
# The dynamical model's initial states are those the two-stage maps use: the
# paired draws saved by `R/two_stage_maps.R prepare`, with the same draws for
# countries and regions without data

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/dynamical_predictions.R")

# the parameter draws and every country's logit relative initial state, as
# saved for the two-stage maps (the same for every type's file)
dynamical <- readRDS(file.path("outputs/two_stage/maps",
                               insecticides_plot_order[1], "dynamical.rds"))
types <- dynamical$parameters$types
design <- dynamical$parameters$options$selection_columns

# an even subset of the draws: the posterior mean of the initial state needs
# few, and every cell is computed for each
n_draws <- 200
keep <- round(seq(1, dynamical$parameters$n_draws, length.out = n_draws))
parameters <- subset_draws(dynamical$parameters, keep)
logit_init <- dynamical$logit_init[keep, , , drop = FALSE]
rm(dynamical)

# the initial state at every cell of the mask: its country's level plus the
# effects of the cell's initial-state covariates (cell_logit_init(), as in
# the maps). Cells without a country in the lookup are NA
mask <- rast("data/clean/raster_mask.tif")
cells <- terra::cells(mask)
cell_country <- as.character(
  terra::extract(rast("data/clean/country_raster.tif"), cells)$country_name)
country_index <- match(cell_country, dimnames(logit_init)[[2]])
x_init <- if (!is.null(parameters$init_coef)) {
  init_covariate_matrix(cells, design)
}
ok <- which(!is.na(country_index))
init_mean <- matrix(NA_real_, length(cells), length(types),
                    dimnames = list(NULL, types))
for (rows in split(ok, ceiling(seq_along(ok) / 1e5))) {
  for (k in seq_along(types)) {
    logit_q0 <- cell_logit_init(
      parameters, k,
      matrix(logit_init[, country_index[rows], k], n_draws),
      x_init[rows, , drop = FALSE])
    init_mean[rows, k] <- colMeans(plogis(logit_q0))
  }
}

dynamical_raster <- rast(mask, nlyrs = length(types))
values_full <- matrix(NA_real_, ncell(mask), length(types))
values_full[cells, ] <- init_mean
values(dynamical_raster) <- values_full
names(dynamical_raster) <- types

# the two-stage model's 1995 map
two_stage_raster <- rast(ir_map_files(insecticides_plot_order, 1995))
names(two_stage_raster) <- insecticides_plot_order

pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")
country_borders <- readRDS("data/clean/country_borders.RDS")

# one panel per insecticide, in class order; susceptibility below 50% is shown
# at the end of the scale
plot_initial <- function(raster, file) {
  raster <- terra::mask(raster[[insecticides_plot_order]], pf_water_mask)
  ggplot() +
    geom_sf(data = country_borders, linewidth = 0, fill = grey(0.75)) +
    geom_spatraster(data = raster) +
    geom_sf(data = country_borders, col = grey(0.4), linewidth = 0.1,
            fill = "transparent") +
    scale_fill_gradient(
      labels = scales::percent,
      high = "palegreen",
      name = "Initial\nmortality",
      limits = c(0.5, 1),
      oob = scales::squish,
      na.value = "transparent") +
    facet_wrap(~lyr, ncol = 3) +
    theme_ir_maps()
  ggsave(file,
         bg = "white",
         width = 8,
         height = 8,
         dpi = 300)
}

plot_initial(dynamical_raster, "figures/initial_susceptibility_map.png")
plot_initial(two_stage_raster,
             "figures/initial_susceptibility_map_two_stage.png")
