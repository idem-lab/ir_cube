# Plot the decline (2000-2025 and projected 2025-2030) in the effectiveness of
# ITNs across Africa due to IR. Plot the decline in both direct entomological
# protection, and in epidemiological outcomes, based on Imperial and MAP models.

# four panel plot:

# A) Declining effectiveness of nets: timeseries with ribbon of uncertainty
# (from Tas' estimates) giving average effectiveness of nets across Africa,
# relative to a fully susceptible population

# B) Declining impact of net distribution:
# plot the 'no resistance' line in the second plot
# add ribbon of uncertainty on Tas'

# C) Spatial variation in effectiveness: plot quantiles of IR (2025) on a map

# D) Significant spatial variation in impact of net distributions:
# plot overall impact of nets as averages in these quantiles


# load packages and functions
source("R/packages.R")
source("R/functions.R")

# data to load:

# spatial:
# - Country borders for plotting
# - pf transmission limits and water mask for plotting
# - IR layers
# - MAP net use (to manually project to 2030)
# - MAP PR 2000 baseline
# - MAP populations by year (projected to 2030)
# other:
# - Symons epi impact estimate
# - Nash/Sherrard-Smith killing effect estimates

# years to plot for
years_plot <- 2000:2030

# load admin borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")

# load mask with limits of transmission and water bodies for plotting
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

# - IR layers (susceptibility to the pyrethroids used in LLINs)
ir_filenames <- sprintf(
  "outputs/ir_maps/llin_effective/ir_%s_susceptibility.tif",
  years_plot
)
ir <- rast(ir_filenames)
names(ir) <- years_plot

# - MAP net use (to manually project to 2030)
nets <- rast("data/clean/net_use_cube.tif")
names(nets) <- str_remove(names(nets), "^nets_")
latest_nets <- max(as.numeric(names(nets)))
end_year <- max(years_plot)
nets <- post_pad_cube(nets, end_year)

# - MAP PR 2000 baseline
pr2000 <- rast("data/clean/pfpr_2000.tif")

# - MAP populations by year (projected to 2030)
pop <- rast("data/clean/pop_cube.tif")
names(pop) <- str_remove(names(pop), "^pop_")
pop_future <- rast("data/clean/pop_cube_future.tif")
names(pop_future) <- str_remove(names(pop_future), "^pop_")
pop_all <- c(pop, pop_future)

# other:

# Encode Symons et al. PMI function for relative effectiveness of ITNs at
# reducing logit-PfPR (reduction in on the log-odds ratio), compared to in the absence
# of resistance.

# when susceptibility is 100%, the ITN effectiveness parameter has full effect
# (multiplier 1), when susceptibility is 0%, the ITN effectiveness parameter has
# 54% (multiplier 0.54).

# compute the relative epi impact from the posterior mean of the parameters
epi_function_mean <- function(susceptibility) {
  1 - 0.46 * (1 - susceptibility)
}

# The log-odds ratio for malaria infection under full ITN coverage (relative to
# no ITN coverage) is -1.68, which corresponds to an odds ratio of 0.186374.
# With a fully-resistant vector population, this is would reduce to 54% of its
# previous value, -0.9072, giving an odds ratio of 0.4036529.

# Alternative interpretation: the odds ratio for a population with 100% coverage
# and 100% resistance is equivalent to that of a population with 54% coverage
# and 0% resistance.

# get the estimate and SD. for the parameters behind this model, to be able to
# propagate variance through the subsequent calculations

# beta_0 = beta_{physical} = 'ITN_physical'
#   beta_0 ~ N(-1.68, 0.11^2)
# beta_1 = beta_{penalty} = 'ITN_insecticidal_penalty'
#   beta_1 ~ N(0.78, 0.12^2)

# the log odds (contribution to the linear predictor) is computed from these
# parameters, with a non-linear function of ITN coverage f(ITN) and
# susceptibility k:

#   log_odds = beta_0 * f(ITN) + beta_1 * (1 - k) * f(ITN)
#            = f(ITN) * (beta_0 + beta_1 * (1 - k))

# we can therefore factorise to get a net effectiveness measure g(k):
#   log_odds = f(ITN) * g(k)
#   g(k) = beta_0 + beta_1 * (1 - k)

# from which we can obtain a *relative* net effectiveness, giving the
# effectiveness relative to a fully susceptible population (k = 1):
#   g*(k) = g(k) / g(1)

# since:
#   g(1) = beta_0 + beta_1 * (1 - k)
#        = beta_0 + beta_1 * 0
#        = beta_0
# we have:
#   g*(k) = g(k) / beta_0
#         = (beta_0 + beta_1 * (1 - k)) / beta_0
#         = 1 + (beta_1 * (1 - k)) / beta_0

# this is equivalent to a ratio of normal distributions, but with non-zero means
# it doesn't have a closed form (that I know of), so I will compute credible
# intervals of effectiveness by Monte Carlo simulation.



# generate a version of the epi function, drawn at random from the posteriors
# given in Symons et al. Note we are assuming that they are independent, since
# the posterior correlation is not given.
make_stochastic_epi_function <- function() {
  beta_0 <- rnorm(1, -1.68, 0.11)
  beta_1 <- rnorm(1, 0.78, 0.12)
  function(susceptibility) {
    1 + (beta_1 * (1 - susceptibility)) / beta_0
  }
}


# make a batch of them
n_sims <- 100
epi_functions_list <- replicate(n_sims,
                                make_stochastic_epi_function(),
                                simplify = FALSE)


# # check this is defined correctly:
# susc <- seq(0, 1, length.out = 100)
# plot(epi_function(susc) ~ susc,
#      ylim = c(0, 1),
#      type = "l")
# for (i in seq_len(n_sims)) {
#   lines(epi_functions_list[[i]](susc) ~ susc,
#         lwd = 0.1)
# }

# perform calculations:

# define the weights for net effectiveness (weight by where nets are in use) and
# population-level impact (weight by population at risk)

# convert a raster layer into spatial weights (pixels sum to 1) to get averages
# over space
normalise <- function(x) {
  total <- global(x, "sum", na.rm = TRUE)[1, 1]
  x / total
}

# cube of net use weights
net_weights <- sapp(nets, normalise)

# cube of population at risk (potential infected population, in absence of
# interventions) weights
infection_weights <- sapp(pop_all * pr2000, normalise)


# First, get the average net use across Africa, weighted by population at risk.
# This gives us the fraction of the maximum possible impact our net
# distributions have, in the absence of resistance. If net use was 100% for all
# populations at risk, this would equal 1. It is lower due to imperfect net
# coverage (but not due to IR).
net_impact_no_ir <- global(
  x = nets * infection_weights,
  fun = "sum",
  na.rm = TRUE
)


# Now average epi impacts, accounting for IR (both per net and at population
# level)

# cube of net relative epi impact of a net in each pixel (using the posterior
# mean function)
relative_epi_impact_mean <- terra::app(ir, epi_function_mean)
names(relative_epi_impact_mean) <- names(ir)

# average the epi impact of the nets in use across Africa to get the temporal
# trends in impact of each net (using the posterior mean function)
average_epi_impact_mean <- global(
  x = relative_epi_impact_mean * net_weights,
  fun = "sum",
  na.rm = TRUE
)

# average the epi impact of the whole net distribution program across Africa to
# get the temporal trends in impact of net distribution (using the posterior
# mean function)
net_impact_epi_impact_mean <- global(
  x = nets * relative_epi_impact_mean * infection_weights,
  fun = "sum",
  na.rm = TRUE
)


# compute population-at-risk-weighted quantiles of 2025 LLIN effectiveness, to
# divide up the map into different regions
vals <- infection_weights[["2025"]][]
cells <- which(!is.na(vals) & vals > 0)

impact_2025 <- relative_epi_impact_mean[["2025"]][cells][, 1]
infection_weights_2025 <- infection_weights[["2025"]][cells][, 1]

# n_bands <- 3
# cuts <- seq(0, 1, length.out = n_bands + 1)
cuts <- c(0, 0.5, 0.9, 1)
quants <- Hmisc::wtd.quantile(x = impact_2025,
                              weights = infection_weights_2025,
                              probs = cuts,
                              normwt = TRUE)

# redefine endpoints, in case there are out-of-bounds effectiveness results
# outside PAR
n_quants <- length(quants)
quants[1] <- 0
quants[n_quants] <- 1

# do reclassification on ir_2025 to map the bands. Do 'right = NA' to include
# both upper and lower bounds
rcl <- cbind(from = quants[-n_quants],
             to = quants[-1],
             becomes = seq_len(n_quants - 1))
bands <- terra::classify(relative_epi_impact_mean[["2025"]],
                         rcl,
                         right = NA)
names <- data.frame(
  from = 3:1,
  to = c("highest", "medium", "lowest")
)
levels(bands) <- names

# create net weights and infection weights subsetted to each of these bands, to
# compute effectiveness over time in each

nets_highest <- nets * (bands == "highest")
net_weights_highest <- sapp(nets_highest, normalise)

nets_medium <- nets * (bands == "medium")
net_weights_medium <- sapp(nets_medium, normalise)

nets_lowest <- nets * (bands == "lowest")
net_weights_lowest <- sapp(nets_lowest, normalise)

# compute posterior means for these
average_epi_impact_mean_highest <- global(
  x = relative_epi_impact_mean * net_weights_highest,
  fun = "sum",
  na.rm = TRUE
)

average_epi_impact_mean_medium <- global(
  x = relative_epi_impact_mean * net_weights_medium,
  fun = "sum",
  na.rm = TRUE
)

average_epi_impact_mean_lowest <- global(
  x = relative_epi_impact_mean * net_weights_lowest,
  fun = "sum",
  na.rm = TRUE
)

# 
# plot(average_epi_impact_mean_lowest$sum ~ years_plot,
#      ylim = c(0, 1),
#      type = "l")
# lines(average_epi_impact_mean_medium$sum ~ years_plot)
# lines(average_epi_impact_mean_highest$sum ~ years_plot)
# abline(v = latest_nets)

# Now do the same under a Monte Carlo simulation, to get credible intervals on
# these
n_years <- length(years_plot)
individual_net_impact_sims <- 
  individual_net_impact_sims_highest <- 
  individual_net_impact_sims_medium <- 
  individual_net_impact_sims_lowest <- 
  net_distribution_impact_sims <- matrix(NA, nrow = n_years, ncol = n_sims)
rownames(individual_net_impact_sims) <-
  rownames(individual_net_impact_sims_highest) <-
  rownames(individual_net_impact_sims_medium) <-
  rownames(individual_net_impact_sims_lowest) <- 
  rownames(net_distribution_impact_sims) <- years_plot

for (i in seq_len(n_sims)) {
  print(i)
  
  relative_epi_impact_sim <- terra::app(ir, epi_functions_list[[i]])
  names(relative_epi_impact_sim) <- names(ir)
  
  average_epi_impact_sim <- global(
    x = relative_epi_impact_sim * net_weights,
    fun = "sum",
    na.rm = TRUE
  )
  
  individual_net_impact_sims[, i] <- average_epi_impact_sim$sum
  
  # same for low, medium, high
  average_epi_impact_sim_highest <- global(
    x = relative_epi_impact_sim * net_weights_highest,
    fun = "sum",
    na.rm = TRUE
  )
  average_epi_impact_sim_medium <- global(
    x = relative_epi_impact_sim * net_weights_medium,
    fun = "sum",
    na.rm = TRUE
  )
  average_epi_impact_sim_lowest <- global(
    x = relative_epi_impact_sim * net_weights_lowest,
    fun = "sum",
    na.rm = TRUE
  )
  individual_net_impact_sims_highest[, i] <- average_epi_impact_sim_highest$sum
  individual_net_impact_sims_medium[, i] <- average_epi_impact_sim_medium$sum
  individual_net_impact_sims_lowest[, i] <- average_epi_impact_sim_lowest$sum
  
  # now (for all bands) do the impact of net distribution
  net_impact_epi_impact_sim <- global(
    x = nets * relative_epi_impact_sim * infection_weights,
    fun = "sum",
    na.rm = TRUE
  )
  
  net_distribution_impact_sims[, i] <- net_impact_epi_impact_sim$sum
  
}


# compute credible intervals of each of these over time
individual_net_impact_cis <- apply(individual_net_impact_sims,
                                   1,
                                   quantile,
                                   c(0.25, 0.75))
individual_net_impact_cis_highest <- apply(individual_net_impact_sims_highest,
                                          1,
                                          quantile,
                                          c(0.25, 0.75))
individual_net_impact_cis_medium <- apply(individual_net_impact_sims_medium,
                                          1,
                                          quantile,
                                          c(0.25, 0.75))
individual_net_impact_cis_lowest <- apply(individual_net_impact_sims_lowest,
                                          1,
                                          quantile,
                                          c(0.25, 0.75))
net_distribution_impact_cis <- apply(net_distribution_impact_sims,
                                     1,
                                     quantile,
                                     c(0.25, 0.75))


# panel A: Declining effectiveness of nets
# mean and CI of "individual_net_impact" over time

# all data
panel_a_data <- tibble(
  year = years_plot,
  mean = average_epi_impact_mean$sum,
  upper = individual_net_impact_cis[1, ],
  lower = individual_net_impact_cis[2, ]
)

# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_a_data_past <- panel_a_data %>%
  filter(
    year <= latest_nets
  )

panel_a_data_future <- panel_a_data %>%
  filter(
    year >= latest_nets
  )

panel_a <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    name = "LLIN effectiveness"
  ) +
  # past data
  geom_ribbon(
    data = panel_a_data_past,
    fill = "#56B1F7"
  ) +
  geom_line(
    data = panel_a_data_past,
    colour = grey(0.2),
    # linewidth = 0.8
  ) + 
  # future projections; 1/3 of the alpha and dashed lines
  geom_ribbon(
    data = panel_a_data_future,
    fill = "#56B1F7",
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = panel_a_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3,
    # linewidth = 0.8
  ) + 
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  xlab("") +
  theme_minimal()


# B) Declining impact of net distribution
# plot the 'no resistance' line in the second plot
# add ribbon of uncertainty on Tas'

# all data
panel_b_data <- tibble(
  year = years_plot,
  mean = net_impact_epi_impact_mean$sum,
  upper = net_distribution_impact_cis[1, ],
  lower = net_distribution_impact_cis[2, ],
  mean_no_ir = net_impact_no_ir$sum  
)

# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_b_data_past <- panel_b_data %>%
  filter(
    year <= latest_nets
  )

panel_b_data_future <- panel_b_data %>%
  filter(
    year >= latest_nets
  )

panel_b_final <- panel_b_data %>%
  filter(
    year == max(year)
  )

panel_b <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.6),
    name = "Impact of LLIN distribution"
  ) +
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  # net impact with no IR, past then future
  geom_line(
    data = panel_b_data_past,
    aes(
      y = mean_no_ir
    ),
    colour = grey(0.5),
    # linewidth = 0.8
  ) +
  # net impact with no IR, past then future
  geom_line(
    data = panel_b_data_future,
    aes(
      y = mean_no_ir
    ),
    colour = grey(0.5),
    # linewidth = 0.8,
    linetype = 3
  ) +
  # past data
  geom_ribbon(
    data = panel_b_data_past,
    fill = "#56B1F7"
  ) +
  geom_line(
    data = panel_b_data_past,
    colour = grey(0.2),
    # linewidth = 0.8
  ) + 
  # future projections; 1/3 of the alpha and dashed lines
  geom_ribbon(
    data = panel_b_data_future,
    fill = "#56B1F7",
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = panel_b_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3,
    # linewidth = 0.8
  ) + 
  xlab("") +
  theme_minimal()

# C) Spatial variation in effectiveness: plot quantiles of LLIN effectiveness
# (2025) on a map

bands_mask <- terra::mask(bands, pf_water_mask)

# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)

band_pal <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(9, "RdPu"))
)
band_cols <- band_pal(5)[2:4]

panel_c <- ggplot() +
  africa_bg +
  geom_spatraster(
    data = bands_mask,
  ) +
  geom_sf(data = borders,
          col = border_col,
          linewidth = 0.1,
          fill = "transparent") +
  scale_fill_manual(
    values = c(
      lowest = band_cols[3],
      medium = band_cols[2],
      highest = band_cols[1]
    ),
    na.translate = FALSE,
    # note we arre using element_markdown() for legend titles so need this html
    # linebreak version
    name = "LLIN effectiveness<br>(2025)",
    na.value = "transparent"
  ) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0),
    plot.margin = unit(rep(0, 4), "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.25),
    legend.key.height = rel(0.8),
    legend.key.width = rel(0.8),
    legend.text.position = "left",
    legend.ticks = element_blank(),
    legend.title = element_markdown(size = rel(0.8),
                                    hjust = 1)
  )

# This displays the continent in classes of LLIN effectiveness, with the third
# of the population at risk with the highest LLIN effectiveness (>25%
# susceptibility), the third with the lowest LLIN effectiveness (>25%
# susceptibility)





# panel D should show effectiveness in these bands, in these colours, like panel
# A


panel_d_data <- tibble(
  year = years_plot,
  mean_highest = average_epi_impact_mean_highest$sum,
  upper_highest = individual_net_impact_cis_highest[1, ],
  lower_highest = individual_net_impact_cis_highest[2, ],
  mean_medium = average_epi_impact_mean_medium$sum,
  upper_medium = individual_net_impact_cis_medium[1, ],
  lower_medium = individual_net_impact_cis_medium[2, ],
  mean_lowest = average_epi_impact_mean_lowest$sum,
  upper_lowest = individual_net_impact_cis_lowest[1, ],
  lower_lowest = individual_net_impact_cis_lowest[2, ],
) %>%
  pivot_longer(
    cols = starts_with(c("mean", "upper", "lower")),
    names_pattern = "(.*)_(.*)",
    names_to = c("statistic", "band"),
    values_to = "effectiveness"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = effectiveness
  ) %>%
  mutate(
    band = factor(band,
                  levels = c("highest",
                             "medium",
                             "lowest"))
  )

# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_d_data_past <- panel_d_data %>%
  filter(
    year <= latest_nets
  )

panel_d_data_future <- panel_d_data %>%
  filter(
    year >= latest_nets
  )

panel_d <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper,
    fill = band
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    name = "LLIN effectiveness"
  ) +
  scale_fill_manual(
    values = c(
      lowest = band_cols[3],
      medium = band_cols[2],
      highest = band_cols[1]
    ),
    guide = "none",
    na.translate = FALSE,
    name = "LLIN effectiveness",
    na.value = "transparent",
  ) +
  # past data
  geom_ribbon(
    data = panel_d_data_past,
  ) +
  geom_line(
    data = panel_d_data_past,
    colour = grey(0.2)
  ) + 
  # future projections; 1/3 of the alpha and dashed lines
  geom_ribbon(
    data = panel_d_data_future,
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = panel_d_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3
  ) + 
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  xlab("") +
  theme_minimal()

plots_all <- panel_a +
  panel_b +
  panel_c +
  panel_d +
  plot_layout(
    ncol = 2,
    nrow = 2,
    widths = rep(1, 4)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  filename = "figures/epi_impact.png",
  plot = plots_all,
  bg = "white",
  width = 7,
  height = 6,
  scale = 1,
  dpi = 300
)


# # maybe plot the killing effects from Nash & Sherrard-Smith
# 
# # - Nash/Sherrard-Smith killing effect estimates.
# 
# # load in and prepare curves from Nash/Sherard-Smith, given in supp to Sherard
# # Smith paper: https://doi.org/10.1038/s41467-022-30700-1
# ss_est <- read_excel("data/raw/41467_2022_30700_MOESM3_ESM.xlsx",
#                      sheet = 4,
#                      range = "A15:Y117") %>%
#   # get required columns
#   select(
#     resistance = "Pyrethroid resistance approximated from susceptibility bioassay",
#     kill_logistic = "dN0_1",  # logistic model
#     kill_log_logistic = "dN0_4" # log-logistic model
#   ) %>%
#   # remove empty row
#   slice(
#     -1
#   ) %>%
#   # split apart by model type
#   pivot_longer(
#     cols = starts_with("kill"),
#     names_to = "model",
#     values_to = "kill",
#     names_prefix = "kill_"
#   ) %>%
#   group_by(
#     model
#   ) %>%
#   mutate(
#     susceptibility = 1 - resistance,
#     relative_killing_effect = kill / max(kill)
#   )
# 
# # # plot to check against paper figures
# # ss_est %>%
# #   ggplot(
# #     aes(
# #       x = resistance,
# #       y = kill
# #     )
# #   ) +
# #   facet_wrap(~model) +
# #   geom_line() +
# #   coord_cartesian(ylim = c(0, 1)) +
# #   theme_minimal()
# 
# ss_est_logistic <- ss_est %>%
#   filter(model == "logistic")
# 
# ss_est_log_logistic <- ss_est %>%
#   filter(model == "log_logistic")
# 
# # produce interpolating functions to get the relative killing effect
# killing_function_logistic <- splinefun(
#   x = ss_est_logistic$susceptibility,
#   y = ss_est_logistic$relative_killing_effect
# )
# 
# killing_function_log_logistic <- splinefun(
#   x = ss_est_log_logistic$susceptibility,
#   y = ss_est_log_logistic$relative_killing_effect
# )
# 
# # check the interpolators behave appropriately
# stopifnot(
#   identical(
#     killing_function_logistic(ss_est_logistic$susceptibility),
#     ss_est_logistic$relative_killing_effect
#   )
# )
# 
# stopifnot(
#   identical(
#     killing_function_log_logistic(ss_est_log_logistic$susceptibility),
#     ss_est_log_logistic$relative_killing_effect
#   )
# )
# 
# 
# # make a cube of the killing effect of ITNs, over time and space
# relative_killing_logistic <- terra::app(ir, killing_function_logistic)
# names(relative_killing_logistic) <- names(ir)
# 
# relative_killing_log_logistic <- terra::app(ir, killing_function_log_logistic)
# names(relative_killing_log_logistic) <- names(ir)
# 
# 
# # average the killing effect and epi impact of nets over the nets in use
# average_killing_logistic <- global(
#   x = relative_killing_logistic * net_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# average_killing_log_logistic <- global(
#   x = relative_killing_log_logistic * net_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# # apply the impact of IR and recalculate average
# net_impact_killing_logistic <- global(
#   x = nets * relative_killing_logistic * infection_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# net_impact_killing_log_logistic <- global(
#   x = nets * relative_killing_log_logistic * infection_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# 
# par(mfrow = c(1, 2),
#     mar = c(2, 5, 3, 2))
# 
# # Effectiveness\n of single-AI nets
# no_change <- rep(1, length(years_plot))
# # impact of each net
# plot(no_change ~ years_plot,
#      type = "l",
#      col = grey(0.3),
#      lwd = 2,
#      xlab = "",
#      ylab = "Relative effectiveness",
#      ylim = c(0, 1))
# lines(average_killing_logistic$sum ~ years_plot,
#       col = "red",
#       lwd = 2)
# lines(average_killing_log_logistic$sum ~ years_plot,
#       col = "red",
#       lty = 2,
#       lwd = 2)
# lines(average_epi_impact$sum ~ years_plot,
#       col = "blue",
#       lwd = 2)
# abline(v = latest_nets, lty = 2)
# # directly label these
# 
# plot(net_impact_no_ir$sum ~ years_plot,
#      type = "l",
#      lwd = 2,
#      col = grey(0.3),
#      ylab = "Relative impact",
#      xlab = "",
#      ylim = c(0, 1))
# lines(net_impact_epi_impact$sum ~ years_plot,
#       col = "blue",
#       lwd = 2)
# abline(v = latest_nets, lty = 2)
# title(main = "Impact of\nnet distribution")
# 
# df_plot <- tibble(
#   year = years_plot,
#   `impact_direct/no resistance` = 1,
#   `impact_direct/vector killing (logistic)` = average_killing_logistic$sum,
#   `impact_direct/vector killing (log-logistic)` = average_killing_log_logistic$sum,
#   `impact_direct/transmission reduction` = average_epi_impact$sum,
#   `impact_overall/no resistance` = net_impact_no_ir$sum,
#   `impact_overall/transmission reduction` = net_impact_epi_impact$sum
# ) %>%
#   pivot_longer(
#     cols = starts_with("impact"),
#     names_to = c("scale", "impact_type"),
#     names_pattern = "impact_(.*)/(.*)",
#     values_to = "impact",
#     # names_prefix = "impact_"
#   )
# 
# plot_direct <- df_plot %>%
#   filter(
#     scale == "direct"
#   ) %>%
#   mutate(
#     phase = case_when(
#       year > latest_nets ~ "future",
#       .default = "past"
#     )
#   ) %>%
#   ggplot(
#     aes(
#       x = year,
#       y = impact,
#       colour = impact_type,
#       linetype = phase
#     )
#   ) +
#   geom_vline(
#     aes(
#       xintercept = latest_nets
#     ),
#     col = grey(0.6),
#     linetype = 2
#   ) +
#   geom_line(
#     linewidth = 0.8
#   ) +
#   scale_y_continuous(
#     labels = scales::percent,
#     limits = c(0, 1)
#   ) +
#   scale_colour_manual(
#     values = c(
#       "transmission reduction" = "blue",
#       "vector killing (log-logistic)" = "coral",
#       "vector killing (logistic)" = "coral3",
#       "no resistance" = grey(0.4)
#     )
#   ) +
#   scale_linetype_manual(
#     values = c(
#       "past" = 1,
#       "future" = 3
#     ),
#     guide = "none"
#   ) +
#   ylab(
#     "Effectiveness"
#   ) +
#   xlab(
#     ""
#   ) +
#   theme_minimal() +
#   theme(
#     legend.title = element_blank(),
#     # legend.position = "bottom"
#   )
# 
# plot_direct
