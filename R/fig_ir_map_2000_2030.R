# figure 3: 2000-2025 maps of changes in net use & LLIN susceptibility

# 6 panel plot, for three time periods: 2000-2010, 2010-2020, 2020-2025. Top row
# showing rate of reduction in susceptibility between periods (in % per year),
# bottom row showing average net use over each period

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load admin borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")

# load mask with limits of transmission and water bodies for plotting
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

# load in and prepare rasters

# First, the LLIN use
nets <- rast("data/clean/net_use_cube.tif")
names(nets) <- str_remove(names(nets), "^nets_")

# Next, Susceptibility to the pyrethroids used in LLINs, of the two-stage model
# (R/two_stage_maps.R):
ir_yrs_all <- 2000:2030
ir <- rast(ir_map_files("llin_effective", ir_yrs_all))
names(ir) <- ir_yrs_all

# compute change in IR, and average net use, over fixed periods

# compute the change from the start to the end of a period. Signed: the
# two-stage model's susceptibility need not fall monotonically, so the range
# over the period would count a rise as a fall
change <- function(x) x[[nlyr(x)]] - x[[1]]

breaks <- c(2000, 2005, 2010, 2015, 2020, 2025)
n_breaks <- length(breaks)
period_start <- breaks[-n_breaks]
period_end <- breaks[-1]
period_name <- paste0(period_start, "-", period_end)
n_periods <- length(period_start)

net_use_list <- list()
ir_change_list <- list()
for(i in seq_len(n_periods)) {
  period_yrs <- seq(period_start[i],
                    period_end[i],
                    by = 1)
  period_ir_change <- change(ir[[names(ir) %in% period_yrs]])
  # alternative: make this annualised to account for different-sized bins?
  # period_ir_change <- period_ir_change / diff(range(period_yrs))
  period_net_use <- mean(nets[[names(nets) %in% period_yrs]])

  names(period_ir_change) <- names(period_net_use) <- period_name[i]
  ir_change_list[[i]] <- period_ir_change
  net_use_list[[i]] <- period_net_use
}

ir_change <- rast(ir_change_list)
net_use <- rast(net_use_list)

# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)

country_borders <- geom_sf(data = borders,
                           col = border_col,
                           linewidth = 0.1,
                           fill = "transparent")

ir_change_mask <- terra::mask(ir_change, pf_water_mask)
ir_change_fig <- ggplot() +
  africa_bg +
  geom_spatraster(data = ir_change_mask) +
  country_borders +
  facet_wrap(~lyr, nrow = 2) +
  # losses in red; any gains in blue
  scale_fill_gradient2(
    labels = scales::percent,
    low = "red",
    mid = grey(0.9),
    high = "#2166ac",
    midpoint = 0,
    na.value = "transparent"
  ) +
  labs(fill = "Change in mortality") +
  theme_ir_maps() +
  theme(
    plot.margin = unit(rep(0, 4), "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.25)
  )

# save the plot
ggsave(
  filename = "figures/ir_map_change.png",
  plot = ir_change_fig,
  bg = "white",
  width = 8,
  height = 6,
  scale = 0.8
)

