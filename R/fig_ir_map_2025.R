# figure 2: 2025 maps of insecticide resistance 

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")

ir_yr <- 2025

# load admin borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")

# load mask with limits of transmission and water bodies for plotting
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

# load in and prepare estimate rasters

# First, the pyrethroids used in LLINs:

llin_filename <- sprintf(
  "outputs/ir_maps/llin_effective/ir_%s_susceptibility.tif",
  ir_yr
)

llin_rast <- rast(llin_filename)
names(llin_rast) <- "A) LLIN pyrethroids*"

# now, all the individual ones
insecticide_names <- c(
  "Alpha-cypermethrin",
  "Deltamethrin",
  "Lambda-cyhalothrin",
  "Permethrin",
  "Fenitrothion", 
  "Malathion",
  "Pirimiphos-methyl", 
  "DDT",
  "Bendiocarb"
)

ir_filenames <- sprintf(
  "outputs/ir_maps/%s/ir_%s_susceptibility.tif",
  insecticide_names,
  ir_yr
)

ir_rasts <- rast(ir_filenames)

llin_pyrethroid_names <- c("Alpha-cypermethrin",
                           "Deltamethrin",
                           "Permethrin")
llin_pyrethroid_star <- ifelse(insecticide_names %in% llin_pyrethroid_names,
                          "*",
                          "")

# make figure names, with panel letters and asterisks
insecticide_figure_names <- paste0(
  LETTERS[2:10],
  ") ",
  insecticide_names,
  llin_pyrethroid_star)

names(ir_rasts) <- insecticide_figure_names

# change the coordinates for plotting all rasters, to ensure the plots are approximately square
new_ratio <- 0.8

ext <- as.vector(ext(llin_rast))
diffs <- diff(ext)[c(1, 3)]

# keep the y limits the same, and adjust the x limits
xmid <- ext[1] + 0.5 * diffs[1]
xdiff <- diffs[2] / new_ratio

# rebuild the limits
xlim <- xmid + xdiff * c(-0.5, 0.5)

llin_rast_mask <- terra::mask(llin_rast, pf_water_mask)
ir_rasts_mask <- terra::mask(ir_rasts, pf_water_mask)

# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)

# make a plot of susceptibility to Pyrethroid insecticides in 2025
pyrethroid_fig <- ggplot() +
  africa_bg +
  geom_spatraster(
    data = llin_rast_mask,
  ) +
  geom_sf(data = borders,
          col = border_col,
          linewidth = 0.1,
          fill = "transparent") +
  facet_wrap(~lyr) +
  scale_fill_gradient(
    labels = scales::percent,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    high = "#56B1F7",
    low = "white",
    na.value = "transparent",
    name = "Susceptibility",
    guide = guide_colorbar(frame.colour = border_col,
                           frame.linewidth = 0.1)
  ) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0),
    plot.margin = unit(rep(0, 4), "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.25),
    legend.text.position = "left",
    legend.ticks = element_blank()
  )

# make 9 different plots, for each of the individual insecticides
n_insecticides <- length(insecticide_figure_names)
all_insecticides_fig_list <- list()
all_insecticides_cols <- rev(scales::hue_pal()(n_insecticides))

for (i in seq_len(n_insecticides)) {
  all_insecticides_fig_list[[i]] <- ggplot() +
    africa_bg +
    geom_spatraster(
      data = ir_rasts_mask[[i]]
    ) +
    geom_sf(data = borders,
            col = border_col,
            linewidth = 0.05,
            fill = "transparent") +
    facet_wrap(~lyr) +
    scale_fill_gradient(
      labels = scales::percent,
      name = "",
      breaks = c(0, 1),
      high = all_insecticides_cols[i],
      low = "white",
      limits = c(0, 1),
      na.value = "transparent",
      guide = guide_colorbar(frame.colour = border_col,
                             frame.linewidth = 0.05)
    ) +
    coord_sf(
      xlim = xlim,
      ylim = c(min(ir_africa$latitude),
               max(ir_africa$latitude))
    ) +
    theme_ir_maps() +
    theme(
      strip.text.x = element_text(hjust = 0),
      plot.margin = unit(rep(0, 4), "cm"),
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.3),
      legend.key.height = rel(0.3),
      legend.key.width = rel(0.5),
      legend.ticks = element_blank(),
      legend.text.position = "left",
      legend.text = element_text(size = rel(0.6),
                                 margin = margin(r = 1))
    )
}

# make a 9-panel plot of susceptibility to all insecticides in 2025
all_insecticides_fig <- patchwork::wrap_plots(all_insecticides_fig_list,
                                              ncol = 3) +
  plot_layout()

# combine them to plot
combined_fig <- pyrethroid_fig + all_insecticides_fig

# save the plot
ggsave(
  filename = "figures/ir_map_all_insecticides_2025.png",
  plot = combined_fig,
  bg = "white",
  width = 12,
  height = 6,
  scale = 0.8
)

