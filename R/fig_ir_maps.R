# plot IR predictions

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load admin borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")

# load mask with limits of transmission and water bodies for plotting
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)
country_borders <- geom_sf(data = borders,
                           col = border_col,
                           linewidth = 0.1,
                           fill = "transparent")

insecticides_plot <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin",
                       "Permethrin",
                       "Fenitrothion", 
                       "Malathion",
                       "Pirimiphos-methyl", 
                       "DDT",
                       "Bendiocarb",
                       "llin_effective")


n_insecticides <- length(insecticides_plot) - 1
insecticides_col <- c(
  # for the base insecticides
  rev(scales::hue_pal()(n_insecticides)),
  # for LLIN effective
  "#56B1F7"
)
  
for (this_insecticide in insecticides_plot) {
  
  this_index <- match(this_insecticide,
                      insecticides_plot)
  plot_years <- c(2000, 2005, 2010, 2015, 2020, 2025, 2030)
  
  directory <- file.path("outputs/ir_maps", this_insecticide)
  files <- file.path(directory,
                     sprintf("ir_%s_susceptibility.tif",
                             plot_years))
  this_raster <- rast(files)
  names(this_raster) <- plot_years
  
  # mask by transmission limits
  this_raster <- terra::mask(this_raster, pf_water_mask)
  
  insecticide_name <- switch(this_insecticide,
                             llin_effective = "LLIN pyrethroids",
                             this_insecticide
  )
  
  years_list <- list()
  for (i in seq_along(plot_years)) {
    years_list[[i]] <-   ggplot() +
      africa_bg +
      geom_spatraster(
        data = this_raster[[i]]
      ) +
      country_borders +
      scale_fill_gradient(
        labels = scales::percent,
        name = "Susceptibility",
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        high = insecticides_col[this_index],
        low = "white",
        na.value = "transparent",
        guide = guide_colorbar(frame.colour = border_col,
                               frame.linewidth = 0.1)) +
      facet_wrap(~lyr, nrow = 1, ncol = 1) +
      theme_ir_maps() +
      theme(
        plot.margin = unit(rep(0, 4), "cm"),
        legend.text.position = "left",
        legend.ticks = element_blank()
      )
  }
  
  plot_list <- c(
    years_list,
    list(patchwork::guide_area())
  )
  
  patchwork::wrap_plots(plot_list) +
    patchwork::plot_layout(
      guides = "collect",
      nrow = 2
    ) +
    patchwork::plot_annotation(
      title = insecticide_name,
      subtitle = "Susceptibility of An. gambiae (s.l./s.s.) in WHO bioassays"
    )
  
  figure_name <- file.path(
    "figures",
    sprintf("%s_ir_map.png",
            this_insecticide)
  )
  
  ggsave(figure_name,
         bg = "white",
         width = 13,
         height = 8,
         scale = 0.8,
         dpi = 300)
  
}
  