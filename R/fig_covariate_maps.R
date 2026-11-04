# plot covariates

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

# plot time-varying covariates: the layers the model uses
# (R/model_covariates.R), at four years from the baseline to the present
plot_years <- c(1995, 2005, 2015, 2025)
source("R/model_covariates.R")

# a cube's layers at plot_years, padded as the model pads them
# (read_padded_cube()): years before the cube's first are its first layer,
# years after its last are its last; the facet label says which years are
# copies
model_years <- function(cube, years_plot = plot_years) {
  years <- as.numeric(str_sub(names(cube), start = -4L))
  used <- pmin(pmax(years_plot, min(years)), max(years))
  out <- cube[[match(used, years)]]
  names(out) <- ifelse(used == years_plot, years_plot,
                       sprintf("%i\n(assumed same as %i)", years_plot, used))
  out
}

# pyrethroid net use, a conventional ITN counting as 0.25 of an LLIN
# (net_use_file(), R/prep_net_use_pyrethroid.R)
itn_plot <- model_years(rast(net_use_file(selection_design())))

# IRS coverage, scaled to its maximum (which is 1, so this is the coverage)
irs_plot <- model_years(rast("data/clean/irs_coverage_scaled_cube.tif"))

# population density, people per km2, with WorldPop's empty pixels at 0
# (pop_zero_filled())
pop <- pop_zero_filled(rast("data/clean/pop_cube.tif"))
pop_dens_plot <- model_years(pop) /
  terra::cellSize(rast("data/clean/raster_mask.tif"), unit = "km")

# plot ITN
itn_plot_mask <- terra::mask(itn_plot, pf_water_mask)
ggplot() +
  africa_bg +
  geom_spatraster(
    data = itn_plot_mask
  ) +
  country_borders +
  scale_fill_gradient(
    name = "Net use",
    labels = scales::percent,
    limits = c(0, 1),
    low = "pink",
    high = "purple",
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "Pyrethroid net use"
  ) +
  theme_ir_maps()

ggsave("figures/itn_map.png",
       bg = "white",
       width = 14,
       height = 4.5,
       dpi = 300)

# plot IRS
irs_plot_mask <- terra::mask(irs_plot, pf_water_mask)
ggplot() +
  africa_bg +
  geom_spatraster(
    data = irs_plot_mask
  ) +
  country_borders +
  scale_fill_gradient(
    labels = scales::percent,
    low = "pink",
    high = "purple",
    name = "Coverage",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "IRS coverage"
  ) +
  theme_ir_maps()

ggsave("figures/irs_map.png",
       bg = "white",
       width = 14,
       height = 4.5,
       dpi = 300)


# plot pop, coloured on the scale the model sees it: the chance of encounter,
# 1 - exp(-d log(2) / d_half) at density d, with d_half
# selection_design()'s pop_d_half (selection_pop_matrix())
d_half <- selection_design()$pop_d_half
encounter_transform <- scales::trans_new(
  "encounter",
  transform = function(d) -expm1(-d * log(2) / d_half),
  inverse = function(s) -d_half * log1p(-s) / log(2)
)
pop_dens_plot_mask <- terra::mask(pop_dens_plot, pf_water_mask)
ggplot() +
  africa_bg +
  geom_spatraster(
    data = pop_dens_plot_mask
  ) +
  country_borders +
  scale_fill_continuous(
    transform = encounter_transform,
    labels = scales::number_format(big.mark = ","),
    low = "lightblue",
    high = "navy",
    breaks = c(10, 100, d_half, 1000),
    name = "People per km<sup>2</sup>",
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "Population density"
  ) +
  theme_ir_maps()

ggsave("figures/pop_map.png",
       bg = "white",
       width = 14,
       height = 4.5,
       dpi = 300)


# collated total production of crop types
crops_group <- rast("data/clean/crop_group_scaled.tif")

# production of individual crops
crops_all <- rast("data/clean/crop_scaled.tif")

# Pull out crop types implicated in risk for IR. refer to this review, crop type
# section:
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1162-4
crops_implicated <- c(
  # "increased resistance at cotton growing sites, a finding subsequently
  # supported in eight other papers from five different African countries", "the
  # cash crop with the highest intensity insecticide use of any crop"
  crops_all$cotton,
  # "In eight studies, vegetable cultivation strongly related to
  # insecticide-resistant field collections", "Vegetable production requires
  # significantly higher quantities and/or more frequent application of
  # pesticides than other food crops"
  crops_all$vegetables,
  # "Seven of the studies reviewed here examined the insecticide susceptibility
  # of vector populations at rice-growing sites, and found low-to-moderate
  # resistance levels in these mosquito populations."
  crops_all$rice)

# combine all temporally-static covariates
crops <- c(crops_group, crops_implicated)
names(crops) <- str_to_sentence(names(crops))

crops_mask <- terra::mask(crops, pf_water_mask)
ggplot() +
  africa_bg +
  geom_spatraster(
    data = crops_mask
  ) +
  country_borders +
  scale_fill_gradient(
    # labels = scales::percent,
    transform = "sqrt",
    low = "#eef8e4",
    high = "#003d16",
    name = "Relative production",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2) +
  ggtitle(
    label = "Agricultural intensities"
  ) +
  theme_ir_maps()

ggsave("figures/crop_map.png",
       bg = "white",
       width = 14,
       height = 6,
       dpi = 300)
