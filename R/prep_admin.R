# prep admin unit/region rasters and shapefiles

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load raster extent
mask <- rast("data/clean/raster_mask.tif")

# load all countries and all regions in Africa (per UN)
africa_countries_un <- country_region_lookup()

# harmonise these with the GADM layers. Need resolution <= 3 to include e.g.
# British Indian Ocean Territory, which is in the UN country list
gadm <- geodata::world(resolution = 3, path = "data/raw") %>%
  st_as_sf() %>%
  mutate(
    country_name = case_when(
      NAME_0 == "Central African Republic" ~ "CAR",
      # different apostrophe!
      NAME_0 == "Côte d'Ivoire" ~ "Côte d’Ivoire",
      NAME_0 == "Democratic Republic of the Congo" ~ "DR Congo",
      NAME_0 == "São Tomé and Príncipe" ~ "Sao Tome & Principe",
      .default = NAME_0
    )
  )

# subset down to the GADM shapefiles in the UN Africa region and the MAP mask
gadm_in_mask <- gadm %>%
  right_join(
    africa_countries_un,
    by = join_by(country_name)
  ) %>%
  mutate(
    mask_val = terra::extract(mask,
                              .,
                              fun = "mean",
                              na.rm = TRUE)[, 2],
  ) %>%
  filter(
    !is.na(mask_val)
  ) %>%
  select(-mask_val)

# Rasterise at cell centres, then give each mask cell that no polygon's
# interior covers (coast, lake shore, islands GADM omits) the country of the
# nearest polygon, if that polygon is within fill_cap_km of the cell centre.
# The cap is about one cell width: every cell a polygon touches is within 3.2
# km of it, so all of those are filled. Cells further out are left NA rather
# than given a distant country's initial condition (the uncapped fill put some
# 1,555 km from their country; see issue #15). None of them carry bioassay
# data; predict.R and two_stage_maps.R leave them unpredicted
fill_cap_km <- 5

gadm_raster_country <- gadm_in_mask %>%
  terra::rasterize(mask, field = "country_name")
gadm_raster_region <- gadm_in_mask %>%
  terra::rasterize(mask, field = "region")

missing_cell_coords <- gadm_raster_country %>%
  extract(cells(mask), xy = TRUE) %>%
  filter(is.na(country_name)) %>%
  select(x, y)

# sf rather than terra::nearest, which on lon/lat polygons measures to the
# wrong point (a cell 7.6 km off the Kenyan coast comes out 508 km from Kenya)
# and so picks the wrong country for many cells
missing_points <- missing_cell_coords %>%
  st_as_sf(coords = c("x", "y"),
           crs = st_crs(gadm_in_mask))
nearest_index <- st_nearest_feature(missing_points, gadm_in_mask)
nearest_km <- as.numeric(st_distance(missing_points,
                                     gadm_in_mask[nearest_index, ],
                                     by_element = TRUE)) / 1000
within_cap <- nearest_km <= fill_cap_km
nearest_country <- gadm_in_mask$country_name[nearest_index]
nearest_region <- gadm_in_mask$region[nearest_index]

cat(sum(within_cap), "of", nrow(missing_cell_coords),
    "cells outside every polygon filled with the nearest country;",
    sum(!within_cap), "more than", fill_cap_km, "km away left NA\n")

gadm_country_levels <- levels(gadm_raster_country)[[1]]
gadm_region_levels <- levels(gadm_raster_region)[[1]]
missing_cell_index <- terra::cellFromXY(mask, as.matrix(missing_cell_coords))
gadm_raster_country[missing_cell_index[within_cap]] <- gadm_country_levels[
  match(nearest_country[within_cap], gadm_country_levels$country_name), 1]
gadm_raster_region[missing_cell_index[within_cap]] <- gadm_region_levels[
  match(nearest_region[within_cap], gadm_region_levels$region), 1]

gadm_raster_country <- terra::mask(gadm_raster_country, mask)
gadm_raster_region <- terra::mask(gadm_raster_region, mask)

writeRaster(gadm_raster_country,
            "data/clean/country_raster.tif",
            overwrite = TRUE)

writeRaster(gadm_raster_region,
            "data/clean/region_raster.tif",
            overwrite = TRUE)

# country borders for plotting are data/clean/country_borders.RDS, written by
# R/prep_country_borders.R from the same GADM layer
