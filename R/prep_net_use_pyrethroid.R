# Pyrethroid-only net use (#26): net use times the share of nets in use that
# select for pyrethroid resistance, for the selection design option
# nets = "pyrethroid" (selection_design() in R/model_covariates.R). PBO and
# dual-AI nets are assumed to add no selection, and a conventional ITN (cITN)
# to expose a fraction w of the vectors an LLIN does:
#
#   nets(x, t) = use(x, t) * (w cITN + LLIN) / total (a(x), t)
#
# with a(x) the admin 1 unit, and the share an annual mean over months,
# weighted by the net crop: the sum over months of w cITN + LLIN over the sum
# over months of all nets. w = 0.25 (main) or 0.47 (sensitivity), from
# R/net_type_weight.R.
#
# Net use from either of two sources (net_use_source in selection_design()):
#   run06   the MITN run 06 use layers (data/raw/itn/net_use_20260929), of the
#           same run as the net crop by type. The 2025 layer is omitted (it has
#           an error), so the model carries 2024 forward. They are missing
#           outside the 44 countries of the net crop (361,587 mask cells: North
#           Africa, South Africa, Lesotho, small islands, and 5 cells of
#           Equatorial Guinea, and up to 39 more cells in some years), and
#           are set to 0 there. The legacy layers are 0 there too but for 166
#           cell-years (150 in Egypt, 12 in Equatorial Guinea, 4 in
#           Guinea-Bissau)
#   legacy  the layers of net_use_cube.tif (prep_rasters.R), whose 2024 layer
#           is a copy of 2023. They differ from run 06's by up to 0.96 at a
#           cell (mean absolute difference 0.07-0.18 over the mask, by year)
#
# Writes, for each w, data/clean/net_use_pyrethroid_cube_w<w>.tif (run06) and
# net_use_pyrethroid_legacy_cube_w<w>.tif (legacy), layers nets_<year>,
# 2000-2024, processed as net_use_cube.tif is in prep_rasters.R;
# data/clean/net_use_run06_cube.tif, all run 06 net use processed the same
# way, and net_use_run06_filled.tif, 1 at the cells set to 0; and data/clean/net_pyrethroid_share_cube_w<w>.tif (layers
# share_<year>, the share at each mask cell). Run after prep_rasters.R.

source("R/packages.R")
source("R/functions.R")

net_type_w <- c(0.25, 0.47)
years <- 2000:2024

mask <- rast("data/clean/raster_mask.tif")


# net crop share by admin 1 unit and year -------------------------------------

# monthly net crop by admin 1 unit and type. admin1_id is the fdef_id of the
# admin 1 raster (all 582 units are in it); area_id matches none of its values
netcrop <- read_csv("data/raw/itn/net_use_20260929/netcrop_multitype_timeseries.csv",
                    show_col_types = FALSE) %>%
  filter(year %in% years)

netcrop_annual <- function(data, ...) {
  data %>%
    group_by(..., year) %>%
    summarise(citn = sum(cITN),
              llin = sum(LLIN),
              total = sum(`Total Nets`),
              .groups = "drop")
}
admin_crop <- netcrop_annual(netcrop, iso, admin1_id)
country_crop <- netcrop_annual(netcrop, iso)

# the share of nets that are pyrethroid-only, weighting cITNs by w. Where an
# admin 1 unit has no nets in a year (the share is then irrelevant where use
# is 0), it takes its country's share, and 1 where the country has none
pyrethroid_share <- function(citn, llin, total, w) {
  ifelse(total > 0, (w * citn + llin) / total, NA_real_)
}


# admin 1 unit of each mask cell ---------------------------------------------

admin1 <- terra::crop(rast("data/raw/admin/admin2023_1_MG_5K.tif"), mask)
admin0 <- terra::crop(rast("data/raw/admin/admin2023_0_MG_5K.tif"), mask)
stopifnot(compareGeom(admin1, mask), compareGeom(admin0, mask))
mask_cells <- which(!is.na(terra::values(mask, mat = FALSE)))
cell_admin <- tibble(cell = mask_cells,
                     admin1_id = terra::values(admin1, mat = FALSE)[mask_cells],
                     admin0_id = terra::values(admin0, mat = FALSE)[mask_cells])

# the country of each iso code in the net crop, as the admin 0 id most of the
# cells of its admin 1 units have
country_ids <- cell_admin %>%
  inner_join(distinct(admin_crop, iso, admin1_id), by = "admin1_id") %>%
  count(iso, admin0_id) %>%
  group_by(iso) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(iso, admin0_id)
stopifnot(!anyDuplicated(country_ids$admin0_id))

raster_admin1 <- setdiff(unique(cell_admin$admin1_id), 0)
cat(sprintf(paste("admin 1 units: %d in the net crop, %d in the raster within",
                  "the mask; %d of the net crop's are not in the raster, %d",
                  "of the raster's are not in the net crop\n"),
            n_distinct(admin_crop$admin1_id), length(raster_admin1),
            length(setdiff(unique(admin_crop$admin1_id), raster_admin1)),
            length(setdiff(raster_admin1, admin_crop$admin1_id))))

# each mask cell takes the share of its admin 1 unit if that is in the net
# crop; else its country's, if the country is; else that of the nearest admin
# 1 unit in the net crop
cell_admin <- cell_admin %>%
  mutate(source = case_when(
    admin1_id %in% admin_crop$admin1_id ~ "admin 1",
    admin0_id %in% country_ids$admin0_id ~ "country",
    .default = "nearest admin 1"
  )) %>%
  left_join(rename(country_ids, country_iso = iso), by = "admin0_id")

matched <- cell_admin$source == "admin 1"
unmatched <- cell_admin$source == "nearest admin 1"
admin_polygons <- terra::as.polygons(
  terra::mask(admin1, admin1 %in% unique(admin_crop$admin1_id),
              maskvalues = FALSE))
nearest <- terra::nearest(
  terra::vect(terra::xyFromCell(mask, cell_admin$cell[unmatched]),
              crs = terra::crs(mask)),
  admin_polygons)
cell_admin$nearest_id <- NA_real_
cell_admin$nearest_id[unmatched] <- terra::values(admin_polygons)[[1]][
  nearest$to_id]

print(count(cell_admin, source))
cat("mask cells with no admin 1 unit (0 or missing in the raster):",
    sum(is.na(cell_admin$admin1_id) | cell_admin$admin1_id == 0), "\n")
cat("countries of the cells without a unit or country in the net crop",
    "(admin 0 fdef_id: cells):\n")
print(cell_admin %>% filter(source == "nearest admin 1") %>% count(admin0_id) %>%
        arrange(desc(n)) %>% as.data.frame())


# share and pyrethroid-only net use cubes -------------------------------------

# as prep_rasters.R does for net use: crop, impute the misaligned coast, and
# mask
process_use <- function(files) {
  use <- terra::extend(terra::crop(rast(files), mask), mask)
  use <- terra::focal(use,
                      w = 9,
                      fun = "mean",
                      na.policy = "only",
                      na.rm = TRUE)
  use <- terra::mask(use, mask)
  names(use) <- paste0("nets_", years)
  use
}

# legacy: net_use_cube.tif, to float precision
use_legacy <- process_use(sprintf("data/raw/itn/net_use/ITN_%d_use_mean.tif",
                                  years))
net_use <- rast("data/clean/net_use_cube.tif")
stopifnot(identical(names(net_use), names(use_legacy)),
          all(terra::global(abs(net_use - use_legacy), "max",
                            na.rm = TRUE)[[1]] < 1e-6))

# run 06, 0 at the cells it lacks
use_run06 <- process_use(
  sprintf("data/raw/itn/net_use_20260929/use_%d_mean.tif", years))
run06_missing <- is.na(use_run06) & !is.na(use_legacy)
legacy_fill <- terra::values(use_legacy, mat = TRUE)[mask_cells, ]
missing <- terra::values(run06_missing, mat = TRUE)[mask_cells, ] == 1
cat("run 06 use missing at mask cells, by year:", colSums(missing), "\n")
cat("of which with legacy use > 0.01, by country (cell-years):\n")
country <- terra::values(rast("data/clean/country_raster.tif"),
                         dataframe = TRUE)[mask_cells, 1]
print(table(rep(as.character(country), length(years))[
  missing & legacy_fill > 0.01]))
use_run06 <- terra::cover(use_run06, use_legacy * 0)
terra::writeRaster(run06_missing, "data/clean/net_use_run06_filled.tif",
                   overwrite = TRUE)
stopifnot(!anyNA(terra::values(use_run06, mat = TRUE)[mask_cells, ]))
terra::writeRaster(use_run06, "data/clean/net_use_run06_cube.tif",
                   overwrite = TRUE)

uses <- list(run06 = use_run06, legacy = use_legacy)
cube_suffix <- c(run06 = "", legacy = "_legacy")

for (w in net_type_w) {
  country_share <- country_crop %>%
    mutate(country_share = pyrethroid_share(citn, llin, total, w)) %>%
    select(iso, year, country_share)
  admin_share <- admin_crop %>%
    mutate(share = pyrethroid_share(citn, llin, total, w)) %>%
    left_join(country_share, by = c("iso", "year")) %>%
    mutate(filled = is.na(share),
           share = coalesce(share, country_share, 1))
  if (w == net_type_w[1]) {
    cat("admin 1 unit-years with no nets:", sum(admin_share$filled),
        "; of which in countries with none:",
        sum(admin_share$filled & is.na(admin_share$country_share)), "\n")
  }
  country_share <- mutate(country_share,
                          country_share = coalesce(country_share, 1))

  # cells x years
  share_by <- function(ids, table, id_column, value_column) {
    wide <- table %>%
      select(all_of(c(id_column, "year", value_column))) %>%
      pivot_wider(names_from = year, values_from = all_of(value_column))
    as.matrix(wide[match(ids, wide[[id_column]]), as.character(years)])
  }
  share <- share_by(cell_admin$admin1_id, admin_share, "admin1_id", "share")
  from_country <- cell_admin$source == "country"
  share[from_country, ] <- share_by(cell_admin$country_iso[from_country],
                                    country_share, "iso", "country_share")
  share[unmatched, ] <- share_by(cell_admin$nearest_id[unmatched],
                                 admin_share, "admin1_id", "share")
  stopifnot(!anyNA(share), all(share >= 0 & share <= 1))

  share_cube <- rast(mask, nlyrs = length(years))
  share_cube[cell_admin$cell] <- share
  names(share_cube) <- paste0("share_", years)
  terra::writeRaster(share_cube,
                     sprintf("data/clean/net_pyrethroid_share_cube_w%.2f.tif",
                             w),
                     overwrite = TRUE)

  for (source in names(uses)) {
    net_use_pyrethroid <- uses[[source]] * share_cube
    names(net_use_pyrethroid) <- paste0("nets_", years)
    terra::writeRaster(net_use_pyrethroid,
                       sprintf("data/clean/net_use_pyrethroid%s_cube_w%.2f.tif",
                               cube_suffix[[source]], w),
                       overwrite = TRUE)
  }
}
