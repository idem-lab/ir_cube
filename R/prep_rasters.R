# prepare rasters

# load packages and functions
source("R/packages.R")
source("R/functions.R")


# ITN cubes

# nets per capita
# from:
# GBD2024/Processing/Stages/13_ITN_Coverage_Africa/Checkpoint_Outputs/20240626/nets_per_capita
itn_pc_files <- list.files("data/raw/itn/nets_per_capita/",
                        pattern = ".tif",
                        full.names = TRUE)
itn_pc_years <- itn_pc_files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_percapita_nets_mean.tif")
nets_per_capita <- rast(itn_pc_files)
names(nets_per_capita) <- paste0("nets_", itn_pc_years)

terra::writeRaster(nets_per_capita,
                   "data/clean/nets_per_capita_cube.tif",
                   overwrite = TRUE)

# and scale it to 0-1
nets_per_capita_max <- max(global(nets_per_capita,
                                  "max",
                                  na.rm = TRUE))
nets_per_capita_scaled <- nets_per_capita / nets_per_capita_max
terra::writeRaster(nets_per_capita_scaled,
                   "data/clean/nets_per_capita_scaled_cube.tif",
                   overwrite = TRUE)

# make a mask for later use
mask <- nets_per_capita[[1]] * 0
names(mask) <- "mask" 
terra::writeRaster(mask,
                   "data/clean/raster_mask.tif",
                   overwrite = TRUE)


# net use
# from:
# GBD2024/Processing/Stages/13_ITN_Coverage_Africa/Checkpoint_Outputs/20240626/
itn_use_files <- list.files("data/raw/itn/net_use/",
                        pattern = ".tif",
                        full.names = TRUE)
itn_use_years <- itn_use_files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_use_mean.tif")
net_use <- rast(itn_use_files)
names(net_use) <- paste0("nets_", itn_use_years)

# crop these as the extent is wrong
net_use <- terra::crop(net_use, mask)

# the coast seems to be misaligned, so impute cells for prediction
net_use_filled <- terra::focal(net_use,
                               w = 9,
                               fun = "mean",
                               na.policy = "only",
                               na.rm = TRUE)
net_use <- terra::mask(net_use_filled, mask)

terra::writeRaster(net_use,
                   "data/clean/net_use_cube.tif",
                   overwrite = TRUE)

# and scale it to 0-1
net_use_max <- max(global(net_use,
                          "max",
                          na.rm = TRUE))
net_use_scaled <- net_use / net_use_max
terra::writeRaster(net_use_scaled,
                   "data/clean/net_use_scaled_cube.tif",
                   overwrite = TRUE)


# IRS cube
# these are admin-level IRS coverage copied from
# GBD2024/Processing/Stages/14_IRS_Coverage_Africa/Intermediary_Outputs/20240620
irs_files <- list.files("data/raw/irs/admin_coverage",
                        pattern = ".tif",
                        full.names = TRUE)
irs_years <- irs_files %>%
  basename() %>%
  str_remove("admin_coverage_") %>%
  str_remove(".tif")
irs_coverage <- rast(irs_files)

# extend these to the ITN mask
irs_coverage <- terra::extend(irs_coverage, mask)

# pad them with zeros (where distribution data missing, including -9999 coding)
irs_coverage[irs_coverage < 0] <- 0
irs_coverage[is.na(irs_coverage)] <- 0

# mask them
irs_coverage <- terra::mask(irs_coverage, mask)

# set their names
names(irs_coverage) <- paste0("irs_", irs_years)

terra::writeRaster(irs_coverage,
                   "data/clean/irs_coverage_cube.tif",
                   overwrite = TRUE)

# scale this
irs_coverage_max <- max(global(irs_coverage,
                               "max",
                               na.rm = TRUE))
irs_coverage_scaled <- irs_coverage / irs_coverage_max
terra::writeRaster(irs_coverage_scaled,
                   "data/clean/irs_coverage_scaled_cube.tif",
                   overwrite = TRUE)

# Human population, combining MAP's WorldPop observed and projected global
# rasters to get 2000-2024

# past population, 2000-2020
pop_past_files <- list.files("data/raw/pop/",
                             pattern = ".tif",
                             full.names = TRUE)
pop_past_years <- pop_past_files %>%
  basename() %>%
  str_remove("WorldPop_UNAdj_v3_DRC_fix.") %>%
  str_remove(".Annual.Data.5km.sum.tif")
pop_past <- rast(pop_past_files)
names(pop_past) <- paste0("pop_", pop_past_years)

# projected future population, 2021-2049
pop_future_files <- list.files("data/raw/pop_future/",
                               pattern = ".tif",
                               full.names = TRUE)
pop_future_years <- pop_future_files %>%
  basename() %>%
  str_remove("WorldPop_UNAdj_v3_DRC_fix_projected.") %>%
  str_remove(".Annual.Data.5km.sum.tif")
pop_future <- rast(pop_future_files)
names(pop_future) <- paste0("pop_", pop_future_years)

# combine, crop, and mask to mastergrids
pop_all <- c(pop_past, pop_future)
pop_all <- crop(pop_all, mask)
pop_all <- terra::mask(pop_all, mask)

# fill in zeros on land with an epsilon, fill in NAs (outside mastergrids) with
# 0, and relevel total populations, just in case I end up using these for population
# stats by mistake.

# Set the epsilon to 1 person per 5km, approximately the 11th percentile of
# non-zero populations, and for Egypt (with most missing population) this
# matches the neighbouring low-population Sahara desert areas well

# annual populations for Africa within the mastergrids file, for relevelling
target_africa_pops <- global(pop_all, "sum", na.rm = TRUE)

# check approximate percentile of the epsilon
epsilon <- 1
epsilon_quantile <- mean(na.omit(c(pop_all[])) <= epsilon)
round(100 * epsilon_quantile)

# fill in the missing populations (and remask the true NA areas)
pop_all[is.na(pop_all)] <- 0
pop_all[pop_all == 0] <- epsilon
pop_all <- terra::mask(pop_all, mask)

# relevel the total continent populations
new_africa_pops <- global(pop_all, "sum", na.rm = TRUE)
ratios <- target_africa_pops / new_africa_pops
for (i in seq_len(terra::nlyr(pop_all))) {
  print(i)
  pop_all[[i]] <- pop_all[[i]] * ratios[[1]][i]
}

# check there's less than a person different in the total populations
final_africa_pops <- global(pop_all, "sum", na.rm = TRUE)
max(abs(final_africa_pops - target_africa_pops)) < 1

# transform and scale to 0-1
# trans_pop_all <- log(pop_all)
trans_pop_all <- identity(pop_all)
trans_pop_all_min <- min(global(trans_pop_all,
                            "min",
                            na.rm = TRUE))
trans_pop_all <- trans_pop_all - trans_pop_all_min
trans_pop_all_max <- max(global(trans_pop_all,
                            "max",
                            na.rm = TRUE))
trans_pop_all <- trans_pop_all / trans_pop_all_max

# clip to modelling years, and save
pop <- pop_all[[paste0("pop_", 2000:2024)]]
trans_pop <- trans_pop_all[[paste0("pop_", 2000:2024)]]

terra::writeRaster(pop,
                   "data/clean/pop_cube.tif",
                   overwrite = TRUE)

terra::writeRaster(trans_pop,
                   "data/clean/pop_scaled_cube.tif",
                   overwrite = TRUE)

# save future projection years, too
pop_future <- pop_all[[paste0("pop_", 2025:2030)]]
trans_pop_future <- trans_pop_all[[paste0("pop_", 2025:2030)]]

terra::writeRaster(pop_future,
                   "data/clean/pop_cube_future.tif",
                   overwrite = TRUE)

terra::writeRaster(trans_pop_future,
                   "data/clean/pop_scaled_cube_future.tif",
                   overwrite = TRUE)

# log population, min-max scaled over the same years, for the selection design
# option pop = "log" (#23, selection_design() in R/model_covariates.R)
log_pop_all <- log(pop_all)
log_pop_range <- range(unlist(global(log_pop_all, "range", na.rm = TRUE)))
log_pop_all <- (log_pop_all - log_pop_range[1]) / diff(log_pop_range)
names(log_pop_all) <- str_replace(names(pop_all), "^pop_", "log_pop_")

terra::writeRaster(log_pop_all[[paste0("log_pop_", 2000:2024)]],
                   "data/clean/pop_log_scaled_cube.tif",
                   overwrite = TRUE)

terra::writeRaster(log_pop_all[[paste0("log_pop_", 2025:2030)]],
                   "data/clean/pop_log_scaled_cube_future.tif",
                   overwrite = TRUE)

# other static ones

# NB I made a symlink from Gerry's dropbox into my into data raw with the
# following mac terminal command (from inside data/raw):
# ln -s ~/Dropbox/sharing/ir_pred_rasters ./

# load crops, both as-is and standardised to 0-1
crops <- rast("data/raw/ir_pred_rasters/crop_2010.tif")


# create layers for total yield of crop groups
cereals <- c("wheat", "rice", "maize", "barley", "pearl millet", "small millet",
             "sorghum", "other cereals")

roots <- c("potato", "sweet potato", "yams", "cassava", "other roots")

pulses <- c("bean", "chickpea",  "cowpea", "pigeonpea", "lentil",
            "other pulses")

oilcrops <- c("soybean", "groundnut", "coconut", "oilpalm", "sunflower",
              "rapeseed", "sesameseed", "other oil crops")

fibrecrops <- c("sugarcane", "sugarbeet", "cotton", "other fibre crops")

othercrops <- c("arabica coffee", "robusta coffee", "cocoa", "tea", "tobacco",
                "banana", "plantain", "tropical fruit", "temperate fruit",
                "vegetables", "rest of crops")

crop_info <- geodata::spamCrops() %>%
  as_tibble() %>%
  mutate(
    code = toupper(code),
    type = case_when(
      crop %in% cereals ~ "cereal crops",
      crop %in% roots ~ "root crops",
      crop %in% pulses ~ "pulse crops",
      crop %in% oilcrops ~ "oil crops",
      crop %in% fibrecrops ~ "fibre crops",
      crop %in% othercrops ~ "other crops"
    )
  )

# compute and standardise the yield of all crops
crops_group <- app(crops, sum)
crops_group <- crops_group / global(crops_group, "max", na.rm = TRUE)[1, 1]
names(crops_group) <- "all crops"

# do the same for groups of these and append them
for (crop_type in unique(crop_info$type)) {
  
  crop_codes_type <- crop_info %>%
    filter(
      type == crop_type
    ) %>%
    pull(code)
  
  crops_type <- app(crops[[crop_codes_type]], sum)
  crops_type <- crops_type / global(crops_type, "max", na.rm = TRUE)[1, 1]
  names(crops_type) <- crop_type
  crops_group <- c(crops_group, crops_type)
  
}

# save these out
terra::writeRaster(crops_group,
                   "data/clean/crop_group_scaled.tif",
                   overwrite = TRUE)


# now rename the individual crops, scale, and save layer for all crops
full_names <- names(crops) %>%
  tibble(code = .) %>%
  left_join(crop_info,
            by = "code") %>% 
  pull(crop)

names(crops) <- full_names

crops_scaled <- crops / global(crops, "max", na.rm = TRUE)$max

terra::writeRaster(crops_scaled,
                   "data/clean/crop_scaled.tif",
                   overwrite = TRUE)

# download and cache (against the risk of the package breaking) MAP's PfPR for
# 2000, and use to to get the limits of transmission for plotting
pr2000_file <- "data/clean/pfpr_2000.tif"

if (!file.exists(pr2000_file)) {
  pr2000_all <- malariaAtlas::getRaster("Malaria__202508_Global_Pf_Parasite_Rate",
                                        year = 2000)
  mask <- rast("data/clean/raster_mask.tif")
  pr2000 <- crop(pr2000_all[[1]], mask)
  pr2000[is.na(pr2000)] <- 0
  pr2000 <- terra::mask(pr2000, mask)
  names(pr2000) <- "PfPR_2000"
  writeRaster(pr2000, pr2000_file, overwrite = TRUE)
}


pf_limits_file <- "data/clean/pfpr_limits.tif"
if (!file.exists(pf_limits_file)) {
  pf_limits <- malariaAtlas::getRaster("Explorer__2010_Pf_Limits_Decompressed")
  pf_limits <- crop(pf_limits, llin_rast)
  pf_limits <- terra::aggregate(pf_limits, 5, fun = "modal")
  pf_limits <- terra::mask(pf_limits, llin_rast)
  writeRaster(pf_limits, pf_limits_file, overwrite = TRUE)
}

# reload the main mask
mask <- rast("data/clean/raster_mask.tif")

# load lakes and rivers layer
water <- rast("data/raw/WorldCover_water_30s.tif")
water <- crop(water, mask)
water <- aggregate(water, 5, fun = "modal")

# convert it into a mask
water_mask <- water
water_mask[is.na(water_mask)] <- 0
water_mask[water_mask == 1] <- NA
water_mask <- terra::mask(water_mask, mask)
water_mask <- water_mask * 0

# combine with the pf limits to amke a mask of limits and water
pf_limits <- rast(pf_limits_file)
pf_limits_mask <- pf_limits
pf_limits_mask[pf_limits_mask != 2] <- NA
pf_limits_mask <- pf_limits_mask * 0
pf_limits_water_mask <- terra::mask(pf_limits_mask, water_mask)

writeRaster(pf_limits_water_mask,
            "data/clean/pfpr_water_mask.tif",
            overwrite = TRUE)

