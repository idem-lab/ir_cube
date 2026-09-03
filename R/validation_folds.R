# Definitions of the cross-validation folds: the data subsetting shared with
# the model fitting scripts, and the three out-of-sample experiments (spatial
# extrapolation by country, spatial interpolation between sampled locations,
# and temporal forecasting of the final years).
#
# Extracted from predictive_validation.R so that the null models and the
# dynamical model are validated against exactly the same splits (#10). Sourcing
# this file defines `df`, `spatial_extrapolation`, `spatial_interpolation` and
# `temporal_forecasting`, along with the indexing the model fitting needs.
#
# The interpolation diagnostic plot is drawn only when `plot_folds` is TRUE, so
# that sourcing this file from a fitting script is silent.

if (!exists("plot_folds")) {
  plot_folds <- FALSE
}

# Run out-of-sample predictive validation experiments

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")


# Note: there is code in fit_model to subset this to some insecticides. Relevant
# code is copied here for now, but move that into the data preparation scripts
# and save only the model-ready version to load in here

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# an Africa polygon for plotting
gadm_polys <- readRDS("data/clean/gadm_polys.RDS")
africa <- gadm_polys %>%
  # st_combine() %>%
  st_union()

baseline_year <- 1995
final_data_year <- 2024

insecticides_keep <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin", 
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")

df <- ir_africa %>%
  filter(
    insecticide_type %in% insecticides_keep
  ) %>%
  group_by(
    insecticide_type
  ) %>%
  # subset to the most common concentration for each insecticide
  filter(
    concentration == sample_mode(concentration)
  ) %>%
  ungroup() %>%
  filter(
    # drop any from before the baseline
    year_start >= baseline_year,
    year_start <= final_data_year
  ) %>%
  mutate(
    # need year id for main dynamical model
    year_id = year_start - baseline_year + 1,
    # add on cell ids corresponding to these observations,
    cell = cellFromXY(mask,
                      as.matrix(select(., longitude, latitude)))
  ) %>%
  # drop a handful of datapoints missing covariates
  filter(
    !is.na(extract(mask, cell)[, 1])
  )


# indexing for main model fitting
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)
regions <- unique(df$region)
countries <- unique(df$country_name)
unique_cells <- unique(df$cell)
years <- baseline_year - 1 + sort(unique(df$year_id))

df <- df %>%
  mutate(
    cell_id = match(cell, unique_cells),
    region_id = match(region, regions),
    country_id = match(country_name, countries),
    class_id = match(insecticide_class, classes),
    type_id = match(insecticide_type, types)
  )

# Define the training and test folds for: spatial extrapolation (country
# dropout), spatial interpolation (multi-area dropout), and temporal forecasting
# (last years dropout)

# Subset the test set to only the recent period, to ensure similarity of
# resistance levels with the present day, whilst retaining a reasonable amount
# of data for testing.
test_min_year <- 2010

# Find the countries with enough data for model validation (at least 10
# bioassays for each insecticide type) and make a table of counts
country_bioassay_counts <- df %>%
  filter(
    year_start >= test_min_year
  ) %>%
  group_by(
    country_name,
    insecticide_type) %>%
  summarise(
    records = n(),
    .groups = "drop"
  ) %>%
  # find country/insecticide combinations with a reasonable number of records in recent years
  filter(
    records > 10
  ) %>%
  # find countries with enough records of all insecticides
  group_by(
    country_name
  ) %>%
  mutate(
    n_insecticides = n()
  ) %>%
  filter(
    n_insecticides == 9
  ) %>%
  pivot_wider(
    names_from = country_name,
    values_from = records
  ) %>%
  select(
    -n_insecticides
  )

# view(country_bioassay_counts)

# pull out the country names
countries_to_validate <- country_bioassay_counts %>%
  select(-insecticide_type) %>%
  colnames()

countries_to_validate

# make a list of training and test datasets for out-of-sample validation

# Subset a dataset based on values in a field. If keep == TRUE (the default),
# return the subset of the dataset where field_name (a character for the column
# name) is one of the elements in field_values (a character vector), if keep =
# FALSE, return everything except for those records. When keep = TRUE, records
# are only kept if they are after 'keep_min_year' - to optionally limit test
# sets to recent years
split_data <- function(field_values,
                       field_name,
                       dataset,
                       keep = TRUE,
                       keep_min_year = 1900) {
  
  if (keep) {
    
    subsetted <- dataset %>%
      filter_at(
        vars(starts_with(field_name)),
        any_vars(. %in% field_values)
      ) %>%
      filter(
        year_start >= keep_min_year
      )
    
  } else {
    
    subsetted <- dataset %>%
      filter_at(
        vars(starts_with(field_name)),
        any_vars(!(. %in% field_values))
      )
    
  }
  
  subsetted
  
}

spatial_extrapolation <- list(
  training = lapply(
    countries_to_validate,
    split_data,
    field_name = "country_name",
    dataset = df,
    keep = FALSE
  ),
  test = lapply(
    countries_to_validate,
    split_data,
    field_name = "country_name",
    dataset = df,
    keep = TRUE,
    keep_min_year = test_min_year
  )
)
names(spatial_extrapolation$training) <- countries_to_validate
names(spatial_extrapolation$test) <- countries_to_validate



# Spatial interpolation

# do k-means clustering to identify region centroids
df_locations <- df %>%
  select(
    longitude,
    latitude
  ) %>%
  as.matrix()

# set the RNG seed, as this is stochastic
set.seed(2025-07-07)
centroids <- df_locations %>%
  unique() %>%
  kmeans(
    centers = 50,
    iter.max = 300,
    nstart = 100
  ) %>%
  `[[`(
    "centers"
  )


if (plot_folds) {
  plot(df_locations,
       asp = 1,
       pch = 16,
       cex = 0.5,
       col = grey(0.8))

  points(centroids,
         pch = 16,
         cex = 0.5,
         col = "red")
}

# identify all points within some distance of these as, as being test data
dists <- fields::rdist.earth(
  x1 = centroids,
  x2 = df_locations,
  miles = FALSE
)

# for each datapoint, get the minimum distance to a centroid
min_dists <- apply(dists, 2, min)

# get the distance in km such that 5% of records are in the test set
test_threshold <- quantile(min_dists, 0.05)

# set half this to be the buffer distance, and define the outer edge of the
# buffer circle
buffer_distance <- 0.5 * test_threshold
buffer_threshold <- test_threshold + buffer_distance

# split into test (in threshold and on or after min test year), training
# (outside buffer) and excluded points (everything else)
df_interp <- df %>%
  mutate(
    fold = case_when(
      min_dists < test_threshold &
        year_start >= test_min_year ~ "test",
      min_dists > buffer_threshold ~ "training",
      .default = "excluded"
    ),
    fold = factor(fold, levels = c("excluded", "training", "test"))
  )

# check the training and test splits contain all the insecticides
df_interp %>%
  filter(
    fold != "excluded"
  ) %>%
  group_by(
    fold,
    insecticide_type
  ) %>%
  summarise(
    records = n()
  ) %>%
  pivot_wider(
    names_from = fold,
              values_from  = records
  )
  

# plot the training and test split, with circles and coloured points.
mask_poly <- mask %>%
  as.polygons() %>%
  simplifyGeom(tolerance = 0.05)

train_test_col <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

interp_plot <- df_interp %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = crs(mask)
  ) %>%
  arrange(fold) %>%
  ggplot(
    aes(
      colour = fold
    )
  ) +
  geom_spatvector(
    data = mask_poly,
    colour = "transparent",
    fill = grey(0.9)
  ) +
  # geom_spatraster(
  #   data = mask
  # ) +
  # scale_fill_gradient(
  #   low = grey(0.95),
  #   high = grey(0.95),
  #   na.value = "transparent",
  #   guide = "none"
  # ) +
  geom_sf() +
  scale_colour_manual(
    values = c(
      "training" = train_test_col[2],
      "test" = train_test_col[1],
      "excluded" = grey(0.8)
    )
  ) +
  coord_sf(
    ylim = range(df$latitude)
  ) +
  theme_minimal()

interp_plot_small <- interp_plot +
  coord_sf(
    xlim = c(0, 5),
    ylim = c(5, 10)
  )

if (plot_folds) {
  print(interp_plot / interp_plot_small)
}

spatial_interpolation <- list(
  training = split_data(
    field_values = "training",
    field_name = "fold",
    dataset = df_interp,
    keep = TRUE
  ),
  test = split_data(
    field_values = "test",
    field_name = "fold",
    dataset = df_interp,
    keep = TRUE
  )
)


# Temporal forecasting

# work out then the latest covariate layers are to get the three test years

nets_cube <- rast("data/clean/net_use_cube.tif")
irs_cube <- rast("data/clean/irs_coverage_scaled_cube.tif")
pop_cube <- rast("data/clean/pop_scaled_cube.tif")
nets_final_year <- nets_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("nets_") %>%
  as.numeric()
irs_final_year <- irs_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("irs_") %>%
  as.numeric()
pop_final_year <- pop_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("pop_") %>%
  as.numeric()
final_year <- min(nets_final_year, irs_final_year, pop_final_year)
validation_years <- final_year + -2:0

# split into a training set (before those years) and a test set
temporal_forecasting <- list(
  training = split_data(
    field_values = validation_years,
    field_name = "year_start",
    dataset = df,
    keep = FALSE
  ),
  test = split_data(
    field_values = validation_years,
    field_name = "year_start",
    dataset = df,
    keep = TRUE
  )
)

# these are the train and test sets
spatial_extrapolation
spatial_interpolation
temporal_forecasting

