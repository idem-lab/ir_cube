# Definitions of the cross-validation folds: the data subsetting shared with
# the model fitting scripts, and the three out-of-sample experiments (spatial
# extrapolation by country, spatial interpolation between sampled locations,
# and temporal forecasting of the final years).
#
# Separated from the scoring so that the null models and the
# dynamical model are validated against exactly the same splits (#10). Sourcing
# this file defines `df`, `spatial_interpolation` and
# `temporal_forecasting_folds`, along with the indexing the model fitting needs.
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
source("R/bioassay_subset.R")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")


# load the mask
mask <- rast("data/clean/raster_mask.tif")

baseline_year <- 1995
final_data_year <- 2024

# the modelled subset, as in fit_model.R (R/bioassay_subset.R)
df <- subset_modelled_bioassays(ir_africa, mask, baseline_year = baseline_year,
                                final_data_year = final_data_year)


# indexing for main model fitting (R/bioassay_subset.R)
list2env(index_bioassays(df), environment())
years <- baseline_year - 1 + sort(unique(df$year_id))

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

# Leave-one-country-out is not defined here. It confounds spatial prediction
# with the country initial condition: a held-out country's init_country_raw has
# no data, reverts to the region prior, and that error is amplified through
# fifteen to twenty-nine years of deterministic selection, correlating with the
# fitted country effect at r = -0.94. The sub-national blocks in
# validation_blocks.R test spatial prediction without that confound.



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
  

# plot the training and test split, with circles and coloured points. Built
# inside the condition: polygonising the mask took a few seconds on every
# source() of this file, including from the fitting scripts, which never draw it
if (plot_folds) {

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

  print(interp_plot / interp_plot_small)

}

# `fold` already encodes the test year cutoff
spatial_interpolation <- list(
  training = df_interp %>%
    filter(fold == "training"),
  test = df_interp %>%
    filter(fold == "test")
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

# Build one forecasting fold from a cut year and a window length.
#
# `training` must be `year_start < cut_year`, not "not in the test years": the
# data run past the test window, and records after it in the training set make
# the forecast an interpolation, which favoured the dynamical model over the
# nearest neighbour null (#12 review).
#
# `before` is the window of equal length immediately before the cut. The
# forecasting experiment is scored on the change in mortality between that
# window and the holdout window, which differences the site level out and leaves
# the local slope — otherwise the test is largely spatial, since most held-out
# site-years have training data a few years earlier (#12 review 5.1). Those
# records are part of the training set; they are named here so the model's
# predictions at them can be requested at fitting time, the only time they can
# be.
forecasting_fold <- function(cut_year, window, data = df) {

  test_years <- cut_year + seq_len(window) - 1
  before_years <- cut_year - rev(seq_len(window))

  training <- data %>%
    filter(year_start < cut_year)

  fold <- list(
    cut_year = cut_year,
    window = window,
    test_years = test_years,
    before_years = before_years,
    training = training,
    test = data %>%
      filter(year_start %in% test_years),
    before = training %>%
      filter(year_start %in% before_years)
  )

  stopifnot(
    max(fold$training$year_start) < cut_year,
    all(fold$test$year_start %in% test_years),
    all(fold$before$year_start %in% before_years),
    nrow(fold$test) > 0,
    nrow(fold$before) > 0
  )

  fold
}

# The rolling origin. Five-year windows, cut at 2014 and 2018.
#
# The original design used a single three-year holdout starting at 2020, the
# latest the covariates allow. That is the worst available window on both
# counts: it is the one pause in the record — observed change at pixels assayed
# in both windows is decisively negative at every origin from 2005 to 2017 and
# flat at 2018-2020 — and it is the thinnest, 280 paired (pixel, insecticide)
# groups against 1,436 at a 2013 origin. Five-year windows give 30-50% more
# paired pixels and a signal roughly 5/3 larger, because the gap between window
# midpoints is the window length, and a five-year window spans the pause as well
# as the decline either side of it. The power analysis behind this, which measures
# this, and the PR #12 discussion.
#
# The two cuts train on 52% and 82% of the data and hold out windows whose true
# rates of decline differ by a factor of two, which is the test: does the model
# track a slowing rate, or carry a fixed slope forward? Earlier origins have a
# stronger signal still but train on too little data to be the model being
# deployed. Covariates end in 2022, so 2018 is the latest feasible five-year
# origin.
forecast_window <- 5
forecast_cuts <- c(2014, 2018)

temporal_forecasting_folds <- lapply(forecast_cuts, forecasting_fold,
                                     window = forecast_window)
names(temporal_forecasting_folds) <- as.character(forecast_cuts)

# The three-year 2020 fold this replaces is gone: it drew its training set with
# `year_start <= cut`, so the cut year itself appeared in both training and
# holdout, and its window landed on the one pause in twenty years of decline
# and was also the thinnest. Do not reinstate it.

# the validation experiments, and the only ones scored anywhere
validation_experiments <- c("spatial_interpolation", "spatial_blocks",
                            "temporal_forecasting")

# these are the train and test sets
spatial_interpolation
temporal_forecasting_folds

