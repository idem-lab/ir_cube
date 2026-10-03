# the subset of the collated bioassay data used in modelling, factored out of
# R/fit_model.R so that the model and the data summaries use the same filter.
# Needs R/packages.R and R/functions.R.

# the insecticide types modelled: those with at least 1000 unique places/times,
# and alpha-cypermethrin (914 unique) because of its use in LLINs (R/functions.R)
modelled_insecticides <- insecticides_plot_order

# numbers of unique location/time records per insecticide, and the mean number
# of mosquitoes tested per location/time
count_location_years <- function(ir_africa) {
  ir_africa %>%
    group_by(insecticide_type, latitude, longitude, year_start) %>%
    summarise(
      mosquito_number = sum(mosquito_number),
      .groups = "drop"
    ) %>%
    group_by(insecticide_type) %>%
    summarise(
      n = n(),
      mosquito_number = mean(mosquito_number),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

# subset to the modelled insecticides, the modal concentration of each, the
# modelled years, and cells inside the mask. Adds year_id (1-indexed from
# baseline_year) and the mask cell of each record.
subset_modelled_bioassays <- function(ir_africa,
                                      mask,
                                      insecticides_keep = modelled_insecticides,
                                      baseline_year = 1995,
                                      final_data_year = 2024) {
  ir_africa %>%
    filter(insecticide_type %in% insecticides_keep) %>%
    group_by(insecticide_type) %>%
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
      # create an index to the simulation year (in 1-indexed integers)
      year_id = year_start - baseline_year + 1,
      # add on cell ids corresponding to these observations,
      cell = terra::cellFromXY(mask,
                               as.matrix(select(., longitude, latitude)))
    ) %>%
    # drop a handful of datapoints missing covariates
    filter(
      !is.na(terra::extract(mask, cell)[, 1])
    )
}

# The indices of the modelled bioassays `df` (subset_modelled_bioassays()), as
# the fits and folds use them: the classes, types, regions, countries and
# unique_cells in order of first appearance, df with their *_id columns added,
# and classes_index, the class of each type
index_bioassays <- function(df) {
  out <- list(classes = unique(df$insecticide_class),
              types = unique(df$insecticide_type),
              regions = unique(df$region),
              countries = unique(df$country_name),
              unique_cells = unique(df$cell))
  out$df <- df %>%
    mutate(
      cell_id = match(cell, out$unique_cells),
      region_id = match(region, out$regions),
      country_id = match(country_name, out$countries),
      class_id = match(insecticide_class, out$classes),
      type_id = match(insecticide_type, out$types)
    )
  out$classes_index <- out$df %>%
    distinct(type_id, class_id) %>%
    arrange(type_id) %>%
    pull(class_id)
  out
}
