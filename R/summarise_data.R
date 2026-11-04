# summarise the bioassay data: numbers of records, mosquitoes tested (total and
# median per bioassay), unique locations, location-years and countries, and
# year range, for the full collated dataset and for the modelled subset,
# overall, by insecticide class and by insecticide type; and the records each
# processing step removes, in R/prep_bioassays.R (per source) and in
# subset_modelled_bioassays()

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")

# the modelled subset, as in R/fit_model.R
df <- subset_modelled_bioassays(ir_africa, mask)

# summary statistics for a (grouped) set of bioassay records
summarise_bioassays <- function(data) {
  data %>%
    summarise(
      records = n(),
      mosquitoes = sum(mosquito_number),
      median_tested = median(mosquito_number),
      locations = n_distinct(latitude, longitude),
      location_years = n_distinct(latitude, longitude, year),
      countries = n_distinct(country_name),
      year_min = min(year),
      year_max = max(year),
      .groups = "drop"
    )
}

# overall, by class, and by type, with the type order used in the figures
summarise_levels <- function(data) {
  bind_rows(
    data %>%
      summarise_bioassays() %>%
      mutate(level = "overall",
             insecticide_class = "All",
             insecticide_type = "All"),
    data %>%
      group_by(insecticide_class) %>%
      summarise_bioassays() %>%
      mutate(level = "class",
             insecticide_type = "All"),
    data %>%
      group_by(insecticide_class, insecticide_type) %>%
      summarise_bioassays() %>%
      mutate(level = "type")
  ) %>%
    mutate(
      level = factor(level, levels = c("overall", "class", "type"))
    ) %>%
    arrange(level, desc(records)) %>%
    relocate(level, insecticide_class, insecticide_type)
}

data_summary <- bind_rows(
  collated = summarise_levels(ir_africa),
  modelled = summarise_levels(df),
  .id = "dataset"
)

write_csv(data_summary, "outputs/data_summary.csv")

# location_years for each type in the collated data are the counts that set the
# inclusion rule for insecticide types (R/fit_model.R)
record_counts <- count_location_years(ir_africa)
stopifnot(
  all.equal(
    record_counts$n,
    data_summary %>%
      filter(dataset == "collated", level == "type") %>%
      pull(location_years, name = insecticide_type) %>%
      `[`(record_counts$insecticide_type),
    check.attributes = FALSE
  )
)

print(data_summary, n = Inf, width = Inf)

# records removed at each step of R/prep_bioassays.R, summed over sources
attrition <- read_csv("outputs/review/bioassay_attrition.csv",
                      show_col_types = FALSE, name_repair = "minimal")
cat("\nrecords removed at each step of R/prep_bioassays.R, all sources:\n")
print(enframe(colSums(select(attrition, -source)), name = "step",
              value = "records"), n = Inf)

# and at each step of subset_modelled_bioassays(), in its order; the last,
# the cells outside the mask, is the difference from its result
modelled_types <- filter(ir_africa, insecticide_type %in% modelled_insecticides)
modal_concentration <- modelled_types %>%
  group_by(insecticide_type) %>%
  filter(concentration == sample_mode(concentration)) %>%
  ungroup()
from_baseline <- filter(modal_concentration, year >= 1995)
subset_steps <- tibble(
  step = c("collated", "not a modelled insecticide",
           "not the insecticide's modal concentration", "before 1995",
           "outside the mask"),
  kept = c(nrow(ir_africa), nrow(modelled_types), nrow(modal_concentration),
           nrow(from_baseline), nrow(df))
) %>%
  mutate(removed = lag(kept) - kept)
stopifnot(all(subset_steps$removed >= 0, na.rm = TRUE))
cat("\nrecords removed by subset_modelled_bioassays():\n")
print(subset_steps)
