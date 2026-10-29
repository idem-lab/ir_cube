# visualise resistance data

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")

# visualise resistance data over time

# only plot pyrethroids
ir_africa %>%
  # fix up some names for plotting
  rename(
    mortality = mortality_adjusted,
    country = country_name,
    insecticide = insecticide_type
  ) %>%
  filter(
    insecticide_class %in% c("Pyrethroids")
  ) %>%
  # for each insecticide, filter to only the most common concentration
  group_by(
    insecticide
  ) %>%
  filter(
    concentration == sample_mode(concentration)
  ) %>%
  group_by(
    country
  ) %>%
  filter(
    n_distinct(year) >= 10,
    n() >= 100
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = year,
      y = mortality,
      group = country,
      colour = insecticide
    )
  ) +
  facet_wrap(~country,
             ncol = 3) +
  geom_point(alpha = 0.5) +
  theme_minimal()
