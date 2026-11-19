# Compute weights for an 'effective resistance to nets' output metric

# Analyse DHS/MIS survey data to approximate the proportion of single-AI LLINs
# that have each of the major pyrethroids (deltamethrin, permethrin,
# alpha-cypermethrin) as the active ingredient. These are the weights of
# llin_effective (R/predict.R, R/two_stage_maps.R): from the nets recorded in
# the household surveys of 2017-2022, deltamethrin 57%, permethrin 22% and
# alpha-cypermethrin 21%.

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the DHS/MIS data,
survey_net_data <- read_csv("data/raw/nets_by_ingredients_20240625.csv",
                            col_types = cols(
                              program = col_character(),
                              region = col_character(),
                              surveyname = col_character(),
                              country = col_character(),
                              year = col_double(),
                              brandid = col_double(),
                              brand = col_character(),
                              brand_reclass = col_character(),
                              net_type = col_character(),
                              active_ingredient = col_character(),
                              count = col_double()
                            ))

# select only those records where the nets unambiguously had a single active
# ingredient (AI). In some cases this was listed, in others, it could be
# determined by the net brand and the year (enabling discrimination between
# model versions), but for others it was ambiguous. Note computing the fraction
# from the unambiguous results is an imperfect approach. This approach assumes
# the information is 'missing at random', i.e. or those where it is not listed,
# it is equally like to follow these proportions
pyrethroids <- c("Deltamethrin", "Permethrin", "Alpha-cypermethrin")

# compute weights for the past 5 years
weights <- survey_net_data %>%
  filter(
    year >= 2022 - 5,
    active_ingredient %in% pyrethroids
  ) %>%
  group_by(
    active_ingredient
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    weight = count / sum(count)
  ) %>%
  select(
    -count
  ) %>%
  pivot_wider(
    names_from = active_ingredient,
    values_from = weight
  ) %>%
  as.list()

# stash these to use in prediction
saveRDS(weights,
        "temporary/ingredient_weights.RDS")

# print the percentages
perc <- round(unlist(weights) * 100)
paste(sprintf("%s: %i%s", names(perc), perc, "%"), collapse = ", ")


# Compute the fraction of treated nets in these surveys which are with unambiguous
# single AI vs unambiguous dual-AI vs ambiguous.
dual_AI <- c("Permethrin, Piperonyl Butoxide (PBO)",
             "Deltamethrin, Piperonyl Butoxide (PBO)",
             "Alpha-cypermethrin, Piperonyl Butoxide (PBO)",
             "Alpha-cypermethrin, Pyriproxyfen",
             "Alpha-cypermethrin, Chlorfenapyr")

survey_net_data %>%
  filter(
    year >= 2022 - 5
  ) %>%
  mutate(
    type = case_when(
      active_ingredient %in% pyrethroids ~ "single-AI",
      active_ingredient %in% dual_AI ~ "dual-AI",
      active_ingredient == "untreated" ~ "untreated",
      .default = "ambiguous"
    )
  ) %>%
  filter(
    type != "untreated"
  ) %>%
  group_by(
    type
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    percent = round(100 * count / sum(count))
  )
