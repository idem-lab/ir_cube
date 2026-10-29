# Are the replicate bioassay groups representative of the data as a whole?
#
# Bioassays carried out in the same 5km pixel, in the same year, with the same
# insecticide, are replicate measurements of a single population susceptible
# fraction, so the variation between them identifies the observation
# overdispersion without reference to any spatial model. That makes the
# resulting rho an external quantity: it can be compared with the rho estimated
# inside the dynamical model, and it sets the noise floor against which
# out-of-sample predictive performance is judged (see idem-lab/ir_cube#10).
#
# rho itself is estimated in fig_illustrate_bioassay_variability.R and read back
# by rho_lookup(). This script asks only whether the groups it rests on look
# like the rest of the data. A quadrature-based marginal likelihood for rho once
# sat here as a second estimator; nothing ever called it, and it is gone (#12
# review).

source("R/validation_functions.R")

suppressMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

mask <- rast("data/clean/raster_mask.tif")
ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")

baseline_year <- 1995

insecticides_keep <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin",
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")

# same subsetting as the model fitting and validation scripts
sample_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

df <- ir_africa %>%
  filter(insecticide_type %in% insecticides_keep) %>%
  group_by(insecticide_type) %>%
  filter(concentration == sample_mode(concentration)) %>%
  ungroup() %>%
  filter(
    year >= baseline_year,
    !is.na(died),
    !is.na(mosquito_number),
    mosquito_number > 1
  ) %>%
  mutate(
    cell = cellFromXY(mask, as.matrix(select(., longitude, latitude)))
  ) %>%
  filter(!is.na(terra::extract(mask, cell)[, 1]))

# identify replicate groups: one pixel, one year, one insecticide
df <- df %>%
  group_by(cell, year, insecticide_type) %>%
  mutate(
    group_id = cur_group_id(),
    group_size = n()
  ) %>%
  ungroup()

replicated <- df %>%
  filter(group_size >= 2)

cat(sprintf("%i assays in %i replicated pixel-year-insecticide groups\n",
            nrow(replicated),
            n_distinct(replicated$group_id)))



# are replicated pixel-years representative? -------------------------------

# The floor rests on these groups, but repeatedly sampled sites are plausibly
# sentinel sites, and may differ systematically from the rest of the data. This
# compares them on the quantities that matter for the floor
representativeness <- df %>%
  mutate(
    replicated = ifelse(group_size >= 2, "replicated", "single")
  ) %>%
  group_by(replicated) %>%
  summarise(
    assays = n(),
    pixel_years = n_distinct(group_id),
    mean_mortality = mean(died / mosquito_number),
    sd_mortality = sd(died / mosquito_number),
    median_mosquito_number = median(mosquito_number),
    median_year = median(year),
    countries = n_distinct(country_name),
    proportion_pyrethroid = mean(insecticide_class == "Pyrethroids"),
    .groups = "drop"
  )

cat("\nreplicated versus singly sampled pixel-years:\n")
print(representativeness %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

write.csv(representativeness,
          "outputs/bioassay_rho_representativeness.csv",
          row.names = FALSE)

# the same comparison by country, to show whether the replicated groups come
# from a narrow set of places
by_country <- df %>%
  mutate(replicated = group_size >= 2) %>%
  group_by(country_name) %>%
  summarise(
    assays = n(),
    proportion_replicated = mean(replicated),
    .groups = "drop"
  ) %>%
  arrange(desc(assays))

cat("\ntop countries by number of assays:\n")
print(head(by_country, 10))

representativeness_plot <- df %>%
  mutate(
    replicated = ifelse(group_size >= 2,
                        "replicated pixel-year",
                        "single assay")
  ) %>%
  ggplot(
    aes(x = died / mosquito_number,
        fill = replicated)
  ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.6,
    position = "identity"
  ) +
  facet_wrap(~ insecticide_class, scales = "free_y") +
  scale_fill_manual(values = c("replicated pixel-year" = "#B2182B",
                               "single assay" = grey(0.5)),
                    name = "") +
  labs(
    x = "observed mortality",
    y = "density",
    title = "Are repeatedly sampled pixel-years representative?",
    subtitle = "the overdispersion estimate, and so the noise floor, rests on the replicated groups"
  ) +
  theme_minimal()

ggsave("figures/bioassay_rho_representativeness.png",
       representativeness_plot,
       bg = "white",
       width = 10,
       height = 6)
