# Make figure of model predictions over time for all of Africa overlaid with
# point-estimates of susceptibility bioassay data, aggregated up at various
# levels.

# Compute the point estimates as *weighted average* susceptibilities (weighted
# sum of number that died over weighted sum of number tested) to represent
# continent-wide sampling for that insecticide. Replace legend title with
# 'Effective sample size'

# Repeat the plots, with multiple lines and points for each regions. (Make
# functions to simplify this?)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/two_stage_predictions.R")

# the modelled data, as R/fit_model.R builds them
baseline_year <- 1995
invisible(list2env(modelled_bioassays(baseline_year), environment()))
final_data_year <- max(df$year)
years_predict <- baseline_year:final_data_year

# set colours
pyrethroid_blue <- "#56B1F7"

# get all the pyrethroids
pyrethroids <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  filter(class == "Pyrethroids") %>%
  pull(insecticide)

# list the pyrethroids used in LLINs (ie. not Lambda-cyhalothrin), these are the
# products in all nets recorded in surveys
llin_pyrethroids <- c("Alpha-cypermethrin",
                      "Deltamethrin",
                      "Permethrin")

# tidy up column names
df_sub <- df %>%
  rename(
    country = country_name,
    insecticide = insecticide_type
  )
  
# subset to pyrethroids, and add on UN geoscheme regions for Africa
df_pyrethroids <- df_sub %>%
  filter(
    insecticide %in% llin_pyrethroids
    # insecticide_class == "Pyrethroids",
  )

# compute predicted population-level susceptibility, weighted over the numbers
# of samples per location, per insecticide type: the two-stage model's
# posterior (R/two_stage_predictions.R), averaged over the locations with any
# bioassays for that insecticide, Africa-wide and per region, for all years

# calculate the fraction of all samples (for each insecticide) that come from
# each cell, and use this to compute a weighted average of predicted
# susceptibility

weights_mat <- df_sub %>%
  group_by(cell_id, type_id) %>%
  # get annual sums 
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  complete(
    cell_id,
    type_id,
    fill = list(mosquito_number = 0)
  ) %>%
  group_by(type_id) %>%
  mutate(weight = mosquito_number / sum(mosquito_number)) %>%
  select(-mosquito_number) %>%
  arrange(cell_id, type_id) %>%
  pivot_wider(
    names_from = type_id,
    values_from = weight
  ) %>%
  select(-cell_id) %>%
  as.matrix() %>%
  `colnames<-`(NULL)

# the region of each cell (that of its first record), and its country in the
# dynamical model
cell_region <- df_sub %>%
  group_by(cell_id) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(cell_id) %>%
  pull(region)
cell_country <- data_cell_country(df)

# for each insecticide, the weights renormalised within each region (all zero
# where a region has no data for it), and draws x (Africa, regions) x years of
# the weighted average susceptibility
regions <- unique(df_sub$region)
type_draws <- lapply(seq_along(types), function(k) {
  regional <- sapply(regions, function(region) {
    w <- weights_mat[, k] * (cell_region == region)
    if (sum(w) > 0) w / sum(w) else w
  })
  weights <- cbind(Africa = weights_mat[, k], regional)
  setup <- two_stage_setup(types[k], years_predict, df)
  two_stage_weighted_draws(setup, unique_cells, cell_country,
                           weights)$two_stage
})

# posterior mean and 95% credible interval of draws x years
summarise_draws <- function(x) {
  tibble(
    year = years_predict,
    susc_pop_mean = colMeans(x),
    susc_pop_lower = apply(x, 2, quantile, 0.025),
    susc_pop_upper = apply(x, 2, quantile, 0.975)
  )
}

# now compute a weighted sum over the pyrethroids (weights given by the numbers
# of mosquitoes tested for each pyrethroid, over the whole timeseries)
pyrethroid_weights <- df_sub %>% 
  group_by(
    type_id,
    insecticide 
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    is_a_pyrethroid = insecticide %in% llin_pyrethroids,
    mosquito_number = mosquito_number * as.numeric(is_a_pyrethroid),
    weight = mosquito_number / sum(mosquito_number)
  ) %>%
  arrange(type_id) %>%
  pull(weight)

# posterior summaries of pyrethroid susceptibility by year; the types' draws
# share their dynamical draws, so are combined draw by draw
pyrethroid_susc <- Reduce(`+`, lapply(which(pyrethroid_weights > 0), function(k) {
  pyrethroid_weights[k] * type_draws[[k]][, "Africa", ]
}))
pop_mort_sry_pyrethroid <- summarise_draws(pyrethroid_susc)

insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  pull(insecticide)

insecticide_type_labels <- sprintf("%s) %s%s",
                                   LETTERS[1 + seq_along(insecticides_plot)],
                                   insecticides_plot,
                                   ifelse(insecticides_plot %in% llin_pyrethroids,
                                          "*",
                                          ""))

insecticides_plot_lookup <- tibble(
  insecticide = insecticides_plot,
  insecticide_type_label = insecticide_type_labels
)

insecticides_plot_lookup <- tibble(
  insecticide = insecticides_plot
) %>%
  mutate(
    idx = row_number(),
    suffix = case_when(
      insecticides_plot %in% llin_pyrethroids ~ "*",
      .default = ""
    ),
    insecticide_type_label = sprintf("%s) %s%s",
                                     LETTERS[1 + idx],
                                     insecticide,
                                     suffix),
    insecticide_type_label_2 = sprintf("%s) %s",
                                     LETTERS[idx],
                                     insecticide)
  ) %>%
  select(-idx, -suffix)

# summarise posteriors of plots against all insecticides
pop_mort_sry <- bind_rows(
  lapply(seq_along(types), function(k) {
    summarise_draws(type_draws[[k]][, "Africa", ]) %>%
      mutate(insecticide = types[k])
  })
) %>%
  filter(
    insecticide %in% insecticides_plot
  ) %>%
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

# now aggregate data for plotting Africa-wide averages
df_overall_plot <- df_sub %>%
  
  # first, group by cells, insecticides, and years
  group_by(cell_id,
           insecticide,
           year) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    bioassays = n(),
    .groups = "drop"
  ) %>%
  
  # now compute the population fraction for each cell, for each insecticide, in
  # each year and in the full dataset, to compute weights reducing the annual
  # variability in the data estimates
  
  # for this insecticide, how many mosquitoes were collected in total
  group_by(insecticide) %>%
  mutate(
    total_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide, how many mosquitoes were collected in total in
  # this cell
  group_by(insecticide, cell_id) %>%
  mutate(
    cell_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and this year, how many mosquitoes were collected in
  # total
  group_by(insecticide, year) %>%
  mutate(
    total_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and this year, how many mosquitoes were collected in
  # this cell
  group_by(insecticide, year, cell_id) %>%
  mutate(
    cell_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # compute the cell's fraction of the total population in each year and
  # overall, and compute a corresponding weight
  ungroup() %>%
  mutate(
    year_fraction = cell_year_mosquito_number / total_year_mosquito_number,
    overall_fraction = cell_overall_mosquito_number / total_overall_mosquito_number,
    weight = overall_fraction / year_fraction
  ) %>%
  
  # now compute Africa-wide weighted susceptibilities for each year and insecticide
  group_by(
    year,
    insecticide
  ) %>%
  mutate(
    # components of the relative variance of the weights
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  summarise(
    bioassays = sum(bioassays),
    died_weighted = sum(died * weight),
    mosquito_number_weighted = sum(mosquito_number * weight),
    # died = sum(died),
    mosquito_number = sum(mosquito_number),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    # Susceptibility_raw = died / mosquito_number,
    Susceptibility = died_weighted / mosquito_number_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_bioassays = bioassays / design_effect,
    effective_samples = mosquito_number / design_effect,
    effective_died = Susceptibility * effective_samples,
    .after = insecticide
  ) %>%
  # add labels for plotting
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

# now aggregate the pyrethroids, ensuring even weighting across pyrethroid types
# between years
df_pyrethroids_plot <- df_overall_plot %>%
  filter(
    insecticide %in% llin_pyrethroids
  ) %>%

  # compute target weights for different pyrethroids over the whole dataset
  mutate(
    total_overall_samples = sum(effective_samples)
  ) %>%
  
  group_by(insecticide) %>%
  mutate(
    insecticide_overall_samples = sum(effective_samples)
  ) %>%

  group_by(year) %>%
  mutate(
    total_year_samples = sum(effective_samples)
  ) %>%
  
  group_by(insecticide, year) %>%
  mutate(
    insecticide_year_samples = sum(effective_samples)
  ) %>%
  
  ungroup() %>%
  mutate(
    year_fraction = insecticide_year_samples / total_year_samples,
    overall_fraction = insecticide_overall_samples / total_overall_samples,
    weight = overall_fraction / year_fraction
  ) %>%
  
  # collapse down to year
  group_by(year) %>%
  
  mutate(
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  # now compute the year-specific population fraction for this cell in  insecticide
  summarise(
    bioassays_weighted = sum(effective_bioassays * weight),
    died_weighted = sum(effective_died * weight),
    samples_weighted = sum(effective_samples * weight),
    died = sum(effective_died),
    samples = sum(effective_samples),
    bioassays = sum(effective_bioassays),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility_raw = died / samples,
    Susceptibility = died_weighted / samples_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_samples = samples / design_effect,
    effective_bioassays = bioassays / design_effect,
    .after = year
  )

# plot(df_pyrethroids_plot$Susceptibility_raw ~ df_pyrethroids_plot$Susceptibility)

# make them share a point size legend
max(c(df_pyrethroids_plot$effective_bioassays,
      df_overall_plot$effective_bioassays))
size_limits <- c(0, 500)

# add lines to prediction ribbon for visibility
line_size <- 0.5

# plot all-pyrethroid figure
pyrethroid_fig <- pop_mort_sry_pyrethroid %>%
  ggplot(
    aes(
      x = year
    )
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    ),
    data = pop_mort_sry_pyrethroid,
    fill = pyrethroid_blue,
    colour = pyrethroid_blue,
    size = line_size,
    lineend = "round"
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = effective_bioassays,
    ),
    shape = 21,
    fill = pyrethroid_blue,
    colour = "black",
    data = df_pyrethroids_plot
  ) +
  guides(
    size = "none"
  ) +
  facet_wrap(~ "A) LLIN pyrethroids*") +
  scale_y_continuous(
    labels = scales::percent) +
  scale_size_area(
    limits = size_limits
  ) +
  xlab("") +
  ylab("Susceptibility") +
  coord_cartesian(xlim = c(1995, 2024),
                  ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

all_insecticides_fig <- pop_mort_sry %>%
  ggplot(
    aes(
      x = year,
      fill = insecticide_type_label,
      colour = insecticide_type_label,
    )
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    ),
    size = line_size
  ) +
  geom_point(
    data = df_overall_plot,
    mapping = aes(
      y = Susceptibility,
      group = "none",
      size = effective_bioassays
    ),
    shape = 21,
    colour = "black"
  ) +
  facet_wrap(~insecticide_type_label,
             ncol = 3) +
  guides(
    size = guide_legend(title = "Effective\nsamples")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = c(0, 1),
    limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2000, 2020)) +
  scale_size_area(
    labels = scales::number_format(
      accuracy = 100,
      big.mark = ","
    ),
    breaks = c(1, 3, 5) * 100,
    limits = size_limits
  ) +
  scale_fill_discrete(
    direction = -1,
    guide = "none") +
  scale_colour_discrete(
    direction = -1,
    guide = "none") +
  xlab("") +
  # suppress ylab, as it is on the other panel
  ylab("") +
  coord_cartesian(xlim = c(1995, 2024)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

# use patchwork to set up the multi-panel plot
pyrethroid_fig + all_insecticides_fig

ggsave("figures/all_africa_pred_data.png",
       bg = "white",
       scale = 0.8,
       width = 14,
       height = 6)

# Annual aggregated bioassay results (circles) and predicted population-level
# resistance across the sites where the samples were collected (bold colour
# band; 95%CI)



# now do the same again, by region

# weighted data summarise first
# now aggregate data for plotting Africa-wide averages
df_region_overall_plot <- df_sub %>%
    
  # first, group by cells, regions, insecticides, and years
  group_by(cell_id,
           region,
           insecticide,
           year) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    bioassays = n(),
    .groups = "drop"
  ) %>%
  
  # now compute the population fraction for each cell, for each insecticide, in
  # each year and in the full dataset, to compute weights reducing the annual
  # variability in the data estimates
  
  # for this insecticide and region, how many mosquitoes were collected in total
  group_by(insecticide, region) %>%
  mutate(
    total_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and region, how many mosquitoes were collected in total in
  # this cell
  group_by(insecticide, region, cell_id) %>%
  mutate(
    cell_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and region this year, how many mosquitoes were collected in
  # total
  group_by(insecticide, region, year) %>%
  mutate(
    total_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and region this year, how many mosquitoes were collected in
  # this cell
  group_by(insecticide, region, year, cell_id) %>%
  mutate(
    cell_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # compute the cell's fraction of the total population in each year and
  # overall, and compute a corresponding weight
  ungroup() %>%
  mutate(
    year_fraction = cell_year_mosquito_number / total_year_mosquito_number,
    overall_fraction = cell_overall_mosquito_number / total_overall_mosquito_number,
    weight = overall_fraction / year_fraction
  ) %>%
  # now compute Africa-wide weighted susceptibilities for each year and insecticide
  group_by(
    year,
    insecticide,
    region
  ) %>%
  mutate(
    # components of the relative variance of the weights
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  summarise(
    bioassays = sum(bioassays),
    died_weighted = sum(died * weight),
    mosquito_number_weighted = sum(mosquito_number * weight),
    # died = sum(died),
    mosquito_number = sum(mosquito_number),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    # Susceptibility_raw = died / mosquito_number,
    Susceptibility = died_weighted / mosquito_number_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_bioassays = bioassays / design_effect,
    effective_samples = mosquito_number / design_effect,
    effective_died = Susceptibility * effective_samples,
    .after = insecticide
  ) %>%
  # add labels for plotting
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

# summarise the regional predictions
region_pop_mort_sry <- bind_rows(
  lapply(seq_along(types), function(k) {
    bind_rows(
      lapply(regions, function(region) {
        summarise_draws(type_draws[[k]][, region, ]) %>%
          mutate(region = region)
      })
    ) %>%
      mutate(insecticide = types[k])
  })
) %>%
  # drop out the predictions that are all zero (zero weights)
  group_by(region,
           insecticide,
           year) %>%
  filter(!all(susc_pop_upper == 0)) %>%
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  ) %>%
  mutate(
    region = factor(region,
                    levels = c("Southern Africa",
                               "Western Africa",
                               "Middle Africa", 
                               "Eastern Africa",
                               "Northern Africa"))
  )


region_all <- region_pop_mort_sry %>%
  ggplot(
    aes(
      x = year
    )
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
      fill = region,
    ),
    colour = grey(0.2),
    size = 0.1
  ) +
  facet_wrap(~insecticide_type_label_2,
             ncol = 3) +
  guides(
    size = guide_legend(title = "Effective\nsamples")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = c(0, 1),
    limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2000, 2020)) +
  scale_size_area(
    labels = scales::number_format(
      accuracy = 100,
      big.mark = ","
    ),
    breaks = c(1, 3, 5) * 100,
    limits = size_limits
  ) +
  scale_fill_brewer(
    palette = "Accent"
  ) +
  xlab("") +
  # suppress ylab, as it is on the other panel
  ylab("") +
  coord_cartesian(xlim = c(1995, 2024)) +
  guides(
    fill = guide_legend(title = "")
  ) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

region_all

ggsave("figures/regional_pred_data.png",
       bg = "white",
       scale = 0.8,
       width = 8,
       height = 6)
