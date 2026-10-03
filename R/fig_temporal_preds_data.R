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
# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

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
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  )
  
# subset to pyrethroids, and add on UN geoscheme regions for Africa
df_pyrethroids <- df_sub %>%
  filter(
    insecticide %in% llin_pyrethroids
    # insecticide_class == "Pyrethroids",
  )

# compute predicted population-level resistance to these pyrethroids, weighted
# over the numbers of samples per location, per insecticide type

# make posterior predictions of Africa-wide average susceptibility (over
# locations with any bioassays) for all years with the average sample size


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

# we already have a greta array for susceptibility in all unique cells in the
# dataset, for all insecticides and all years: dynamic_cells$all_states

# expand out into an array, replicating by year
weights_array <- array(rep(c(weights_mat), n_times),
                       dim = c(n_unique_cells, n_types, n_times))

# check dimensions are the right way around
# identical(weights_array[, , 1],  weights_mat)
# all(weights_array[1, 1, ] == weights_array[1, 1, 1])

# multiply through and sum to get average susceptibilities
weighted_susc_array <- dynamic_cells$all_states * weights_array
overall_susc <- apply(weighted_susc_array, 2:3, "sum")

# expand these out into a vector by insecticides and years, to later compile
# sims
indices <- expand_grid(
  type_id = seq_len(n_types),
  time_id = seq_len(n_times)
)

overall_susc_vec <- overall_susc[as.matrix(indices)]
overall_pop_mort_vec <- overall_susc_vec

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

pyrethroid_susc <- as_data(t(pyrethroid_weights)) %*% overall_susc


# get weighted average pyrethroid susceptibility for all cells and times
pyrethroid_weights_array <- array(pyrethroid_weights,
                                  dim = c(n_types, n_times, n_unique_cells))
pyrethroid_weights_array <- aperm(pyrethroid_weights_array, c(3, 1, 2))
pyrethroid_all_cells <- apply(pyrethroid_weights_array * dynamic_cells$all_states, c(1, 3), FUN = "sum")

# vector version of all cells and times for the weighted pyrethroid estimate
indices_all_pyrethroid <- expand_grid(
  cell_id = seq_len(n_unique_cells),
  time_id = seq_len(n_times)
)

pyrethroid_all_cells_vec <- pyrethroid_all_cells[as.matrix(indices_all_pyrethroid)]

# get posterior samples of these
sims <- calculate(overall_pop_mort_vec,
                  pyrethroid_susc,
                  pyrethroid_all_cells_vec,
                  values = draws,
                  nsim = 1e3)

# posterior summaries of pyrethroid susceptibility by year
pop_mort_sry_pyrethroid <- sims$pyrethroid_susc[,  1, ] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    time_id = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  mutate(
    year = baseline_year + time_id - 1
  )

# posterior means for all cells
pop_mort_sry_pyrethroid_all <- sims$pyrethroid_all_cells_vec[,  , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    indices_all_pyrethroid,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(cell_id, time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    .groups = "drop"
  ) %>%
  mutate(
    year = baseline_year + time_id - 1
  )

type_years <- indices %>%
  mutate(
    year = baseline_year + time_id - 1,
    insecticide = types[type_id],
    row = row_number()
  )  %>%
  select(row, year, insecticide)

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
pop_mort_sry <- sims$overall_pop_mort_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    row = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(row) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  left_join(type_years,
            by = "row") %>%
  select(-row) %>%
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

# ggplot() +
#   geom_point(
#     data = df_overall_plot,
#     mapping = aes(
#       y = Susceptibility,
#       group = "none",
#       size = effective_bioassays
#     ),
#     shape = 21,
#     colour = "black"
#   ) +
#   
# make a matrix mapping from cell IDs to regions, to aggregate by region
# for each insecticide.
regions <- unique(df_sub$region)
n_regions <- length(regions)
cell_region_mask <- df_sub %>%
  group_by(cell_id, type_id) %>%
  # need to collapse down to weights
  group_by(cell_id, ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    region_id = match(region, regions)
  ) %>%
  select(
    cell_id,
    region_id
  ) %>%
  mutate(
    value = 1
  ) %>%
  complete(
    cell_id,
    region_id,
    fill = list(value = 0)
  ) %>%
  arrange(
    cell_id,
    region_id
  ) %>%
  pivot_wider(
    names_from = region_id,
    values_from = value
  ) %>%
  select(-cell_id) %>%
  as.matrix() %>%
  `colnames<-`(NULL)

# make a matrix of weights, for the cells, that aggregate up over regions
regional_weights_mat <- df_sub %>%
  mutate(
    region_id = match(region, regions)
  ) %>%
  group_by(cell_id, type_id, region_id) %>%
  # collapse totals by cell and type 
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  # pad out empty combinations with 0s (nesting cell within region)
  complete(
    type_id,
    nesting(cell_id, region_id),
    fill = list(mosquito_number = 0)
  ) %>%
  # normalise within type and region
  group_by(type_id, region_id) %>%
  mutate(
    weight = mosquito_number / sum(mosquito_number),
    # remove divide-by-zero errors (not all regions represented for each
    # insecticide)
    weight = replace_na(weight, 0)
  ) %>%
  ungroup() %>%
  # then drop region info and mosquito numbers (we pull region info back later
  # from cell_region_mask)
  select(cell_id, type_id, weight) %>%
  # turn into a cell-by-type matrix of weights
  arrange(cell_id, type_id) %>%
  pivot_wider(
    names_from = type_id,
    values_from = weight
  ) %>%
  select(-cell_id) %>%
  as.matrix() %>%
  `colnames<-`(NULL)

# colSums(regional_weights_mat)
regional_weights_array <- array(rep(c(regional_weights_mat), n_times),
                                dim = c(n_unique_cells, n_types, n_times))

# multiply through and sum to get average susceptibilities for all insecticides, for all regions

regional_weighted_susc_mat <- dynamic_cells$all_states * regional_weights_array
dim(regional_weighted_susc_mat) <- c(dim(regional_weighted_susc_mat)[1], prod(dim(regional_weighted_susc_mat)[2:3]))
weighted_susc_region_array <- t(cell_region_mask) %*% regional_weighted_susc_mat
dim(weighted_susc_region_array) <- c(dim(cell_region_mask)[2], dim(weighted_susc_array)[2:3])

region_type_time_index <- expand_grid(
  region_id = seq_len(n_regions),
  type_id = seq_len(n_types),
  time_id = seq_len(n_times)
)
weighted_susc_region_vec <- weighted_susc_region_array[as.matrix(region_type_time_index)]

# now calculate these values and summarise for plotting
sims_regional <- calculate(
  weighted_susc_region_vec,
  values = draws,
  nsim = 1e3
)

region_pop_mort_sry <- sims_regional$weighted_susc_region_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    region_type_time_index,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(region_id, time_id, type_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    year = baseline_year + time_id - 1,
    insecticide = types[type_id],
    region = regions[region_id],
    .before = everything()
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
