# plot modelled and data trends in pyrethroid resistance in high and low net
# coverage places

# load packages and functions
# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# the covariates at the data cells, on their own scales (R/model_covariates.R)
source("R/model_covariates.R")
all_extract <- covariate_extract(unique_cells, baseline_year, final_data_year,
                                 model_options$selection_columns)

# load time-varying net use data and flatten it
nets_cube <- rast("data/clean/net_use_cube.tif")
years_sub <- 2010:2024
nets_cube_sub <- nets_cube[[paste0("nets_", years_sub)]]
nets_flat <- app(nets_cube_sub, "mean")

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

# tidy up column names
df_sub <- df %>%
  rename(
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  )


# add the low, medium, high LLINs underneath
# find how much data by different coverages

net_use_cell_lookup <- df %>%
  group_by(
    cell_id, cell
  ) %>%
  # get one record per observed cell
  slice(1) %>%
  ungroup() %>%
  mutate(
    net_use = terra::extract(nets_flat, pull(., cell))[, 1],
    net_use_class = case_when(
      net_use < 0.3 ~ "A) Low use",
      .default = "B) High use")
  ) %>%
  select(
    cell,
    cell_id,
    net_use_class
  )

# extract net use over time at all observation locations
df_net_use <- all_extract %>%
  left_join(
    net_use_cell_lookup,
    by = "cell_id"
  ) %>%
  group_by(
    year_id,
    net_use_class
  ) %>%
  summarise(
    lower = quantile(nets, 0.25),
    upper = quantile(nets, 0.75),
    mean = mean(nets),
    .groups = "drop"
  ) %>%
  mutate(
    year = year_id + baseline_year - 1
  )

# aggregate the same data across these coverage classes to plot each against time
df_net_class_points <- df_sub %>%
  # first, group by cells, insecticides, and years
  group_by(year,
           insecticide,
           cell_id) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    bioassays = n(),
    .groups = "drop"
  ) %>%
  
  # add net coverage classes for these cells
  left_join(
    net_use_cell_lookup,
    by = "cell_id"
  ) %>%
  
  # now compute survey weights across the cells within these insecticides and net
  # coverage classes, to reduce inter-year variability in representativeness
  
  # for this insecticide and coverage class, how many mosquitoes were collected in total
  group_by(insecticide, net_use_class) %>%
  mutate(
    total_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and coverage class, how many mosquitoes were collected in total in
  # this cell
  group_by(insecticide, net_use_class, cell_id) %>%
  mutate(
    cell_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and coverage class and this year, how many mosquitoes
  # were collected in total
  group_by(insecticide, net_use_class, year) %>%
  mutate(
    total_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and coverage class and year, how many mosquitoes were
  # collected in this cell
  group_by(insecticide, net_use_class, year, cell_id) %>%
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
  
  # now compute Africa-wide weighted susceptibilities for each year and
  # insecticide and net coverage class
  group_by(
    year,
    insecticide,
    net_use_class
  ) %>%
  mutate(
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  # now compute the year-specific population fraction for this cell in  insecticide
  summarise(
    died_weighted = sum(died * weight),
    mosquito_number_weighted = sum(mosquito_number * weight),
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    bioassays = sum(bioassays),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility_raw = died / mosquito_number,
    Susceptibility = died_weighted / mosquito_number_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_samples = mosquito_number / design_effect,
    effective_bioassays = bioassays / design_effect,
    effective_died = Susceptibility * effective_samples,
    .after = insecticide
  )

# plot(df_net_class_points$Susceptibility ~ df_net_class_points$Susceptibility_raw)


# compute predicted population-level susceptibility to the insecticides
# average net coverage over time in each of those places

# now compute a weighted sum over the pyrethroids (weights given by the numbers
# of mosquitoes tested for each pyrethroid, over the whole timeseries)
pyrethroid_net_class_points <- df_net_class_points %>%
  filter(
    insecticide %in% pyrethroids
  ) %>%
  
  # compute target weights for different pyrethroids over the whole dataset for each coverage class
  group_by(net_use_class) %>%
  mutate(
    total_overall_samples = sum(effective_samples)
  ) %>%
  group_by(insecticide, net_use_class) %>%
  mutate(
    insecticide_overall_samples = sum(effective_samples)
  ) %>%
  group_by(year, net_use_class) %>%
  mutate(
    total_year_samples = sum(effective_samples)
  ) %>%
  group_by(insecticide, year, net_use_class) %>%
  mutate(
    insecticide_year_samples = sum(effective_samples)
  ) %>%
  
  ungroup() %>%
  mutate(
    year_fraction = insecticide_year_samples / total_year_samples,
    overall_fraction = insecticide_overall_samples / total_overall_samples,
    weight = overall_fraction / year_fraction
  ) %>%
  
  # collapse down to year and coverage class
  group_by(year, net_use_class) %>%
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

# now compute a weighted sum over the pyrethroids (weights given by the numbers
# of mosquitoes tested for each pyrethroid, over the whole timeseries)
pyrethroid_net_class_weights <- df_sub %>%
  left_join(
    net_use_cell_lookup,
    by = "cell_id"
  ) %>%
  group_by(
    type_id,
    insecticide,
    net_use_class
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  group_by(net_use_class) %>%
  mutate(
    is_a_pyrethroid = insecticide %in% pyrethroids,
    mosquito_number = mosquito_number * as.numeric(is_a_pyrethroid),
    weight = mosquito_number / sum(mosquito_number)
  ) %>%
  ungroup() %>%
  select(
    type_id,
    weight,
    net_use_class
  ) %>%
  pivot_wider(
    names_from = net_use_class,
    values_from = "weight"
  ) %>%
  arrange(type_id) %>%
  select(-type_id) %>%
  as.matrix()

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

# AAARGH, this code is horrible, but best I can do for now on this plane.

# compute a weighted average of the insecticides in low, medium, high net
# coverage locations

# just subset the data and repeat what happens in the overall plot

low_idx <- net_use_cell_lookup %>%
  filter(
    net_use_class == "A) Low use"
  ) %>%
  pull(cell_id)

# subset and renormalise weights
weights_mat_low <- weights_mat[low_idx, ]
weights_mat_low <- sweep(weights_mat_low,
                         2,
                         colSums(weights_mat_low),
                         "/")

# expand out into an array, replicating by year
weights_array_low <- array(rep(c(weights_mat_low), n_times),
                           dim = c(nrow(weights_mat_low), n_types, n_times))

# check dimensions are the right way around
# identical(weights_array_low[, , 1],  weights_mat_low)
# all(weights_array_low[1, 1, ] == weights_array_low[1, 1, 1])

# multiply through and sum to get average susceptibilities
weighted_susc_array_low <- dynamic_cells$all_states[low_idx, , ] * weights_array_low
overall_susc_low <- apply(weighted_susc_array_low, 2:3, "sum")
pyrethroid_net_class_susc_low <- as_data(t(pyrethroid_net_class_weights[, "A) Low use"])) %*% overall_susc_low

# and repeat for high
high_idx <- net_use_cell_lookup %>%
  filter(
    net_use_class == "B) High use"
  ) %>%
  pull(cell_id)


# subset and renormalise weights
weights_mat_high <- weights_mat[high_idx, ]
weights_mat_high <- sweep(weights_mat_high,
                         2,
                         colSums(weights_mat_high),
                         "/")

# expand out into an array, replicating by year
weights_array_high <- array(rep(c(weights_mat_high), n_times),
                            dim = c(nrow(weights_mat_high), n_types, n_times))

# multiply through and sum to get average susceptibilities
weighted_susc_array_high <- dynamic_cells$all_states[high_idx, , ] * weights_array_high
overall_susc_high <- apply(weighted_susc_array_high, 2:3, "sum")
pyrethroid_net_class_susc_high <- as_data(t(pyrethroid_net_class_weights[, "B) High use"])) %*% overall_susc_high

overall_susc_index <- expand_grid(
  type_id = seq_len(n_types),
  time_id = seq_len(n_times)
)

overall_susc_low_vec <- overall_susc_low[as.matrix(overall_susc_index)]
overall_susc_high_vec <- overall_susc_high[as.matrix(overall_susc_index)]


# get posterior samples of these
sims_net_class <- calculate(pyrethroid_net_class_susc_low,
                            pyrethroid_net_class_susc_high,
                            overall_susc_low_vec,
                            overall_susc_high_vec,
                            values = draws,
                            nsim = 1e3)


overall_mort_sry_low <- sims_net_class$overall_susc_low_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    overall_susc_index,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id, type_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    year = baseline_year + time_id - 1,
    insecticide = types[type_id],
    .before = everything()
  )

overall_mort_sry_high <- sims_net_class$overall_susc_high_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    overall_susc_index,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id, type_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    year = baseline_year + time_id - 1,
    insecticide = types[type_id],
    .before = everything()
  )

pop_mort_sry_low <- sims_net_class$pyrethroid_net_class_susc_low[, 1, ] %>%
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

pop_mort_sry_high <- sims_net_class$pyrethroid_net_class_susc_high[, 1, ] %>%
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

pop_mort_sry_net_class <- bind_rows(
  `A) Low use` = pop_mort_sry_low,
  `B) High use` = pop_mort_sry_high,
  .id = "net_use_class"
)

overall_mort_sry_net_class <- bind_rows(
  `A) Low use` = overall_mort_sry_low,
  `B) High use` = overall_mort_sry_high,
  .id = "net_use_class"
)

pyrethroid_mort_sry_net_class <- overall_mort_sry_net_class %>%
  filter(insecticide %in% pyrethroids)

net_class_fig <- pyrethroid_net_class_points %>%
  ggplot(
    aes(
      x = year
    )
  ) +
  geom_ribbon(
    aes(
      x = year,
      ymax = upper,
      ymin = lower,
      group = net_use_class
    ),
    data = df_net_use,
    fill = grey(0.9),
    colour = grey(0.7)
  ) +
  # geom_ribbon(
  #   aes(
  #     ymax = susc_pop_upper,
  #     ymin = susc_pop_lower,
  #     fill = insecticide
  #   ),
  #   data = pyrethroid_mort_sry_net_class,
  #   # fill = pyrethroid_blue
  # ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
      fill = insecticide
    ),
    data = pop_mort_sry_net_class,
    fill = pyrethroid_blue
  ) +
  # geom_point(
  #   aes(
  #     y = Susceptibility,
  #     size = effective_bioassays,
  #     fill = insecticide
  #   ),
  #   data = filter(df_net_class_points, insecticide %in% pyrethroids),
  #   shape = 21
  # ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = effective_bioassays
    ),
    shape = 21,
    fill = pyrethroid_blue
  ) +
  facet_wrap(~net_use_class) +
  guides(
    size = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    sec.axis = sec_axis(
      ~.,
      name = "LLIN use (grey)",
      labels = scales::percent,
    )
  ) +
  scale_size_area() +
  xlab("") +
  ylab("Susceptibility to pyrethroids") +
  coord_cartesian(xlim = c(2000, 2024),
                  ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

net_class_fig

ggsave("figures/all_africa_net_use_pred_data.png",
       bg = "white",
       scale = 0.8,
       width = 8,
       height = 4)
