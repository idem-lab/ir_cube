# make a figure showing the fitted relationships between covariates and
# selection for each insecticide

# load packages and functions
# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# get posterior means and standard deviations of these
sims <- calculate(effect_type, values = draws, nsim = 1e4)

# make into tibbles for plotting
mean_tibble <- sims$effect_type %>%
  apply(2:3, mean) %>%
  `colnames<-`(types) %>%
  as_tibble() %>%
  mutate(covariate = colnames(x_cell_years)) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "insecticide",
    values_to = "mean"
  )

sd_tibble <- sims$effect_type %>%
  apply(2:3, sd) %>%
  `colnames<-`(types) %>%
  as_tibble() %>%
  mutate(covariate = colnames(x_cell_years)) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "insecticide",
    values_to = "sd"
  )

# get the mean and sd of the scaled covariates in the data
data_sry <- x_cell_years %>%
  as_tibble() %>%
  pivot_longer(
    cols = everything(),
    names_to = "covariate",
    values_to = "value"
  ) %>%
  group_by(covariate) %>%
  summarise(
    data_mean = mean(value),
    data_sd = sd(value),
    .groups = "drop"
  )

stats <- mean_tibble %>%
  left_join(
    sd_tibble,
    by = c("covariate", "insecticide")
  ) %>%
  left_join(
    data_sry,
    by = "covariate"
  ) %>%
  mutate(
    mean_scaled = mean * data_mean,
    # is this right foraffine transform of a normal? Can't check as I'm on a
    # plane
    sd_scaled = sd * data_mean,
    precision = 1 / sd ^ 2,
    precision_scaled = 1 / sd_scaled ^ 2,
    inv_coef_variation = mean / sd ^ 2,
    inv_coef_variation_scaled = mean_scaled / sd_scaled ^ 2
  ) %>%
  # rename some things for plotting
  mutate(
    covariate = case_when(
      # the population column carries its transform and trend, e.g.
      # pop_enc:g_dom, and the crops their trend, e.g. rice:g_ag
      # (selection_column_names())
      str_detect(covariate, "^pop") ~ "Human population",
      covariate == "nets" ~ "LLIN use",
      covariate == "irs" ~ "IRS coverage",
      .default = str_to_sentence(str_remove(covariate, ":g_ag$"))
    )
  )

# order the covariates by the mean effect
covariate_order_scaled <- stats %>%
  group_by(covariate) %>%
  summarise(mean_scaled = mean(mean_scaled)) %>%
  arrange(mean_scaled) %>%
  pull(covariate)

stats <- stats %>%
  mutate(
    covariate = factor(covariate,
                       levels = covariate_order_scaled)
  )

# stats for combined effects
combined_stats <- stats %>%
  mutate(
    covariate_type = case_when(
      covariate %in% c("LLIN use", "IRS coverage") ~ "Vector control",
      covariate == "Human population" ~ "Human population",
      .default = "Agriculture"
    ),
    covariate_type = factor(covariate_type,
                            levels = c("Human population",
                                       "Agriculture",
                                       "Vector control")),
    .before = everything()
  ) %>%
  group_by(covariate_type, insecticide) %>%
  summarise(
    mean_scaled = sum(mean_scaled),
    .groups = "drop"
  )


overall_stats <- stats %>%
  group_by(insecticide) %>%
  summarise(
    mean_scaled = sum(mean_scaled),
    .groups = "drop"
  ) %>%
  mutate(
    covariate = "Overall"
  )

# plot the overall selection pressure
overall_stats_effect_grid <- overall_stats %>%
  rename(
    `Average\nselection\npressure` = mean_scaled
  ) %>%
  ggplot(
    aes(
      x = insecticide,
      y = covariate,
      fill = `Average\nselection\npressure`,
    )
  ) +
  geom_tile(
    colour = grey(0.5)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1)
  ) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0)
  )

# combine the different types of covariates
combined_effect_grid <- combined_stats %>%
  rename(
    `Average\nselection\npressure` = mean_scaled
  ) %>%
  ggplot(
    aes(
      x = insecticide,
      y = covariate_type,
      fill = `Average\nselection\npressure`,
    )
  ) +
  geom_tile(
    colour = grey(0.5)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1)
  ) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0)
  )

effect_grid <- stats %>%
  rename(
    `Average\nselection\npressure` = mean_scaled
  ) %>%
  ggplot(
    aes(
      x = insecticide,
      y = covariate,
      fill = `Average\nselection\npressure`
    )
  ) +
  geom_tile(
    colour = grey(0.5)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1)
  ) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0)
  )

# The average contribution of each covariate to the selection coefficient for
# resistance to each insecticide. Computed as the posteiror mean of the partial
# selection coefficient, multiplied by the mean over the dataset of the
# indicator variable for the covariate. This captures both the relative strength
# of the covariate in selecting for resistance, and its prevalence.


# do the same on the coefficients themselves
regression_grid <- stats %>%
  rename(
    `Effect\nsize` = mean
  ) %>%
  ggplot(
    aes(
      x = insecticide,
      y = covariate,
      fill = `Effect\nsize`
    )
  ) +
  geom_tile(
    colour = grey(0.5)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1)
  ) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0)
  )

# plot the combined and specific selection pressures

# suppress some plot elements in each
max_pressure <- max(c(overall_stats$mean_scaled,
                      combined_stats$mean_scaled,
                      stats$mean_scaled))
fill_scale <- scale_fill_gradient(low = "white",
                                  high = "red",
                                  limits = c(0, max_pressure)) 
(overall_stats_effect_grid +
    fill_scale +
    theme(legend.position = "none",
          plot.title.position = "plot")) +
(combined_effect_grid +
    fill_scale +
   theme(axis.text.x = element_blank(),
         legend.position = "none",
         plot.title.position = "plot")) +  
  (effect_grid +
     fill_scale +
     theme(axis.text.x = element_blank(),
           plot.title.position = "plot")) +
  plot_layout(nrow = 3, heights = c(1, 3, 9)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12))

ggsave(
  "figures/covariate_loadings.png",
  bg = "white",
  scale = 0.8,
  width = 8,
  height = 8
)  


regression_grid

ggsave(
  "figures/covariate_loadings_effect_size.png",
  bg = "white",
  scale = 0.8,
  width = 8,
  height = 6
)  
