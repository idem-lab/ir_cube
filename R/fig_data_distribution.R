# make plots of the available data by time, space, and insecticide type

source("R/packages.R")
source("R/functions.R")

source("R/bioassay_subset.R")

# the modelled bioassays and their indices, as in R/fit_model.R
list2env(modelled_bioassays(), environment())

insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  pull(insecticide)

colour_types <- scales::hue_pal(direction = -1)(9)

df_plot <- df %>%
  rename(
    country = country_name,
    insecticide = insecticide_type
  ) %>%
  group_by(
    country, year, insecticide
  ) %>%
  filter(
    year >= 1995
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    country = factor(country, levels = rev(sort(unique(country))))
  )

for (this_insecticide in insecticides_plot) {
  df_plot %>%
    filter(
      insecticide == this_insecticide
    ) %>%
    ggplot(
      aes(
        x = year,
        y = country,
        size = mosquito_number,
      )
    ) +
    geom_point(
      fill = colour_types[match(this_insecticide, insecticides_plot)],
      shape = 21,
      alpha = 1
    ) +
    ylab("") +
    xlab("") +
    scale_size_continuous(
      labels = scales::number_format(accuracy = 1000, big.mark = ","),
      limits = range(df_plot$mosquito_number)
    ) +
    scale_x_continuous(
      limits = range(df_plot$year)
    ) +
    scale_y_discrete(
      # labels = rev(sort(unique(df_plot$country))),
      drop = FALSE
    ) +
    guides(
      size = guide_legend(title = "No. tested")
    ) +
    theme_minimal() +
    ggtitle(
      sprintf("Data availability for %s",
              this_insecticide)
    )
  
  ggsave(sprintf("figures/data_time_%s.png",
                 this_insecticide),
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
}
