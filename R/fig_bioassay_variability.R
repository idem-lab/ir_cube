# Figure of the variability between replicate bioassays at the six most
# sampled pixel-years of the LLIN pyrethroids, against the 95% sampling
# intervals at 100 mosquitoes: beta-binomial at the type's overdispersion rho
# between bioassays of one study, and binomial (as if the mosquitoes were
# independent). Points are coloured by study.
#
#   Rscript R/fig_bioassay_variability.R
#
# Reads rho from outputs/bioassay_rho_hierarchical.csv, so run after
# R/fig_illustrate_bioassay_variability.R, which fits it.
#
# Multiple discriminating concentration bioassay results at single locations
# (5 km pixels) and years, as an example of the inherent variability in this
# type of data. Dark grey horizontal lines give the (sample size weighted) mean
# of the bioassay mortality, the bands the 95% sampling intervals with this
# mean and 100 mosquitoes.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")

mask <- rast("data/clean/raster_mask.tif")
df <- subset_modelled_bioassays(
  readRDS("data/clean/all_gambiae_complex_data.RDS"), mask)
rho <- read.csv("outputs/bioassay_rho_hierarchical.csv")
rho <- setNames(rho$rho, rho$insecticide_type)

insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

# the six pixel-year-types with the most bioassays, labelled by their site
# (the most common site name among their records) and country, and year; and
# by the insecticide only if they are not all the same one
most_sampled <- df %>%
  filter(insecticide_type %in% insecticides_plot_small) %>%
  group_by(cell, year, insecticide_type) %>%
  summarise(
    count = n(),
    site = names(which.max(table(site))),
    country_name = names(which.max(table(country_name))),
    .groups = "drop"
  ) %>%
  arrange(desc(count)) %>%
  slice(1:6)
one_type <- n_distinct(most_sampled$insecticide_type) == 1
most_sampled <- most_sampled %>%
  mutate(
    label = sprintf("%s\n%s\n%i%s", site, country_name, year,
                    if (one_type) "" else paste0("\n", insecticide_type))
  )

df_most_sampled <- most_sampled %>%
  select(cell, year, insecticide_type, label) %>%
  left_join(df, by = c("cell", "year", "insecticide_type")) %>%
  mutate(
    Susceptibility = died / mosquito_number
  )

# the pooled mortality and its 95% sampling intervals at 100 mosquitoes
interval_labels <- c(betabinomial = "95% CI",
                     binomial = "95% CI if independent")
df_most_sampled_stats <- df_most_sampled %>%
  group_by(label, insecticide_type) %>%
  summarise(
    Susceptibility = sum(died) / sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    binom_lower_100 = qbinom(0.025, 100, Susceptibility) / 100,
    binom_upper_100 = qbinom(0.975, 100, Susceptibility) / 100,
    rho = rho[insecticide_type],
    alpha = Susceptibility * (1 / rho - 1),
    beta = alpha * (1 - Susceptibility) / Susceptibility,
    betabinom_lower_100 = qbbinom(0.025, 100, alpha, beta) / 100,
    betabinom_upper_100 = qbbinom(0.975, 100, alpha, beta) / 100
  )

# the fill scale for the bands
fill_values <- setNames(c(grey(0.9), grey(0.7)), interval_labels)

# each point coloured by its study, labelled by these names
study_names <- c(
  "mbogo 2016, Kenya" = "Mbogo (2016)",
  "Kenya Medical Research Institute, Kenya" = "KEMRI",
  "githinji 2020, Kenya" = "Githinji (2020)",
  "Project IIR: Implications of insecticide resistance for vector control, Kenya" =
    "Project IIR",
  "cook 2018, Kenya" = "Cook (2018)",
  "cook 2018, Benin" = "Cook (2018)",
  "hessoudjossou 2024, Benin" = "Hessoudjossou (2024)",
  "ismail 2018, Sudan" = "Ismail (2018)",
  "National Malaria Programme, Sudan" = "Sudan NMCP",
  "yahouedo 2016, Benin" = "Yahouedo (2016)")
stopifnot(all(df_most_sampled$study %in% names(study_names)))

set.seed(1)
df_points <- df_most_sampled %>%
  # shuffle the x axes so the points don't overlap too much
  mutate(
    x_random = sample(seq(0.2, 1.8, length.out = n())),
    study_name = unname(study_names[study]),
    # the panels, in the order they are plotted
    label = factor(label)
  )

# the studies in the order they appear from left to right: by panel, then by
# their leftmost point. The colours are Paul Tol's muted palette, the most
# distinct first, in order of appearance
study_order <- df_points %>%
  group_by(study_name) %>%
  summarise(panel = min(as.integer(label)),
            x = min(x_random[as.integer(label) == panel]),
            .groups = "drop") %>%
  arrange(panel, x) %>%
  pull(study_name)
study_colours <- setNames(
  c("#332288", "#DDCC77", "#117733", "#CC6677", "#88CCEE", "#882255",
    "#44AA99", "#AA4499", "#999933")[seq_along(study_order)],
  study_order)
stopifnot(!anyNA(study_colours))

df_points %>%
  ggplot(
    aes(
      xmin = 0,
      xmax = 2,
      group = label
    )
  ) +
  geom_rect(
    aes(
      ymax = betabinom_upper_100,
      ymin = betabinom_lower_100,
      fill = interval_labels[["betabinomial"]]
    ),
    data = df_most_sampled_stats
  ) +
  geom_rect(
    aes(
      ymax = binom_upper_100,
      ymin = binom_lower_100,
      fill = interval_labels[["binomial"]]
    ),
    data = df_most_sampled_stats
  ) +
  geom_rect(
    aes(
      ymax = Susceptibility,
      ymin = Susceptibility
    ),
    data = df_most_sampled_stats,
    linewidth = 1,
    colour = grey(0.4)
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = pmin(mosquito_number, 500),
      x = x_random,
      colour = factor(study_name, levels = study_order)
    ),
    shape = 16
  ) +
  # the points' black outlines: a second layer, because the fill scale is the
  # bands' and the colour scale the studies'
  geom_point(
    aes(
      y = Susceptibility,
      size = pmin(mosquito_number, 500),
      x = x_random
    ),
    shape = 1,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = study_colours,
    name = "Study"
  ) +
  # point areas proportional to the number tested, capped at 500 so that the
  # size legend spans the common numbers tested
  scale_size_area(
    breaks = c(100, 250, 500),
    labels = c("100", "250", "500+"),
    limits = c(0, 500)
  ) +
  facet_wrap(~label,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    limits = c(0, 2)
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = unname(interval_labels),
    name = NULL
  ) +
  labs(x = NULL, y = "Mortality") +
  theme_minimal() +
  guides(
    size = guide_legend(title = "No. tested", order = 2),
    colour = guide_legend(order = 3, override.aes = list(size = 3)),
    fill = guide_legend(order = 1)
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggtitle(
    "Significant variability in replicated bioassay results"
  )

ggsave(
  "figures/bioassay_variability.png",
  bg = "white",
  scale = 0.9,
  width = 9.5,
  height = 6
)
