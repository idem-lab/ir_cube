# The two-stage model's mortality over time at the most heavily sampled
# pixels of the LLIN pyrethroids, with the bioassays there: the population
# mortality and its 95% interval, and 95% posterior predictive intervals for a
# bioassay of 100 mosquitoes, under binomial sampling and under the model's own
# observation process.
#
#   Rscript R/fig_fit_subset.R
#
# Run after R/two_stage_maps.R. The predictions are cached in
# temporary/fig_fit_subset.rds (R/figure_cache.R), so changes to the plot
# alone rerun in seconds.

source("R/packages.R")
source("R/functions.R")
source("R/two_stage_predictions.R")
source("R/validation_functions.R")
source("R/figure_cache.R")

# the modelled bioassays and their cells, as the full fit has them
fit <- new.env()
load("temporary/fitted_model.RData", envir = fit)
df <- fit$df
unique_cells <- fit$unique_cells
years_predict <- fit$baseline_year:fit$final_data_year
rm(fit)
invisible(gc())
cell_country <- data_cell_country(df)

insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

# the pixels with the most pyrethroid bioassay data: at least 27 bioassays
# over at least 8 years, of all three pyrethroids (8 pixels), each named by
# the most common site and country of its records
locations_plot <- df %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  group_by(cell, cell_id) %>%
  filter(
    n() >= 27,
    n_distinct(year) >= 8,
    n_distinct(insecticide_type) >= 3
  ) %>%
  summarise(
    site = names(which.max(table(site))),
    country_name = names(which.max(table(country_name))),
    .groups = "drop"
  ) %>%
  mutate(
    place = paste(site, country_name, sep = "\n")
  ) %>%
  select(
    place,
    cell_id
  )

# predict these with the two-stage model: the population-level mortality,
# ilogit(m + omega + xi), and the mortality in new bioassays of
# `sample_size` mosquitoes under binomial sampling from it, and under the
# model's own observation process: beta-binomial sampling, at the per-type
# overdispersion, from the mortality with fresh pixel-year, pixel and
# study-year noise, ilogit(m + omega + xi + u + p + sy). The posterior mean
# and 95% intervals per type, cell and year
predict_locations <- function(types, years, df, unique_cells, cell_country,
                              cell_ids, sample_size) {
  bind_rows(lapply(types, function(type) {
    setup <- two_stage_setup(type, years, df)
    identity <- diag(length(cell_ids))
    cells <- unique_cells[cell_ids]
    country <- cell_country[cell_ids]
    population <- two_stage_weighted_draws(setup, cells, country,
                                           identity)$two_stage
    assay <- two_stage_weighted_draws(setup, cells, country, identity,
                                      noise = TRUE)$two_stage
    rho <- rho_for_record(tibble(insecticide_type = type), rho_lookup())
    binomial_mortality <- array(
      rbinom(length(population), sample_size, population),
      dim(population)) / sample_size
    betabinomial_mortality <- array(
      rbetabinom(length(assay), sample_size, assay, rho),
      dim(assay)) / sample_size
    # a quantile over the draws, for every cell and year
    quantiles <- function(x, prob) as.vector(apply(x, 2:3, quantile, prob))

    # posterior mean mortality rate, and intervals for the population value
    # (posterior uncertainty), and posterior predictive intervals (posterior
    # uncertainty and sampling error) for observed mortality from binomial
    # sampling, and from the model's observation process
    expand_grid(year = years, cell_id = cell_ids) %>%
      mutate(
        insecticide_type = type,
        Mortality = as.vector(apply(population, 2:3, mean)),
        pop_lower = quantiles(population, 0.025),
        pop_upper = quantiles(population, 0.975),
        binomial_lower = quantiles(binomial_mortality, 0.025),
        binomial_upper = quantiles(binomial_mortality, 0.975),
        betabinomial_lower = quantiles(betabinomial_mortality, 0.025),
        betabinomial_upper = quantiles(betabinomial_mortality, 0.975)
      )
  }))
}

preds_plot <- cached_figure_data(
  "temporary/fig_fit_subset.rds", predict_locations,
  args = list(types = insecticides_plot_small, years = years_predict, df = df,
              unique_cells = unique_cells, cell_country = cell_country,
              cell_ids = locations_plot$cell_id, sample_size = 100),
  sources = c("R/packages.R", "R/functions.R", "R/two_stage_predictions.R",
              "R/validation_functions.R"),
  inputs = two_stage_figure_inputs(insecticides_plot_small)
) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  ) %>%
  mutate(
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_small)
  )

points_plot <- df %>%
  filter(
    cell_id %in% locations_plot$cell_id,
    insecticide_type %in% insecticides_plot_small
  ) %>%
  mutate(
    Mortality = died / mosquito_number
  ) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  ) %>%
  mutate(
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_small)
  )
  

# set the colours for these insecticides
colours_plot <- insecticide_colours()[insecticides_plot_small]

# plot these, then add cell data over the top
preds_plot %>%
  ggplot(
    aes(
      x = year,
      y = Mortality,
      group = place,
      fill = insecticide_type
    )
  ) +
  # betabinomial posterior predictive intervals on bioassay data (captures
  # sample size and non-independence effect)
  geom_ribbon(
    aes(
      ymax = betabinomial_lower,
      ymin = betabinomial_upper,
    ),
    colour = grey(0.4),
    linewidth = 0.25,
    linetype = 2,
    alpha = 0.1
  ) +
  # binomial posterior predictive intervals on bioassay data (captures sample
  # size but assumes independence, so underestimates variance)
  geom_ribbon(
    aes(
      ymax = binomial_lower,
      ymin = binomial_upper,
    ),
    alpha = 0.4
  ) +
  # credible intervals on population-level proportion (our best guess at the
  # 'truth')
  geom_ribbon(
    aes(
      ymax = pop_lower,
      ymin = pop_upper,
    ),
  ) +
  geom_point(
    aes(
      x = year,
      y = Mortality,
      group = place,
      size = mosquito_number
    ),
    shape = 21,
    data = points_plot,
  ) +
  scale_fill_manual(
    values = colours_plot,
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_size_binned(
    range = c(1, 4),
    n.breaks = 6
  ) +
  facet_grid(insecticide_type ~ place) +
  # decades, so the labels of neighbouring panels don't run together
  scale_x_continuous(breaks = c(2000, 2010, 2020)) +
  guides(
    size = guide_legend(title = "No. tested")
  ) +
  labs(x = NULL, y = "Mortality") +
  coord_cartesian(
    xlim = c(1998, 2024),
    ylim = c(0, 1)
  ) +
  theme_minimal() +
  ggtitle(
    "Modelled population-level mortality and bioassay results at heavily sampled locations"
  )

ggsave("figures/fit_subset.png",
       bg = "white",
       scale = 0.8,
       # 18 for the four pixels this showed before
       width = 4.5 + 3.4 * nrow(locations_plot),
       height = 7)
