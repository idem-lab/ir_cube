# summarise model: the dynamical model's selection effect sizes and its
# response to net use at three exemplar places, and the fit of the two-stage
# model (R/two_stage_predictions.R) to the data

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

# the two-stage model's predictions, after the fit is loaded, so that the
# current code replaces the functions saved with it
source("R/two_stage_predictions.R")
years_predict <- baseline_year:final_data_year
cell_country <- data_cell_country(df)

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# summarise the covariate effect sizes (at insecticide class level)
effect_sizes <- summary(calculate(exp(beta_class[, 1]), values = draws))$statistics[, c("Mean", "SD")]
rownames(effect_sizes) <- colnames(x_cell_years)
round(effect_sizes, 2)

# posterior predictive draws of the mortality of new assays at the observed
# ones: for the two-stage model, with fresh pixel-year, pixel and study-year
# noise, at the external per-type overdispersion; and the dynamical model's, at
# its own
sims <- assay_draws(df, types)

# get RMSE and MAE for observed data and posterior mean within-sample
# predictions, of both models
obs <- df$died / df$mosquito_number
error <- function(model) colMeans(sims[[model]]) - obs
within_sample_error <- tibble(
  model = c("two_stage", "dynamical"),
  rmse = sapply(model, function(m) sqrt(mean(error(m) ^ 2))),
  mae = sapply(model, function(m) mean(abs(error(m))))
)
within_sample_error

# get posterior predictive simulations of observations, from the two-stage
# model
set.seed(2024)
died <- ppd_simulate(df$mosquito_number, sims$two_stage, sims$rho_two_stage)

# create a dharma object to compute randomised quantile residuals and
# corresponding residual z scores
dharma <- DHARMa::createDHARMa(
  simulatedResponse = t(died),
  observedResponse = df$died,
  integerResponse = TRUE
)

# some statistical evidence of misspecification, but it's a very large sample
# size and a very small deviation
plot(dharma)

png("figures/internal_validation_dharma.png",
    width = 10,
    height = 5,
    units = "in",
    res = 300)
plot(dharma)
dev.off()

# looks relatively uniform for each insecticide
tibble(
  scaled_residual = dharma$scaledResiduals,
  insecticide_type = factor(df$insecticide_type,
                            levels = insecticides_plot_order)
) %>%
  ggplot(
    aes(
      x = scaled_residual,
      fill = insecticide_type
    )
  ) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.02)
  ) +
  geom_hline(
    aes(yintercept = expected),
    # the count per bin expected if the residuals are uniform
    data = function(d) {
      d %>%
        count(insecticide_type) %>%
        mutate(expected = n / 50)
    },
    linetype = 2
  ) +
  scale_fill_manual(
    values = insecticide_colours(),
    guide = "none"
  ) +
  facet_wrap(~insecticide_type,
             scales = "free_y") +
  xlab("Scaled (quantile) residual") +
  ylab("Number of bioassays") +
  theme_minimal()

ggsave("figures/internal_validation_dharma_histograms.png",
       bg = "white",
       width = 9,
       height = 7)

# DDT the most obviously skewed
not_ddt_index <- df$insecticide_type != "DDT"
par(mfrow = c(2, 1))
hist(dharma$scaledResiduals,
     breaks = 100,
     main = "all")
hist(dharma$scaledResiduals[not_ddt_index],
     breaks = 100,
     main = "not DDT")

# is there excessive dispersion in DDT and underdispersion in Deltamethrin?
# make rho vary by insecticide class with a hierarchy in logit_rho?

dharma_not_ddt <- DHARMa::createDHARMa(
  simulatedResponse = t(died[, not_ddt_index]),
  observedResponse = df$died[not_ddt_index],
  integerResponse = TRUE
)

par(mfrow = c(1, 1))
plot(dharma_not_ddt)

# no evidence of temporal variation in model misspecification
z <- qnorm(dharma$scaledResiduals)
plot(z ~ jitter(df$year),
     cex = 0.5)
abline(h = 0)

# plot predicted trends and numbers of nets at a few unique locations

# look up the net coverage for different cells, find those with low, medium, and
# high net coverage
set.seed(2)
year_ids_keep <- which(years >= 2010 & years <= 2024)
exemplar_cells <- all_extract %>%
  filter(
    year_id %in% year_ids_keep,
    irs == 0,
    pop < 0.01,
    `all crops` < 0.01,
    `cereal crops` < 0.01,
    `root crops` < 0.01,
  ) %>%
  group_by(cell_id) %>%
  summarise(
    net_use = mean(nets),
    .groups = "drop"
  ) %>%
  # then pull the quartile threshold values
  mutate(
    which = case_when(
      net_use > 0 & net_use < 0.1 ~ "low",
      net_use >= 0.2 & net_use < 0.4 ~ "medium",
      net_use >= 0.5 ~ "high",
      .default = NA
    )
  ) %>%
  filter(
    !is.na(which)
  ) %>%
  group_by(
    which
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(
    which = factor(which, levels = c("low", "medium", "high"))
  ) %>%
  # find the coordinates and reverse geocode them
  mutate(
    cell = unique_cells[cell_id],
  ) %>%
  bind_cols(
    xyFromCell(mask, .$cell)
  ) %>%
  reverse_geocode(
    lat = y,
    long = x,
    method = 'osm',
    full_results = TRUE,
    # place names in English, not the local language
    custom_query = list("accept-language" = "en")
  )

# find a short name for these places
place_lookup <- exemplar_cells %>%
  mutate(
    precise_place = case_when(
      str_detect(address, "Dano") ~ "Dano",
      str_detect(address, "Tadjoura") ~ "Tadjoura",
      str_detect(address, "Ali Sabieh") ~ "Ali Sabieh",
      .default = str_split_i(address, ",", 1),
    ),
    place = paste(precise_place, country, sep = ", ")
  ) %>%
  select(
    cell_id,
    which,
    place)

place_order <- place_lookup %>%
  arrange(which) %>%
  pull(place)

data_plot <- all_extract %>%
  filter(
    cell_id %in% place_lookup$cell_id
  ) %>%
  left_join(place_lookup,
            by = "cell_id") %>%
  mutate(
    year = baseline_year + year_id - 1,
    place = factor(place,
                   levels = place_order)
  )

pred_lookup <- data_plot %>%
  select(
    cell_id, which, place
  ) %>%
  group_by(which) %>%
  slice(1) %>%
  ungroup()

# pull the IR and ITN timeseries for these, and plot

# do predictions of the effective resistance to ITNS
ingredient_weights <- readRDS("temporary/ingredient_weights.RDS")

# the predicted susceptibility to each LLIN insecticide at these cells,
# combined draw by draw with the ingredient weights: draws x cells x years.
# These cells were chosen to show the response to net use alone, which is the
# dynamical model's, so this figure shows the dynamical model, not the
# two-stage model (whose correction is not driven by the covariates)
exemplar_cell_ids <- pred_lookup$cell_id
effective_susc <- Reduce(`+`, lapply(names(ingredient_weights), function(type) {
  setup <- two_stage_setup(type, years_predict, df)
  ingredient_weights[[type]] *
    two_stage_weighted_draws(setup,
                             unique_cells[exemplar_cell_ids],
                             cell_country[exemplar_cell_ids],
                             diag(length(exemplar_cell_ids)))$dynamical
}))

preds <- expand_grid(
  year_id = seq_along(years_predict),
  cell_id = exemplar_cell_ids
) %>%
  mutate(
    post_mean = as.vector(apply(effective_susc, 2:3, mean)),
    post_lower = as.vector(apply(effective_susc, 2:3, quantile, 0.025)),
    post_upper = as.vector(apply(effective_susc, 2:3, quantile, 0.975))
  )

ir_plot <- data_plot %>%
  left_join(
    preds,
    by = c("cell_id", "year_id")
  ) %>%
  ggplot(
    aes(
      x = year,
      y = post_mean,
      ymax = post_upper,
      ymin = post_lower,
      group = place
    )
  ) +
  facet_wrap(~place) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  ylab("Susceptibility to LLIN insecticides") +
  xlab("") +
  geom_ribbon(fill = "#56B1F7") +
  theme_minimal()

itns_plot <- data_plot %>%
  ggplot(
    aes(
      x = year,
      y = nets,
      group = place
    )
  ) +
  facet_wrap(~place) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  ylab("LLIN use") +
  xlab("") +
  geom_line(
    colour = grey(0.5),
    linewidth = 1.2
  ) + 
  theme_minimal() +
  # suppress facet labels for bottom row
  theme(
    strip.text.x = element_blank()
  )

ir_plot / itns_plot

ggsave("figures/exemplar_itn_susc.png",
       bg = "white",
       width = 8,
       height = 5)


insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

# find some locations with lots of bioassay data for the pyrethroids and plot
# predictions and data for these
locations_plot <- df %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  group_by(cell) %>%
  filter(
    n() >= 30,
    n_distinct(year) >= 8,
    n_distinct(insecticide_type) >= 3
  ) %>%
  select(
    country_name,
    cell_id,
    cell
  ) %>%
  distinct() %>%
  bind_cols(
    xyFromCell(mask, .$cell)
  ) %>%
  rename(
    longitude = x,
    latitude = y
  ) %>%
  # geocode then find more interpretable names (google to see if this is how
  # they are referred to in IR papers)
  reverse_geocode(
    lat = latitude,
    long = longitude,
    method = 'osm',
    full_results = TRUE,
    # place names in English, not the local language
    custom_query = list("accept-language" = "en")
  ) %>%
  mutate(
    precise_place = str_split_i(address, ",", 1),
    # tidy up some of these where possible
    precise_place = case_when(
      grepl("Tiassalé", address) ~ "Tiassalé",
      grepl("Homa Bay", address) ~ "Homa Bay",
      grepl("Garoua", address) ~ "Garoua",
      grepl("Bandiagara", address) ~ "Bandiagara",
      grepl("Yaoundé", address) ~ "Yaoundé",
      grepl("Pitoa", address) ~ "Pitoa",
      grepl("Houet", address) ~ "Houet",
      grepl("Dakar", address) ~ "Dakar",
      grepl("Cotonou", address) ~ "Cotonou",
      grepl("Kéréwane", address) ~ "Kéréwane, Kolda",
      grepl("Soumousso", address) ~ "Soumousso",
      grepl("Busia", address) ~ "Busia",
      .default = precise_place
    ),
    place = paste(precise_place, country_name, sep = ", ")
  ) %>%
  select(
    address,
    place,
    latitude,
    longitude,
    cell_id
  )

# predict these with the two-stage model: the population-level susceptibility,
# ilogit(m + omega + xi), and the mortality in new bioassays of 100 mosquitoes
# under binomial sampling from it, and under the model's own observation
# process: beta-binomial sampling, at the per-type overdispersion, from the
# mortality with fresh pixel-year, pixel and study-year noise, ilogit(m +
# omega + xi + u + p + sy). Draws x cells x years
sample_size_plot <- 100
location_cell_ids <- locations_plot$cell_id
preds_plot <- bind_rows(
  lapply(insecticides_plot_small, function(type) {
    setup <- two_stage_setup(type, years_predict, df)
    identity <- diag(length(location_cell_ids))
    cells <- unique_cells[location_cell_ids]
    country <- cell_country[location_cell_ids]
    population <- two_stage_weighted_draws(setup, cells, country,
                                           identity)$two_stage
    assay <- two_stage_weighted_draws(setup, cells, country, identity,
                                      noise = TRUE)$two_stage
    rho <- rho_for_record(tibble(insecticide_type = type), rho_lookup())
    binomial_mortality <- array(
      rbinom(length(population), sample_size_plot, population),
      dim(population)) / sample_size_plot
    betabinomial_mortality <- array(
      rbetabinom(length(assay), sample_size_plot, assay, rho),
      dim(assay)) / sample_size_plot
    # a quantile over the draws, for every cell and year
    quantiles <- function(x, prob) as.vector(apply(x, 2:3, quantile, prob))

    # posterior mean mortality rate, and intervals for the population value
    # (posterior uncertainty), and posterior predictive intervals (posterior
    # uncertainty and sampling error) for observed mortality from binomial
    # sampling, and from the model's observation process
    expand_grid(
      year = years_predict,
      cell_id = location_cell_ids
    ) %>%
      mutate(
        insecticide_type = type,
        Susceptibility = as.vector(apply(population, 2:3, mean)),
        pop_lower = quantiles(population, 0.025),
        pop_upper = quantiles(population, 0.975),
        binomial_lower = quantiles(binomial_mortality, 0.025),
        binomial_upper = quantiles(binomial_mortality, 0.975),
        betabinomial_lower = quantiles(betabinomial_mortality, 0.025),
        betabinomial_upper = quantiles(betabinomial_mortality, 0.975)
      )
  })
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
    Susceptibility = died / mosquito_number
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
      y = Susceptibility,
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
      y = Susceptibility,
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
  xlab("") +
  coord_cartesian(
    xlim = c(1998, 2024),
    ylim = c(0, 1)
  ) +
  theme_minimal() +
  ggtitle(
    "Modelled population-level susceptibility and bioassay results at heavily sampled locations"
  )

ggsave("figures/fit_subset.png",
       bg = "white",
       scale = 0.8,
       width = 18,
       height = 7)

# these are the only sites with at least 30 pyrethroid bioassay datapoints,
# spanning at least 8 years, and with data for all three major pyrethroids

