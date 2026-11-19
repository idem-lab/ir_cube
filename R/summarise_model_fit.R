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
