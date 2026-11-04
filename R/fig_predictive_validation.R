# Figures summarising out-of-sample posterior predictive performance (#10).
#
# Reads only the tables written by validation_metrics.R and
# in_sample_calibration.R, so figures can be redrawn without refitting or
# rescoring anything.
#
# The main figure leads with the two questions a reader can check without any
# distributional vocabulary: when the model said there was a 95% chance, how
# often was it right, and when it predicts 60% mortality, is the average
# outcome 60%. Scoring rules have no natural zero and so are left to the table
# and the supplement, where they are reported against the null model and the
# bioassay noise floor.

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

summaries <- read.csv("outputs/cv_summary.csv")
coverage <- read.csv("outputs/cv_coverage.csv")
reliability <- read.csv("outputs/cv_reliability.csv")
scores <- read.csv("outputs/cv_scores.csv")
aggregated <- read.csv("outputs/cv_aggregate_summary.csv")
rho_comparison <- read.csv("outputs/cv_rho_comparison.csv")
by_fold <- read.csv("outputs/cv_by_fold.csv", encoding = "UTF-8")

# The nearest neighbour null is the practice baseline a person would actually
# apply: one neighbour, the most recent two available years.
model_labels <- c(dynamical = "Dynamical model",
                  two_stage = "Two-stage model",
                  nearest_neighbour = "Nearest neighbour",
                  intercept = "Insecticide mean")
model_colours <- c("Dynamical model" = "#2166AC",
                   "Two-stage model" = "#C51B7D",
                   "Nearest neighbour" = "#B2182B",
                   "Insecticide mean" = grey(0.55))

# The forecasting experiment is one per origin, since pooling a 2014 forecast
# with a 2018 one would average over holdout windows whose true rates of decline
# differ by a factor of two. Origins are labelled by their cut year and ordered
# after the spatial experiments, earliest first.
experiment_labels <- c(spatial_interpolation = "Spatial interpolation",
                       spatial_blocks = "Spatial extrapolation")

label_experiments <- function(experiment) {
  cut_year <- sub("^temporal_forecasting_", "", experiment)
  is_forecast <- cut_year != experiment
  out <- ifelse(is_forecast,
                paste("Forecast from", cut_year),
                experiment_labels[experiment])
  forecast_levels <- unique(out[is_forecast])
  forecast_levels <- forecast_levels[order(forecast_levels)]
  factor(out, levels = c(unname(experiment_labels), forecast_levels))
}

tidy_labels <- function(data) {
  data %>%
    mutate(
      model = factor(model_labels[model], levels = model_labels),
      experiment = label_experiments(experiment)
    )
}


# main figure --------------------------------------------------------------

# The models the figure shows: the insecticide mean is left to the table. The
# two-stage model is drawn last and widest, so it is the one read first; the
# legend lists it first
figure_models <- c("Two-stage model", "Dynamical model", "Nearest neighbour")
figure_data <- function(data) {
  data %>%
    tidy_labels() %>%
    filter(model %in% figure_models) %>%
    mutate(model = factor(model, levels = rev(figure_models)))
}
model_widths <- c("Two-stage model" = 1.2, "Dynamical model" = 0.6,
                  "Nearest neighbour" = 0.6)

# in-sample coverage (R/in_sample_calibration.R), at the modelled bioassays
# of the full fit, beside each held-out curve
in_sample <- read.csv("outputs/in_sample_coverage.csv") %>%
  select(-experiment) %>%
  mutate(model = factor(model_labels[model], levels = rev(figure_models)))
coverage_lines <- bind_rows(
  figure_data(coverage) %>% mutate(bioassays = "Held out"),
  tidyr::crossing(experiment = unique(figure_data(coverage)$experiment),
                  in_sample) %>%
    mutate(bioassays = "In sample (full fit)")
)

coverage_plot <- coverage_lines %>%
  ggplot(
    aes(x = nominal,
        y = empirical,
        colour = model,
        linewidth = model,
        linetype = bioassays)
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_line() +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, breaks = figure_models,
                      name = "") +
  scale_linewidth_manual(values = model_widths, breaks = figure_models,
                         name = "") +
  scale_linetype_manual(values = c("Held out" = 1,
                                   "In sample (full fit)" = 3),
                        name = "") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Nominal coverage of the central predictive interval",
    y = "Observed coverage",
    subtitle = "Coverage of central predictive intervals"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

reliability_plot <- reliability %>%
  figure_data() %>%
  ggplot(
    aes(x = predicted,
        y = observed,
        colour = model)
  ) +
  # the scatter a perfect model would still show: the 95% range of the gap
  # between bin mean observation and bin mean prediction over replicate datasets
  # drawn from the model's own posterior predictive distribution, so bioassay
  # noise, overdispersion, posterior uncertainty and the bin size are all in it
  geom_ribbon(
    aes(ymin = predicted + ppc_lower,
        ymax = predicted + ppc_upper,
        group = model),
    fill = grey(0.85),
    colour = NA,
    alpha = 0.6
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_point(size = 1.6) +
  geom_line(aes(linewidth = model), alpha = 0.7) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, breaks = figure_models,
                      name = "") +
  scale_linewidth_manual(values = model_widths * 0.7, breaks = figure_models,
                         name = "") +
  coord_equal() +
  guides(colour = "none", linewidth = "none") +
  labs(
    x = "Predicted mortality (map), mean in each tenth of predictions",
    y = "Observed mortality, mean",
    subtitle = "Observed against predicted mortality of held-out bioassays"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

main_figure <- coverage_plot / reliability_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave("figures/CV_predictive_calibration.png",
       main_figure,
       bg = "white",
       width = 11,
       height = 8)


# the table ----------------------------------------------------------------

# one row per experiment and model, with the column glosses that belong in the
# caption rather than the header. MSE scores the map (#32)
table_out <- summaries %>%
  tidy_labels() %>%
  transmute(
    experiment,
    model,
    `held-out bioassays` = n,
    `95% interval coverage` = round(coverage_95, 3),
    `mean PIT` = round(mean_pit, 3),
    `CRPS (mortality)` = round(crps, 4),
    `MSE` = round(mse, 4),
    `bioassay noise floor` = round(mse_floor, 4),
    `MSE above the floor` = round(excess, 4),
    `RMS error in the fraction` = round(rms_p, 3),
    `Cramer-von Mises` = round(cvm, 2)
  ) %>%
  arrange(experiment, model)

write.csv(table_out, "outputs/cv_table.csv", row.names = FALSE)
print(as.data.frame(table_out))


# supplementary figures ----------------------------------------------------

# The dynamical model's fitted overdispersion against the external,
# replicate-based estimate. A fitted value above the external one means the
# model is absorbing process misfit into the observation process.
#
# The dynamical model only. The null models no longer fit an overdispersion of
# their own: they earn their place on point prediction - mean squared error and
# variance explained - and fitting one per null made the comparison a contest in
# vagueness rather than in prediction (#12 review). So this is a diagnostic of
# the one model, not a comparison across models.
rho_plot <- rho_comparison %>%
  filter(model == "dynamical") %>%
  tidy_labels() %>%
  ggplot(
    aes(x = rho_external,
        y = rho_fitted,
        shape = insecticide_class)
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_point(size = 3, colour = model_colours[["Dynamical model"]]) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_shape_discrete(name = "") +
  coord_equal() +
  labs(
    x = "overdispersion estimated from replicate bioassays",
    y = "overdispersion fitted by the dynamical model",
    title = "Is the dynamical model treating its own error as bioassay noise?",
    subtitle = "points above the line indicate process misfit absorbed into the observation model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/CV_rho_comparison.png", rho_plot, bg = "white",
       width = 11, height = 5)

# calibration of pooled groups of assays. A single bioassay is a noisy measure
# of the population fraction; pooling brings the comparison to bear on that
# quantity
aggregate_plot <- aggregated %>%
  tidy_labels() %>%
  mutate(
    grouping = recode(grouping,
                      country_year = "pooled within country, year and insecticide")
  ) %>%
  ggplot(
    aes(x = mean_assays,
        y = rmse,
        colour = model)
  ) +
  geom_point(size = 3) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "") +
  labs(
    x = "mean bioassays pooled per group",
    y = "root mean squared error of pooled mortality",
    title = "Predictive error against the population quantity",
    subtitle = "pooling assays reduces measurement noise but not model error"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/CV_aggregated.png", aggregate_plot, bg = "white",
       width = 11, height = 4)


# where the error sits: per fold, and against separation from the training data.
# The pooled numbers hide which folds carry the result, and the geometry panel
# is the practical question — at what separation from the data should the
# mechanistic model be preferred to local interpolation (#12 review)
fold_plot <- by_fold %>%
  filter(experiment == "spatial_blocks") %>%
  tidy_labels() %>%
  ggplot(
    aes(x = reorder(fold, excess),
        y = excess,
        fill = model)
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = model_colours, name = "") +
  labs(
    x = "",
    y = "mean squared error above the bioassay noise floor",
    title = "Which held-out block carries the extrapolation result?",
    subtitle = "lower is better; the insecticide mean is the no-information baseline"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/CV_by_fold.png", fold_plot, bg = "white",
       width = 9.5, height = 5)

