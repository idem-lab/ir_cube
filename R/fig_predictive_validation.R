# Figures summarising out-of-sample posterior predictive performance (#10).
#
# Reads only the tables written by validation_metrics.R, so figures can be
# redrawn without refitting or rescoring anything.
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

# The nearest neighbour null is reported twice: as the practice baseline a
# person would actually apply (one neighbour, the most recent two available
# years), and as an oracle bound at whichever neighbour count minimises its own
# error on the held-out records - hindsight the dynamical model is not given.
model_labels <- c(dynamical = "dynamical model",
                  two_stage = "two-stage model",
                  nearest_neighbour = "nearest neighbour",
                  nearest_neighbour_oracle = "nearest neighbour (best k)",
                  intercept = "insecticide mean")
model_colours <- c("dynamical model" = "#2166AC",
                   "two-stage model" = "#C51B7D",
                   "nearest neighbour" = "#B2182B",
                   "nearest neighbour (best k)" = "#E08214",
                   "insecticide mean" = grey(0.55))

# The forecasting experiment is one per origin, since pooling a 2014 forecast
# with a 2018 one would average over holdout windows whose true rates of decline
# differ by a factor of two. Origins are labelled by their cut year and ordered
# after the spatial experiments, earliest first.
experiment_labels <- c(spatial_interpolation = "spatial interpolation",
                       spatial_blocks = "spatial extrapolation")

label_experiments <- function(experiment) {
  cut_year <- sub("^temporal_forecasting_", "", experiment)
  is_forecast <- cut_year != experiment
  out <- ifelse(is_forecast,
                paste("forecast from", cut_year),
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

coverage_plot <- coverage %>%
  tidy_labels() %>%
  ggplot(
    aes(x = nominal,
        y = empirical,
        colour = model)
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "stated chance of covering the result",
    y = "how often it did",
    subtitle = "Are the uncertainty ranges the right width? A well calibrated model follows the dashed line"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

reliability_plot <- reliability %>%
  tidy_labels() %>%
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
        ymax = predicted + ppc_upper),
    fill = grey(0.85),
    colour = NA,
    alpha = 0.6
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_point(size = 1.6) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "") +
  coord_equal() +
  guides(colour = "none") +
  labs(
    x = "predicted mortality (map)",
    y = "observed mortality",
    subtitle = "When the model predicts a mortality, is that the average outcome? Grey band is the scatter a correct model would still show"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

main_figure <- coverage_plot / reliability_plot +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Out-of-sample predictive performance",
    subtitle = "held-out bioassays, three cross-validation experiments",
    tag_levels = "a"
  ) &
  theme(legend.position = "bottom")

ggsave("figures/CV_predictive_calibration.png",
       main_figure,
       bg = "white",
       width = 11,
       height = 8)


# the table ----------------------------------------------------------------

# one row per experiment and model, with the column glosses that belong in the
# caption rather than the header. MSE and bias score the map (#32). The
# two-stage model's floor (u + p) is the part of its excess it attributes to
# its site noise, and its coverage without u and p is NA until its folds carry
# map draws
table_out <- summaries %>%
  tidy_labels() %>%
  transmute(
    experiment,
    model,
    `held-out bioassays` = n,
    `95% interval coverage` = round(coverage_95, 3),
    `95% coverage without u, p` = round(coverage_95_map, 3),
    `mean PIT` = round(mean_pit, 3),
    `CRPS (mortality)` = round(crps, 4),
    `MSE` = round(mse, 4),
    `bioassay noise floor` = round(mse_floor, 4),
    `MSE above the floor` = round(excess, 4),
    `model floor (u + p)` = round(model_floor, 4),
    `bias` = round(bias, 3),
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
  geom_point(size = 3, colour = model_colours[["dynamical model"]]) +
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

