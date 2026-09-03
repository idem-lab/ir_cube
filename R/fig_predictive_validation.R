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
who_thresholds <- read.csv("outputs/cv_who_thresholds.csv")

model_labels <- c(dynamical = "dynamical model",
                  nearest_neighbour = "nearest neighbour",
                  intercept = "insecticide mean")
model_colours <- c("dynamical model" = "#2166AC",
                   "nearest neighbour" = "#B2182B",
                   "insecticide mean" = grey(0.55))

experiment_labels <- c(spatial_interpolation = "spatial interpolation",
                       spatial_extrapolation = "spatial extrapolation",
                       temporal_forecasting = "temporal forecasting")

tidy_labels <- function(data) {
  data %>%
    mutate(
      model = factor(model_labels[model], levels = model_labels),
      experiment = factor(experiment_labels[experiment],
                          levels = experiment_labels)
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
  # the scatter a perfect model would still show, given bioassay noise and the
  # number of assays in each bin
  geom_ribbon(
    aes(ymin = predicted - 1.96 * envelope,
        ymax = predicted + 1.96 * envelope),
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
    x = "predicted mortality",
    y = "observed mortality",
    subtitle = "When the model predicts a mortality, is that the average outcome? Grey band is the scatter bioassay noise alone would produce"
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
# caption rather than the header
table_out <- summaries %>%
  tidy_labels() %>%
  transmute(
    experiment,
    model,
    `held-out bioassays` = n,
    `95% interval coverage` = round(coverage_95, 3),
    `mean PIT` = round(mean_pit, 3),
    `CRPS (mortality)` = round(crps, 4),
    `skill vs nearest neighbour` = round(skill, 3),
    `Cramer-von Mises` = round(cvm, 2),
    `CvM null upper` = round(cvm_null_upper, 2)
  ) %>%
  arrange(experiment, model)

write.csv(table_out, "outputs/cv_table.csv", row.names = FALSE)
print(as.data.frame(table_out))


# supplementary figures ----------------------------------------------------

# departure of the randomised PIT values from uniformity, which is the same
# check as the coverage curve at full resolution
pit_ecdf <- scores %>%
  tidy_labels() %>%
  group_by(model, experiment) %>%
  arrange(pit, .by_group = TRUE) %>%
  mutate(
    empirical = row_number() / n(),
    difference = empirical - pit
  ) %>%
  ungroup()

# a simultaneous band, from the same statistic under uniform draws
band_width <- scores %>%
  count(model, experiment) %>%
  summarise(width = 1.36 / sqrt(min(n))) %>%
  pull(width)

pit_plot <- pit_ecdf %>%
  ggplot(
    aes(x = pit,
        y = difference,
        colour = model)
  ) +
  annotate("rect", xmin = 0, xmax = 1,
           ymin = -band_width, ymax = band_width,
           fill = grey(0.9)) +
  geom_hline(yintercept = 0, linetype = 2, colour = grey(0.6)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "") +
  labs(
    x = "probability integral transform",
    y = "empirical minus expected",
    title = "Departure of held-out data from the predictive distribution",
    subtitle = "grey band is the range expected by chance for a calibrated model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/CV_pit.png", pit_plot, bg = "white", width = 11, height = 4)

# fitted overdispersion against the external, replicate-based estimate. A
# fitted value above the external one means the model is absorbing process
# misfit into the observation process
rho_plot <- rho_comparison %>%
  tidy_labels() %>%
  ggplot(
    aes(x = rho_external,
        y = rho_fitted,
        colour = model,
        shape = insecticide_class)
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = grey(0.6)) +
  geom_point(size = 3) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "") +
  scale_shape_discrete(name = "") +
  coord_equal() +
  labs(
    x = "overdispersion estimated from replicate bioassays",
    y = "overdispersion fitted by the model",
    title = "Is the model treating its own error as bioassay noise?",
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
                      pixel_year = "pooled within pixel, year and insecticide",
                      country_year = "pooled within country, year and insecticide")
  ) %>%
  ggplot(
    aes(x = mean_assays,
        y = rmse,
        colour = model)
  ) +
  geom_point(size = 3) +
  facet_grid(grouping ~ experiment) +
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
       width = 11, height = 6)

# the same comparison in WHO resistance categories
who_plot <- who_thresholds %>%
  tidy_labels() %>%
  group_by(model, experiment, category) %>%
  summarise(across(c(observed, predicted, lower, upper), mean),
            .groups = "drop") %>%
  mutate(
    category = factor(category,
                      levels = c("confirmed", "possible", "susceptible"),
                      labels = c("confirmed resistance\n(< 90% mortality)",
                                 "possible resistance\n(90 - 98%)",
                                 "susceptible\n(> 98%)"))
  ) %>%
  ggplot(
    aes(x = category,
        colour = model)
  ) +
  geom_linerange(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.5),
    linewidth = 1
  ) +
  geom_point(
    aes(y = observed),
    position = position_dodge(width = 0.5),
    shape = 4,
    size = 3,
    colour = "black"
  ) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_colour_manual(values = model_colours, name = "predicted") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "",
    y = "proportion of held-out bioassays",
    title = "Does the model reproduce the observed pattern of resistance?",
    subtitle = "crosses are the observed proportions, bars the 95% predictive range"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/CV_who_thresholds.png", who_plot, bg = "white",
       width = 11, height = 5)
