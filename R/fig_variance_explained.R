# The main cross-validation figure (#36): should someone trust bioassay-based
# estimates or the model, with local data, without it, and for change?
#
#   top row, level: % of variance explained in single held-out bioassays
#     A  in sample       the interpolation test assays, full-data fit
#     B  interpolation   the same assays, held out
#     C  extrapolation   the held-out block assays
#   bottom row, change: skill against no change, pairs of single bioassays
#     D  in sample       the forecast pairs, full-data fit
#     E  forecast        the same pairs, 2014 and 2018 origins pooled
#
# Bioassay-based bars first (local, best K, nearest; green), then the models
# (dynamical, two-stage). The bar encoding and the key are in
# R/fig_variance_bars.R; the scores in R/validation_level.R (A-C) and
# R/validation_change.R (D, E); definitions, n* and the distances in
# doc/cv_figure_captions.md.
#
# Also the per-insecticide breakdown of the level scores
# (outputs/cv_variance_explained_by_insecticide.csv, R/variance_explained.R),
# in the same encoding.
source("R/packages.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/fig_variance_bars.R")

level <- read.csv("outputs/cv_level_skill.csv")
change <- read.csv("outputs/cv_change_skill.csv")

panel_order <- c("A", "B", "C", "D", "E")
bar_order <- c("local", "best_k", "nearest", "dynamical", "two_stage")
# the gap between the bioassay-based and the model-based bars, in bar widths
group_gap <- 0.35
placeholder_note <- "not yet\navailable"

km <- function(panel) {
  round(unique(level$median_km[level$panel == panel & !is.na(level$median_km)]))
}
panel_titles <- c(
  A = "A  in sample",
  B = sprintf("B  interpolation, median %i km", km("B")),
  C = sprintf("C  extrapolation, median %i km", km("C")),
  D = "D  in sample",
  E = "E  forecast")

# one row per bar, with its panel's ceiling, laid out within its panel
bars <- bind_rows(level, change) %>%
  group_by(panel) %>%
  mutate(ceiling = estimate[bar == "ceiling"]) %>%
  ungroup() %>%
  filter(bar %in% bar_order) %>%
  mutate(placeholder = source == "placeholder" | is.na(estimate),
         note = placeholder_note,
         panel = factor(panel, levels = panel_order),
         bar = factor(bar, levels = bar_order)) %>%
  arrange(panel, bar) %>%
  group_by(panel) %>%
  # each panel on its own stretch of the axis, so that one axis can label
  # every panel's bars
  mutate(position = 10 * as.integer(panel) + row_number() +
           group_gap * (as.character(bar) %in% c("dynamical", "two_stage")) *
           any(as.character(bar) %in% bioassay_bars)) %>%
  ungroup() %>%
  mutate(bar = as.character(bar),
         title = factor(panel_titles[as.character(panel)],
                        levels = panel_titles))
stopifnot(!anyNA(bars$ceiling))

# a panel's width along the axis, in bar positions, so that bars are the same
# width in every panel and in both rows
expand_x <- 0.12
panel_span <- bars %>%
  group_by(panel) %>%
  summarise(span = max(position) - min(position) + 0.8 + 2 * expand_x,
            .groups = "drop")
span <- setNames(panel_span$span, panel_span$panel)

skill_row <- function(data, y_title, y_limits, y_step = 20) {
  ggplot(data) +
    skill_bar_layers(data) +
    geom_hline(yintercept = 0, colour = grey(0.3), linewidth = 0.3) +
    facet_grid(~ title, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = data$position,
                       labels = function(x) {
                         unname(bar_labels[data$bar[match(x, data$position)]])
                       },
                       expand = expansion(add = 0.4 + expand_x)) +
    scale_y_continuous(breaks = seq(y_step * ceiling(y_limits[1] / y_step),
                                    100, y_step),
                       expand = expansion(add = 0)) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = y_title) +
    base_theme +
    theme(panel.spacing.x = unit(14, "pt"))
}

# change scores can fall below zero: the axis runs down to the next 10 below
# the lowest interval
change_floor <- 10 * floor(min(c(bars$lower[bars$panel %in% c("D", "E")], 0),
                               na.rm = TRUE) / 10)

level_row <- skill_row(filter(bars, panel %in% c("A", "B", "C")),
                       "level: % of variance explained", c(0, 100))
change_row <- skill_row(filter(bars, panel %in% c("D", "E")),
                        "change: skill against no change (%)",
                        c(change_floor, 100))

# The key takes the space to the right of D and E, so that bars keep one width
# across the rows: the bottom panels take their share of the top row's span
top_span <- sum(span[c("A", "B", "C")])
bottom_span <- sum(span[c("D", "E")])
main <- level_row /
  (change_row + skill_key() +
     plot_layout(widths = c(bottom_span + 0.3, top_span - bottom_span - 0.3))) +
  plot_layout(heights = c(100, 100 - change_floor))

ggsave("figures/CV_variance_explained.png", main, width = 11,
       height = 3.6 + 3.6 * (100 - change_floor) / 100, dpi = 300,
       bg = "white")


# per insecticide, level only ------------------------------------------------------

# One panel per experiment, the insecticides grouped by class as in
# fig_ir_maps.R. Only cells where the comparison can be read are drawn (the
# `shown` rule in variance_explained.R: ceiling >= 50% and 95% intervals
# narrower than 100 points); the rest, and why, are in the csv.
by_insecticide_all <- read.csv("outputs/cv_variance_explained_by_insecticide.csv")
by_insecticide <- by_insecticide_all %>% filter(shown)

experiment_order <- c("spatial interpolation", "spatial extrapolation",
                      "temporal change")
insecticide_class_order <- c(
  "Alpha-cypermethrin", "Deltamethrin", "Lambda-cyhalothrin", "Permethrin",
  "Fenitrothion", "Malathion", "Pirimiphos-methyl",
  "DDT", "Bendiocarb")
present <- unique(by_insecticide$stratum)
stopifnot(setequal(present, intersect(insecticide_class_order, present)))
insecticide_levels <- intersect(insecticide_class_order, present)
insecticide_bars <- c("nearest recent survey" = "nearest",
                      "dynamical model" = "dynamical",
                      "two-stage model" = "two_stage")

insecticide_data <- by_insecticide %>%
  group_by(experiment, stratum) %>%
  mutate(ceiling = 100 - estimate[kind == "noise"]) %>%
  ungroup() %>%
  filter(kind == "model", quantity %in% names(insecticide_bars)) %>%
  mutate(bar = unname(insecticide_bars[quantity]),
         offset = match(bar, insecticide_bars),
         position = match(stratum, insecticide_levels) + (offset - 2) * 0.28,
         placeholder = FALSE, note = NA_character_,
         experiment = factor(experiment, levels = experiment_order)) %>%
  # A bar whose estimate is below zero gets no colour: on this axis its
  # interval would show above zero with no rule in it, as though the estimate
  # lay inside it. Plain grey says what is meant, that it explains none of the
  # variance; the values are in the csv
  mutate(across(c(lower, upper),
                ~ ifelse(estimate > 0, .x, 0)),
         estimate = pmax(estimate, 0))

by_type_figure <- ggplot(insecticide_data) +
  skill_bar_layers(insecticide_data, width = 0.24) +
  facet_wrap(~ experiment, ncol = 1) +
  scale_x_continuous(breaks = seq_along(insecticide_levels),
                     labels = gsub("-", "-\n", insecticide_levels),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     expand = expansion(add = 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "% of variance explained") +
  base_theme +
  theme(axis.text.x = element_text(size = 7.5),
        panel.spacing.y = unit(14, "pt"))

# a key for the bar colours as well as the parts, since the bars are not
# labelled one by one
colour_key <- ggplot(tibble(bar = unname(insecticide_bars),
                            x = seq_along(insecticide_bars))) +
  geom_tile(aes(x = x, y = 1, fill = bar), width = 0.12, height = 0.6) +
  geom_text(aes(x = x + 0.1, y = 1, label = bar_labels[bar]), hjust = 0,
            size = 3.3) +
  scale_fill_manual(values = bar_colours, guide = "none") +
  coord_cartesian(xlim = c(0.8, 4), ylim = c(0, 2)) +
  theme_void()

by_type_with_key <- by_type_figure /
  (colour_key + skill_key(swatch_bar = "dynamical") +
     plot_layout(widths = c(1, 1.3))) +
  plot_layout(heights = c(10, 2))

ggsave("figures/CV_variance_explained_by_insecticide.png", by_type_with_key,
       width = 10, height = 10.5, dpi = 300, bg = "white", limitsize = FALSE)

shown_cells <- by_insecticide_all %>% distinct(experiment, stratum, shown)
cat("\nper-insecticide panels: ", sum(shown_cells$shown), " cells shown of ",
    nrow(shown_cells), "\n", sep = "")

saveRDS(list(main = main, by_insecticide = by_type_with_key, bars = bars,
             by_insecticide_data = by_insecticide_all),
        "outputs/figure_variance_explained.RDS")

cat("written:\n",
    " figures/CV_variance_explained.png\n",
    " figures/CV_variance_explained_by_insecticide.png\n", sep = "")
