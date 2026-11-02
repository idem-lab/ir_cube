# The main cross-validation figure (#36): should someone trust bioassay-based
# estimates or the model, with local data, without it, and ahead of the data?
# % of variance explained in single held-out bioassays:
#
#   top row, in space
#     A  in sample       the interpolation test assays, full-data fit
#     B  interpolation   the same assays, held out
#     C  extrapolation   the held-out block assays
#   bottom row, in time
#     D  forecast        the assays of 2014-2018, from the data before 2014
#     E  forecast        the assays of 2018-2022, from the data before 2018
#
# Bioassay-based bars first (local in A, nearest; green), then the models
# (dynamical, two-stage). The bar encoding and the key are in
# R/fig_variance_bars.R; the scores in R/validation_level.R; definitions, n*
# and the distances in doc/cv_figure_captions.md.
#
# Also the pooled in-sample figure: level (A's assays) and change (skill
# against no change, for pairs of single bioassays at one pixel,
# R/validation_change.R), predicted by the full-data fit and by 1-3
# independent bioassays or pairs at the location.
#
# Also the per-insecticide breakdown of B and C
# (outputs/cv_level_skill_by_insecticide.csv, R/validation_level.R), in the
# same encoding.
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
bar_order <- c("local", "nearest", "dynamical", "two_stage")

km <- function(panel) {
  round(unique(level$median_km[level$panel == panel & !is.na(level$median_km)]))
}
# Two-line titles, the panel's name in bold and, below it, how far its
# held-out bioassays sit from the data: the median distance (and, for the
# forecasts, years) to the nearest bioassay the nearest bioassay bar draws on
years <- function(panel) {
  unique(level$median_years[level$panel == panel &
                              !is.na(level$median_years)])
}
title <- function(name, detail) sprintf("**%s**<br>%s", name, detail)
distance <- function(panel) {
  sprintf("Median %i km from nearest bioassay", km(panel))
}
distance_years <- function(panel) {
  sprintf("Median %i km %g years from nearest bioassay", km(panel),
          years(panel))
}
panel_titles <- c(
  A = title("A) In sample", "Prediction at location with a bioassay"),
  B = title("B) Interpolation", distance("B")),
  C = title("C) Extrapolation", distance("C")),
  D = title("D) Forecast 2014-2018", distance_years("D")),
  E = title("E) Forecast 2018-2022", distance_years("E")))

# one row per bar, with its panel's ceiling, laid out within its panel
bars <- level %>%
  group_by(panel) %>%
  mutate(ceiling = estimate[bar == "ceiling"]) %>%
  ungroup() %>%
  filter(bar %in% bar_order) %>%
  mutate(panel = factor(panel, levels = panel_order),
         bar = factor(bar, levels = bar_order)) %>%
  arrange(panel, bar) %>%
  group_by(panel) %>%
  # each panel on its own stretch of the axis, so that one axis can label
  # every panel's bars
  mutate(position = 10 * as.integer(panel) + row_number()) %>%
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
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(x = NULL, y = y_title) +
    base_theme +
    # titles in markdown (bold first line), and may run past a narrow panel
    theme(panel.spacing.x = unit(14, "pt"), strip.clip = "off",
          strip.text = element_markdown(face = "plain", hjust = 0, size = 10,
                                        lineheight = 1.15,
                                        margin = margin(b = 4, l = 0)))
}

space_panels <- c("A", "B", "C")
time_panels <- c("D", "E")

# scores can fall below zero: the axis runs down to the next 10 below the
# lowest interval
axis_floor <- function(lower) {
  10 * floor(min(c(lower, 0), na.rm = TRUE) / 10)
}
level_floor <- axis_floor(bars$lower[bars$panel %in%
                                       c(space_panels, time_panels)])
level_title <- "% of variance explained"

space_row <- skill_row(filter(bars, panel %in% space_panels), level_title,
                       c(level_floor, 100))
time_row <- skill_row(filter(bars, panel %in% time_panels), level_title,
                      c(level_floor, 100))

# The key labels the parts of E's full-model bar, as a plot of its own on the
# same y axis, in the space to the right of D and E: set apart from the bars,
# as it applies to all of them. So that bars keep one width across the rows,
# the bottom panels take their share of the top row's span
key_plot <- function(reference, y_limits, gap = 1.6) {
  ggplot() +
    region_key(reference, x = 0.05, label_gap = 0.03, gap = gap) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 1), expand = expansion(add = 0)) +
    scale_y_continuous(expand = expansion(add = 0)) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_void()
}
top_span <- sum(span[space_panels])
bottom_span <- sum(span[time_panels])
main <- space_row /
  (time_row +
     key_plot(filter(bars, panel == "E", bar == "two_stage"),
              c(level_floor, 100)) +
     plot_layout(widths = c(bottom_span, top_span - bottom_span)))

ggsave("figures/CV_variance_explained.png", main, width = 12,
       height = 7.6 * (100 - level_floor) / 100, dpi = 300, bg = "white")


# in sample, pooled ------------------------------------------------------------------

# The in-sample bars of level (panel A's assays) and change (pairs of single
# bioassays, R/validation_change.R) with the local bar for 1-3 independent
# bioassays, or pairs, at the location, pooled: lighter to darker with more
# data, up to the main figure's local bioassay colour. After
# the main figure, so its key is a reminder, drawn against the full-model bar
# with the parts' lines as close to their extents as can be told apart
pool_sizes <- 1:3
local_bars <- c("local", paste0("local_", pool_sizes[-1]))
pool_colours <- colorRampPalette(c(bar_colours[["nearest"]],
                                  bar_colours[["local"]]))(length(pool_sizes))
pool_units <- c(level = "bioassay", change = "pair")
for (kind in names(pool_units)) {
  names_k <- paste0(kind, "_", pool_sizes)
  bar_colours[names_k] <- pool_colours
  light_bars <- c(light_bars, names_k[1])
  bar_labels[names_k] <- sprintf("%i\n%s%s", pool_sizes, pool_units[[kind]],
                                 ifelse(pool_sizes == 1, "", "s"))
}

pooled_panel <- function(data, kind, panel_title) {
  order <- c(paste0(kind, "_", pool_sizes), "dynamical", "two_stage")
  data %>%
    mutate(ceiling = estimate[bar == "ceiling"]) %>%
    filter(bar %in% c(local_bars, "dynamical", "two_stage")) %>%
    mutate(bar = ifelse(bar %in% local_bars,
                        paste0(kind, "_", match(bar, local_bars)), bar)) %>%
    arrange(match(bar, order)) %>%
    mutate(position = row_number(), title = panel_title)
}
pooled_level <- pooled_panel(
  filter(level, panel == "A"), "level",
  title("A) Level, in sample", "Independent bioassays at the location"))
pooled_change <- pooled_panel(
  filter(change, panel == "in_sample"), "change",
  title("B) Change, in sample", sprintf(paste(
    "Independent pairs of bioassays at the location, median %g years",
    "apart"), change$median_gap[change$panel == "in_sample"][1])))
stopifnot(nrow(pooled_level) == 5, nrow(pooled_change) == 5)

pooled_row <- function(data, y_title) {
  limits <- c(axis_floor(data$lower), 100)
  skill_row(data, y_title, limits) +
    key_plot(filter(data, bar == "two_stage"), limits, gap = 0.4) +
    plot_layout(widths = c(1, 0.45))
}
pooled_floors <- sapply(list(pooled_level, pooled_change), function(data) {
  axis_floor(data$lower)
})
pooled_figure <- pooled_row(pooled_level, level_title) /
  pooled_row(pooled_change, "skill against no change (%)") +
  plot_layout(heights = 100 - pooled_floors)

ggsave("figures/CV_in_sample_pooled.png", pooled_figure, width = 7.5,
       height = 3.6 * sum(100 - pooled_floors) / 100, dpi = 300, bg = "white")


# per insecticide, B and C -----------------------------------------------------------

# One panel per experiment, the insecticides grouped by class as in
# fig_ir_maps.R. Only cells where the comparison can be read are drawn (the
# `shown` rule in validation_level.R: ceiling >= 50% and 95% intervals
# narrower than 100 points); the rest, and why, are in the csv.
by_insecticide_all <- read.csv("outputs/cv_level_skill_by_insecticide.csv")
by_insecticide <- by_insecticide_all %>% filter(shown)

experiment_labels <- c(B = "B) Interpolation", C = "C) Extrapolation")
insecticide_class_order <- c(
  "Alpha-cypermethrin", "Deltamethrin", "Lambda-cyhalothrin", "Permethrin",
  "Fenitrothion", "Malathion", "Pirimiphos-methyl",
  "DDT", "Bendiocarb")
present <- unique(by_insecticide$insecticide_type)
stopifnot(setequal(present, intersect(insecticide_class_order, present)))
insecticide_levels <- intersect(insecticide_class_order, present)
insecticide_bars <- c("nearest", "dynamical", "two_stage")

insecticide_data <- by_insecticide %>%
  filter(bar %in% insecticide_bars) %>%
  mutate(offset = match(bar, insecticide_bars),
         position = match(insecticide_type, insecticide_levels) +
           (offset - 2) * 0.28,
         experiment = factor(experiment_labels[panel],
                             levels = experiment_labels)) %>%
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
colour_key <- ggplot(tibble(bar = insecticide_bars,
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
  plot_layout(heights = c(7, 2))

ggsave("figures/CV_variance_explained_by_insecticide.png", by_type_with_key,
       width = 10, height = 7.7, dpi = 300, bg = "white", limitsize = FALSE)

shown_cells <- by_insecticide_all %>% distinct(panel, insecticide_type, shown)
cat("\nper-insecticide panels: ", sum(shown_cells$shown), " cells shown of ",
    nrow(shown_cells), "\n", sep = "")

saveRDS(list(main = main, in_sample_pooled = pooled_figure,
             by_insecticide = by_type_with_key, bars = bars,
             by_insecticide_data = by_insecticide_all),
        "outputs/figure_variance_explained.RDS")

cat("written:\n",
    " figures/CV_variance_explained.png\n",
    " figures/CV_in_sample_pooled.png\n",
    " figures/CV_variance_explained_by_insecticide.png\n", sep = "")
