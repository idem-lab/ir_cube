# The bar encoding shared by the cross-validation skill figures (#36).
#
# Every bar has the same stack, from the bottom:
#
#   colour   the score: the share of the target predicted. Solid to the lower
#            bound of its 95% interval, translucent to the upper, with a rule
#            at the estimate. Where the lower bound is below zero (change
#            scores only), only the translucent interval and the rule are
#            drawn, so no solid stub reads as the estimate
#   grey     the remainder up to the ceiling: population mortality, or its
#            change, not predicted
#   dotted   the ceiling: the score of a perfect prediction of population
#            mortality, short of 100 by the target's own bioassay noise
#   white    the bioassay noise, above the ceiling: cannot be predicted
#
# Reading the interval as a translucent extension of the bar rather than as an
# error bar keeps the quantities on one additive scale, so the grey is always
# the shortfall.
#
# Sourcing this file defines the colours, skill_bar_layers(), skill_key() and
# base_theme; it draws nothing on its own.
suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# The bioassay-based bars take one green hue (OKLCH h = 150) at three
# saturations, lightest for the most local. The models keep the colours of the
# earlier figures: blue for the dynamical model, magenta for the two-stage
# model.
bar_colours <- c(
  local = "#B6DDBD",
  best_k = "#6FB07D",
  nearest = "#137738",
  dynamical = "#2166AC",
  two_stage = "#C51B7D")
bar_labels <- c(local = "local", best_k = "best K", nearest = "nearest",
                dynamical = "dynamical", two_stage = "two-stage")
bioassay_bars <- c("local", "best_k", "nearest")

noise_fill <- "#F7F7F7"
remainder_fill <- "#DEDEDE"
outline_colour <- grey(0.7)
ceiling_colour <- grey(0.25)
interval_alpha <- 0.45

# The rule at the estimate: white on every bar but the lightest, where it
# would vanish
rule_colour <- function(bar) {
  ifelse(bar == "local", bar_colours[["nearest"]], "white")
}

# One bar per row of `data`: position, bar (a name of bar_colours), estimate,
# lower, upper, ceiling, and placeholder (TRUE where the estimate is not yet
# available: drawn as a dashed outline up to the ceiling, labelled `note`)
skill_bar_layers <- function(data, width = 0.8, note_size = 2.6) {

  data <- data %>%
    mutate(xmin = position - width / 2, xmax = position + width / 2)
  scored <- data %>% filter(!placeholder) %>%
    mutate(solid_top = pmax(lower, 0),
           rule = rule_colour(bar))
  waiting <- data %>% filter(placeholder)

  list(
    # the noise, above the ceiling
    geom_rect(data = data,
              aes(xmin = xmin, xmax = xmax, ymin = ceiling, ymax = 100),
              fill = noise_fill, colour = NA),
    # the remainder up to the ceiling
    geom_rect(data = scored,
              aes(xmin = xmin, xmax = xmax,
                  ymin = pmin(pmax(upper, 0), ceiling), ymax = ceiling),
              fill = remainder_fill, colour = NA),
    geom_rect(data = waiting,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ceiling),
              fill = remainder_fill, colour = NA),
    # the score: its interval translucent, solid to its nearer bound
    geom_rect(data = scored,
              aes(xmin = xmin, xmax = xmax, ymin = lower, ymax = upper,
                  fill = bar),
              alpha = interval_alpha, colour = NA),
    geom_rect(data = scored,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = solid_top,
                  fill = bar),
              colour = NA),
    geom_segment(data = scored,
                 aes(x = xmin, xend = xmax, y = estimate, yend = estimate,
                     colour = rule),
                 linewidth = 0.7),
    # the ceiling
    geom_segment(data = data,
                 aes(x = xmin, xend = xmax, y = ceiling, yend = ceiling),
                 colour = ceiling_colour, linetype = "dotted",
                 linewidth = 0.45),
    # the outline last, over the full extent of the bar
    geom_rect(data = scored,
              aes(xmin = xmin, xmax = xmax, ymin = pmin(lower, 0),
                  ymax = 100),
              fill = NA, colour = outline_colour, linewidth = 0.3),
    geom_rect(data = waiting,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 100),
              fill = NA, colour = grey(0.4), linewidth = 0.4,
              linetype = "22"),
    geom_text(data = waiting,
              aes(x = position, y = ceiling / 2, label = note),
              angle = 90, size = note_size, colour = grey(0.25),
              lineheight = 0.9),
    scale_fill_manual(values = bar_colours, guide = "none"),
    scale_colour_identity()
  )
}

# The one key: the bar parts only, as a ggplot of its own to set beside or
# below the panels. `swatch_bar` gives the colour of the score's swatch
skill_key <- function(swatch_bar = "two_stage", text_size = 3.1) {
  swatch <- function(y, ...) {
    annotate("rect", xmin = 0, xmax = 0.45, ymin = y - 0.3, ymax = y + 0.3,
             ...)
  }
  label <- function(y, text) {
    annotate("text", x = 0.7, y = y, label = text, hjust = 0,
             size = text_size, lineheight = 0.9)
  }
  ggplot() +
    swatch(4, fill = noise_fill, colour = outline_colour, linewidth = 0.3) +
    label(4, "bioassay noise: cannot be predicted") +
    annotate("segment", x = 0, xend = 0.45, y = 3, yend = 3,
             colour = ceiling_colour, linetype = "dotted", linewidth = 0.55) +
    label(3, "population mortality predicted perfectly") +
    swatch(2, fill = remainder_fill, colour = outline_colour,
           linewidth = 0.3) +
    label(2, "population mortality not predicted") +
    annotate("rect", xmin = 0, xmax = 0.45, ymin = 0.7, ymax = 1,
             fill = bar_colours[[swatch_bar]]) +
    annotate("rect", xmin = 0, xmax = 0.45, ymin = 1, ymax = 1.3,
             fill = bar_colours[[swatch_bar]], alpha = interval_alpha) +
    annotate("segment", x = 0, xend = 0.45, y = 1, yend = 1,
             colour = "white", linewidth = 0.8) +
    label(1, "population mortality predicted\n(estimate and 95% interval)") +
    coord_cartesian(xlim = c(-0.3, 7), ylim = c(-0.6, 5.4), clip = "off") +
    theme_void()
}

base_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = grey(0.93)),
        axis.text.x = element_text(size = 7.8),
        legend.position = "none",
        strip.text = element_text(face = "bold", hjust = 0, size = 10,
                                  margin = margin(b = 4, l = 0)),
        plot.margin = margin(4, 6, 4, 6))
