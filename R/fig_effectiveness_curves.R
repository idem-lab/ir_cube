# The relative effectiveness of a pyrethroid-only net against bioassay
# mortality under the two models of the epi impact figure
# (R/fig_ento_epi_impact.R): Symons et al., and Churcher et al. 2024's MINT
# simulations, the median over the settings and each setting. With the
# break-even effectiveness at the mean delivery cost, and each model's
# break-even susceptibility k*. For the supplement, alongside
# figures/epi_impact.png and figures/epi_impact_churcher.png

source("R/packages.R")
source("R/net_effectiveness.R")

# the two models in black and grey, not the net type colours of the epi
# impact figure; the settings in light grey
model_labels <- c(symons = "Symons et al.",
                  churcher = "Churcher et al., median",
                  setting = "Churcher et al., settings")
model_cols <- c(symons = "black",
                churcher = grey(0.45),
                setting = grey(0.85))

curves <- bind_rows(effectiveness_points, .id = "model")

effectiveness_curves <- ggplot(
  mapping = aes(
    x = susceptibility,
    y = effectiveness,
    colour = model
  )
) +
  geom_line(
    data = mutate(mint_settings, model = "setting"),
    aes(group = setting),
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = break_even,
    colour = grey(0.4),
    linetype = 2
  ) +
  annotate(
    "text",
    x = 0,
    y = break_even,
    label = "Break-even",
    hjust = 0,
    vjust = -0.6,
    size = 3,
    colour = grey(0.4)
  ) +
  geom_line(
    data = curves,
    linewidth = 0.9
  ) +
  # the simulated levels of resistance
  geom_point(
    data = filter(curves, model == "churcher"),
    size = 1.5
  ) +
  # each model's break-even susceptibility, marked down to the axis
  annotate(
    "segment",
    x = break_even_susceptibility,
    xend = break_even_susceptibility,
    y = 0,
    yend = break_even,
    linetype = 3,
    colour = model_cols[names(break_even_susceptibility)]
  ) +
  annotate(
    "point",
    x = break_even_susceptibility,
    y = break_even,
    shape = 21,
    fill = "white",
    size = 2
  ) +
  annotate(
    "text",
    x = break_even_susceptibility,
    y = 0,
    label = sprintf("%.1f%%", 100 * break_even_susceptibility),
    hjust = if_else(names(break_even_susceptibility) == "symons", -0.1, 1.1),
    vjust = 0,
    size = 3,
    colour = model_cols[names(break_even_susceptibility)]
  ) +
  scale_x_continuous(
    labels = scales::percent,
    name = "Bioassay mortality"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = 0.02),
    name = "LLIN effectiveness"
  ) +
  # capped at 100%: the few settings simulated above it are clipped
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  scale_colour_manual(
    values = model_cols,
    labels = model_labels,
    breaks = names(model_labels),
    name = NULL
  ) +
  theme_minimal() +
  # the legend below, clear of the curves and the break-even marks
  theme(
    legend.position = "bottom"
  )

# drawn at three quarters of the size, as the epi impact figure
ggsave(
  filename = "figures/effectiveness_curves.png",
  plot = effectiveness_curves,
  bg = "white",
  width = 7,
  height = 5,
  scale = 0.75,
  dpi = 400
)
