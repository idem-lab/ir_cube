# Supplementary figures for the two-stage model (#21): how the final model's
# components build a prediction at a pixel, for Deltamethrin.
#
#   figures/two_stage/supp_components_realisations.png (and .pdf)
#     joint random realisations of every component through time, stacked on a
#     shared time axis: net use, the dynamical model's logit m, the annual
#     anomalies eta, their accumulation xi, the static omega, the target
#     m + omega + xi, and mortality with the data
#   figures/two_stage/supp_components_intervals.png (and .pdf)
#     posterior mean and 95% intervals of mortality: the dynamical model alone
#     against the two-stage model, with the data
#
#   Rscript R/fig_two_stage_components.R
#
# Run after R/two_stage_maps.R (prepare, and the type's fit), whose saved fit
# and paired dynamical draws this reads; nothing is refitted. Run with
# OpenBLAS, e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=3 nice -n 10 Rscript ...
#
# The target, on the logit scale, is m(s, t) + omega(s) + xi(s, t): m the
# dynamical model's prediction, xi(s, t) the sum of the annual anomalies eta up
# to t, eta AR(1) in time with Matern innovations; beyond the last data year T
# eta is simulated forward. u and p are observation-level noise and appear in
# neither figure.

type <- "Deltamethrin"
figure_dir <- "figures/two_stage"
maps_dir <- file.path("outputs/two_stage/maps", type)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

source("R/two_stage_helpers.R")
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/dynamical_predictions.R")
source("R/two_stage_correction.R")
source("R/two_stage_map_functions.R")

end_year <- 2030
years_all <- baseline_year:end_year
k <- match(type, types)
rows_k <- which(df$type_id == k)

# the fit and the 2000 paired dynamical draws saved by R/two_stage_maps.R
fit <- readRDS(file.path(maps_dir, "fit.rds"))
dynamical <- readRDS(file.path(maps_dir, "dynamical.rds"))
train_k <- tibble(lon = df$longitude[rows_k], lat = df$latitude[rows_k],
                  year = df$year_start[rows_k], cell = df$cell[rows_k],
                  died = df$died[rows_k],
                  mosquito_number = df$mosquito_number[rows_k])
stopifnot(fit$n_obs == nrow(train_k),
          ncol(dynamical$logit_train) == nrow(train_k),
          max(abs(fit$m_ref - colMeans(dynamical$logit_train))) < 1e-12)
T_k <- fit$T


# 2. example pixels ---------------------------------------------------------------

# The data per pixel, and the fitted smooth correction at the mode
# (omega + xi at T) at each sampled pixel
coords_train <- coords_km(train_k)
pixel_summary <- train_k %>%
  mutate(x_km = coords_train[, 1], y_km = coords_train[, 2]) %>%
  group_by(cell) %>%
  summarise(lon = mean(lon), lat = mean(lat),
            x_km = mean(x_km), y_km = mean(y_km),
            n_assays = n(), n_years = n_distinct(year),
            first = min(year), last = max(year),
            n_before = sum(year <= 2008), n_after = sum(year >= 2014),
            .groups = "drop")
pixel_coords <- as.matrix(pixel_summary[, c("x_km", "y_km")])
pixel_summary$smooth_T <- as.vector(project_correction(
  fit, correction_node_draws(fit, T_k, mean = TRUE),
  mutate(pixel_summary, year = T_k)))

# (i) well sampled: data both before the net scale-up (to 2008) and after it
# (from 2014), and the most distinct years of data among such pixels
well <- pixel_summary %>%
  filter(n_before > 0, n_after > 0) %>%
  arrange(desc(n_years), desc(n_assays)) %>%
  slice(1)

# (iii) large correction: among pixels with at least three years of data, the
# largest |omega + xi| at T at the mode, i.e. where the data pull the
# prediction furthest from the dynamical model
large <- pixel_summary %>%
  filter(n_years >= 3, cell != well$cell) %>%
  arrange(desc(abs(smooth_T))) %>%
  slice(1)

# (ii) no data: a cell of the prediction mask inside the limits of Pf
# transmission, 400-600 km from the nearest Deltamethrin assay, in a country
# with Deltamethrin data: far beyond the range of omega (~43 km), so omega and
# p revert to their priors and the prediction to the dynamical model, but
# within the reach of the large-scale xi (range ~1000 km). The candidate with
# the highest net use in 2020 is taken, so that the dynamical model has a
# trajectory to follow
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")
country_raster <- rast("data/clean/country_raster.tif")
mask_cells <- terra::cells(terra::mask(mask, pf_water_mask))
set.seed(1)
candidate_cells <- sort(sample(mask_cells, min(length(mask_cells), 40000)))
xy_candidates <- terra::xyFromCell(mask, candidate_cells)
coords_candidates <- project_km(xy_candidates[, 1], xy_candidates[, 2])
nearest_km <- apply(coords_candidates, 1, function(xy) {
  sqrt(min((pixel_coords[, 1] - xy[1]) ^ 2 + (pixel_coords[, 2] - xy[2]) ^ 2))
})
candidate_country <- as.character(
  terra::extract(country_raster, candidate_cells)$country_name)
data_countries <- unique(df$country_name[rows_k])
nets_2020 <- terra::extract(rast("data/clean/net_use_cube.tif")[["nets_2020"]],
                            candidate_cells)[, 1]
ok <- nearest_km > 400 & nearest_km < 600 &
  candidate_country %in% data_countries & !is.na(nets_2020)
none_cell <- candidate_cells[ok][which.max(nets_2020[ok])]
none <- tibble(cell = none_cell,
               lon = xy_candidates[match(none_cell, candidate_cells), 1],
               lat = xy_candidates[match(none_cell, candidate_cells), 2],
               n_assays = 0L, n_years = 0L,
               nearest_km = nearest_km[match(none_cell, candidate_cells)])

pixels <- bind_rows(
  mutate(well, role = "well"),
  mutate(none, role = "none"),
  mutate(large, role = "large")
) %>%
  mutate(country = as.character(
    terra::extract(country_raster, cell)$country_name))
report("pixels: %s", paste(sprintf("%s cell %i (%s, %.2f, %.2f; %i assays, %i years)",
                                  pixels$role, pixels$cell, pixels$country,
                                  pixels$lon, pixels$lat, pixels$n_assays,
                                  pixels$n_years), collapse = "; "))


# 3. draws at the pixels ------------------------------------------------------------

# Every component at each pixel and year 1995-2030, for the 2000 paired draws:
# draw d pairs dynamical draw d with a latent draw from N(mode, H^-1) shifted by
# the cut-posterior formula for it, and the AR(1) forecast of eta beyond T
# (correction_node_draws(), project_correction())
n_draws <- nrow(dynamical$logit_train)
n_pix <- nrow(pixels)
n_years_all <- length(years_all)

# dynamical logit at the pixels, as a years x pixels x draws array
covariates <- map_covariates(pixels$cell, baseline_year, end_year,
                             dynamical$parameters$options$selection_columns)
pixel_country <- match(pixels$country, dimnames(dynamical$logit_init)[[2]])
stopifnot(!anyNA(pixel_country))
m <- dynamical_logit_cells(dynamical$parameters, k,
                           matrix(dynamical$logit_init[, pixel_country, k],
                                  n_draws),
                           map_x(covariates, seq_len(n_pix), n_years_all),
                           seq_len(n_years_all), x_init = covariates$init)
m <- aperm(simplify2array(m), c(3, 2, 1))

# check: at the sampled pixels' data years, the same draws as at the assays
for (i in which(pixels$n_assays > 0)) {
  rows_i <- which(train_k$cell == pixels$cell[i])
  j <- match(train_k$year[rows_i], years_all)
  difference <- max(abs(t(matrix(m[j, i, ], length(j))) -
                          dynamical$logit_train[, rows_i]))
  report("pixel %i: max |m - m at the assays| = %.1e", pixels$cell[i],
         difference)
  stopifnot(difference < 1e-6)
}

# omega + xi at every pixel-year, and omega alone (xi is 0 in year t0)
pixel_xy <- terra::xyFromCell(mask, pixels$cell)
coords_pix <- project_km(pixel_xy[, 1], pixel_xy[, 2])
new <- tibble(x_km = rep(coords_pix[, 1], each = n_years_all),
              y_km = rep(coords_pix[, 2], each = n_years_all),
              year = rep(years_all, n_pix))
new_omega <- tibble(x_km = coords_pix[, 1], y_km = coords_pix[, 2],
                    year = fit$t0)
correction <- array(NA_real_, c(n_years_all, n_pix, n_draws))
omega <- matrix(NA_real_, n_pix, n_draws)
set.seed(4026 + k)
for (batch in split(seq_len(n_draws), ceiling(seq_len(n_draws) / 250))) {
  fields <- correction_node_draws(
    fit, years_all, length(batch),
    m_draws_train = dynamical$logit_train[batch, , drop = FALSE])
  correction[, , batch] <- project_correction(fit, fields, new)
  omega[, batch] <- project_correction(fit, fields, new_omega)
}
omega_array <- array(rep(omega, each = n_years_all), dim(correction))
xi <- correction - omega_array
eta <- xi - xi[c(1, seq_len(n_years_all - 1)), , , drop = FALSE]
lambda <- m + correction
nets <- covariates$time_varying[, , "nets"]
rm(dynamical, covariates, correction, omega_array)


# 4. figures ------------------------------------------------------------------------

years <- years_all

# panel titles, in the order of `pixels`
role_title <- c(well = "A) Well sampled",
                none = "B) No data",
                large = "C) Large correction")
pixels <- pixels %>%
  mutate(title = sprintf("%s: %s\n%s; %s",
                         role_title[role], country,
                         sprintf("%.1f°%s, %.1f°%s",
                                 abs(lat), ifelse(lat >= 0, "N", "S"),
                                 abs(lon), ifelse(lon >= 0, "E", "W")),
                         ifelse(n_assays > 0,
                                sprintf("%i assays, %i years", n_assays,
                                        n_years),
                                sprintf("nearest assay %.0f km",
                                        nearest_km))),
         title = factor(title, levels = title))

# the style of the dynamical-model time-series figures (R/summarise_model_fit.R,
# R/fig_temporal_preds_net_use.R): theme_minimal, the pyrethroid blue, net use
# as a thick grey line, percentages, no x label
pyrethroid_blue <- "#56B1F7"
net_grey <- grey(0.5)
realisation_cols <- c("#E69F00", "#009E73", "#CC79A7")
year_breaks <- seq(2000, 2030, by = 10)

projection_shading <- list(
  annotate("rect", xmin = T_k + 0.5, xmax = 2030.5, ymin = -Inf, ymax = Inf,
           fill = grey(0.93)),
  geom_vline(xintercept = T_k + 0.5, colour = grey(0.5), linewidth = 0.3,
             linetype = "dashed")
)
base_theme <- theme_minimal(base_size = 9) +
  theme(strip.text.x = element_text(hjust = 0, size = 8.5),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 8.5),
        legend.position = "none",
        plot.margin = margin(2, 4, 2, 4))
x_scale <- scale_x_continuous(breaks = year_breaks, limits = c(1994.5, 2030.5),
                              expand = c(0, 0))

# long data frames: one row per pixel x year (x draw)
pixel_year <- function(values, name) {
  tibble(title = rep(pixels$title, each = length(years)),
         year = rep(years, nrow(pixels)),
         !!name := as.vector(values))
}
# realisation draws: years x pixels x draws array (or pixels x draws matrix
# for the static terms) at the chosen draws
realisations <- function(values, which) {
  if (length(dim(values)) == 2) {
    values <- array(rep(values, each = length(years)),
                    c(length(years), dim(values)))
  }
  bind_rows(lapply(seq_along(which), function(r) {
    pixel_year(values[, , which[r]], "value") %>%
      mutate(realisation = factor(r))
  }))
}

bioassays <- train_k %>%
  filter(cell %in% pixels$cell) %>%
  mutate(title = pixels$title[match(cell, pixels$cell)],
         mortality = died / mosquito_number)

# 4a. realisations ------------------------------------------------------------------

# three joint draws, spread over the posterior of the dynamical model: draw d
# pairs the dynamical draw d with the correction's latent draw d
set.seed(5)
which_draws <- sort(sample(n_draws, 3))

row_plot <- function(data, ylab, geoms, y_scale = NULL, strip = FALSE,
                     bottom = FALSE) {
  plot <- ggplot(data, aes(x = year)) +
    projection_shading +
    geoms +
    facet_wrap(~title, nrow = 1) +
    x_scale +
    ylab(ylab) +
    xlab(NULL) +
    base_theme
  if (!is.null(y_scale)) plot <- plot + y_scale
  if (!strip) plot <- plot + theme(strip.text.x = element_blank())
  if (!bottom) plot <- plot + theme(axis.text.x = element_blank())
  plot
}
realisation_colour <- scale_colour_manual(values = realisation_cols)
zero_line <- geom_hline(yintercept = 0, colour = grey(0.6), linewidth = 0.3)

m_mean <- pixel_year(apply(m, 1:2, mean), "value")

p_nets <- row_plot(
  pixel_year(t(nets), "value"), "LLIN use",
  geom_line(aes(y = value), colour = net_grey, linewidth = 1),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  strip = TRUE)
p_m <- row_plot(
  realisations(m, which_draws), "dynamical\nlogit m",
  list(geom_line(aes(y = value), data = m_mean, colour = "black",
                 linewidth = 0.6),
       geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour))
p_eta <- row_plot(
  realisations(eta, which_draws), "annual\nanomaly η",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.3,
                 alpha = 0.6),
       geom_point(aes(y = value, colour = realisation), size = 0.6),
       realisation_colour))
p_xi <- row_plot(
  realisations(xi, which_draws), "cumulative\nξ = Ση",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.5),
       realisation_colour))
p_omega <- row_plot(
  realisations(omega, which_draws), "static\nfield ω",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.5),
       realisation_colour))
p_lambda <- row_plot(
  realisations(lambda, which_draws), "logit\nm+ω+ξ",
  list(geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour))
p_mort <- row_plot(
  realisations(plogis(lambda), which_draws), "mortality",
  list(geom_point(aes(y = mortality, size = mosquito_number),
                  data = bioassays, shape = 21, fill = grey(0.8),
                  colour = grey(0.35), stroke = 0.2, alpha = 0.7),
       geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour,
       scale_size_area(max_size = 2.5)),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  bottom = TRUE)

fig_realisations <- patchwork::wrap_plots(p_nets, p_m, p_eta, p_xi, p_omega,
                                          p_lambda, p_mort, ncol = 1,
                                          heights = c(1, 1.1, 1, 1.1, 0.8, 1.1,
                                                      1.2)) &
  theme(plot.margin = margin(1, 4, 1, 4))
for (ext in c("png", "pdf")) {
  ggsave(file.path(figure_dir, sprintf("supp_components_realisations.%s", ext)),
         plot = fig_realisations,
         bg = "white", width = 7.5, height = 8.4, dpi = 300,
         device = if (ext == "pdf") cairo_pdf else NULL)
}


# 4b. intervals ---------------------------------------------------------------------

# mortality intervals for the population fraction at each pixel-year: the
# dynamical model alone (plogis(m)) and the two-stage model
# (plogis(m + omega + xi), with the cut-posterior shift). u and p are
# observation noise, like the assay-level (beta-binomial) noise, so none is
# included
summarise_draws <- function(values, model) {
  pixel_year(apply(values, 1:2, mean), "mean") %>%
    mutate(lower = as.vector(apply(values, 1:2, quantile, 0.025)),
           upper = as.vector(apply(values, 1:2, quantile, 0.975)),
           model = model)
}
intervals <- bind_rows(
  summarise_draws(plogis(m), "dynamical model"),
  summarise_draws(plogis(lambda), "two-stage model")
) %>%
  mutate(model = factor(model, c("dynamical model", "two-stage model")))
model_cols <- c(`dynamical model` = grey(0.45),
                `two-stage model` = pyrethroid_blue)

p_int <- ggplot(intervals, aes(x = year)) +
  projection_shading +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.35) +
  geom_line(aes(y = mean, colour = model), linewidth = 0.6) +
  geom_point(aes(y = mortality, size = mosquito_number), data = bioassays,
             shape = 21, fill = "white", colour = "black", stroke = 0.3,
             alpha = 0.8) +
  facet_wrap(~title, nrow = 1) +
  scale_fill_manual(values = model_cols, name = NULL) +
  scale_colour_manual(values = model_cols, name = NULL) +
  scale_size_area(max_size = 2.5, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  x_scale +
  xlab(NULL) +
  ylab("Susceptibility to deltamethrin") +
  base_theme +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.text.x = element_blank())
p_int_nets <- row_plot(
  pixel_year(t(nets), "value"), "LLIN use",
  geom_line(aes(y = value), colour = net_grey, linewidth = 1),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  bottom = TRUE)

fig_intervals <- patchwork::wrap_plots(p_int, p_int_nets, ncol = 1, heights = c(3, 1))
for (ext in c("png", "pdf")) {
  ggsave(file.path(figure_dir, sprintf("supp_components_intervals.%s", ext)),
         plot = fig_intervals,
         bg = "white", width = 8, height = 4.5, dpi = 300,
         device = if (ext == "pdf") cairo_pdf else NULL)
}
report("figures written")
