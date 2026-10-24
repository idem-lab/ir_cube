# figures to validate internal consistency of model predictions, along important
# gradients
#
#   Rscript R/fig_internal_validation.R [two_stage|dynamical]
#
# for the two-stage model (the default), or the dynamical model alone, whose
# figures and tables carry the suffix _dynamical

model <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(model)) model <- "two_stage"
stopifnot(model %in% c("two_stage", "dynamical"))
suffix <- if (model == "dynamical") "_dynamical" else ""

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/two_stage_predictions.R")

# the modelled data, as R/fit_model.R builds them
baseline_year <- 1995
invisible(list2env(modelled_bioassays(baseline_year), environment()))
final_data_year <- max(df$year)

# the covariates at the data cells, on their own scales (R/model_covariates.R)
all_extract <- covariate_extract(unique_cells, baseline_year, final_data_year,
                                 two_stage_design(types[1]))

mask <- rast("data/clean/raster_mask.tif")

# country borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.85),
                     inherit.aes = FALSE)
country_borders <- geom_sf(data = borders,
                           col = grey(0.4),
                           linewidth = 0.1,
                           fill = "transparent",
                           inherit.aes = FALSE)

insecticides_col <- insecticide_colours()

# draws of the posterior predictive distribution at each observation
# (R/two_stage_predictions.R): for the two-stage model, the mortality of a new
# assay, with fresh pixel-year and pixel noise, at the external per-type
# overdispersion; for the dynamical model, its predicted mortality at its own
# overdispersion, each paired with the same posterior draw
sims <- assay_draws(df, types)
rho_observations <- sims[[paste0("rho_", model)]]

# randomised quantile residuals, computed from the analytic beta-binomial
# mixture rather than by simulation, so the residuals carry no Monte Carlo
# noise and the out-of-sample residuals in validation_metrics.R are computed
# the same way (#10)
set.seed(2024)
ppd <- ppd_summary(df$died,
                   df$mosquito_number,
                   sims[[model]],
                   rho_observations)

# scale to a normal distribution for easier checking
df_validate <- df %>%
  mutate(
    # one randomisation replicate, not the mean of them: averaging converges
    # to the mid-P value, which is not uniform under calibration for discrete
    # data and so would show a spurious bulge here (#12 review)
    z_resid = pit_to_z(ppd_pit(ppd, n_rep = 1)[, 1]),
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_order)
  ) %>%
  # add on covariate values
  left_join(
    all_extract,
    by = c("cell_id", "year_id")
  ) %>%
  # make spatial clusters
  mutate(
    cluster = kmeans(
      x = as.matrix(select(., latitude, longitude)),
      centers = 15
    )$cluster
  )

# map the residuals, largest on top
df_validate %>%
  arrange(abs(z_resid)) %>%
  ggplot(
    aes(
      x = longitude,
      y = latitude,
      colour = z_resid
    )
  ) +
  africa_bg +
  country_borders +
  geom_point(
    size = 0.6
  ) +
  scale_colour_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Residual<br>z-score"
  ) +
  facet_wrap(~insecticide_type) +
  coord_sf(xlim = c(-18, 52), ylim = c(-35, 38)) +
  theme_ir_maps()

ggsave(sprintf("figures/internal_validation_residual_map%s.png", suffix),
       bg = "white",
       width = 9,
       height = 9)

space_fit <- df_validate %>%
  filter(insecticide_class == "Pyrethroids") %>%
  mgcv::gam(
    z_resid ~ s(latitude, longitude, bs = "gp", k = 200),
    method = "REML",
    data = .
  )

# make a raster of the GAM residual fit
mask_lores <- mask %>%
  terra::aggregate(10)
coords_pred <- mask_lores %>%
  terra::xyFromCell(cells(.)) %>%
  as_tibble() %>%
  rename(
    latitude = y,
    longitude = x
  )
space_pred <- predict(space_fit, coords_pred, se.fit = TRUE)

z_resid_raster <- c(mask_lores, mask_lores)
names(z_resid_raster) <- c("mean", "sd")
z_resid_raster$mean[cells(z_resid_raster)] <- as.matrix(space_pred$fit)
z_resid_raster$sd[cells(z_resid_raster)] <- as.matrix(space_pred$se.fit)

# drop the smooth where it extrapolates beyond the pyrethroid data: where its
# standard error exceeds that at 95% of data locations
data_pred <- predict(space_fit, se.fit = TRUE)
z_resid_raster$mean[z_resid_raster$sd > quantile(data_pred$se.fit, 0.95)] <- NA

# symmetric colour limits about zero, from the smooth at the data
smooth_limit <- max(abs(data_pred$fit))

ggplot() +
  geom_spatraster(
    aes(
      fill = mean,
    ),
    data = z_resid_raster
  ) +
  country_borders +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    na.value = "transparent",
    limits = c(-1, 1) * smooth_limit,
    oob = scales::squish,
    name = "Smoothed<br>residual<br>z-score"
  ) +
  coord_sf(xlim = c(-18, 52), ylim = c(-35, 38)) +
  theme_ir_maps()

ggsave(sprintf("figures/internal_validation_residual_smooth%s.png", suffix),
       bg = "white",
       width = 6,
       height = 6)

# plot against covariates, on a square-root scale except for net use
cov_names <- colnames(all_extract)[-(1:2)]
df_validate %>%
  select(
    all_of(cov_names),
    insecticide_type,
    z_resid
  ) %>%
  pivot_longer(
    cols = all_of(cov_names),
    names_to = "covariate_name",
    values_to = "covariate_value"
  ) %>%
  mutate(
    covariate_value = case_when(
      covariate_name != "nets" ~ sqrt(covariate_value),
      .default = covariate_value
    ),
    covariate_name = factor(covariate_name, levels = cov_names)
  ) %>%
  group_by(
    covariate_name,
    insecticide_type
  ) %>%
  # remove some outlying data skewing the smooths
  filter(
    covariate_value <= quantile(covariate_value, 0.95),
    covariate_value >= quantile(covariate_value, 0.05)
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = covariate_value,
      y = z_resid,
      colour = insecticide_type
    )
  ) +
  geom_point(
    alpha = 0.1,
    size = 0.3,
    colour = grey(0.4)
  ) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5)
  ) +
  scale_colour_manual(
    values = insecticides_col,
    guide = "none"
  ) +
  facet_grid(insecticide_type ~ covariate_name,
             scales = "free_x") +
  xlab("Covariate value (square root, except nets)") +
  ylab("Residual z-score") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 6)
  )

ggsave(sprintf("figures/internal_validation_residual_covariates%s.png",
               suffix),
       bg = "white",
       width = 18,
       height = 12)

# plot against time, in countries with many records
df_validate %>%
  group_by(country_name) %>%
  filter(n() >= 500) %>%
  group_by(
    country_name,
    insecticide_type
  ) %>%
  # remove some outlying data skewing the smooths
  filter(
    year <= quantile(year, 0.95),
    year >= quantile(year, 0.05)
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = year,
      y = z_resid,
      colour = insecticide_type
    )
  ) +
  geom_point(
    alpha = 0.2,
    size = 0.3,
    colour = grey(0.4)
  ) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  geom_smooth(
    method = "loess",
    formula = y ~ x
  ) +
  scale_colour_manual(
    values = insecticides_col,
    guide = "none"
  ) +
  facet_grid(country_name ~ insecticide_type) +
  xlab("") +
  ylab("Residual z-score") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1)
  )

ggsave(sprintf("figures/internal_validation_residual_year_country%s.png",
               suffix),
       bg = "white",
       width = 14,
       height = 14)


# group distributions in different ways, and compute Kolmogorov-Smirnov
# D statistics (and p values) for each
ks_summary <- function(z) {
  test <- ks.test(z, pnorm)
  tibble(n = length(z),
         D = unname(test$statistic),
         p = test$p.value)
}

ks_type <- df_validate %>%
  group_by(insecticide_class, insecticide_type) %>%
  reframe(ks_summary(z_resid)) %>%
  arrange(desc(D))

ks_year_class <- df_validate %>%
  group_by(year, insecticide_class) %>%
  reframe(ks_summary(z_resid)) %>%
  arrange(desc(D))

ks_cluster_class <- df_validate %>%
  filter(insecticide_type != "DDT") %>%
  group_by(cluster, insecticide_class) %>%
  reframe(ks_summary(z_resid)) %>%
  arrange(desc(D))

write_csv(ks_type,
          sprintf("outputs/internal_validation_ks_type%s.csv", suffix))
write_csv(ks_year_class,
          sprintf("outputs/internal_validation_ks_year_class%s.csv", suffix))
write_csv(ks_cluster_class,
          sprintf("outputs/internal_validation_ks_cluster_class%s.csv", suffix))

ks_type
ks_year_class
ks_cluster_class
