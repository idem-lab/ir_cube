# figures to validate internal consistency of model predictions, along important
# gradients

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/validation_functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

mask <- rast("data/clean/raster_mask.tif")

# draw the posterior predictive distribution at each observation. The
# predicted fraction and the overdispersion are drawn in one call so that each
# pair comes from the same posterior sample
rho_observations <- rho_classes[df$class_id]
sims <- calculate(population_mortality_vec,
                  rho_observations,
                  values = draws,
                  nsim = 1e3)

# randomised quantile residuals, computed from the analytic beta-binomial
# mixture rather than by simulation, so the residuals carry no Monte Carlo
# noise and the out-of-sample residuals in validation_metrics.R are computed
# the same way (#10)
ppd <- ppd_summary(df$died,
                   df$mosquito_number,
                   sims$population_mortality_vec[, , 1],
                   sims$rho_observations[, , 1])

# scale to a normal distribution for easier checking
df_validate <- df %>%
  mutate(
    z_resid = pit_to_z(rowMeans(ppd_pit(ppd, n_rep = 100))),
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

# plot and fit models comparing this with covariate values

df_validate %>%
  ggplot(
    aes(
      x = longitude,
      y = latitude,
      colour = z_resid
    )
  ) +
  geom_spatraster(
    data = mask
  ) +
  geom_point(
    alpha = 0.5,
    size = 2
  ) +
  scale_fill_gradient(
    low = grey(0.9),
    high = grey(0.9),
    na.value = "transparent",
    guide = "none"
  ) +
  scale_colour_viridis_c() +
  theme_minimal() +
  theme_ir_maps() +
  ggtitle("Residuals")

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
  ) %>%
  mutate(
    z_resid_pred = predict(space_fit, ., se.fit = TRUE)$fit,
    z_resid_sd = predict(space_fit, ., se.fit = TRUE)$se.fit
  )

z_resid_raster <- c(mask_lores, mask_lores)
names(z_resid_raster) <- c("mean", "sd") 
z_resid_raster$mean[cells(z_resid_raster)] <- as.matrix(coords_pred$z_resid_pred)
z_resid_raster$sd[cells(z_resid_raster)] <- as.matrix(coords_pred$z_resid_sd)

ggplot() +
  geom_spatraster(
    aes(
      fill = mean,
    ),
    data = z_resid_raster
  ) +
  scale_fill_distiller(
    type = "div",
    na.value = "transparent",
    limits = c(-4, 4)
  ) +
  theme_ir_maps() +
  ggtitle(
    "Residual spatial correlation in model"
  )

# plot against covariates
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
      covariate_name != "net_coverage" ~ sqrt(covariate_value),
      .default = covariate_value
    )
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
      group = covariate_name
    )
  ) +
  geom_point(
    alpha = 0.1
  ) +
  geom_hline(yintercept = 0,
             colour = "red",
             linetype = 2) +
  geom_smooth() +
  facet_grid(insecticide_type ~ covariate_name,
             scales = "free_x") +
  theme_minimal()

# plot against time, in different regions
df_validate %>%
  group_by(country_name) %>%
  filter(n() >= 500) %>%
  group_by(
    country_name,
    insecticide_type
  ) %>%
  # remove some outlying data skewing the smooths
  filter(
    year_start <= quantile(year_start, 0.95),
    year_start >= quantile(year_start, 0.05)
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = year_start,
      y = z_resid
    )
  ) +
  geom_point(
    alpha = 0.1
  ) +
  geom_hline(yintercept = 0,
             colour = "red",
             linetype = 2) +
  geom_smooth() +
  facet_grid(country_name ~ insecticide_type,
             scales = "free_x") +
  theme_minimal()


# group distributions in different ways, and compute Kolmogorov-Smirnov
# D statistics (and p values) for each



df_validate %>%
  group_by(year_start, insecticide_class) %>%
  summarise(
    n = n(),
    D = ks.test(z_resid, pnorm)$statistic,
    p = ks.test(z_resid, pnorm)$p.value
  ) %>%
  arrange(desc(D))


df_validate %>%
  filter(insecticide_type != "DDT") %>%
  group_by(cluster, insecticide_class) %>%
  summarise(
    n = n(),
    D = ks.test(z_resid, pnorm)$statistic,
    p = ks.test(z_resid, pnorm)$p.value
  ) %>%
  arrange(desc(D))

df_validate %>%
  arrange(z_resid) %>%
  ggplot(
    aes(x = longitude,
        y = latitude,
        fill = abs(z_resid))#factor(cluster))# == 7)
  ) +
  geom_point(
    shape = 21
  ) +
  facet_wrap(~insecticide_type) +
  coord_equal() +
  theme_minimal()

