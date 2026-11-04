# The trend in agricultural insecticide use intensity, from FAOSTAT: whether a
# linear trend common to all regions is adequate for the trend g_ag(t) on the
# crop-layer selection covariates (R/model_covariates.R), which is linear and
# the same in every region.
#
# Response: log kg insecticide (FAOSTAT item 1309, active ingredient) per ha
# cropland, 1995-2024, for the modelled countries. Cropland is implied by
# FAOSTAT's own pesticides total and pesticides total per ha of cropland.
# Imputed (I) values, which are carry-forwards and straight-line
# interpolations, are dropped, as are zeros.
#
# Two data sets: the values reported by the countries (flag A), the main
# analysis, and those with FAO's estimates (flags E and X) added. FAO's
# estimates for non-reporting countries are mostly net trade (imports less
# exports) converted to active ingredient with global factors (FAO Statistics
# Working Paper 24-45), and are much smoother than the reported values.
#
# Linear mixed models, fitted by maximum likelihood (nlme):
#
#   y_ct = alpha + b t + a_c + d_c t + e_ct
#   a_c, d_c   country intercept and slope deviations, independent normals
#   e_ct       continuous-time AR(1) within each country (corCAR1), with its
#              own variance for estimated values when they are included
#
# with t the year centred on 2009.5, in decades. Two tests, each of nested
# fixed effects with the same random effects and correlation, so a likelihood
# ratio test (and AIC) of the ML fits is valid:
#   non-linearity       + q t^2 against the linear model;
#   regional variation  + region x t against + region (regional intercepts in
#                       both, so only the slopes are tested).
#
#   Rscript R/faostat_trend_model.R
# Needs data/clean/faostat_pesticide_use.RDS from R/prep_faostat.R.

source("R/packages.R")
source("R/functions.R")
library(nlme)

first_year <- 1995
last_year <- 2024
centre_year <- (first_year + last_year) / 2

regions <- sort(unique(country_region_lookup()$region))

# data ------------------------------------------------------------------------

use <- readRDS("data/clean/faostat_pesticide_use.RDS") %>%
  filter(!aggregate, modelled | country_name == "Sudan (former)") %>%
  mutate(country_name = if_else(country_name == "Sudan (former)",
                                "Sudan", country_name))

# cropland (ha) implied by the pesticides total and its per-ha rate (kg/ha)
cropland <- use %>%
  filter(item_code == "1357",
         element %in% c("Agricultural Use", "Use per area of cropland")) %>%
  select(country_name, year, element, value) %>%
  pivot_wider(names_from = element, values_from = value) %>%
  transmute(country_name, year,
            cropland_ha = `Agricultural Use` * 1000 /
              `Use per area of cropland`)

insecticide <- use %>%
  filter(item_code == "1309", element == "Agricultural Use",
         between(year, first_year, last_year)) %>%
  left_join(cropland, by = join_by(country_name, year)) %>%
  mutate(kg_per_ha = value * 1000 / cropland_ha,
         official = flag == "A")

n_zero <- sum(insecticide$flag %in% c("A", "E", "X") &
                insecticide$kg_per_ha <= 0, na.rm = TRUE)

trend_data <- insecticide %>%
  filter(flag %in% c("A", "E", "X"),
         is.finite(kg_per_ha), kg_per_ha > 0) %>%
  mutate(y = log(kg_per_ha),
         t = (year - centre_year) / 10) %>%
  select(country_name, region, year, t, y, flag, official)

# models ----------------------------------------------------------------------

datasets <- list(reported = filter(trend_data, official),
                 all = trend_data)

fit_trend <- function(fixed, data) {
  data <- as.data.frame(data)
  weights <- if (n_distinct(data$official) > 1) {
    varIdent(form = ~ 1 | official)
  }
  lme(fixed, random = list(country_name = pdDiag(~ t)), data = data,
      correlation = corCAR1(form = ~ year | country_name),
      weights = weights, method = "ML",
      control = lmeControl(maxIter = 500, msMaxIter = 500, opt = "optim"))
}

fits <- lapply(datasets, function(data) {
  list(linear = fit_trend(y ~ t, data),
       quadratic = fit_trend(y ~ t + I(t ^ 2), data),
       regional_level = fit_trend(y ~ t + region, data),
       regional_slope = fit_trend(y ~ t * region, data))
})

# lag-1 autocorrelation of a fit's normalised residuals within countries, in
# consecutive years: near 0 if the AR(1) term captures the serial correlation
residual_lag1 <- function(fit, data) {
  tibble(country_name = data$country_name, year = data$year,
         r = residuals(fit, type = "normalized")) %>%
    arrange(country_name, year) %>%
    group_by(country_name) %>%
    mutate(r_lag = if_else(lag(year) == year - 1, lag(r), NA_real_)) %>%
    ungroup() %>%
    summarise(lag1 = cor(r, r_lag, use = "complete.obs")) %>%
    pull(lag1)
}

options(width = 200)
cat("zero values dropped:", n_zero, "\n")
for (name in names(fits)) {
  data <- datasets[[name]]
  fit <- fits[[name]]
  slope <- summary(fit$linear)$tTable["t", ]
  quadratic <- summary(fit$quadratic)$tTable["I(t^2)", ]
  phi <- coef(fit$linear$modelStruct$corStruct, unconstrained = FALSE)
  cat(sprintf(paste("\n%s: %i values (%i reported), %i countries;",
                    "AR(1) correlation between consecutive years %.2f;",
                    "lag-1 autocorrelation of the normalised residuals %.2f\n"),
              name, nrow(data), sum(data$official),
              n_distinct(data$country_name), phi,
              residual_lag1(fit$linear, data)))
  cat(sprintf(paste("linear trend %.3f (95%% CI %.3f to %.3f) per year on",
                    "the log scale: 2025 use %.2f times 1995's\n"),
              slope[["Value"]] / 10,
              (slope[["Value"]] - 1.96 * slope[["Std.Error"]]) / 10,
              (slope[["Value"]] + 1.96 * slope[["Std.Error"]]) / 10,
              exp(slope[["Value"]] * (2025 - 1995) / 10)))
  cat(sprintf("quadratic term %.3f (SE %.3f) per decade^2\n",
              quadratic[["Value"]], quadratic[["Std.Error"]]))
  cat("non-linearity: quadratic against linear\n")
  print(anova(fit$linear, fit$quadratic))
  cat("regional variation: regional slopes against regional intercepts\n")
  print(anova(fit$regional_level, fit$regional_slope))
  cat("countries per region:\n")
  print(count(distinct(data, region, country_name), region))
}

# figures ---------------------------------------------------------------------

# Caption: FAOSTAT insecticide use (item 1309, kg active ingredient) per ha of
# cropland, per country and year, by region; each country's values are joined
# by a line. Black: the common linear trend fitted to the values reported by
# the countries (the fixed effects, with their 95% interval), at the mean level
# of the region's reporting countries.
value_labels <- c(official = "Reported by the country (FAOSTAT flag A)",
                  estimated = "Estimated by FAO (flags E and X)")
value_type <- function(d) {
  mutate(d, value = if_else(official, value_labels[["official"]],
                            value_labels[["estimated"]]))
}

main <- fits$reported$linear
year_grid <- first_year:last_year
t_grid <- (year_grid - centre_year) / 10
design <- cbind(1, t_grid)
fixed <- design %*% fixef(main)
fixed_se <- sqrt(rowSums((design %*% vcov(main)) * design))
country_levels <- ranef(main) %>%
  as_tibble(rownames = "country_name") %>%
  rename(level = `(Intercept)`, slope = t) %>%
  left_join(distinct(datasets$reported, country_name, region),
            by = "country_name")
region_lines <- country_levels %>%
  group_by(region) %>%
  summarise(level = mean(level), slope = mean(slope), .groups = "drop") %>%
  cross_join(tibble(year = year_grid, t = t_grid, fixed = fixed[, 1],
                    se = fixed_se)) %>%
  mutate(median = fixed + level + slope * t,
         lower = median - 1.96 * se,
         upper = median + 1.96 * se)

trend_plot <- ggplot(value_type(trend_data), aes(year, exp(y))) +
  geom_line(aes(group = country_name), colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = value, shape = value), size = 1.6) +
  geom_ribbon(aes(year, ymin = exp(lower), ymax = exp(upper)),
              data = region_lines, inherit.aes = FALSE,
              fill = "black", alpha = 0.15) +
  geom_line(aes(year, exp(median)), data = region_lines,
            inherit.aes = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = setNames(c("#1f5fa8", "#9a9a9a"),
                                        value_labels),
                      name = NULL) +
  scale_shape_manual(values = setNames(c(16, 1), value_labels), name = NULL) +
  scale_y_log10() +
  facet_wrap(~ region, nrow = 1) +
  labs(x = NULL,
       y = "Insecticide use (kg active ingredient per ha of cropland)") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("figures/faostat_trend_regional_fits.png", trend_plot,
       width = 12, height = 4, dpi = 150, bg = "white")

# Caption: the trend per region (slope of log insecticide use per ha of
# cropland, per year) in the model with regional slopes, with 95% intervals,
# for the reported values and for all values; the dashed line is the common
# trend of the linear model.
region_slopes <- imap(fits, function(fit, name) {
  model <- fit$regional_slope
  contrast <- t(sapply(regions, function(r) {
    as.numeric(names(fixef(model)) %in%
                 c("t", paste0("t:region", r)))
  }))
  estimate <- contrast %*% fixef(model)
  se <- sqrt(rowSums((contrast %*% vcov(model)) * contrast))
  tibble(data = name, region = regions, slope = estimate[, 1] / 10,
         lower = (estimate[, 1] - 1.96 * se) / 10,
         upper = (estimate[, 1] + 1.96 * se) / 10,
         common = fixef(fit$linear)[["t"]] / 10)
}) %>%
  bind_rows() %>%
  mutate(data = factor(c(reported = "Reported values",
                         all = "All values")[data],
                       levels = c("Reported values", "All values")),
         region = factor(region, levels = rev(regions)))

slope_plot <- ggplot(region_slopes, aes(slope, region)) +
  geom_vline(aes(xintercept = common), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.3) +
  facet_wrap(~ data, scales = "free_x") +
  labs(x = "Slope of log insecticide use per ha of cropland, per year",
       y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/faostat_trend_regional_slopes.png", slope_plot,
       width = 7, height = 3.5, dpi = 150, bg = "white")
