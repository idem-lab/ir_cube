# Hierarchical linear trend in agricultural insecticide use intensity, from
# FAOSTAT, to set the relative regional scales of the time trend g_ag(t) on the
# crop-layer selection covariates.
#
# Response: log kg insecticide (FAOSTAT item 1309, a.i.) per ha cropland,
# 1995-2024, for the modelled countries. Cropland is implied by FAOSTAT's own
# pesticides total and pesticides total per ha of cropland. Official (A) and
# estimated (E, X) values are used; imputed (I) values, which are carry-forwards
# and straight-line interpolations, are dropped, as are zeros.
#
#   y_i = a_c + (b + beta_r + delta_c) * t_i + e_i
#   a_c     ~ N(alpha, sigma_a)        country intercepts
#   beta_r  ~ N(0, sigma_beta)         regional slope deviations
#   delta_c ~ N(0, sigma_delta)        country slope deviations
#   e_i     ~ Student-t(nu, 0, s_official or s_estimated)
#
# with t the year centred on 2009.5, in decades. The regional slope is
# b + beta_r.
#
# Seven countries (Benin, DR Congo, Eswatini, Gabon, Liberia, Nigeria, Sierra
# Leone) have estimated insecticide use that is their net-trade pesticides
# total times one common insecticide share, so their E values carry no
# country information on the insecticide share. They are kept in the main fit,
# and dropped in a sensitivity fit, along with a fit to official values only.
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md, section 1):
#   R_LIBS=~/R/greta06-lib \
#   RETICULATE_PYTHON=~/.local/share/r-miniconda/envs/greta06-env/bin/python \
#   Rscript R/faostat_trend_model.R
# Needs data/clean/faostat_pesticide_use.RDS from R/prep_faostat.R.

source("R/greta_setup.R")
start_greta(threads = 2)
source("R/packages.R")
source("R/functions.R")

set.seed(2025)

first_year <- 1995
last_year <- 2024
centre_year <- (first_year + last_year) / 2
n_samples <- 3000
warmup <- 3000
n_chains <- 4

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

# countries whose estimated insecticide share is the common one: their E share
# series matches Nigeria's (to 0.005) in at least 90% of years
shares <- use %>%
  filter(item_code %in% c("1309", "1357"), element == "Agricultural Use") %>%
  select(country_name, year, item_code, value, flag) %>%
  pivot_wider(names_from = item_code, values_from = c(value, flag)) %>%
  filter(flag_1309 == "E") %>%
  mutate(share = value_1309 / value_1357)
common_share_countries <- shares %>%
  left_join(shares %>%
              filter(country_name == "Nigeria") %>%
              select(year, share_common = share),
            by = join_by(year)) %>%
  group_by(country_name) %>%
  summarise(match = mean(abs(share - share_common) < 0.005, na.rm = TRUE)) %>%
  filter(match >= 0.9) %>%
  pull(country_name)

insecticide <- use %>%
  filter(item_code == "1309", element == "Agricultural Use",
         between(year, first_year, last_year)) %>%
  left_join(cropland, by = join_by(country_name, year)) %>%
  mutate(kg_per_ha = value * 1000 / cropland_ha,
         official = flag == "A",
         common_share = country_name %in% common_share_countries &
           flag == "E")

n_zero <- sum(insecticide$flag %in% c("A", "E", "X") &
                insecticide$kg_per_ha <= 0, na.rm = TRUE)

trend_data <- insecticide %>%
  filter(flag %in% c("A", "E", "X"),
         is.finite(kg_per_ha), kg_per_ha > 0) %>%
  mutate(y = log(kg_per_ha),
         t = (year - centre_year) / 10) %>%
  select(country_name, region, year, t, y, flag, official, common_share)

# model -----------------------------------------------------------------------

# build the model for one data set; returns the greta model and the arrays to
# summarise
build_trend_model <- function(data) {
  countries <- sort(unique(data$country_name))
  country_region <- data %>%
    distinct(country_name, region) %>%
    arrange(match(country_name, countries))
  stopifnot(identical(country_region$country_name, countries))
  n_country <- length(countries)
  n_region <- length(regions)
  country_id <- match(data$country_name, countries)
  region_of_country <- match(country_region$region, regions)

  alpha <- normal(-3, 2)
  b <- normal(0, 1)
  sigma_a <- normal(0, 2, truncation = c(0, Inf))
  sigma_beta <- normal(0, 0.5, truncation = c(0, Inf))
  sigma_delta <- normal(0, 0.5, truncation = c(0, Inf))
  # country effects centred (each country has many values, and the
  # non-centred form mixed poorly for alpha and b); regional deviations
  # non-centred (five regions)
  a <- normal(alpha, sigma_a, dim = n_country)
  delta <- normal(0, sigma_delta, dim = n_country)
  beta_raw <- normal(0, 1, dim = n_region)
  beta <- sigma_beta * beta_raw
  region_slope <- b + beta

  # one scale per value type, the estimated one only if there are E values
  s_official <- normal(0, 1, truncation = c(0, Inf))
  if (all(data$official)) {
    s_obs <- s_official[rep(1, nrow(data))]
    s_estimated <- s_official
  } else {
    s_estimated <- normal(0, 1, truncation = c(0, Inf))
    s <- c(s_official, s_estimated)
    s_obs <- s[ifelse(data$official, 1, 2)]
  }
  nu <- gamma(2, 0.1)

  country_slope <- region_slope[region_of_country] + delta
  mu <- a[country_id] + country_slope[country_id] * data$t
  distribution(data$y) <- student(nu, mu, s_obs)

  list(
    model = model(alpha, b, sigma_a, sigma_beta, sigma_delta,
                  s_official, s_estimated, nu),
    arrays = list(b = b, region_slope = region_slope, a = a,
                  country_slope = country_slope, mu = mu, s_obs = s_obs,
                  nu = nu, s_official = s_official, s_estimated = s_estimated,
                  sigma_beta = sigma_beta, sigma_delta = sigma_delta),
    countries = countries,
    region_of_country = region_of_country
  )
}

fit_trend_model <- function(data) {
  built <- build_trend_model(data)
  draws <- mcmc(built$model, n_samples = n_samples, warmup = warmup,
                chains = n_chains, verbose = FALSE)
  rhat <- coda::gelman.diag(draws, multivariate = FALSE)$psrf[, 1]
  ess <- coda::effectiveSize(draws)
  # posterior draws of each array, as a draws x elements matrix
  sims <- lapply(built$arrays, function(x) {
    as.matrix(calculate(x, values = draws))
  })
  list(data = data, built = built, draws = draws, sims = sims,
       max_rhat = max(rhat), min_ess = min(ess),
       worst = names(sort(rhat, decreasing = TRUE))[1:3])
}

fits <- list(
  main = trend_data,
  official_only = filter(trend_data, official),
  no_common_share = filter(trend_data, !common_share)
) %>%
  map(fit_trend_model)

cat("zero values dropped:", n_zero, "\n")
cat("common-share countries:", paste(common_share_countries, collapse = ", "),
    "\n")
for (name in names(fits)) {
  cat(sprintf("%s: %d obs (%d official), %d countries, max rhat %.3f, min ess %.0f (worst: %s)\n",
              name, nrow(fits[[name]]$data), sum(fits[[name]]$data$official),
              length(fits[[name]]$built$countries),
              fits[[name]]$max_rhat, fits[[name]]$min_ess,
              paste(fits[[name]]$worst, collapse = ", ")))
}

# regional summary ------------------------------------------------------------

# slopes on the log scale per year, from per decade
slope_draws <- function(fit) {
  s <- cbind(fit$sims$b, fit$sims$region_slope) / 10
  colnames(s) <- c("Africa", regions)
  s
}

span_years <- 2025 - 1995

summarise_slopes <- function(fit, name) {
  s <- slope_draws(fit)
  ratio <- exp(s * span_years)
  # g_ag scale relative to the continent, two ways: ratio of log-scale slopes,
  # and ratio of the proportional increase in intensity since 1995
  scale_log <- s / s[, "Africa"]
  scale_increase <- (ratio - 1) / (ratio[, "Africa"] - 1)
  # probability that the region's slope is above the continental one
  p_above <- colMeans(s > s[, "Africa"])
  q <- function(x, p) apply(x, 2, quantile, p)
  data <- fit$data
  tibble(
    fit = name,
    region = colnames(s),
    n_countries = c(n_distinct(data$country_name),
                    map_int(regions, ~ n_distinct(data$country_name[data$region == .x]))),
    n_obs = c(nrow(data), map_int(regions, ~ sum(data$region == .x))),
    n_official = c(sum(data$official),
                   map_int(regions, ~ sum(data$official[data$region == .x]))),
    slope_per_year = colMeans(s),
    slope_lower = q(s, 0.025),
    slope_upper = q(s, 0.975),
    ratio_2025_1995 = apply(ratio, 2, median),
    ratio_lower = q(ratio, 0.025),
    ratio_upper = q(ratio, 0.975),
    g_scale_log_slope = apply(scale_log, 2, median),
    g_scale_increase = apply(scale_increase, 2, median),
    g_scale_increase_lower = q(scale_increase, 0.025),
    g_scale_increase_upper = q(scale_increase, 0.975),
    p_slope_above_africa = if_else(region == "Africa", NA_real_, p_above)
  )
}

regional_trend <- imap(fits, summarise_slopes) %>% bind_rows()

write_csv(filter(regional_trend, fit == "main") %>% select(-fit),
          "data/clean/faostat_regional_trend.csv")

options(width = 200)
print(regional_trend %>%
        select(fit, region, n_obs, n_official, slope_per_year, slope_lower,
               slope_upper, ratio_2025_1995, ratio_lower, ratio_upper,
               g_scale_increase, p_slope_above_africa) %>%
        mutate(across(where(is.double), ~ signif(.x, 3))),
      n = Inf)

main <- fits$main
cat(sprintf("main: sigma_beta %.3f [%.3f, %.3f]; sigma_delta %.3f; s_official %.2f; s_estimated %.2f; nu %.1f\n",
            mean(main$sims$sigma_beta),
            quantile(main$sims$sigma_beta, 0.025),
            quantile(main$sims$sigma_beta, 0.975),
            mean(main$sims$sigma_delta), mean(main$sims$s_official),
            mean(main$sims$s_estimated), mean(main$sims$nu)))

# posterior predictive checks --------------------------------------------------

ppc <- function(fit) {
  data <- fit$data
  mu <- fit$sims$mu
  s <- fit$sims$s_obs
  nu <- as.vector(fit$sims$nu)
  n_draws <- nrow(mu)
  y_rep <- mu + s * matrix(rt(length(mu), df = rep(nu, ncol(mu))),
                           nrow = n_draws)
  lower <- apply(y_rep, 2, quantile, 0.025)
  upper <- apply(y_rep, 2, quantile, 0.975)
  inside <- data$y >= lower & data$y <= upper

  # residual lag-1 autocorrelation within country series, observed vs
  # replicated, and the sd of residuals, by value type
  lag1 <- function(resid) {
    d <- tibble(country_name = data$country_name, year = data$year,
                official = data$official, r = resid) %>%
      arrange(country_name, year) %>%
      group_by(country_name) %>%
      mutate(r_lag = if_else(lag(year) == year - 1, lag(r), NA_real_)) %>%
      ungroup()
    d %>%
      group_by(official) %>%
      summarise(acf = cor(r, r_lag, use = "complete.obs"), .groups = "drop")
  }
  mu_mean <- colMeans(mu)
  observed <- lag1(data$y - mu_mean)
  replicated <- map(sample(n_draws, 200),
                    ~ lag1(y_rep[.x, ] - mu[.x, ])) %>%
    bind_rows() %>%
    group_by(official) %>%
    summarise(acf_rep_lower = quantile(acf, 0.025),
              acf_rep_upper = quantile(acf, 0.975), .groups = "drop")

  tibble(official = data$official, inside = inside,
         resid = data$y - mu_mean,
         rep_sd = apply(y_rep - mu, 2, sd)) %>%
    group_by(official) %>%
    summarise(n = n(), coverage_95 = mean(inside),
              resid_sd = sd(resid), resid_mad = mad(resid),
              .groups = "drop") %>%
    left_join(observed, by = join_by(official)) %>%
    left_join(replicated, by = join_by(official))
}

cat("posterior predictive checks (by value type: official TRUE/FALSE)\n")
imap(fits, ~ mutate(ppc(.x), fit = .y)) %>%
  bind_rows() %>%
  relocate(fit) %>%
  mutate(across(where(is.double), ~ signif(.x, 3))) %>%
  print(width = 200)

# figures ---------------------------------------------------------------------

# fitted regional trend: the region's mean country intercept plus the regional
# slope, over the country series
year_grid <- first_year:last_year
region_lines <- map(regions, function(r) {
  in_region <- main$built$region_of_country == match(r, regions)
  intercept <- rowMeans(main$sims$a[, in_region, drop = FALSE])
  slope <- main$sims$region_slope[, match(r, regions)]
  pred <- outer(intercept, rep(1, length(year_grid))) +
    outer(slope, (year_grid - centre_year) / 10)
  tibble(region = r, year = year_grid,
         median = apply(pred, 2, median),
         lower = apply(pred, 2, quantile, 0.025),
         upper = apply(pred, 2, quantile, 0.975))
}) %>%
  bind_rows()

country_lines <- map(seq_along(main$built$countries), function(i) {
  pred <- outer(main$sims$a[, i], rep(1, length(year_grid))) +
    outer(main$sims$country_slope[, i], (year_grid - centre_year) / 10)
  tibble(country_name = main$built$countries[i],
         region = regions[main$built$region_of_country[i]],
         year = year_grid, median = apply(pred, 2, median))
}) %>%
  bind_rows()

value_type <- function(d) {
  mutate(d, value = if_else(official, "official (A)", "estimated (E, X)"))
}

trend_plot <- ggplot(value_type(trend_data), aes(year, exp(y))) +
  geom_line(aes(year, exp(median), group = country_name),
            data = country_lines,
            colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = value, shape = value), size = 0.9) +
  geom_ribbon(aes(year, ymin = exp(lower), ymax = exp(upper)),
              data = region_lines, inherit.aes = FALSE,
              fill = "black", alpha = 0.15) +
  geom_line(aes(year, exp(median)), data = region_lines,
            inherit.aes = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = c("official (A)" = "#1f5fa8",
                                 "estimated (E, X)" = "#9a9a9a"),
                      name = NULL) +
  scale_shape_manual(values = c("official (A)" = 16,
                                "estimated (E, X)" = 1), name = NULL) +
  scale_y_log10() +
  facet_wrap(~ region, nrow = 1) +
  labs(x = NULL, y = "insecticide, kg a.i. per ha cropland",
       caption = paste("Black: regional trend (mean country intercept, regional slope), 95% interval.",
                       "Grey lines: fitted country trends. FAOSTAT RP, item 1309.")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/faostat_trend_regional_fits.png", trend_plot,
       width = 12, height = 4, dpi = 150, bg = "white")

slope_long <- imap(fits, function(fit, name) {
  s <- slope_draws(fit)
  as_tibble(s) %>%
    pivot_longer(everything(), names_to = "region",
                 values_to = "slope") %>%
    mutate(fit = name)
}) %>%
  bind_rows() %>%
  mutate(region = factor(region, levels = rev(c("Africa", regions))),
         fit = factor(fit, levels = names(fits)))

slope_plot <- ggplot(slope_long, aes(slope, region, fill = fit)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.7,
              position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  scale_fill_manual(values = c(main = "#1f5fa8",
                               official_only = "#e08a1e",
                               no_common_share = "#9a9a9a"),
                    name = NULL) +
  labs(x = "slope of log insecticide per ha cropland, per year",
       y = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/faostat_trend_regional_slopes.png", slope_plot,
       width = 7, height = 5, dpi = 150, bg = "white")
