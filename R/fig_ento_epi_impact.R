# Plot the decline (2000-2025 and projected 2025-2030) in the effectiveness of
# pyrethroid-only ITNs across Africa due to IR, its cost to the effective
# coverage of the nets in use, and where pyrethroid-only or next-generation
# nets are the better buy (#75). Effectiveness is the relative epidemiological
# effectiveness g*(k) of Symons et al. (below).

# four panel plot:

# A) Declining effectiveness of nets: timeseries with ribbon of uncertainty
# (from Tas' estimates) giving the effectiveness, relative to a fully
# susceptible population, that a pyrethroid-only LLIN would have for each
# person at risk, whether or not nets are in use there: averaged across Africa
# over the population at risk (population times PfPR, both in 2000, fixed over
# time)

# B) Effective coverage: net use times effectiveness, summed over three groups
# of net types (R/prep_net_use_pyrethroid.R) and averaged across Africa over
# the population at risk of each year (the year's population times PfPR in
# 2000), as stacked areas: dual-AI nets, PBO nets, then pyrethroid-only nets
# (LLINs, and conventional ITNs counting as 0.25 of an LLIN). PBO and dual-AI
# nets are assumed to be fully effective (effectiveness 1), losing none of it
# to resistance. The line above the stack is the counterfactual without
# resistance, all nets in use fully effective; the gap between them is the
# coverage lost to resistance

# C) Spatial variation in effectiveness: classes of 2025 susceptibility on a
# map, by whether pyrethroid-only or dual-AI nets give more effective coverage
# for the same budget (break_even_effectiveness()), with the delivery cost per
# net from two mass campaigns in Mali (Scates et al. 2020), under two models of
# effectiveness (Symons et al. and Churcher et al. 2024); uncertain where the
# models and delivery costs disagree

# D) Effectiveness over time within each class, as A

# After the last year of net use data (2025), use and the mix of net types are
# held at that year's values. The ribbons are 50% intervals of Symons et al.'s
# parameter uncertainty, conditional on the posterior mean susceptibility
# maps. The numbers for the text are printed after the panels' data.



# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/model_covariates.R")
source("R/figure_cache.R")

# The panels' data are cached in temporary/fig_ento_epi_impact.rds
# (R/figure_cache.R), so changes to the plots alone rerun in seconds

# data to load:

# spatial:
# - Country borders for plotting
# - pf transmission limits and water mask for plotting
# - IR layers
# - MAP net use by type (to manually project to 2030)
# - MAP PR 2000 baseline
# - MAP populations by year (projected to 2030)
# other:
# - Symons epi impact estimate
# - net prices and delivery costs (data/raw/costs/net_costs.csv)

# years to plot for
years_plot <- 2000:2030

# load admin borders for plotting
borders <- readRDS("data/clean/country_borders.RDS")


# The panels' data: effectiveness and effective coverage over time,
# Africa-wide and in the bands of 2025 effectiveness, the bands themselves
# (wrapped, to be saved), and the numbers for the text, from the two-stage IR
# maps, net use by type (net_type_files), PfPR 2000 and population. n_sims
# draws of the Symons et al. coefficients, from `seed`; break_even, the lowest
# and highest break-even susceptibility k* (below)
compute_epi_impact <- function(years_plot, n_sims, seed, net_type_files,
                               break_even) {
  # write the cubes terra computes to temporary files: held in memory, they
  # peaked at 26 GB; on disk, at 8 GB, for a minute and a half longer
  terraOptions(todisk = TRUE)

  # load mask with limits of transmission and water bodies for plotting
  pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

  # - IR layers (susceptibility to the pyrethroids used in LLINs), of the
  # two-stage model (R/two_stage_maps.R)
  ir <- rast(ir_map_files("llin_effective", years_plot))
  names(ir) <- years_plot

  # - net use by type (net_type_files), carried forward from its last year to
  # 2030, which holds both use and the mix of types
  nets_by_type <- lapply(net_type_files, function(file) {
    nets <- rast(file)
    names(nets) <- str_remove(names(nets), "^nets_")
    nets
  })
  latest_nets <- max(as.numeric(names(nets_by_type$pyrethroid)))
  end_year <- max(years_plot)
  nets_by_type <- lapply(nets_by_type, post_pad_cube, end_year)
  stopifnot(all(vapply(nets_by_type, function(nets) {
    identical(names(nets), names(ir))
  }, logical(1))))

  # - MAP PR 2000 baseline
  pr2000 <- rast("data/clean/pfpr_2000.tif")

  # - MAP populations by year (projected to 2030)
  pop <- rast("data/clean/pop_cube.tif")
  names(pop) <- str_remove(names(pop), "^pop_")
  pop_future <- rast("data/clean/pop_cube_future.tif")
  names(pop_future) <- str_remove(names(pop_future), "^pop_")
  pop_all <- c(pop, pop_future)

  # other:

  # Encode Symons et al. PMI function for relative effectiveness of ITNs at
  # reducing logit-PfPR (reduction in on the log-odds ratio), compared to in the absence
  # of resistance.

  # when susceptibility is 100%, the ITN effectiveness parameter has full effect
  # (multiplier 1), when susceptibility is 0%, the ITN effectiveness parameter has
  # 54% (multiplier 0.54).

  # compute the relative epi impact from the posterior mean of the parameters
  epi_function_mean <- function(susceptibility) {
    1 - 0.46 * (1 - susceptibility)
  }

  # The log-odds ratio for malaria infection under full ITN coverage (relative to
  # no ITN coverage) is -1.68, which corresponds to an odds ratio of 0.186374.
  # With a fully-resistant vector population, this is would reduce to 54% of its
  # previous value, -0.9072, giving an odds ratio of 0.4036529.

  # Alternative interpretation: the odds ratio for a population with 100% coverage
  # and 100% resistance is equivalent to that of a population with 54% coverage
  # and 0% resistance.

  # get the estimate and SD. for the parameters behind this model, to be able to
  # propagate variance through the subsequent calculations

  # beta_0 = beta_{physical} = 'ITN_physical'
  #   beta_0 ~ N(-1.68, 0.11^2)
  # beta_1 = beta_{penalty} = 'ITN_insecticidal_penalty'
  #   beta_1 ~ N(0.78, 0.12^2)

  # the log odds (contribution to the linear predictor) is computed from these
  # parameters, with a non-linear function of ITN coverage f(ITN) and
  # susceptibility k:

  #   log_odds = beta_0 * f(ITN) + beta_1 * (1 - k) * f(ITN)
  #            = f(ITN) * (beta_0 + beta_1 * (1 - k))

  # we can therefore factorise to get a net effectiveness measure g(k):
  #   log_odds = f(ITN) * g(k)
  #   g(k) = beta_0 + beta_1 * (1 - k)

  # from which we can obtain a *relative* net effectiveness, giving the
  # effectiveness relative to a fully susceptible population (k = 1):
  #   g*(k) = g(k) / g(1)

  # since:
  #   g(1) = beta_0 + beta_1 * (1 - k)
  #        = beta_0 + beta_1 * 0
  #        = beta_0
  # we have:
  #   g*(k) = g(k) / beta_0
  #         = (beta_0 + beta_1 * (1 - k)) / beta_0
  #         = 1 + (beta_1 * (1 - k)) / beta_0

  # this is equivalent to a ratio of normal distributions, but with non-zero means
  # it doesn't have a closed form (that I know of), so I will compute credible
  # intervals of effectiveness by Monte Carlo simulation.



  # draws of the coefficients from the posteriors given in Symons et al., taken
  # as independent, since the posterior correlation is not given. The relative
  # effectiveness is linear in susceptibility k, g*(k) = 1 + r (1 - k), with
  # r = beta_1 / beta_0, so each draw needs only its r

  set.seed(seed)
  beta_0_sims <- rnorm(n_sims, -1.68, 0.11)
  beta_1_sims <- rnorm(n_sims, 0.78, 0.12)
  ratio_sims <- beta_1_sims / beta_0_sims

  # perform calculations:

  # define the weights for net effectiveness and effective coverage, from the
  # population at risk (potential infected population, in absence of
  # interventions): population times PfPR in 2000, at the cells with a
  # susceptibility map (all but 1e-5 of the population at risk)

  # convert a raster layer into spatial weights (pixels sum to 1) to get averages
  # over space
  normalise <- function(x) {
    total <- global(x, "sum", na.rm = TRUE)[1, 1]
    x / total
  }

  population_at_risk <- terra::mask(pop_all * pr2000, ir)

  # cube of population at risk weights, by year, for effective coverage
  infection_weights <- sapp(population_at_risk, normalise)

  # cube of population at risk weights fixed at 2000, for net effectiveness:
  # the effectiveness a pyrethroid-only LLIN would have for each person at
  # risk, whether or not nets are in use there, so its trend is that of
  # resistance alone
  risk_weights <- rep(normalise(population_at_risk[["2000"]]), nlyr(ir))
  names(risk_weights) <- names(ir)


  # Effective coverage: net use times effectiveness, averaged across Africa
  # over the population at risk, for each group of net types. PBO and dual-AI
  # nets are assumed to be fully effective (effectiveness 1), losing none of it
  # to resistance; pyrethroid-only nets have g*(k). Without resistance, all
  # nets would be fully effective, and effective coverage would be the average
  # net use. If net use was 100% for all populations at risk, this would equal
  # 1. It is lower due to imperfect net coverage (but not due to IR).
  coverage <- function(x) {
    global(x * infection_weights, fun = "sum", na.rm = TRUE)$sum
  }
  coverage_pbo <- coverage(nets_by_type$pbo)
  coverage_dual_ai <- coverage(nets_by_type$dual_ai)
  coverage_no_ir <- coverage(nets_by_type$pyrethroid) + coverage_pbo +
    coverage_dual_ai


  # Now average epi impacts, accounting for IR (both per net and at population
  # level)

  # cube of net relative epi impact of a net in each pixel (using the posterior
  # mean function)
  relative_epi_impact_mean <- terra::app(ir, epi_function_mean)
  names(relative_epi_impact_mean) <- names(ir)

  # average the epi impact of a net across Africa, over the population at
  # risk, to get the temporal trends in impact of each net (using the posterior
  # mean function)
  average_epi_impact_mean <- global(
    x = relative_epi_impact_mean * risk_weights,
    fun = "sum",
    na.rm = TRUE
  )

  # the effective coverage of pyrethroid-only nets across Africa (using the
  # posterior mean function)
  coverage_pyrethroid_mean <- coverage(
    nets_by_type$pyrethroid * relative_epi_impact_mean
  )


  # classes of 2025 susceptibility, to divide up the map into regions by
  # which nets give more effective coverage for a fixed budget: pyrethroid-only
  # nets above the highest break-even susceptibility k*, next-generation nets
  # below the lowest, and uncertain between (break_even)
  susceptibility_2025 <- ir[["2025"]]
  bands <- 1 + (susceptibility_2025 >= break_even[["low"]]) +
    (susceptibility_2025 > break_even[["high"]])
  levels(bands) <- data.frame(
    id = 1:3,
    band = c("next_generation", "uncertain", "pyrethroid")
  )
  band_names <- c("pyrethroid", "uncertain", "next_generation")

  # create population at risk weights subsetted to each of these bands, to
  # compute effectiveness over time in each
  risk_weights_bands <- lapply(band_names, function(band) {
    sapp(risk_weights * (bands == band), normalise)
  })
  names(risk_weights_bands) <- band_names

  # compute posterior means for these
  average_epi_impact_mean_bands <- lapply(risk_weights_bands, function(weights) {
    global(
      x = relative_epi_impact_mean * weights,
      fun = "sum",
      na.rm = TRUE
    )$sum
  })

  # the population at risk in 2025 by susceptibility, in bins of 1e-4, for the
  # classes' shares of it and their mean effectiveness under each model
  people_by_susceptibility <- terra::zonal(
    population_at_risk[["2025"]],
    round(susceptibility_2025 * 1e4),
    fun = "sum",
    na.rm = TRUE
  ) %>%
    as_tibble(.name_repair = ~ c("susceptibility", "people")) %>%
    mutate(susceptibility = susceptibility / 1e4)
  weighted_total <- function(x) {
    global(x, "sum", na.rm = TRUE)$sum
  }

  # each type's share of net use in 2024, Africa-wide: the people using a net
  # of the type over the people using any
  use_2024 <- vapply(nets_by_type, function(nets) {
    weighted_total(nets[["2024"]] * pop_all[["2024"]])
  }, numeric(1))
  type_shares_2024 <- use_2024 / sum(use_2024)

  # Now do the same under a Monte Carlo simulation, to get credible intervals on
  # these. As g*(k) is linear in k, a weighted sum of it over pixels is
  # S_w + r (S_w - S_wk), with S_w the sum of the weights and S_wk the sum of
  # the weights times k, over the pixels with a susceptibility. So the sums are
  # computed once per year, and each draw is arithmetic on them: years x draws
  weighted_sums <- function(weights) {
    weights <- terra::mask(weights, ir)
    list(w = global(weights, "sum", na.rm = TRUE)$sum,
         wk = global(weights * ir, "sum", na.rm = TRUE)$sum)
  }
  simulate_impact <- function(weights) {
    sums <- weighted_sums(weights)
    impact <- sums$w + outer(sums$w - sums$wk, ratio_sims)
    rownames(impact) <- years_plot
    impact
  }
  individual_net_impact_sims <- simulate_impact(risk_weights)
  individual_net_impact_sims_bands <- lapply(risk_weights_bands,
                                             simulate_impact)
  # all nets' effective coverage, the PBO and dual-AI nets' being certain
  coverage_sims <- simulate_impact(nets_by_type$pyrethroid * infection_weights) +
    coverage_pbo + coverage_dual_ai

  # check: at the posterior mean coefficients (r = -0.46), the sums give what
  # the map does
  mean_check <- weighted_sums(risk_weights)
  stopifnot(isTRUE(all.equal(
    mean_check$w - 0.46 * (mean_check$w - mean_check$wk),
    average_epi_impact_mean$sum, tolerance = 1e-6)))


  # compute credible intervals of each of these over time: the 50% intervals
  # the panels show, and 95% intervals for the text
  intervals <- function(sims) {
    cis <- apply(sims, 1, quantile, c(0.25, 0.75, 0.025, 0.975))
    tibble(
      lower = cis[1, ],
      upper = cis[2, ],
      lower_95 = cis[3, ],
      upper_95 = cis[4, ]
    )
  }

  # panel A: the effectiveness of a net over time
  panel_a_data <- tibble(
    year = years_plot,
    mean = average_epi_impact_mean$sum
  ) %>%
    bind_cols(intervals(individual_net_impact_sims))

  # panel B: effective coverage by net type, its total and interval, and the
  # counterfactual without IR
  panel_b_data <- tibble(
    year = years_plot,
    pyrethroid = coverage_pyrethroid_mean,
    pbo = coverage_pbo,
    dual_ai = coverage_dual_ai,
    mean = pyrethroid + pbo + dual_ai,
    mean_no_ir = coverage_no_ir
  ) %>%
    bind_cols(intervals(coverage_sims))

  # panel C: the bands, within the limits of transmission, in memory to be
  # saved (not in a temporary file)
  bands_mask <- terra::toMemory(terra::mask(bands, pf_water_mask))

  # panel D: the effectiveness of a net in each band
  panel_d_data <- bind_rows(lapply(band_names, function(band_name) {
    tibble(
      year = years_plot,
      band = band_name,
      mean = average_epi_impact_mean_bands[[band_name]]
    ) %>%
      bind_cols(intervals(individual_net_impact_sims_bands[[band_name]]))
  })) %>%
    mutate(
      band = factor(band,
                    levels = band_names)
    )

  list(panel_a_data = panel_a_data,
       panel_b_data = panel_b_data,
       panel_d_data = panel_d_data,
       bands = terra::wrap(bands_mask),
       people_by_susceptibility = people_by_susceptibility,
       type_shares_2024 = type_shares_2024,
       latest_nets = latest_nets)
}

# The break-even effectiveness of a pyrethroid-only net (#75). For a fixed
# budget, a net of price P with delivery cost D (the cost of delivering a net,
# without the net itself) buys budget / (P + D) nets. A pyrethroid-only net of
# relative effectiveness e then gives more effective coverage than a dual-AI
# net where e / (P_pyr + D) > 1 / (P_dual + D), that is where
#   e > e*(D) = (P_pyr + D) / (P_dual + D).
# This assumes a fixed budget, the same delivery cost for every net type,
# dual-AI nets fully effective, and effective coverage as the outcome
break_even_effectiveness <- function(delivery, price_pyrethroid, price_dual_ai) {
  (price_pyrethroid + delivery) / (price_dual_ai + delivery)
}

# net prices, campaign costs and CPI, with their sources
net_costs <- read_csv("data/raw/costs/net_costs.csv", show_col_types = FALSE)
cost_value <- function(quantity, item) {
  value <- net_costs$value[net_costs$quantity == quantity &
                             net_costs$item == item]
  stopifnot(length(value) == 1)
  value
}

# The delivery cost per net in the mass campaigns of Scates et al. 2020, in
# Segou and Bamako, Mali. They report the total economic cost in USD and in
# international dollars (Int$), not the nets' part of it. Nets are tradable
# goods, so cost the same in both; all other costs are converted at purchasing
# power parity, k times higher in Int$. With the same net cost N per net at
# both sites, a site's cost per net is N + D in USD and N + k D in Int$, so k
# is the ratio of the sites' differences in Int$ and in USD, and
# D = (Int$ - USD) / (k - 1). Inflated from 2016 to 2026 USD with US CPI-U
campaigns <- tibble(site = c("segou", "bamako")) %>%
  mutate(
    usd = map_dbl(site, ~ cost_value("campaign_cost_usd", .x)),
    int = map_dbl(site, ~ cost_value("campaign_cost_int", .x)),
    nets = map_dbl(site, ~ cost_value("nets_distributed", .x)),
    usd = usd / nets,
    int = int / nets
  )
ppp_ratio <- diff(campaigns$int) / diff(campaigns$usd)
inflation <- cost_value("cpi", "2026") / cost_value("cpi", "2016")
campaigns <- campaigns %>%
  mutate(
    delivery_2016 = (int - usd) / (ppp_ratio - 1),
    net_2016 = usd - delivery_2016,
    delivery = delivery_2016 * inflation
  )
delivery_range <- c(low = min(campaigns$delivery),
                    high = max(campaigns$delivery))
break_even <- break_even_effectiveness(
  delivery = delivery_range,
  price_pyrethroid = cost_value("net_price", "pyrethroid_only"),
  price_dual_ai = cost_value("net_price", "dual_ai")
)

# The relative effectiveness of pyrethroid-only nets in the simulations of
# the MINT transmission model behind Churcher et al. 2024 (Lancet Global
# Health 12:e1973-e1983; MINT-data-v20230208,
# https://doi.org/10.5281/zenodo.8344312): a mass distribution of
# pyrethroid-only nets ("llin") at 60% use against none, with no nets or IRS
# before, in each of 56 settings (prevalence in children under 5 of 5-60%,
# human and indoor biting high or low, seasonal or perennial), at bioassay
# survival 0, 20, ..., 100%. The nets' effect is the difference in the
# log-odds of prevalence in children under 5 between none and nets, averaged
# over the three years after distribution, as Symons et al. measure it,
# relative to that in the same setting at k = 1; the median over the settings
mint_dir <- "data/raw/mint_zenodo_8344312"
mint_index <- readRDS(file.path(mint_dir, "index.rds")) %>%
  filter(itnUsage == "0%", sprayInput == "0%")
stopifnot(nrow(mint_index) == 56 * 6)
mint_files <- paste0("table/", mint_index$index, ".rds")
unzip(file.path(mint_dir, "table.zip"), files = mint_files,
      exdir = tempdir())
mint_prevalence <- function(table) {
  c(table$prevYear1, table$prevYear2, table$prevYear3)
}
mint_effectiveness <- mint_index %>%
  mutate(
    susceptibility = 1 - parse_number(levelOfResistance) / 100,
    effect = map_dbl(file.path(tempdir(), mint_files), function(file) {
      table <- readRDS(file)
      none <- filter(table, intervention == "none")
      nets <- filter(table, intervention == "llin", near(netUse, 0.6))
      stopifnot(nrow(none) == 1, nrow(nets) == 1)
      mean(qlogis(mint_prevalence(none)) - qlogis(mint_prevalence(nets)))
    })
  ) %>%
  group_by(currentPrevalence, bitingPeople, bitingIndoors,
           seasonalityOfTransmission) %>%
  mutate(effectiveness = effect / effect[susceptibility == 1]) %>%
  group_by(susceptibility) %>%
  summarise(effectiveness = median(effectiveness))

# Two models of the relative effectiveness e(k) of a pyrethroid-only net, as a
# function of susceptibility k (bioassay mortality), 1 at k = 1, both the
# nets' effect on the log-odds of prevalence:
#   symons    Symons et al., at the posterior mean of their coefficients
#             (epi_function_mean() in compute_epi_impact())
#   churcher  Churcher et al. 2024's MINT simulations, above, interpolated
#             linearly between the simulated k
effectiveness_models <- list(
  symons = function(k) {
    1 - 0.46 * (1 - k)
  },
  churcher = approxfun(
    x = mint_effectiveness$susceptibility,
    y = mint_effectiveness$effectiveness
  )
)

# the susceptibility at which a model's effectiveness is e: 0 where it is
# above e at any susceptibility
susceptibility_at <- function(model, e) {
  if (model(0) >= e) {
    return(0)
  }
  uniroot(function(k) model(k) - e, c(0, 1), tol = 1e-8)$root
}

# The break-even susceptibility k*, where e(k*) = e*(D), under each model at
# the lowest and highest delivery cost D. e*(D) rises with D, so the lowest
# k* is that of the more optimistic model at the lowest D, and the highest
# that of the other model at the highest D. Pyrethroid-only nets are the
# better buy above the highest k*, under either model at any D in the range,
# and next-generation nets below the lowest; between, it is uncertain
break_even_thresholds <- expand_grid(
  model = names(effectiveness_models),
  delivery = c("low", "high")
) %>%
  mutate(
    effectiveness = break_even[delivery],
    susceptibility = map2_dbl(model, effectiveness, function(model, e) {
      susceptibility_at(effectiveness_models[[model]], e)
    })
  )
break_even_susceptibility <- c(
  low = min(break_even_thresholds$susceptibility),
  high = max(break_even_thresholds$susceptibility)
)

# net use by type (R/prep_net_use_pyrethroid.R): MITN run 06 use times each
# group's share of the net crop in the admin 1 unit. The pyrethroid-only nets
# are the model's (net_use_file(), R/model_covariates.R): LLINs, and
# conventional ITNs counting as 0.25 of an LLIN
net_type_files <- c(pyrethroid = net_use_file(selection_design()),
                    pbo = "data/clean/net_use_pbo_cube.tif",
                    dual_ai = "data/clean/net_use_dual_ai_cube.tif")

epi <- cached_figure_data(
  "temporary/fig_ento_epi_impact.rds", compute_epi_impact,
  args = list(years_plot = years_plot, n_sims = 1000, seed = 2026,
              net_type_files = net_type_files,
              break_even = break_even_susceptibility),
  sources = c("R/packages.R", "R/functions.R", "R/model_covariates.R"),
  inputs = c(ir_map_files("llin_effective", years_plot), net_type_files,
             "data/clean/pfpr_water_mask.tif", "data/clean/pfpr_2000.tif",
             "data/clean/pop_cube.tif", "data/clean/pop_cube_future.tif")
)
panel_a_data <- epi$panel_a_data
panel_b_data <- epi$panel_b_data
panel_d_data <- epi$panel_d_data
bands_mask <- terra::unwrap(epi$bands)
latest_nets <- epi$latest_nets


# numbers for the text --------------------------------------------------------

# percentages: a mean with its 50% and 95% intervals
with_intervals <- function(mean, lower, upper, lower_95, upper_95) {
  sprintf("%.1f [50%%: %.1f, %.1f; 95%%: %.1f, %.1f]", 100 * mean,
          100 * lower, 100 * upper, 100 * lower_95, 100 * upper_95)
}

cat("\nAfrica-wide effective coverage (% of the population at risk, as panel B;",
    "pyrethroid-only nets at the posterior mean of Symons et al.):\n")
panel_b_data %>%
  filter(year %in% c(2017, 2024)) %>%
  mutate(lost = mean_no_ir - mean) %>%
  rowwise() %>%
  mutate(
    text = paste0(
      year, ": pyrethroid-only ", sprintf("%.1f", 100 * pyrethroid),
      ", PBO ", sprintf("%.1f", 100 * pbo),
      ", dual-AI ", sprintf("%.1f", 100 * dual_ai),
      "\n  total ", with_intervals(mean, lower, upper, lower_95, upper_95),
      "\n  no resistance ", sprintf("%.1f", 100 * mean_no_ir),
      "\n  lost to resistance (percentage points) ",
      with_intervals(lost, mean_no_ir - upper, mean_no_ir - lower,
                     mean_no_ir - upper_95, mean_no_ir - lower_95),
      "\n  lost to resistance (% of no resistance) ",
      with_intervals(lost / mean_no_ir, 1 - upper / mean_no_ir,
                     1 - lower / mean_no_ir, 1 - upper_95 / mean_no_ir,
                     1 - lower_95 / mean_no_ir)
    )
  ) %>%
  pull(text) %>%
  cat(sep = "\n")

cat("\nEffectiveness of a pyrethroid-only LLIN, weighted by the population at",
    "risk in 2000 (%, panel A):\n")
panel_a_data %>%
  filter(year %in% c(2017, 2024, 2030)) %>%
  mutate(text = paste0(year, ": ", with_intervals(mean, lower, upper,
                                                  lower_95, upper_95))) %>%
  pull(text) %>%
  cat(sep = "\n")

cat("\nShares of net use in 2024, Africa-wide (% of the people using a net,",
    "conventional ITNs counting 0.25):\n")
cat(sprintf("%s %.1f", names(epi$type_shares_2024),
            100 * epi$type_shares_2024), sep = "\n")

cat("\nDelivery cost per net (USD per net distributed; Scates et al. 2020):",
    sprintf("\n  PPP ratio k %.3f, net cost %s (2016 USD)", ppp_ratio,
            paste(sprintf("%.3f", campaigns$net_2016), collapse = ", ")),
    sprintf("\n  %s: %.2f (2016 USD), %.2f (2026 USD)", campaigns$site,
            campaigns$delivery_2016, campaigns$delivery),
    sprintf("\n  CPI 2026 / 2016: %.4f", inflation), "\n")
cat(sprintf("D_low %.2f, D_high %.2f (2026 USD); e*(D_low) %.4f, e*(D_high) %.4f\n",
            delivery_range[["low"]], delivery_range[["high"]],
            break_even[["low"]], break_even[["high"]]))

cat("\nBreak-even susceptibility k* (bioassay mortality, %), where e(k*) =",
    "e*(D), by effectiveness model and delivery cost:\n")
break_even_thresholds %>%
  mutate(text = sprintf("%s, D_%s: e* %.1f, k* %.1f", model, delivery,
                        100 * effectiveness, 100 * susceptibility)) %>%
  pull(text) %>%
  cat(sep = "\n")

# the classes' shares of the population at risk in 2025, and their mean 2025
# effectiveness over it under each model
band_summary <- epi$people_by_susceptibility %>%
  mutate(
    band = case_when(
      susceptibility > break_even_susceptibility[["high"]] ~ "pyrethroid",
      susceptibility >= break_even_susceptibility[["low"]] ~ "uncertain",
      .default = "next_generation"
    )
  ) %>%
  group_by(band) %>%
  summarise(
    share = sum(people),
    symons = weighted.mean(effectiveness_models$symons(susceptibility),
                           people),
    churcher = weighted.mean(effectiveness_models$churcher(susceptibility),
                             people)
  ) %>%
  mutate(share = share / sum(share))

cat("\nClasses of 2025 susceptibility (share: % of the population at risk in",
    "2025; effectiveness: its mean over them under each model, %):\n")
band_summary %>%
  mutate(text = sprintf(paste("%s: share %.1f, effectiveness %.1f (Symons),",
                              "%.1f (Churcher)"),
                        band, 100 * share, 100 * symons,
                        100 * churcher)) %>%
  pull(text) %>%
  cat(sep = "\n")

# the gain in effective coverage from buying pyrethroid-only instead of dual-AI
# nets with the same budget, where pyrethroid-only nets are the better buy:
# (P_dual + D) / (P_pyr + D) e - 1, linear in e, so its mean over the
# population at risk is that at the class's mean effectiveness
pyrethroid_better <- band_summary %>%
  filter(band == "pyrethroid")
cat("\nGain in effective coverage from pyrethroid-only nets, where they are",
    "the better buy (% of that from dual-AI nets):\n")
expand_grid(model = names(effectiveness_models),
            delivery = c("low", "high")) %>%
  mutate(
    gain = (cost_value("net_price", "dual_ai") + delivery_range[delivery]) /
      (cost_value("net_price", "pyrethroid_only") +
         delivery_range[delivery]) *
      map_dbl(model, ~ pyrethroid_better[[.x]]) - 1,
    text = sprintf("%s, D_%s: %.1f", model, delivery, 100 * gain)
  ) %>%
  pull(text) %>%
  cat(sep = "\n")


# panel A: Declining effectiveness of nets
# mean and CI of "individual_net_impact" over time


# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_a_data_past <- panel_a_data %>%
  filter(
    year <= latest_nets
  )

panel_a_data_future <- panel_a_data %>%
  filter(
    year >= latest_nets
  )

panel_a <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    name = "LLIN effectiveness"
  ) +
  # past data
  geom_ribbon(
    data = panel_a_data_past,
    fill = "#56B1F7"
  ) +
  geom_line(
    data = panel_a_data_past,
    colour = grey(0.2),
    # linewidth = 0.8
  ) + 
  # future projections; 1/3 of the alpha and dashed lines
  geom_ribbon(
    data = panel_a_data_future,
    fill = "#56B1F7",
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = panel_a_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3,
    # linewidth = 0.8
  ) + 
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  xlab("") +
  theme_minimal()


# B) Effective coverage by net type: stacked areas of dual-AI, PBO, then
# pyrethroid-only nets, with the interval of the total, and the line of the
# counterfactual without resistance above. The gap between the line and the
# top of the stack is the coverage lost to resistance

net_type_labels <- c(pyrethroid = "Pyrethroid-only",
                     pbo = "PBO",
                     dual_ai = "Dual-AI")
# pyrethroid-only nets blue, dual-AI nets maroon, and PBO nets, a less good
# next-generation net, the maroon greyed. None is close to the insecticide
# colours of the other figures (insecticide_colours(), R/functions.R)
dual_ai_col <- "#AE017E"
net_type_cols <- c(pyrethroid = "#56B1F7",
                   pbo = colorRampPalette(c(dual_ai_col, grey(0.8)))(101)[66],
                   dual_ai = dual_ai_col)

# the bands of the stack, from the bottom
panel_b_stack <- panel_b_data %>%
  select(year, dual_ai, pbo, pyrethroid) %>%
  pivot_longer(
    cols = -year,
    names_to = "type",
    values_to = "coverage"
  ) %>%
  mutate(
    type = factor(type,
                  levels = c("dual_ai", "pbo", "pyrethroid"))
  ) %>%
  arrange(year, type) %>%
  group_by(year) %>%
  mutate(
    top = cumsum(coverage),
    bottom = top - coverage
  ) %>%
  ungroup()

# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_b_data_past <- panel_b_data %>%
  filter(
    year <= latest_nets
  )

panel_b_data_future <- panel_b_data %>%
  filter(
    year >= latest_nets
  )

panel_b_stack_past <- panel_b_stack %>%
  filter(
    year <= latest_nets
  )

panel_b_stack_future <- panel_b_stack %>%
  filter(
    year >= latest_nets
  )

panel_b_final <- panel_b_data %>%
  filter(
    year == max(year)
  )

panel_b <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.6),
    name = "Effective coverage"
  ) +
  scale_fill_manual(
    values = net_type_cols,
    labels = net_type_labels,
    breaks = c("pyrethroid", "pbo", "dual_ai"),
    name = NULL
  ) +
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  # the stack, past then future; future at half the alpha
  geom_ribbon(
    data = panel_b_stack_past,
    aes(
      y = NULL,
      ymin = bottom,
      ymax = top,
      fill = type
    )
  ) +
  geom_ribbon(
    data = panel_b_stack_future,
    aes(
      y = NULL,
      ymin = bottom,
      ymax = top,
      fill = type
    ),
    alpha = 0.5
  ) +
  # the interval of the total
  geom_ribbon(
    data = panel_b_data_past,
    fill = grey(0.2),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = panel_b_data_future,
    fill = grey(0.2),
    alpha = 0.15
  ) +
  geom_line(
    data = panel_b_data_past,
    colour = grey(0.2)
  ) +
  geom_line(
    data = panel_b_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3
  ) +
  # dotted lines between the bands in the projection
  geom_line(
    data = filter(panel_b_stack_future, type != "pyrethroid"),
    aes(
      y = top,
      ymin = NULL,
      ymax = NULL,
      group = type
    ),
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3
  ) +
  # effective coverage with no IR, past then future
  geom_line(
    data = panel_b_data_past,
    aes(
      y = mean_no_ir,
      colour = "No resistance"
    )
  ) +
  geom_line(
    data = panel_b_data_future,
    aes(
      y = mean_no_ir,
      colour = "No resistance"
    ),
    linetype = 3
  ) +
  scale_colour_manual(
    values = c("No resistance" = grey(0.5)),
    name = NULL
  ) +
  xlab("") +
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = rel(0.7)),
    legend.spacing.y = unit(0, "cm"),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    colour = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

# C) Spatial variation in effectiveness: plot the classes of susceptibility
# (2025) by the better buy, from the break-even susceptibility, on a map


# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)

# the colours of panel B's pyrethroid-only and dual-AI nets for the classes
# where they are the most cost-effective, and for the uncertain class a colour
# halfway between them, lightened halfway to white. The thresholds are in the
# log, for the caption
pyrethroid_col <- net_type_cols[["pyrethroid"]]
next_generation_col <- net_type_cols[["dual_ai"]]
uncertain_col <- colorRampPalette(
  c(colorRampPalette(c(pyrethroid_col, next_generation_col))(3)[2], "white")
)(3)[2]
band_cols <- c(pyrethroid = pyrethroid_col,
               uncertain = uncertain_col,
               next_generation = next_generation_col)
band_labels <- c(pyrethroid = "Pyrethroid-only",
                 uncertain = "Uncertain",
                 next_generation = "Next-generation")

panel_c <- ggplot() +
  africa_bg +
  geom_spatraster(
    data = bands_mask,
  ) +
  geom_sf(data = borders,
          col = border_col,
          linewidth = 0.1,
          fill = "transparent") +
  scale_fill_manual(
    values = band_cols,
    labels = band_labels,
    breaks = names(band_cols),
    na.translate = FALSE,
    # note we arre using element_markdown() for legend titles so need this html
    # linebreak version
    name = "Most cost-effective<br>(2025)",
    na.value = "transparent"
  ) +
  # no padding around Africa, so it fills the panel
  coord_sf(expand = FALSE) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0),
    plot.margin = unit(rep(0, 4), "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.22),
    legend.key.height = rel(0.8),
    legend.key.width = rel(0.8),
    legend.text.position = "left",
    legend.text = element_text(size = rel(0.7),
                               hjust = 1),
    legend.ticks = element_blank(),
    legend.title = element_markdown(size = rel(0.8),
                                    hjust = 1)
  )

# This displays the continent in classes of susceptibility: where
# pyrethroid-only nets give more effective coverage than dual-AI nets for the
# same budget, under either effectiveness model and any delivery cost in the
# range; where dual-AI nets do; and where it is uncertain





# panel D should show effectiveness in these bands, in these colours, like panel
# A



# split into past and future, both including the latest year of data so there is
# no break in lines or ribbons
panel_d_data_past <- panel_d_data %>%
  filter(
    year <= latest_nets
  )

panel_d_data_future <- panel_d_data %>%
  filter(
    year >= latest_nets
  )

panel_d <- ggplot(
  mapping = aes(
    x = year,
    y = mean,
    ymin = lower,
    ymax = upper,
    fill = band
  )
) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1),
    name = "LLIN effectiveness"
  ) +
  scale_fill_manual(
    values = band_cols,
    guide = "none",
    na.translate = FALSE,
    name = "LLIN effectiveness",
    na.value = "transparent",
  ) +
  # past data
  geom_ribbon(
    data = panel_d_data_past,
  ) +
  geom_line(
    data = panel_d_data_past,
    colour = grey(0.2)
  ) + 
  # future projections; 1/3 of the alpha and dashed lines
  geom_ribbon(
    data = panel_d_data_future,
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = panel_d_data_future,
    colour = grey(0.2),
    alpha = 0.5,
    linetype = 3
  ) + 
  # vertical dividing line for past and future
  geom_vline(
    xintercept = latest_nets,
    color = grey(0.6),
    linetype = 2
  ) +
  xlab("") +
  theme_minimal()

plots_all <- panel_a +
  panel_b +
  panel_c +
  panel_d +
  plot_layout(
    ncol = 2,
    nrow = 2,
    widths = rep(1, 4)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  filename = "figures/epi_impact.png",
  plot = plots_all,
  bg = "white",
  width = 7,
  height = 6,
  scale = 1,
  dpi = 300
)


# # maybe plot the killing effects from Nash & Sherrard-Smith
# 
# # - Nash/Sherrard-Smith killing effect estimates.
# 
# # load in and prepare curves from Nash/Sherard-Smith, given in supp to Sherard
# # Smith paper: https://doi.org/10.1038/s41467-022-30700-1
# ss_est <- read_excel("data/raw/41467_2022_30700_MOESM3_ESM.xlsx",
#                      sheet = 4,
#                      range = "A15:Y117") %>%
#   # get required columns
#   select(
#     resistance = "Pyrethroid resistance approximated from susceptibility bioassay",
#     kill_logistic = "dN0_1",  # logistic model
#     kill_log_logistic = "dN0_4" # log-logistic model
#   ) %>%
#   # remove empty row
#   slice(
#     -1
#   ) %>%
#   # split apart by model type
#   pivot_longer(
#     cols = starts_with("kill"),
#     names_to = "model",
#     values_to = "kill",
#     names_prefix = "kill_"
#   ) %>%
#   group_by(
#     model
#   ) %>%
#   mutate(
#     susceptibility = 1 - resistance,
#     relative_killing_effect = kill / max(kill)
#   )
# 
# # # plot to check against paper figures
# # ss_est %>%
# #   ggplot(
# #     aes(
# #       x = resistance,
# #       y = kill
# #     )
# #   ) +
# #   facet_wrap(~model) +
# #   geom_line() +
# #   coord_cartesian(ylim = c(0, 1)) +
# #   theme_minimal()
# 
# ss_est_logistic <- ss_est %>%
#   filter(model == "logistic")
# 
# ss_est_log_logistic <- ss_est %>%
#   filter(model == "log_logistic")
# 
# # produce interpolating functions to get the relative killing effect
# killing_function_logistic <- splinefun(
#   x = ss_est_logistic$susceptibility,
#   y = ss_est_logistic$relative_killing_effect
# )
# 
# killing_function_log_logistic <- splinefun(
#   x = ss_est_log_logistic$susceptibility,
#   y = ss_est_log_logistic$relative_killing_effect
# )
# 
# # check the interpolators behave appropriately
# stopifnot(
#   identical(
#     killing_function_logistic(ss_est_logistic$susceptibility),
#     ss_est_logistic$relative_killing_effect
#   )
# )
# 
# stopifnot(
#   identical(
#     killing_function_log_logistic(ss_est_log_logistic$susceptibility),
#     ss_est_log_logistic$relative_killing_effect
#   )
# )
# 
# 
# # make a cube of the killing effect of ITNs, over time and space
# relative_killing_logistic <- terra::app(ir, killing_function_logistic)
# names(relative_killing_logistic) <- names(ir)
# 
# relative_killing_log_logistic <- terra::app(ir, killing_function_log_logistic)
# names(relative_killing_log_logistic) <- names(ir)
# 
# 
# # average the killing effect and epi impact of nets over the nets in use
# average_killing_logistic <- global(
#   x = relative_killing_logistic * net_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# average_killing_log_logistic <- global(
#   x = relative_killing_log_logistic * net_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# # apply the impact of IR and recalculate average
# net_impact_killing_logistic <- global(
#   x = nets * relative_killing_logistic * infection_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# net_impact_killing_log_logistic <- global(
#   x = nets * relative_killing_log_logistic * infection_weights,
#   fun = "sum",
#   na.rm = TRUE
# )
# 
# 
# par(mfrow = c(1, 2),
#     mar = c(2, 5, 3, 2))
# 
# # Effectiveness\n of single-AI nets
# no_change <- rep(1, length(years_plot))
# # impact of each net
# plot(no_change ~ years_plot,
#      type = "l",
#      col = grey(0.3),
#      lwd = 2,
#      xlab = "",
#      ylab = "Relative effectiveness",
#      ylim = c(0, 1))
# lines(average_killing_logistic$sum ~ years_plot,
#       col = "red",
#       lwd = 2)
# lines(average_killing_log_logistic$sum ~ years_plot,
#       col = "red",
#       lty = 2,
#       lwd = 2)
# lines(average_epi_impact$sum ~ years_plot,
#       col = "blue",
#       lwd = 2)
# abline(v = latest_nets, lty = 2)
# # directly label these
# 
# plot(net_impact_no_ir$sum ~ years_plot,
#      type = "l",
#      lwd = 2,
#      col = grey(0.3),
#      ylab = "Relative impact",
#      xlab = "",
#      ylim = c(0, 1))
# lines(net_impact_epi_impact$sum ~ years_plot,
#       col = "blue",
#       lwd = 2)
# abline(v = latest_nets, lty = 2)
# title(main = "Impact of\nnet distribution")
# 
# df_plot <- tibble(
#   year = years_plot,
#   `impact_direct/no resistance` = 1,
#   `impact_direct/vector killing (logistic)` = average_killing_logistic$sum,
#   `impact_direct/vector killing (log-logistic)` = average_killing_log_logistic$sum,
#   `impact_direct/transmission reduction` = average_epi_impact$sum,
#   `impact_overall/no resistance` = net_impact_no_ir$sum,
#   `impact_overall/transmission reduction` = net_impact_epi_impact$sum
# ) %>%
#   pivot_longer(
#     cols = starts_with("impact"),
#     names_to = c("scale", "impact_type"),
#     names_pattern = "impact_(.*)/(.*)",
#     values_to = "impact",
#     # names_prefix = "impact_"
#   )
# 
# plot_direct <- df_plot %>%
#   filter(
#     scale == "direct"
#   ) %>%
#   mutate(
#     phase = case_when(
#       year > latest_nets ~ "future",
#       .default = "past"
#     )
#   ) %>%
#   ggplot(
#     aes(
#       x = year,
#       y = impact,
#       colour = impact_type,
#       linetype = phase
#     )
#   ) +
#   geom_vline(
#     aes(
#       xintercept = latest_nets
#     ),
#     col = grey(0.6),
#     linetype = 2
#   ) +
#   geom_line(
#     linewidth = 0.8
#   ) +
#   scale_y_continuous(
#     labels = scales::percent,
#     limits = c(0, 1)
#   ) +
#   scale_colour_manual(
#     values = c(
#       "transmission reduction" = "blue",
#       "vector killing (log-logistic)" = "coral",
#       "vector killing (logistic)" = "coral3",
#       "no resistance" = grey(0.4)
#     )
#   ) +
#   scale_linetype_manual(
#     values = c(
#       "past" = 1,
#       "future" = 3
#     ),
#     guide = "none"
#   ) +
#   ylab(
#     "Effectiveness"
#   ) +
#   xlab(
#     ""
#   ) +
#   theme_minimal() +
#   theme(
#     legend.title = element_blank(),
#     # legend.position = "bottom"
#   )
# 
# plot_direct
