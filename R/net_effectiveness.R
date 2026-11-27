# The effectiveness of pyrethroid-only nets against resistance, and the
# break-even effectiveness and susceptibility at which they and dual-AI nets
# give the same effective coverage for the same budget (#75). Used by
# R/fig_ento_epi_impact.R and R/fig_effectiveness_curves.R; needs the
# tidyverse (R/packages.R)

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
# D = (Int$ - USD) / (k - 1). Inflated from 2016 to 2026 USD with US CPI-U;
# the break-even is at the mean of the two sites
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
delivery <- mean(campaigns$delivery)
break_even <- break_even_effectiveness(
  delivery = delivery,
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
# relative to that in the same setting at k = 1 (mint_settings); the median
# over the settings (mint_effectiveness)
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
mint_settings <- mint_index %>%
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
  mutate(effectiveness = effect / effect[susceptibility == 1],
         setting = cur_group_id()) %>%
  ungroup()
mint_effectiveness <- mint_settings %>%
  group_by(susceptibility) %>%
  summarise(effectiveness = median(effectiveness))

# Encode Symons et al. PMI function for relative effectiveness of ITNs at
# reducing logit-PfPR (reduction in on the log-odds ratio), compared to in the absence
# of resistance.

# when susceptibility is 100%, the ITN effectiveness parameter has full effect
# (multiplier 1), when susceptibility is 0%, the ITN effectiveness parameter has
# 54% (multiplier 0.54).

# The log-odds ratio for malaria infection under full ITN coverage (relative to
# no ITN coverage) is -1.68, which corresponds to an odds ratio of 0.186374.
# With a fully-resistant vector population, this is would reduce to 54% of its
# previous value, -0.9072, giving an odds ratio of 0.4036529.

# Alternative interpretation: the odds ratio for a population with 100% coverage
# and 100% resistance is equivalent to that of a population with 54% coverage
# and 0% resistance.

# the posteriors of the parameters behind this model:

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

# which is used at the posterior means, r = beta_1 / beta_0 = -0.46. The
# uncertainty shown is that of the susceptibility maps (compute_epi_draws()).
# g*(k) is linear in k, so the effectiveness at a weighted mean of k is the
# weighted mean of the effectiveness

# Two models of the relative effectiveness e(k) of a pyrethroid-only net, as a
# function of susceptibility k (bioassay mortality), 1 at k = 1, both the
# nets' effect on the log-odds of prevalence, as the points of piecewise
# linear curves:
#   symons    Symons et al., at the posterior mean of their coefficients
#             (above): linear, 1 - 0.46 (1 - k)
#   churcher  Churcher et al. 2024's MINT simulations, above, interpolated
#             linearly between the simulated k
effectiveness_points <- list(
  symons = tibble(susceptibility = c(0, 1),
                  effectiveness = 1 - 0.46 * (1 - c(0, 1))),
  churcher = mint_effectiveness
)
effectiveness_models <- lapply(effectiveness_points, function(points) {
  approxfun(x = points$susceptibility, y = points$effectiveness)
})

# the susceptibility at which a model's effectiveness is e: 0 where it is
# above e at any susceptibility
susceptibility_at <- function(model, e) {
  if (model(0) >= e) {
    return(0)
  }
  uniroot(function(k) model(k) - e, c(0, 1), tol = 1e-8)$root
}

# The break-even susceptibility k*, where e(k*) = e*(D), under each model.
# Pyrethroid-only nets are the better buy above it, dual-AI nets below. The
# figure uses Symons et al.'s model, the supplement's version Churcher et
# al.'s
break_even_susceptibility <- vapply(effectiveness_models, susceptibility_at,
                                    numeric(1), e = break_even)
