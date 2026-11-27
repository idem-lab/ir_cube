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
source("R/bioassay_subset.R")
source("R/two_stage_predictions.R")
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


# The deterministic parts of the panels: effective coverage of the PBO and
# dual-AI nets (assumed fully effective), the pyrethroid-only nets' use, and
# all nets' effective coverage without resistance, averaged over the
# population at risk of each year; each type's share of net use in 2024; and
# the year of the last net use data. From net use by type (net_type_files),
# PfPR 2000 and population, at the cells with a susceptibility map (all but
# 1e-5 of the population at risk)
compute_epi_impact <- function(years_plot, net_type_files) {
  # write the cubes terra computes to temporary files: held in memory, they
  # peaked at 26 GB; on disk, at 8 GB, for a minute and a half longer
  terraOptions(todisk = TRUE)

  # - IR layers (susceptibility to the pyrethroids used in LLINs), of the
  # two-stage model (R/two_stage_maps.R), for the cells they cover
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

  # convert a raster layer into spatial weights (pixels sum to 1) to get averages
  # over space
  normalise <- function(x) {
    total <- global(x, "sum", na.rm = TRUE)[1, 1]
    x / total
  }

  # cube of population at risk (potential infected population, in absence of
  # interventions) weights, by year: population times PfPR in 2000
  infection_weights <- sapp(terra::mask(pop_all * pr2000, ir), normalise)

  # Effective coverage: net use times effectiveness, averaged across Africa
  # over the population at risk, for each group of net types. PBO and dual-AI
  # nets are assumed to be fully effective (effectiveness 1), losing none of it
  # to resistance; pyrethroid-only nets' effectiveness is computed from the
  # posterior (compute_epi_draws()). Without resistance, all nets would be
  # fully effective, and effective coverage would be the average net use. If
  # net use was 100% for all populations at risk, this would equal 1. It is
  # lower due to imperfect net coverage (but not due to IR).
  coverage <- function(x) {
    global(x * infection_weights, fun = "sum", na.rm = TRUE)$sum
  }
  coverage_pyrethroid_use <- coverage(nets_by_type$pyrethroid)
  coverage_pbo <- coverage(nets_by_type$pbo)
  coverage_dual_ai <- coverage(nets_by_type$dual_ai)

  # each type's share of net use in 2024, Africa-wide: the people using a net
  # of the type over the people using any
  use_2024 <- vapply(nets_by_type, function(nets) {
    global(nets[["2024"]] * pop_all[["2024"]], "sum", na.rm = TRUE)$sum
  }, numeric(1))

  list(coverage = tibble(
         year = years_plot,
         pyrethroid_use = coverage_pyrethroid_use,
         pbo = coverage_pbo,
         dual_ai = coverage_dual_ai,
         no_ir = coverage_pyrethroid_use + coverage_pbo + coverage_dual_ai
       ),
       type_shares_2024 = use_2024 / sum(use_2024),
       latest_nets = latest_nets)
}

# The two-stage posterior of llin_effective susceptibility k as the panels
# weight it. The draws are those `map` of R/two_stage_maps.R summarises: 1000
# paired draws of each LLIN pyrethroid's two-stage model
# (R/two_stage_predictions.R), combined draw by draw with the llin_effective
# weights (llin_weights), at every cell with a population at risk in any year
# or within the limits of transmission. For each draw and year:
#   a        the mean of k over the population at risk in 2000 (panel A)
#   b        the sum of k times the pyrethroid-only nets' use and the year's
#            population at risk weights (panel B); b_weight the sum of the
#            weights
#   d        the sum of k times the population at risk in 2000 weights in each
#            class of 2025 susceptibility, under each rule (panel D);
#            d_weight the sum of the weights. The classes are those above
#            the highest break-even susceptibility k* (break_even), below the
#            lowest, and between; the rule "mean" puts a cell in a class by its
#            posterior mean, and rule p by the posterior probability p (in
#            `confidence`) of being above the highest or below the lowest
#            k*, uncertain otherwise
# and for each cell, its population at risk in 2025 and the posterior mean,
# quantiles (probs) and probabilities of being above the highest k* and below
# the lowest of its 2025 susceptibility (panel C). The cells are processed in
# forked workers (IR_CUBE_MAP_WORKERS; see below), in chunks of chunk_size
# cells, each holding all the draws of its cells; the results do not depend
# on either
compute_epi_draws <- function(years_plot, net_type_files, llin_weights,
                              break_even, confidence, probs,
                              chunk_size = 2500) {
  terraOptions(todisk = TRUE)
  baseline_year <- 1995
  years_all <- baseline_year:max(years_plot)
  year_index <- match(years_plot, years_all)
  index_2025 <- match(2025, years_plot)

  # the cells, and their weights: population at risk as compute_epi_impact()
  ir <- rast(ir_map_files("llin_effective", years_plot))
  names(ir) <- years_plot
  pr2000 <- rast("data/clean/pfpr_2000.tif")
  pop_all <- c(rast("data/clean/pop_cube.tif"),
               rast("data/clean/pop_cube_future.tif"))
  names(pop_all) <- str_remove(names(pop_all), "^pop_")
  population_at_risk <- terra::mask(pop_all[[as.character(years_plot)]] *
                                      pr2000, ir)
  nets <- rast(net_type_files[["pyrethroid"]])
  names(nets) <- str_remove(names(nets), "^nets_")
  nets <- post_pad_cube(nets, max(years_plot))[[as.character(years_plot)]]
  selected <- terra::app(population_at_risk, sum, na.rm = TRUE) > 0 |
    !is.na(rast("data/clean/pfpr_water_mask.tif"))
  selected <- terra::mask(selected, ir[["2025"]])
  cells <- which(terra::values(selected, mat = FALSE) == 1)
  country <- as.character(terra::extract(
    rast("data/clean/country_raster.tif"), cells)$country_name)
  cells <- cells[!is.na(country)]
  country <- country[!is.na(country)]
  at_cells <- function(x) {
    values <- as.matrix(terra::extract(x, cells))
    values[is.na(values)] <- 0
    values
  }
  people <- at_cells(population_at_risk)
  stopifnot(isTRUE(all.equal(
    colSums(people),
    global(population_at_risk, "sum", na.rm = TRUE)$sum,
    check.attributes = FALSE)))
  weight_a <- people[, "2000"] / sum(people[, "2000"])
  weight_b <- sweep(people, 2, colSums(people), "/") * at_cells(nets)
  people_2025 <- people[, "2025"]
  map_2025 <- at_cells(ir[["2025"]])[, 1]
  rm(people, population_at_risk, nets, selected, ir, pop_all, pr2000)

  # the draws of each pyrethroid (the map's: same draws, batches and seeds)
  data <- modelled_bioassays(baseline_year)$df
  setups <- lapply(names(llin_weights), function(type) {
    setup <- two_stage_setup(type, years = years_all, data = data,
                             n_draws = 1000, batch_size = 100,
                             baseline_year = baseline_year)
    message(format(Sys.time(), "%H:%M:%S"), " set up ", type)
    setup
  })
  rm(data)
  stopifnot(length(unique(lapply(setups, `[[`, "design"))) == 1,
            length(unique(lapply(setups, `[[`, "batches"))) == 1)
  n_draws <- setups[[1]]$n_draws

  # the classes of each cell under each rule, from the draws of its 2025
  # susceptibility (draws x cells): 1 pyrethroid, 2 uncertain, 3
  # next-generation
  class_names <- c("pyrethroid", "uncertain", "next_generation")
  rules <- c("mean", as.character(confidence))
  classify <- function(k_2025) {
    mean <- colMeans(k_2025)
    above <- colMeans(k_2025 > break_even[["high"]])
    below <- colMeans(k_2025 < break_even[["low"]])
    by_mean <- 3 - (mean >= break_even[["low"]]) -
      (mean > break_even[["high"]])
    by_probability <- sapply(confidence, function(p) {
      if_else(above >= p, 1, if_else(below >= p, 3, 2))
    })
    classes <- cbind(by_mean, by_probability)
    colnames(classes) <- paste0("class_", rules)
    list(classes = classes,
         stats = cbind(mean = mean, above = above, below = below,
                       t(apply(k_2025, 2, quantile, probs)), classes))
  }

  # the sums of a group of cells (rows of the cells above): the draws of all
  # their cells and years, a chunk of cells at a time, then weighted and summed
  group_sums <- function(group) {
    group_cells <- two_stage_cells(setups[[1]], cells[group], country[group])
    chunks <- split(seq_along(group), ceiling(seq_along(group) / chunk_size))
    sums <- list(
      a = matrix(0, n_draws, length(years_plot)),
      b = matrix(0, n_draws, length(years_plot)),
      d = array(0, c(n_draws, length(rules) * 3, length(years_plot))),
      d_weight = numeric(length(rules) * 3)
    )
    stats <- vector("list", length(chunks))
    for (i in seq_along(chunks)) {
      rows <- group[chunks[[i]]]
      chunk <- two_stage_chunk(setups[[1]], group_cells, chunks[[i]])
      k <- array(0, c(n_draws, length(rows), length(years_plot)))
      start <- 0
      for (b in seq_along(setups[[1]]$batches)) {
        index <- start + seq_along(setups[[1]]$batches[[b]])
        for (t in seq_along(setups)) {
          logit <- two_stage_logit_batch(setups[[t]], b, chunk)$two_stage
          for (j in seq_along(years_plot)) {
            k[index, , j] <- k[index, , j] +
              llin_weights[[t]] * plogis(logit[[year_index[j]]])
          }
        }
        start <- max(index)
      }
      # check: these are the map's draws
      stopifnot(max(abs(colMeans(k[, , index_2025]) - map_2025[rows])) <
                  1e-5)
      classified <- classify(k[, , index_2025])
      class_weights <- do.call(cbind, lapply(seq_along(rules), function(r) {
        sapply(1:3, function(c) weight_a[rows] * (classified$classes[, r] == c))
      }))
      for (j in seq_along(years_plot)) {
        sums$a[, j] <- sums$a[, j] + k[, , j] %*% weight_a[rows]
        sums$b[, j] <- sums$b[, j] + k[, , j] %*% weight_b[rows, j]
        sums$d[, , j] <- sums$d[, , j] + k[, , j] %*% class_weights
      }
      sums$d_weight <- sums$d_weight + colSums(class_weights)
      stats[[i]] <- cbind(cell = cells[rows], people_2025 = people_2025[rows],
                          classified$stats)
      rm(k, chunk)
      if (i %% 5 == 0) {
        message(format(Sys.time(), "%H:%M:%S"), " cells ", min(group), "-",
                max(group), ": ", i, " of ", length(chunks), " chunks")
      }
    }
    sums$stats <- do.call(rbind, stats)
    sums
  }

  # Forked workers deadlocked on the Linux box (each waiting on a lock held
  # at the fork), so this ran in one process (IR_CUBE_MAP_WORKERS = 1): 3 h
  # with OpenBLAS (8 threads), 21 GB at peak. At 4 workers, memory also ran
  # short: each added a few GB to the parent's ~11 GB
  invisible(gc())
  n_workers <- as.integer(Sys.getenv("IR_CUBE_MAP_WORKERS", "1"))
  message(format(Sys.time(), "%H:%M:%S"), " drawing ", length(cells),
          " cells in ", n_workers, " workers")
  groups <- split(seq_along(cells),
                  ceiling(seq_along(cells) /
                            ceiling(length(cells) / n_workers)))
  done <- parallel::mclapply(groups, group_sums, mc.cores = n_workers,
                             mc.preschedule = FALSE)
  failed <- !vapply(done, is.list, logical(1))
  if (any(failed)) {
    stop("drawing failed for ", sum(failed), " groups of cells: ",
         paste(unique(unlist(done[failed])), collapse = "; "))
  }
  add <- function(name) Reduce(`+`, lapply(done, `[[`, name))
  d_names <- paste(rep(rules, each = 3), class_names, sep = ":")
  list(
    years = years_plot,
    a = `colnames<-`(add("a"), years_plot),
    b = `colnames<-`(add("b"), years_plot),
    b_weight = colSums(weight_b),
    d = `dimnames<-`(add("d"), list(NULL, d_names, years_plot)),
    d_weight = setNames(add("d_weight"), d_names),
    cells = as_tibble(do.call(rbind, lapply(done, `[[`, "stats")))
  )
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
# nets' effect on the log-odds of prevalence:
#   symons    Symons et al., at the posterior mean of their coefficients
#             (above)
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
  args = list(years_plot = years_plot, net_type_files = net_type_files),
  sources = c("R/packages.R", "R/functions.R", "R/model_covariates.R"),
  inputs = c(ir_map_files("llin_effective", years_plot), net_type_files,
             "data/clean/pfpr_2000.tif", "data/clean/pop_cube.tif",
             "data/clean/pop_cube_future.tif")
)

# the posterior of susceptibility, as the panels weight it (3 h in one
# process; see compute_epi_draws())
llin_weights <- unlist(readRDS("temporary/ingredient_weights.RDS"))
draws <- cached_figure_data(
  "temporary/fig_ento_epi_impact_draws.rds", compute_epi_draws,
  args = list(years_plot = years_plot, net_type_files = net_type_files,
              llin_weights = llin_weights,
              break_even = break_even_susceptibility,
              confidence = c(0.8, 0.9, 0.95),
              probs = c(0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975)),
  sources = c("R/packages.R", "R/functions.R", "R/model_covariates.R",
              "R/bioassay_subset.R", "R/two_stage_predictions.R"),
  inputs = c(ir_map_files("llin_effective", years_plot), net_type_files,
             two_stage_figure_inputs(names(llin_weights)))
)
latest_nets <- epi$latest_nets

# The rule for the classes of panels C and D: "mean", by the posterior mean
# 2025 susceptibility, or a confidence (0.8, 0.9 or 0.95), by the posterior
# probability of being above the highest break-even susceptibility or below
# the lowest (compute_epi_draws()); and the credible intervals of the ribbons
classification <- "mean"
interval <- 0.95

# the posterior mean and credible interval of draws x years
posterior_summary <- function(x) {
  tibble(
    year = years_plot,
    mean = colMeans(x),
    lower = apply(x, 2, quantile, (1 - interval) / 2),
    upper = apply(x, 2, quantile, 1 - (1 - interval) / 2)
  )
}
symons <- effectiveness_models$symons

# panel A: the effectiveness of a pyrethroid-only net at the mean
# susceptibility over the population at risk in 2000
panel_a_data <- posterior_summary(symons(draws$a))

# panel B: effective coverage by net type, its total and interval, and the
# counterfactual without IR. The pyrethroid-only nets' is the sum of their
# weights times the effectiveness at the weighted mean susceptibility
stopifnot(isTRUE(all.equal(draws$b_weight, epi$coverage$pyrethroid_use,
                           check.attributes = FALSE, tolerance = 1e-6)))
coverage_pyrethroid <- sweep(
  symons(sweep(draws$b, 2, draws$b_weight, "/")), 2, draws$b_weight, "*"
)
coverage_total <- sweep(coverage_pyrethroid, 2,
                        epi$coverage$pbo + epi$coverage$dual_ai, "+")
panel_b_data <- epi$coverage %>%
  mutate(
    pyrethroid = colMeans(coverage_pyrethroid),
    mean = pyrethroid + pbo + dual_ai,
    mean_no_ir = no_ir
  ) %>%
  bind_cols(select(posterior_summary(coverage_total), lower, upper))

# panel C: the classes of each cell under the rule, within the limits of
# transmission
class_names <- c("pyrethroid", "uncertain", "next_generation")
bands <- rast(rast("data/clean/raster_mask.tif"))
bands[draws$cells$cell] <- draws$cells[[paste0("class_", classification)]]
levels(bands) <- data.frame(id = 1:3, band = class_names)
bands_mask <- terra::mask(bands, rast("data/clean/pfpr_water_mask.tif"))

# panel D: the effectiveness of a pyrethroid-only net in each class, at the
# mean susceptibility over the class's population at risk in 2000
panel_d_data <- bind_rows(lapply(class_names, function(band_name) {
  column <- paste(classification, band_name, sep = ":")
  posterior_summary(symons(draws$d[, column, ] / draws$d_weight[[column]])) %>%
    mutate(band = band_name)
})) %>%
  mutate(band = factor(band, levels = class_names))


# numbers for the text --------------------------------------------------------

# percentages: a mean with its credible interval
with_interval <- function(mean, lower, upper) {
  sprintf("%.1f [%.1f, %.1f]", 100 * mean, 100 * lower, 100 * upper)
}

cat(sprintf(paste0("\nAfrica-wide effective coverage (%% of the population at",
                   " risk, as panel B; posterior mean [%.0f%% credible",
                   " interval]):\n"), 100 * interval))
panel_b_data %>%
  filter(year %in% c(2017, 2024)) %>%
  mutate(lost = mean_no_ir - mean) %>%
  rowwise() %>%
  mutate(
    text = paste0(
      year, ": pyrethroid-only ", sprintf("%.1f", 100 * pyrethroid),
      ", PBO ", sprintf("%.1f", 100 * pbo),
      ", dual-AI ", sprintf("%.1f", 100 * dual_ai),
      "\n  total ", with_interval(mean, lower, upper),
      "\n  no resistance ", sprintf("%.1f", 100 * mean_no_ir),
      "\n  lost to resistance (percentage points) ",
      with_interval(lost, mean_no_ir - upper, mean_no_ir - lower),
      "\n  lost to resistance (% of no resistance) ",
      with_interval(lost / mean_no_ir, 1 - upper / mean_no_ir,
                    1 - lower / mean_no_ir)
    )
  ) %>%
  pull(text) %>%
  cat(sep = "\n")

cat("\nEffectiveness of a pyrethroid-only LLIN, weighted by the population at",
    "risk in 2000 (%, panel A):\n")
panel_a_data %>%
  filter(year %in% c(2017, 2024, 2030)) %>%
  mutate(text = paste0(year, ": ", with_interval(mean, lower, upper))) %>%
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

# the classes' shares of the population at risk in 2025 under each rule, and
# their mean 2025 effectiveness over it under each model, at the posterior
# mean susceptibility
band_summary <- draws$cells %>%
  select(people_2025, mean, starts_with("class_")) %>%
  pivot_longer(starts_with("class_"), names_to = "rule",
               names_prefix = "class_", values_to = "class") %>%
  mutate(band = class_names[class]) %>%
  group_by(rule, band) %>%
  summarise(
    people = sum(people_2025),
    symons = weighted.mean(effectiveness_models$symons(mean), people_2025),
    churcher = weighted.mean(effectiveness_models$churcher(mean),
                             people_2025),
    .groups = "drop"
  ) %>%
  group_by(rule) %>%
  mutate(share = people / sum(people)) %>%
  ungroup()

cat("\nClasses of 2025 susceptibility by rule (share: % of the population at",
    "risk in 2025; effectiveness: its mean over them under each model at the",
    "posterior mean susceptibility, %):\n")
band_summary %>%
  mutate(text = sprintf(paste("%s, %s: share %.1f, effectiveness %.1f",
                              "(Symons), %.1f (Churcher)"),
                        rule, band, 100 * share, 100 * symons,
                        100 * churcher)) %>%
  pull(text) %>%
  cat(sep = "\n")

# the gain in effective coverage from buying pyrethroid-only instead of dual-AI
# nets with the same budget, where pyrethroid-only nets are the better buy
# under the rule: (P_dual + D) / (P_pyr + D) e - 1, linear in e, so its mean
# over the population at risk is that at the class's mean effectiveness
pyrethroid_better <- band_summary %>%
  filter(rule == classification, band == "pyrethroid")
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
# halfway between them, lightened three quarters of the way to white. The
# thresholds are in the log, for the caption
pyrethroid_col <- net_type_cols[["pyrethroid"]]
next_generation_col <- net_type_cols[["dual_ai"]]
uncertain_col <- colorRampPalette(
  c(colorRampPalette(c(pyrethroid_col, next_generation_col))(3)[2], "white")
)(5)[4]
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
