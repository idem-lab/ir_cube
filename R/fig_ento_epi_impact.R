# Plot the decline (2000-2025 and projected 2025-2030) in the effectiveness of
# pyrethroid-only ITNs across Africa due to IR, its cost to the effective
# coverage of the nets in use, and where pyrethroid-only or dual-AI nets are
# the better buy (#75). Effectiveness is the relative epidemiological
# effectiveness g*(k) of Symons et al. (below).

# three panel plot:

# A) Effective coverage, 2000-2025 (the years of net use data): net use times
# effectiveness, summed over three groups of net types
# (R/prep_net_use_pyrethroid.R) and averaged across Africa over the population
# at risk of each year (the year's population times PfPR in
# 2000), as stacked areas: dual-AI nets, PBO nets, then pyrethroid-only nets
# (LLINs, and conventional ITNs counting as 0.25 of an LLIN). PBO and dual-AI
# nets are assumed to be fully effective (effectiveness 1), losing none of it
# to resistance. The dotted lines are the counterfactual without resistance,
# all nets in use fully effective; the gap between the top line and the stack
# is the coverage lost to resistance

# B) Declining effectiveness of nets: timeseries with ribbon of uncertainty
# giving the effectiveness, relative to a fully susceptible population, that a
# pyrethroid-only LLIN would have for each person at risk, whether or not nets
# are in use there: averaged across Africa over the population at risk
# (population times PfPR, both in 2000, fixed over time)

# C) Spatial variation in effectiveness: where pyrethroid-only or dual-AI nets
# give more effective coverage for the same budget in 2025
# (break_even_effectiveness()), at the mean delivery cost per net of two mass
# campaigns in Mali (Scates et al. 2020): pyrethroid-only nets where the
# posterior probability that they do is at least 0.9, dual-AI nets where it is
# at most 0.1, uncertain between. The supplement's version uses the
# effectiveness model of Churcher et al. 2024 instead

# B is projected to 2030 from the two-stage model's forecasts
# (compute_epi_impact() also holds net use at its 2025 values after 2025,
# unplotted). The ribbons are credible intervals of the
# two-stage posterior (compute_epi_draws()), with Symons et al.'s coefficients
# at their posterior means. The numbers for the text are printed after the
# panels' data.



# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/model_covariates.R")
source("R/bioassay_subset.R")
source("R/two_stage_predictions.R")
source("R/figure_cache.R")
# the effectiveness models and the break-even (net prices and delivery costs)
source("R/net_effectiveness.R")

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
# or within the limits of transmission. For each effectiveness model e(k)
# (effectiveness: the points of a piecewise linear curve), draw and year:
#   effectiveness  the mean of e(k) over the population at risk in 2000
#                  (panel B)
#   coverage       the effective coverage of the pyrethroid-only nets: the
#                  sum of e(k) times their use and the year's population at
#                  risk weights (panel A); use_weight the sum of the weights
# and for each cell, its population at risk in 2025, and the posterior mean
# of its 2025 susceptibility and the posterior probability that it is above
# each break-even susceptibility k* (break_even, one per effectiveness model;
# panel C). The cells are processed in forked workers (IR_CUBE_MAP_WORKERS;
# see below), in chunks of chunk_size cells, each holding all the draws of its
# cells; the results do not depend on either
compute_epi_draws <- function(years_plot, net_type_files, llin_weights,
                              effectiveness, break_even, chunk_size = 2500) {
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
  weight_2000 <- people[, "2000"] / sum(people[, "2000"])
  weight_use <- sweep(people, 2, colSums(people), "/") * at_cells(nets)
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

  # the posterior mean of each cell's 2025 susceptibility, and its posterior
  # probability of being above each k*, from its draws (draws x cells)
  exceedance <- function(k_2025) {
    above <- sapply(break_even, function(k_star) colMeans(k_2025 > k_star))
    colnames(above) <- paste0("above_", names(break_even))
    cbind(mean = colMeans(k_2025), above)
  }

  # the sums of a group of cells (rows of the cells above): the draws of all
  # their cells and years, a chunk of cells at a time, then weighted and summed
  group_sums <- function(group) {
    group_cells <- two_stage_cells(setups[[1]], cells[group], country[group])
    chunks <- split(seq_along(group), ceiling(seq_along(group) / chunk_size))
    empty <- lapply(effectiveness, function(points) {
      matrix(0, n_draws, length(years_plot))
    })
    sums <- list(effectiveness = empty, coverage = empty)
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
      for (j in seq_along(years_plot)) {
        for (m in names(effectiveness)) {
          e <- matrix(approx(effectiveness[[m]]$susceptibility,
                             effectiveness[[m]]$effectiveness,
                             xout = k[, , j])$y, n_draws)
          sums$effectiveness[[m]][, j] <- sums$effectiveness[[m]][, j] +
            e %*% weight_2000[rows]
          sums$coverage[[m]][, j] <- sums$coverage[[m]][, j] +
            e %*% weight_use[rows, j]
        }
      }
      stats[[i]] <- cbind(cell = cells[rows], people_2025 = people_2025[rows],
                          exceedance(k[, , index_2025]))
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
  add <- function(name) {
    lapply(setNames(nm = names(effectiveness)), function(m) {
      total <- Reduce(`+`, lapply(done, function(sums) sums[[name]][[m]]))
      `colnames<-`(total, years_plot)
    })
  }
  list(
    years = years_plot,
    effectiveness = add("effectiveness"),
    coverage = add("coverage"),
    use_weight = colSums(weight_use),
    cells = as_tibble(do.call(rbind, lapply(done, `[[`, "stats")))
  )
}

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

# the posterior of effectiveness and susceptibility, as the panels weight it
# (3 h in one process; see compute_epi_draws())
llin_weights <- unlist(readRDS("temporary/ingredient_weights.RDS"))
draws <- cached_figure_data(
  "temporary/fig_ento_epi_impact_draws.rds", compute_epi_draws,
  args = list(years_plot = years_plot, net_type_files = net_type_files,
              llin_weights = llin_weights,
              effectiveness = effectiveness_points,
              break_even = break_even_susceptibility),
  sources = c("R/packages.R", "R/functions.R", "R/model_covariates.R",
              "R/bioassay_subset.R", "R/two_stage_predictions.R"),
  inputs = c(ir_map_files("llin_effective", years_plot), net_type_files,
             two_stage_figure_inputs(names(llin_weights)))
)
latest_nets <- epi$latest_nets

# the credible intervals of the ribbons
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
# the data of panels A and B under each effectiveness model:
#   a  effective coverage by net type, its total and interval, and the
#      counterfactual without IR
#   b  the effectiveness of a pyrethroid-only net over the population at risk
#      in 2000
stopifnot(isTRUE(all.equal(draws$use_weight, epi$coverage$pyrethroid_use,
                           check.attributes = FALSE, tolerance = 1e-6)))
panel_data <- lapply(setNames(nm = names(effectiveness_models)), function(model) {
  coverage_pyrethroid <- draws$coverage[[model]]
  coverage_total <- sweep(coverage_pyrethroid, 2,
                          epi$coverage$pbo + epi$coverage$dual_ai, "+")
  list(
    a = epi$coverage %>%
      mutate(
        pyrethroid = colMeans(coverage_pyrethroid),
        mean = pyrethroid + pbo + dual_ai,
        mean_no_ir = no_ir
      ) %>%
      bind_cols(select(posterior_summary(coverage_total), lower, upper)),
    b = posterior_summary(draws$effectiveness[[model]])
  )
})

# panel C: the classes of 2025 susceptibility by the posterior probability p
# that pyrethroid-only nets are the better buy (susceptibility above the
# model's k*): pyrethroid-only where p is at least `confidence`, dual-AI where
# it is at most 1 - confidence, uncertain between; within the limits of
# transmission. Symons et al.'s for panel C, Churcher et al.'s for the
# supplement
confidence <- 0.9
class_names <- c("pyrethroid", "uncertain", "dual_ai")
better_buy_class <- function(model) {
  p <- draws$cells[[paste0("above_", model)]]
  if_else(p >= confidence, 1, if_else(p <= 1 - confidence, 3, 2))
}
better_buy_map <- function(model) {
  classes <- rast(rast("data/clean/raster_mask.tif"))
  classes[draws$cells$cell] <- better_buy_class(model)
  levels(classes) <- data.frame(id = 1:3, band = class_names)
  terra::mask(classes, rast("data/clean/pfpr_water_mask.tif"))
}


# numbers for the text --------------------------------------------------------

# percentages: a mean with its credible interval
with_interval <- function(mean, lower, upper) {
  sprintf("%.1f [%.1f, %.1f]", 100 * mean, 100 * lower, 100 * upper)
}

for (model in names(panel_data)) {
  cat(sprintf(paste0("\nAfrica-wide effective coverage under %s's model (%% of",
                     " the population at risk, as panel A; posterior mean [%.0f%% credible",
                     " interval]):\n"), model, 100 * interval))
  panel_data[[model]]$a %>%
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

  cat("\nEffectiveness of a pyrethroid-only LLIN under", model, "model,",
      "weighted by the population at risk in 2000 (%, panel B):\n")
  panel_data[[model]]$b %>%
    filter(year %in% c(2017, 2024, 2030)) %>%
    mutate(text = paste0(year, ": ", with_interval(mean, lower, upper))) %>%
    pull(text) %>%
    cat(sep = "\n")
}

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
cat(sprintf(paste0("Mean of the two sites, D: %.2f (2026 USD); break-even",
                   " effectiveness e*(D): %.4f\n"), delivery, break_even))

cat("\nBreak-even susceptibility k* (bioassay mortality, %), where e(k*) =",
    "e*(D), by effectiveness model:\n")
cat(sprintf("%s: %.1f", names(break_even_susceptibility),
            100 * break_even_susceptibility), sep = "\n")

# the classes' shares of the population at risk in 2025, under each model
class_shares <- map_dfr(names(effectiveness_models), function(model) {
  tibble(model = model,
         class = factor(class_names[better_buy_class(model)],
                        levels = class_names),
         people = draws$cells$people_2025)
}) %>%
  group_by(model, class) %>%
  summarise(people = sum(people), .groups = "drop_last") %>%
  mutate(share = people / sum(people)) %>%
  ungroup()

cat(sprintf(paste0("\nClasses of panel C, at posterior probability %.2f",
                   " (share: %% of the population at risk in 2025):\n"),
            confidence))
class_shares %>%
  mutate(text = sprintf("%s, %s: share %.1f", model, class, 100 * share)) %>%
  pull(text) %>%
  cat(sep = "\n")

# the gain in effective coverage from buying pyrethroid-only instead of dual-AI
# nets with the same budget, where pyrethroid-only nets are the better buy
# (their class in panel C): (P_dual + D) / (P_pyr + D) e - 1
# = e / e*(D) - 1, linear in e, so its mean over the population at risk is
# that at their mean effectiveness, here at each cell's posterior mean
# susceptibility
gain <- vapply(names(effectiveness_models), function(model) {
  better <- better_buy_class(model) == 1
  effectiveness <- weighted.mean(
    effectiveness_models[[model]](draws$cells$mean[better]),
    draws$cells$people_2025[better]
  )
  effectiveness / break_even - 1
}, numeric(1))
cat("\nGain in effective coverage from pyrethroid-only nets, where they are",
    "the better buy (% of that from dual-AI nets):\n")
cat(sprintf("%s: %.1f", names(gain), 100 * gain), sep = "\n")

# panel B: Declining effectiveness of nets, under an effectiveness model
effectiveness_panel <- function(model) {
  panel_b_data <- panel_data[[model]]$b

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
      limits = c(0.6, 1),
      name = "LLIN effectiveness"
    ) +
    # vertical dividing line for past and future, light grey and solid,
    # unlike the break-even line
    geom_vline(
      xintercept = latest_nets,
      colour = grey(0.85)
    ) +
    # past data
    geom_ribbon(
      data = panel_b_data_past,
      fill = "#56B1F7"
    ) +
    geom_line(
      data = panel_b_data_past,
      colour = grey(0.2),
      # linewidth = 0.8
    ) + 
    # future projections, at half the alpha; the line solid, as dotted lines
    # are counterfactuals
    geom_ribbon(
      data = panel_b_data_future,
      fill = "#56B1F7",
      alpha = 0.5
    ) +
    geom_line(
      data = panel_b_data_future,
      colour = grey(0.2),
      alpha = 0.5,
      # linewidth = 0.8
    ) + 
    # the break-even effectiveness e*(D) at the mean delivery cost: above it,
    # pyrethroid-only nets give more effective coverage than dual-AI nets for
    # the same budget
    geom_hline(
      yintercept = break_even,
      colour = grey(0.4),
      linetype = 2
    ) +
    # which nets are more cost-effective on either side of the line, in the
    # corners the effectiveness leaves clear under either model: top right
    # above it, bottom left below it; on three lines, so the top one stays
    # clear of Churcher's curve, which is still high in 2010
    annotate(
      "text",
      x = c(max(years_plot), min(years_plot)),
      y = c(1, 0.6),
      label = c("Pyrethroid-only\nmore\ncost-effective",
                "Dual-AI\nmore\ncost-effective"),
      hjust = c(1, 0),
      vjust = c(1, 0),
      size = 3,
      lineheight = 0.9,
      colour = grey(0.4)
    ) +
    xlab("") +
    theme_minimal()

  panel_b
}


# A) Effective coverage by net type: stacked areas of dual-AI, PBO, then
# pyrethroid-only nets, at their posterior means, and the line of the
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

# panel A, under an effectiveness model
coverage_panel <- function(model) {
  # the years of net use data only
  panel_a_data <- panel_data[[model]]$a %>%
    filter(
      year <= latest_nets
    )

  # the bands of the stack, from the bottom
  panel_a_stack <- panel_a_data %>%
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

  panel_a <- ggplot(
    mapping = aes(
      x = year,
      y = mean,
      ymin = lower,
      ymax = upper
    )
  ) +
    scale_x_continuous(
      breaks = c(2000, 2010, 2020)
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 0.55),
      name = "Effective coverage"
    ) +
    scale_fill_manual(
      values = net_type_cols,
      labels = net_type_labels,
      breaks = c("pyrethroid", "pbo", "dual_ai"),
      name = NULL
    ) +
    # the stack
    geom_ribbon(
      data = panel_a_stack,
      aes(
        y = NULL,
        ymin = bottom,
        ymax = top,
        fill = type
      )
    ) +
    # the counterfactual without resistance, dotted: the tops of the dual-AI
    # and PBO bands (where those nets are in use), which resistance leaves
    # unchanged, and the total for all nets
    geom_line(
      data = filter(panel_a_stack, type != "pyrethroid", top > 0),
      aes(
        y = top,
        ymin = NULL,
        ymax = NULL,
        group = type,
        colour = "No resistance"
      ),
      linetype = 3
    ) +
    geom_line(
      data = panel_a_data,
      aes(
        y = mean_no_ir,
        colour = "No resistance"
      ),
      linetype = 3
    ) +
    scale_colour_manual(
      values = c("No resistance" = "black"),
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

  panel_a
}

# C) Spatial variation in effectiveness: the classes of 2025 susceptibility by
# the better buy, on a map


# grey background for Africa
africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))

border_col <- grey(0.4)

# the colours of panel A's pyrethroid-only and dual-AI nets for the classes
# where they are the better buy, and for the uncertain class a colour halfway
# between them, lightened three quarters of the way to white. The break-even
# susceptibilities and the cut are in the log, for the caption
pyrethroid_col <- net_type_cols[["pyrethroid"]]
uncertain_col <- colorRampPalette(
  c(colorRampPalette(c(pyrethroid_col, dual_ai_col))(3)[2], "white")
)(5)[4]
band_cols <- c(pyrethroid = pyrethroid_col,
               uncertain = uncertain_col,
               dual_ai = dual_ai_col)
band_labels <- c(pyrethroid = "Pyrethroid-only",
                 uncertain = "Uncertain",
                 dual_ai = "Dual-AI")

better_buy_panel <- function(model) {
  ggplot() +
    africa_bg +
    geom_spatraster(
      data = better_buy_map(model)
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
      # note we arre using element_markdown() for legend titles so need this
      # html linebreak version
      name = "Most cost-effective<br>(2025)",
      na.value = "transparent"
    ) +
    # no padding around Africa, so it fills the panel
    coord_sf(expand = FALSE) +
    theme_ir_maps() +
    theme(
      plot.margin = unit(rep(0, 4), "cm"),
      # in the Atlantic, at the bottom left
      legend.position = "inside",
      legend.position.inside = c(0.05, 0.15),
      legend.justification.inside = c(0, 0),
      legend.key.height = rel(0.8),
      legend.key.width = rel(0.8),
      legend.text.position = "left",
      legend.text = element_text(size = rel(0.7),
                                 hjust = 1),
      legend.ticks = element_blank(),
      legend.title = element_markdown(size = rel(0.8),
                                      hjust = 1)
    )
}

# the figure under an effectiveness model: A and B side by side above, each
# 3.5 by 2.7, and C below, 7 by 7
epi_impact_figure <- function(model) {
  coverage_panel(model) +
    effectiveness_panel(model) +
    better_buy_panel(model) +
    plot_layout(
      design = "AB\nCC",
      heights = c(2.7, 7)
    ) +
    plot_annotation(
      tag_levels = "A"
    )
}

# drawn at three quarters of the size, so the text is larger
ggsave(
  filename = "figures/epi_impact.png",
  plot = epi_impact_figure("symons"),
  bg = "white",
  width = 7,
  height = 9.7,
  scale = 0.75,
  dpi = 400
)

# the supplement's version, under the effectiveness model of Churcher et al.
# 2024
ggsave(
  filename = "figures/epi_impact_churcher.png",
  plot = epi_impact_figure("churcher"),
  bg = "white",
  width = 7,
  height = 9.7,
  scale = 0.75,
  dpi = 400
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
