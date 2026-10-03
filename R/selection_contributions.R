# The contribution of each selection term (nets, IRS, population, crops) to
# the cumulative selection in fits of the dynamical model, per insecticide
# type and pooled over the pyrethroids:
#
#   r2_main            round-2 main fit, chain 4 dropped
#   r2_dhalf5, r2_dhalf200  round 2 with d_half 5 and 200
#
#   Rscript R/selection_contributions.R
#
# Each fit is loaded once, in its own process, for draws of the selection
# effects (dynamical_parameter_draws(), with the fit's own options), cached in
# outputs/selection_contributions/. The covariates are rebuilt with
# map_covariates() and the fit's own design, and checked against the fit's
# x_cell_years at its bioassay pixels.
#
# The cumulative selection from 1995 to year Y on the logit scale is
# sum_t log w_t, log w_t = log(1 + sum_j x_jt e_j). It is not additive over
# the terms, so within each cell-year log w_t is split over them in proportion
# to x_jt e_j, which sums to log w_t exactly. Reversion adds t kappa (<= 0)
# and is reported separately. Summaries, per draw:
#   bioassay  mean over the pixels with bioassays of the type in the r2_main
#             data (the same pixels for every fit), weighted by the number of
#             assays of the type there; the pyrethroid pool weights the four
#             types' pixels by their assays
#   all       unweighted mean over 20,000 random mask pixels; the pyrethroid
#             pool is the mean over the four types
# Shares are ratios of these means: term / (nets + IRS + population + crops).
# Writes outputs/selection_contributions.csv and
# figures/selection_contributions.png.

suppressMessages({
  library(tidyverse)
  library(terra)
})
source("R/functions.R")
source("R/model_covariates.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

fits <- tribble(
  ~fit, ~file, ~label,
  "r2_main", "outputs/fits/r2_main.RData", "Round 2 (main)",
  "r2_dhalf5", "outputs/fits/r2_sens_dhalf5.RData", "Round 2, d_half 5",
  "r2_dhalf200", "outputs/fits/r2_sens_dhalf200.RData", "Round 2, d_half 200")
cache_dir <- "outputs/selection_contributions"
n_draws_keep <- 400
report_years <- c(2010, 2020, 2024)
pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Lambda-cyhalothrin",
                 "Permethrin")
terms <- c("Nets", "IRS", "Population", "Crops")


# draws of the selection effects from one fit, in its own process -------------

extract_fit <- function(fit_name, file) {
  e <- new.env()
  load(file, envir = e)
  draws <- e$draws
  # draws has already lost any chains drop_stuck_chains.R dropped
  fold <- list(draws = draws, options = e$model_options)
  index <- paired_draw_index(fold)
  index <- index[round(seq(1, length(index), length.out = n_draws_keep))]
  par <- dynamical_parameter_draws(fold, df = e$df,
                                   classes_index = e$classes_index,
                                   types = e$types, draw_index = index)
  design <- complete_selection_design(fold$options$selection_columns)
  columns <- selection_column_names(design)
  stopifnot(identical(columns, colnames(e$x_cell_years)))

  # the rebuilt covariates against the fit's own, at its bioassay pixels
  n_years <- nrow(e$x_cell_years) / length(e$unique_cells)
  rebuilt <- selection_design_matrix(e$unique_cells, e$baseline_year,
                                     e$baseline_year + n_years - 1, design)
  difference <- max(abs(rebuilt$x_cell_years - e$x_cell_years))
  cat(fit_name, ": chains", length(draws), ", draws", par$n_draws,
      ", max covariate difference", signif(difference, 3), "\n")

  out <- list(effect = par$effect_type, kappa = par$kappa_type,
              floor = par$mortality_floor, types = e$types,
              classes = e$classes[e$classes_index], columns = columns,
              design = design, baseline_year = e$baseline_year,
              unique_cells = e$unique_cells,
              df = select(e$df, cell_id, type_id),
              covariate_difference = difference)
  saveRDS(out, file.path(cache_dir, paste0(fit_name, ".RDS")))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3 && args[1] == "extract") {
  i <- match(args[2], fits$fit)
  extract_fit(fits$fit[i], fits$file[i])
  quit(save = "no")
}

dir.create(cache_dir, showWarnings = FALSE)
for (i in seq_len(nrow(fits))) {
  cache <- file.path(cache_dir, paste0(fits$fit[i], ".RDS"))
  if (!file.exists(cache)) {
    status <- system2("Rscript", c("R/selection_contributions.R", "extract",
                                   fits$fit[i], "x"))
    stopifnot(status == 0)
  }
}
fit_draws <- lapply(setNames(nm = fits$fit),
                    function(f) readRDS(file.path(cache_dir,
                                                  paste0(f, ".RDS"))))
stopifnot(all(vapply(fit_draws, function(f) {
  f$covariate_difference < 1e-6
}, logical(1))))


# the pixels -------------------------------------------------------------------

# bioassay pixels: those of the r2_main data, with their assay counts by type
main <- fit_draws$r2_main
types <- main$types
stopifnot(all(vapply(fit_draws, function(f) identical(f$types, types),
                     logical(1))))
assays <- main$df %>%
  mutate(cell = main$unique_cells[cell_id], type = types[type_id]) %>%
  count(type, cell)
bioassay_cells <- sort(unique(assays$cell))

mask <- rast("data/clean/raster_mask.tif")
mask_cells <- which(!is.na(values(mask, mat = FALSE)))
set.seed(1)
sample_cells <- sort(sample(mask_cells, 20000))

years <- 1995:max(report_years)
stopifnot(all(vapply(fit_draws, function(f) f$baseline_year == min(years),
                     logical(1))))


# cumulative selection by term -------------------------------------------------

# the term of each design column
term_of <- function(columns) {
  case_when(grepl("^nets", columns) ~ "Nets",
            grepl("^irs", columns) ~ "IRS",
            grepl("pop", columns) ~ "Population",
            TRUE ~ "Crops")
}

# For covariates `covariates` (map_covariates()) at cells with weights `w`
# (summing to 1) and effects `effect` (draws x covariates): the weighted mean
# over cells of the cumulative log w by term in each of report_years, as a
# tibble of year, term, draw and value
cumulative_terms <- function(covariates, w, effect, columns) {
  term <- term_of(columns)
  cumulative <- matrix(0, length(terms), nrow(effect),
                       dimnames = list(terms, NULL))
  out <- list()
  for (t in seq_along(years)) {
    x_t <- cbind(covariates$time_varying[, t, ], covariates$flat)
    total <- x_t %*% t(effect)
    share <- ifelse(total > 0, log1p(total) / total, 0)
    for (name in terms) {
      j <- term == name
      if (!any(j)) next
      part <- x_t[, j, drop = FALSE] %*% t(effect[, j, drop = FALSE])
      cumulative[name, ] <- cumulative[name, ] + colSums(w * part * share)
    }
    if (years[t] %in% report_years) {
      out[[length(out) + 1]] <- tibble(
        year = years[t], term = rep(terms, nrow(effect)),
        draw = rep(seq_len(nrow(effect)), each = length(terms)),
        value = c(cumulative))
    }
  }
  bind_rows(out)
}

results <- list()
for (f in fits$fit) {
  fd <- fit_draws[[f]]
  cov_bioassay <- map_covariates(bioassay_cells, min(years), max(years),
                                 fd$design)
  cov_all <- map_covariates(sample_cells, min(years), max(years), fd$design)
  stopifnot(!anyNA(cov_bioassay$time_varying))
  # sample pixels with any covariate missing are left out
  ok <- !apply(is.na(cov_all$time_varying), 1, any) &
    !apply(is.na(cov_all$flat), 1, any)
  for (k in seq_along(types)) {
    effect <- matrix(fd$effect[, , k], nrow = dim(fd$effect)[1])
    a <- filter(assays, type == types[k])
    rows <- match(a$cell, bioassay_cells)
    subset_cells <- function(cv, r) {
      list(time_varying = cv$time_varying[r, , , drop = FALSE],
           flat = cv$flat[r, , drop = FALSE])
    }
    kappa <- if (!is.null(fd$kappa)) fd$kappa[, k] else rep(0, nrow(effect))
    for (measure in c("bioassay", "all")) {
      if (measure == "bioassay") {
        cv <- subset_cells(cov_bioassay, rows)
        w <- a$n / sum(a$n)
      } else {
        cv <- subset_cells(cov_all, which(ok))
        w <- rep(1 / sum(ok), sum(ok))
      }
      ct <- cumulative_terms(cv, w, effect, fd$columns)
      reversion <- tibble(year = rep(report_years, each = nrow(effect)),
                          term = "Reversion",
                          draw = rep(seq_len(nrow(effect)),
                                     length(report_years)),
                          value = (year - min(years) + 1) *
                            kappa[draw])
      results[[length(results) + 1]] <- bind_rows(ct, reversion) %>%
        mutate(fit = f, type = types[k], measure = measure,
               weight = sum(a$n))
    }
  }
  rm(cov_bioassay, cov_all)
  invisible(gc())
}
results <- bind_rows(results)

# the pyrethroid pool, per draw: weighted by assays at the bioassay pixels,
# the mean of the types over all pixels
pooled <- results %>%
  filter(type %in% pyrethroids) %>%
  group_by(fit, measure, year, term, draw) %>%
  summarise(value = if (measure[1] == "bioassay") {
    weighted.mean(value, weight)
  } else {
    mean(value)
  }, .groups = "drop") %>%
  mutate(type = "Pyrethroids")
results <- bind_rows(select(results, -weight), pooled)


# summaries --------------------------------------------------------------------

per_draw <- results %>%
  pivot_wider(names_from = term, values_from = value) %>%
  mutate(total = Nets + IRS + Population + Crops,
         across(all_of(terms), ~ .x / total, .names = "share_{.col}"))

summary <- per_draw %>%
  select(fit, type, measure, year, draw, total, Reversion,
         starts_with("share_")) %>%
  pivot_longer(c(total, Reversion, starts_with("share_")),
               names_to = "quantity") %>%
  group_by(fit, type, measure, year, quantity) %>%
  summarise(median = median(value), lower = quantile(value, 0.025),
            upper = quantile(value, 0.975), .groups = "drop") %>%
  mutate(quantity = recode(quantity, total = "total_log_selection",
                           Reversion = "reversion"),
         quantity = str_to_lower(quantity),
         fit = factor(fit, levels = fits$fit),
         type = factor(type, levels = c("Pyrethroids",
                                        insecticides_plot_order)),
         measure = factor(measure, levels = c("bioassay", "all"))) %>%
  arrange(fit, type, measure, year, quantity)
write_csv(summary, "outputs/selection_contributions.csv")

table_2020 <- summary %>%
  filter(year == 2020, type %in% c("Pyrethroids", "DDT", "Bendiocarb",
                                   "Pirimiphos-methyl")) %>%
  mutate(value = if_else(grepl("^share", quantity),
                         sprintf("%.0f (%.0f-%.0f)", 100 * median,
                                 100 * lower, 100 * upper),
                         sprintf("%.2f (%.2f-%.2f)", median, lower,
                                 upper))) %>%
  select(fit, type, measure, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value) %>%
  select(fit, type, measure, share_nets, share_population, share_crops,
         share_irs, total_log_selection, reversion) %>%
  arrange(type, measure, fit)
print(table_2020, n = Inf, width = Inf)


# figure -----------------------------------------------------------------------

term_colours <- c(Nets = "#1b9e77", IRS = "#d95f02", Population = "#7570b3",
                  Crops = "#e6ab02")
plot_data <- summary %>%
  filter(year == 2020, type %in% c("Pyrethroids", "DDT"),
         grepl("^share", quantity)) %>%
  mutate(term = factor(c(share_nets = "Nets", share_irs = "IRS",
                         share_population = "Population",
                         share_crops = "Crops")[quantity], levels = terms),
         fit = factor(fits$label[match(fit, fits$fit)], levels = fits$label),
         measure = factor(c(bioassay = "Bioassay pixels (assay-weighted)",
                            all = "All pixels (unweighted)")[
                              as.character(measure)],
                          levels = c("Bioassay pixels (assay-weighted)",
                                     "All pixels (unweighted)")))
p <- ggplot(plot_data, aes(term, median, fill = fit)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.8), width = 0.25,
                linewidth = 0.3) +
  facet_grid(type ~ measure) +
  scale_fill_manual(values = c("#08306b", "#e6550d", "#fdae6b"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share of cumulative selection, 1995-2020",
       caption = paste(
         "Bars: posterior median; lines: 95% interval (400 draws).",
         "log w = log(1 + sum_j x_j e_j) split over the terms in proportion",
         "to x_j e_j in each cell-year,\nsummed over 1995-2020 and averaged",
         "over the pixels; share = term / sum of the four terms (reversion",
         "excluded). Pyrethroids: the four types pooled.")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom")
ggsave("figures/selection_contributions.png", p, width = 9, height = 6,
       dpi = 150, bg = "white")
cat("done\n")
