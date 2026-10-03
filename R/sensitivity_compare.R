# Compare sensitivity fits with the main fit: predicted mortality at the
# bioassay pixels and over the whole map, and the key parameters.
#
#   Rscript R/sensitivity_compare.R predict <fitted_model.RData> <out.rds>
#     [n_px] [n_all]
#   Rscript R/sensitivity_compare.R compare <output dir> <main.rds>
#     <name>=<out.rds> [<name>=<out.rds> ...]
#
# predict: posterior mean predicted mortality of one fit, by predict.R's path
# (map_covariates(), map_logit_init(), dynamical_logit_cells(), with the fit's
# model options), for every type,
#   - at the bioassay pixels (the cells of the fit's data), every year from the
#     baseline year to 2025, over n_px draws (default 500, the draws predict.R
#     maps);
#   - at every 4th mask cell (SENS_MAP_THIN in the environment) with a
#     country, in 2025, mean and SD over n_all draws (default 100);
# and the posterior mean and SD of the key parameters: the mortality floor,
# reversion per class, and the initial-state and selection coefficients
# (init_coef, beta_overall) averaged over types, per covariate. Draws are
# paired_draw_index()'s, so chains dropped by R/drop_stuck_chains.R stay out.
# map_covariates() needs more than 3 GB at any number of cells. For testing,
# SENS_COVARIATES=fit uses the fit's own covariates instead (x_cell_years,
# carrying the last data year forward), at the fit's cells only, which needs
# under 2 GB.
#
# compare: for each sensitivity fit against the main fit, per type,
#   - mean |difference| in 2025 mortality at the bioassay pixels, and in the
#     change from 2014 to 2024;
#   - over the map cells, mean and 95th percentile of |difference| in 2025, and
#     the Monte Carlo noise in that mean, from the two fits' posterior SDs
#     (E|N(0, s^2/n_a + s^2/n_b)|);
# and flags the cross-validation trigger: either bioassay-pixel difference
# above 0.03 for any type. Also the in-sample fit of each (correlation and
# RMSE of predicted against observed mortality at the modelled assays, per
# type). Writes sensitivity_*.csv to the output dir.

arguments <- commandArgs(trailingOnly = TRUE)
mode <- arguments[1]
trigger <- 0.03

if (identical(mode, "predict")) {

  fit_file <- arguments[2]
  out_file <- arguments[3]
  n_px <- if (length(arguments) > 3) as.integer(arguments[4]) else 500L
  n_all <- if (length(arguments) > 4) as.integer(arguments[5]) else 100L
  map_thin <- as.integer(Sys.getenv("SENS_MAP_THIN", "4"))
  fit_covariates <- Sys.getenv("SENS_COVARIATES") == "fit"

  source("R/packages.R")
  source("R/functions.R")
  source("R/dynamical_predictions.R")
  source("R/two_stage_map_functions.R")
  # report()
  source("R/two_stage_helpers.R")
  end_year <- 2025

  fit_env <- new.env()
  load(fit_file, envir = fit_env)
  baseline_year <- fit_env$baseline_year
  types <- fit_env$types
  classes <- fit_env$classes
  classes_index <- fit_env$classes_index
  countries <- fit_env$countries
  regions <- fit_env$regions
  df <- fit_env$df
  covariate_names <- colnames(fit_env$x_cell_years)
  init_names <- colnames(fit_env$x_cells_init)
  if (fit_covariates) {
    x_cell_years <- fit_env$x_cell_years
    x_cells_init <- fit_env$x_cells_init
    unique_cells <- fit_env$unique_cells
    final_data_year <- fit_env$final_data_year
  }
  fold <- list(draws = fit_env$draws, options = fit_env$model_options,
               x_cells_init = fit_env$x_cells_init)
  options <- fold$options
  rm(fit_env)
  invisible(gc())

  draw_index <- paired_draw_index(fold)
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]

  # key parameters: posterior mean and SD, averaging over the types
  parameter_rows <- list()
  add_parameter <- function(label, values) {
    parameter_rows[[length(parameter_rows) + 1]] <<- data.frame(
      parameter = label, mean = mean(values), sd = sd(values))
  }
  columns <- colnames(draws_matrix)
  pick <- function(pattern) columns[grepl(pattern, columns)]
  if (length(pick("^mortality_floor"))) {
    add_parameter("mortality_floor", draws_matrix[, pick("^mortality_floor")[1]])
  }
  reversion <- pick("^reversion_rate\\[")
  for (j in seq_along(reversion)) {
    add_parameter(sprintf("reversion_rate[%s]",
                          if (length(classes) == length(reversion))
                            classes[j] else j),
                  draws_matrix[, reversion[j]])
  }
  for (j in seq_along(init_names)) {
    cols <- pick(sprintf("^init_coef\\[%i,", j))
    if (length(cols)) {
      add_parameter(sprintf("init_coef[%s]", init_names[j]),
                    rowMeans(draws_matrix[, cols, drop = FALSE]))
    }
  }
  for (j in seq_along(covariate_names)) {
    cols <- pick(sprintf("^beta_overall\\[%i,", j))
    if (length(cols)) {
      add_parameter(sprintf("beta_overall[%s]", covariate_names[j]),
                    rowMeans(draws_matrix[, cols, drop = FALSE]))
    }
  }
  parameters_table <- do.call(rbind, parameter_rows)
  print(parameters_table, digits = 3)

  parameters <- dynamical_parameter_draws(fold, classes_index, types, df,
                                          draw_index = draw_index,
                                          options = options)
  set.seed(1)
  logit_init_all <- map_logit_init(parameters, countries, regions, df)
  select_draws <- function(n) {
    d <- round(seq(1, length(draw_index), length.out = n))
    list(parameters = subset_draws(parameters, d),
         logit_init = logit_init_all[d, , , drop = FALSE])
  }
  draws_px <- select_draws(n_px)
  draws_all <- select_draws(n_all)
  chains_kept <- length(fold$draws) - length(fold$stuck_chains)
  rm(parameters, logit_init_all, draws_matrix, fold)
  invisible(gc())
  report("parameters done, %i draws at the pixels, %i over the map",
         n_px, n_all)

  years_predict <- baseline_year:end_year
  if (fit_covariates) {
    cells <- unique_cells
    n_years_fit <- final_data_year - baseline_year + 1
    stopifnot(nrow(x_cell_years) == length(cells) * n_years_fit)
    time_varying <- aperm(array(x_cell_years, c(n_years_fit, length(cells),
                                                ncol(x_cell_years))),
                          c(2, 1, 3))
    carried <- c(seq_len(n_years_fit),
                 rep(n_years_fit, max(0, end_year - final_data_year)))
    covariates <- list(time_varying = time_varying[, carried, , drop = FALSE],
                       flat = matrix(0, length(cells), 0),
                       init = x_cells_init)
    rm(time_varying, x_cell_years)
  } else {
    mask <- rast("data/clean/raster_mask.tif")
    cells <- terra::cells(mask)
    cells <- sort(union(cells[seq(1, length(cells), by = map_thin)],
                        unique(df$cell)))
    covariates <- map_covariates(cells, baseline_year, end_year,
                                 options$selection_columns)
  }
  invisible(gc())
  report("covariates at %i cells", length(cells))
  country_raster <- rast("data/clean/country_raster.tif")
  cell_country <- as.character(terra::extract(country_raster,
                                              cells)$country_name)
  country_index <- match(cell_country, dimnames(draws_px$logit_init)[[2]])
  stopifnot(identical(dimnames(draws_px$logit_init)[[2]],
                      dimnames(draws_all$logit_init)[[2]]))

  # posterior mean (and SD) of each type at rows `rows` of the covariates, in
  # years `keep`
  predict_rows <- function(rows, d, keep, sd = FALSE) {
    mean_out <- array(NA_real_, c(length(rows), length(types), length(keep)),
                      dimnames = list(NULL, types, keep))
    sd_out <- if (sd) mean_out
    for (chunk in split(seq_along(rows), ceiling(seq_along(rows) / 1000))) {
      ok <- chunk[!is.na(country_index[rows[chunk]])]
      if (length(ok) == 0) next
      x <- map_x(covariates, rows[ok], max(keep) - baseline_year + 1)
      for (k in seq_along(types)) {
        dyn <- dynamical_logit_cells(
          d$parameters, k,
          matrix(d$logit_init[, country_index[rows[ok]], k],
                 d$parameters$n_draws),
          x, keep - baseline_year + 1,
          x_init = covariates$init[rows[ok], , drop = FALSE])
        for (j in seq_along(keep)) {
          p <- plogis(dyn[[j]])
          mean_out[ok, k, j] <- colMeans(p)
          if (sd) sd_out[ok, k, j] <- col_sds(p)
        }
      }
    }
    list(mean = mean_out, sd = sd_out)
  }

  time_start <- Sys.time()
  px_cells <- sort(unique(df$cell))
  px_rows <- match(px_cells, cells)
  stopifnot(!anyNA(px_rows))
  px <- predict_rows(px_rows, draws_px, years_predict)$mean
  report("%i bioassay pixels done in %.1f min", length(px_cells),
         as.numeric(difftime(Sys.time(), time_start, units = "mins")))

  all_mean <- all_sd <- matrix(NA_real_, length(cells), length(types),
                               dimnames = list(NULL, types))
  blocks <- split(seq_along(cells), ceiling(seq_along(cells) / 100000))
  for (block in blocks) {
    result <- predict_rows(block, draws_all, end_year, sd = TRUE)
    all_mean[block, ] <- result$mean[, , 1]
    all_sd[block, ] <- result$sd[, , 1]
    report("map cells %i of %i", max(block), length(cells))
  }

  observed <- df[, c("cell", "insecticide_type", "year_start", "died",
                     "mosquito_number")]
  saveRDS(list(fit = fit_file, chains_kept = chains_kept,
               parameters = parameters_table, px_cells = px_cells, px = px,
               cells = cells, all_mean = all_mean, all_sd = all_sd,
               n_px = n_px, n_all = n_all, observed = observed),
          out_file)
  report("done in %.1f min",
         as.numeric(difftime(Sys.time(), time_start, units = "mins")))

} else if (identical(mode, "compare")) {

  output_dir <- arguments[2]
  main_file <- arguments[3]
  specs <- arguments[-(1:3)]
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fits <- c(main = main_file,
            setNames(sub("^[^=]*=", "", specs), sub("=.*$", "", specs)))
  fits <- fits[file.exists(fits)]
  main <- readRDS(fits[["main"]])
  types <- dimnames(main$px)[[2]]
  year <- function(x, y) x$px[, , as.character(y)]
  change <- function(x) year(x, 2024) - year(x, 2014)

  # in-sample fit at the modelled assays, per fit and type
  observed <- main$observed
  observed$obs <- observed$died / observed$mosquito_number
  in_sample <- list()
  for (name in names(fits)) {
    x <- if (name == "main") main else readRDS(fits[[name]])
    index <- cbind(match(observed$cell, x$px_cells),
                   match(observed$insecticide_type, types),
                   match(as.character(observed$year_start),
                         dimnames(x$px)[[3]]))
    predicted <- x$px[index]
    for (type in types) {
      i <- observed$insecticide_type == type & !is.na(predicted)
      in_sample[[length(in_sample) + 1]] <- data.frame(
        fit = name, type = type, n = sum(i),
        r = cor(predicted[i], observed$obs[i]),
        r_weighted = cov.wt(cbind(predicted[i], observed$obs[i]),
                            wt = observed$mosquito_number[i],
                            cor = TRUE)$cor[1, 2],
        rmse = sqrt(mean((predicted[i] - observed$obs[i]) ^ 2)),
        mean_predicted = mean(predicted[i]),
        mean_observed = mean(observed$obs[i]))
    }
    rm(x)
  }
  in_sample <- do.call(rbind, in_sample)
  write.csv(in_sample, file.path(output_dir, "sensitivity_in_sample.csv"),
            row.names = FALSE)

  by_type <- list()
  parameters <- list(cbind(fit = "main", main$parameters))
  for (name in setdiff(names(fits), "main")) {
    x <- readRDS(fits[[name]])
    stopifnot(identical(x$px_cells, main$px_cells),
              identical(dimnames(x$px)[[2]], types),
              identical(x$cells, main$cells))
    parameters[[name]] <- cbind(fit = name, x$parameters)
    difference_2025 <- year(x, 2025) - year(main, 2025)
    difference_change <- change(x) - change(main)
    difference_all <- x$all_mean - main$all_mean
    noise <- sqrt(2 / pi) * sqrt(x$all_sd ^ 2 / x$n_all +
                                   main$all_sd ^ 2 / main$n_all)
    for (k in seq_along(types)) {
      by_type[[length(by_type) + 1]] <- data.frame(
        fit = name, type = types[k],
        px_abs_2025 = mean(abs(difference_2025[, k]), na.rm = TRUE),
        px_mean_2025 = mean(difference_2025[, k], na.rm = TRUE),
        px_abs_change = mean(abs(difference_change[, k]), na.rm = TRUE),
        px_mean_change = mean(difference_change[, k], na.rm = TRUE),
        map_abs_2025 = mean(abs(difference_all[, k]), na.rm = TRUE),
        map_mean_2025 = mean(difference_all[, k], na.rm = TRUE),
        map_p95_abs_2025 = unname(quantile(abs(difference_all[, k]), 0.95,
                                           na.rm = TRUE)),
        map_share_above_0.05 = mean(abs(difference_all[, k]) > 0.05,
                                    na.rm = TRUE),
        map_noise = mean(noise[, k], na.rm = TRUE),
        map_correlation = cor(x$all_mean[, k], main$all_mean[, k],
                              use = "complete.obs"))
    }
    rm(x)
    invisible(gc())
  }
  by_type <- do.call(rbind, by_type)
  write.csv(by_type, file.path(output_dir, "sensitivity_by_type.csv"),
            row.names = FALSE)

  summary_table <- do.call(rbind, lapply(split(by_type, by_type$fit),
                                         function(x) {
    data.frame(
      fit = x$fit[1],
      px_abs_2025_max = max(x$px_abs_2025),
      px_abs_2025_type = x$type[which.max(x$px_abs_2025)],
      px_abs_2025_mean = mean(x$px_abs_2025),
      px_abs_change_max = max(x$px_abs_change),
      px_abs_change_type = x$type[which.max(x$px_abs_change)],
      px_abs_change_mean = mean(x$px_abs_change),
      map_abs_2025_mean = mean(x$map_abs_2025),
      map_abs_2025_max = max(x$map_abs_2025),
      map_p95_abs_2025_max = max(x$map_p95_abs_2025),
      map_noise_mean = mean(x$map_noise),
      map_correlation_min = min(x$map_correlation),
      cv_trigger = max(x$px_abs_2025) > trigger ||
        max(x$px_abs_change) > trigger)
  }))
  summary_table <- summary_table[match(unique(by_type$fit),
                                       summary_table$fit), ]
  write.csv(summary_table, file.path(output_dir, "sensitivity_summary.csv"),
            row.names = FALSE)

  parameters <- do.call(rbind, parameters)
  write.csv(parameters, file.path(output_dir, "sensitivity_parameters.csv"),
            row.names = FALSE)

  options(width = 200)
  print(summary_table, digits = 3)
  print(reshape(parameters[, c("fit", "parameter", "mean")],
                idvar = "parameter", timevar = "fit", direction = "wide"),
        digits = 3)

} else {
  stop("usage: sensitivity_compare.R predict <fit> <out.rds> [n_px] [n_all] ",
       "| compare <output dir> <main.rds> <name>=<out.rds> ...")
}
