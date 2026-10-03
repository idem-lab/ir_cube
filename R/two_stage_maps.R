# Maps of the two-stage model (#21): the final correction model
# (R/two_stage_correction.R) fitted per insecticide type to ALL the bioassay
# data, on top of the full dynamical fit (R/fit_model.R), on the full
# prediction grid for the panel years of the dynamical-model maps
# (R/fig_ir_maps.R).
#
#   Rscript R/two_stage_maps.R prepare
#   Rscript R/two_stage_maps.R <type index>     (1-9, one process each)
#   Rscript R/two_stage_maps.R figures
#
# Run with OpenBLAS (reference BLAS is ~10x slower), e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=3 nice -n 10 Rscript ...
#
#   prepare  load the full dynamical fit (temporary/fitted_model.RData), check
#            that its data and covariates are the ones the current scripts
#            build, and recompute its 2000 paired posterior logit draws at
#            every assay and its initial conditions in every country
#            (R/dynamical_predictions.R), with the fit's model options. Save,
#            per type, the draws at its assays and the parameters the grid
#            recursion needs to outputs/two_stage/maps/<type>/dynamical.rds
#            (also read by R/fig_two_stage_components.R);
#   <k>      fit the final model for type k (fit_correction(), m_ref = the
#            posterior mean logit at the assays, per-type rho, t0 = 1995,
#            T = the type's last data year) and save it without its TMB object
#            to outputs/two_stage/maps/<type>/fit.rds. Then, on every mask cell
#            and map year, from n_map_draws paired draws (dynamical draw d with
#            a latent draw shifted for it, the cut posterior), write rasters of
#              two_stage_mortality  posterior mean of ilogit(m + omega + xi)
#              two_stage_sd_pp      its posterior SD, in percentage points:
#                                   the uncertainty of the dynamical model and
#                                   of the correction together
#              dynamical_mortality  posterior mean of ilogit(m)
#              difference_pp        two-stage minus dynamical, percentage points
#              correction_mean      posterior mean of omega + xi (logit)
#            The target is m + omega + xi. u and p are observation-level noise
#            and are not mapped. Beyond T, xi is the AR(1) forecast;
#   figures  maps in figures/two_stage/, in the layout of R/fig_ir_maps.R.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1)
step <- arguments[1]

source("R/two_stage_helpers.R")
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

output_dir <- "outputs/two_stage/maps"
figure_dir <- "figures/two_stage"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
type_dir <- function(type) file.path(output_dir, type)
raster_file <- function(type, quantity) {
  file.path(type_dir(type), sprintf("%s.tif", quantity))
}

# the panel years of R/fig_ir_maps.R; covariate year index 1 is baseline_year
map_years <- c(2000, 2005, 2010, 2015, 2020, 2025, 2030)
end_year <- max(map_years)
year_index <- map_years - baseline_year + 1

# 1000 paired draws (every other one of the 2000), in batches of 100: the sums
# behind the means and SD are accumulated batch by batch, so memory is set by
# chunk_size x map_batch_size x years, not by the number of draws. The Monte
# Carlo SE of a mean mortality, or of the difference map, is the posterior SD
# / sqrt(1000): ~0.3 percentage points where the SD is 10, ~1.3 where it is 40
# (far from data, late years). The draws are joint fields, so this error is
# spatially smooth
n_map_draws <- 1000
map_batch_size <- 100
chunk_size <- 25000
quantities <- c("two_stage_mortality", "two_stage_sd_pp",
                "dynamical_mortality", "difference_pp", "correction_mean")

# the grid cells and each one's country (NA outside the UNSD lookup)
grid_cells <- function() {
  cells <- terra::cells(mask)
  country <- as.character(terra::extract(rast("data/clean/country_raster.tif"),
                                         cells)$country_name)
  list(cells = cells, country = country)
}


# prepare: the dynamical draws ---------------------------------------------------

if (step == "prepare") {

  # save.image() of R/fit_model.R, November 2025. Its data must be exactly what
  # the current scripts build, or its draws cannot be paired with the assays
  fit_env <- new.env()
  load("temporary/fitted_model.RData", envir = fit_env)
  stopifnot(
    isTRUE(all.equal(fit_env$df, df)),
    identical(fit_env$types, types),
    identical(fit_env$classes, classes),
    identical(fit_env$countries, countries),
    identical(fit_env$regions, regions),
    identical(fit_env$unique_cells, unique_cells),
    identical(fit_env$classes_index, classes_index),
    isTRUE(all.equal(fit_env$x_cell_years, x_cell_years)),
    isTRUE(all.equal(fit_env$cell_years_index, cell_years_index))
  )
  fold <- list(draws = fit_env$draws, options = fit_env$model_options,
               x_cells_init = fit_env$x_cells_init)
  rm(fit_env)
  parameters <- dynamical_parameter_draws(fold, classes_index, types, df)
  design <- parameters$options$selection_columns
  rm(fold)
  invisible(gc())
  logit_assays <- dynamical_logit(parameters, df, df, x_cell_years,
                                  cell_years_index)
  set.seed(21)
  logit_init <- map_logit_init(parameters, countries, regions, df)
  # countries with data keep the fit's initial states
  stopifnot(isTRUE(all.equal(logit_init[, countries, , drop = FALSE],
                             parameters$logit_init_relative,
                             check.attributes = FALSE)))
  report("dynamical draws: %i assays x %i draws", ncol(logit_assays),
         nrow(logit_assays))

  grid <- grid_cells()

  # the grid's covariates must be the model's at the data cells, 1995-2024
  data_cells <- match(unique_cells, grid$cells)
  stopifnot(!anyNA(data_cells))
  n_fit_years <- max(cell_years_index$year_id)
  x_data <- map_x(map_covariates(grid$cells[data_cells], baseline_year,
                                 end_year, design),
                  seq_along(data_cells), n_fit_years)
  x_grid <- sapply(seq_len(ncol(x_cell_years)), function(j) {
    x_data[cbind(cell_years_index$cell_id, cell_years_index$year_id, j)]
  })
  stopifnot(max(abs(x_grid - x_cell_years)) < 1e-12)
  rm(x_data, x_grid)

  for (k in seq_along(types)) {
    dir.create(type_dir(types[k]), showWarnings = FALSE)
    saveRDS(list(logit_train = logit_assays[, df$type_id == k, drop = FALSE],
                 parameters = parameters,
                 logit_init = logit_init),
            file.path(type_dir(types[k]), "dynamical.rds"))
  }
  report("saved; peak memory %.1f GB", peak_memory_gb())
  quit(save = "no")
}


# one type: fit and map -------------------------------------------------------------

if (step != "figures") {

  source("R/two_stage_correction.R")
  k <- as.integer(step)
  stopifnot(!is.na(k), k >= 1, k <= length(types))
  type <- types[k]
  dynamical <- readRDS(file.path(type_dir(type), "dynamical.rds"))
  rows_k <- which(df$type_id == k)
  train_k <- tibble(lon = df$longitude[rows_k], lat = df$latitude[rows_k],
                    year = df$year_start[rows_k], cell = df$cell[rows_k],
                    died = df$died[rows_k],
                    mosquito_number = df$mosquito_number[rows_k],
                    m = colMeans(dynamical$logit_train),
                    rho = rho_for_record(tibble(insecticide_type = type),
                                         rho_lookup()))

  set.seed(string_seed(paste("maps", type, sep = "__")))
  meshes <- build_correction_meshes(coords_km(train_k),
                                    prediction_mask_coords())
  time_fit <- system.time(fit <- fit_correction(train_k, t0 = baseline_year,
                                                meshes = meshes))
  stopifnot(fit$opt$convergence == 0, isTRUE(fit$stage_b$converged_first),
            !isFALSE(fit$stage_b$converged_second))
  fit$obj <- NULL
  saveRDS(fit, file.path(type_dir(type), "fit.rds"))
  report("%s fitted in %.0f s and saved; peak memory %.1f GB", type,
         time_fit[["elapsed"]], peak_memory_gb())

  time_map <- system.time({
    grid <- grid_cells()
    covariates <- map_covariates(grid$cells, baseline_year, end_year,
                                 dynamical$parameters$options$selection_columns)
    xy <- terra::xyFromCell(mask, grid$cells)
    coords <- project_km(xy[, 1], xy[, 2])
    cell_country <- match(grid$country, dimnames(dynamical$logit_init)[[2]])
    rm(xy)

    # the map draws are an even subset of the 2000, by the scoring's rule.
    # The node fields of each batch of draws (omega and xi at the map years,
    # without the full latent vector) are small, so they are drawn up front
    map_draws <- thin_draws(matrix(seq_len(nrow(dynamical$logit_train))),
                            n_map_draws)[, 1]
    batches <- split(map_draws, ceiling(seq_along(map_draws) / map_batch_size))
    batch_parameters <- lapply(batches, subset_draws,
                               parameters = dynamical$parameters)
    logit_init <- dynamical$logit_init
    set.seed(string_seed(paste("maps draws", type, sep = "__")))
    fields <- lapply(batches, function(draws) {
      f <- correction_node_draws(
        fit, map_years, length(draws),
        m_draws_train = dynamical$logit_train[draws, , drop = FALSE])
      f$theta <- NULL
      f
    })
    fields_mean <- correction_node_draws(fit, map_years, mean = TRUE)
    rm(dynamical)
    invisible(gc())

    # project_correction() warns about points outside the mesh; count them
    n_outside <- 0
    # (counted once per chunk, from the draws' projection)
    project <- function(fields, new, count = TRUE) {
      withCallingHandlers(
        project_correction(fit, fields, new),
        warning = function(w) {
          if (grepl("outside the mesh", conditionMessage(w))) {
            if (count) n_outside <<- n_outside +
              as.numeric(sub(" .*", "", conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        })
    }

    results <- lapply(setNames(quantities, quantities), function(q) {
      matrix(NA_real_, length(grid$cells), length(map_years))
    })
    chunks <- split(seq_along(grid$cells),
                    ceiling(seq_along(grid$cells) / chunk_size))
    for (chunk in chunks) {
      ok <- chunk[!is.na(cell_country[chunk])]
      if (length(ok) == 0) next
      x_chunk <- map_x(covariates, ok, max(year_index))
      new <- tibble(x_km = rep(coords[ok, 1], length(map_years)),
                    y_km = rep(coords[ok, 2], length(map_years)),
                    year = rep(map_years, each = length(ok)))
      correction_mean <- project(fields_mean, new, count = FALSE)

      # sums over the draws, cells x years, accumulated batch by batch
      sum_two_stage <- matrix(0, length(ok), length(map_years))
      sum_sq_two_stage <- sum_two_stage
      sum_dynamical <- sum_two_stage
      for (b in seq_along(batches)) {
        draws <- batches[[b]]
        m <- dynamical_logit_cells(batch_parameters[[b]], k,
                                   matrix(logit_init[draws, cell_country[ok],
                                                     k], length(draws)),
                                   x_chunk, year_index,
                                   x_init = covariates$init[ok, , drop = FALSE])
        correction <- project(fields[[b]], new, count = b == 1)
        for (j in seq_along(map_years)) {
          rows <- (j - 1) * length(ok) + seq_along(ok)
          m_j <- t(m[[as.character(year_index[j])]])
          p_two_stage <- plogis(m_j + correction[rows, , drop = FALSE])
          sum_two_stage[, j] <- sum_two_stage[, j] + rowSums(p_two_stage)
          sum_sq_two_stage[, j] <- sum_sq_two_stage[, j] +
            rowSums(p_two_stage ^ 2)
          sum_dynamical[, j] <- sum_dynamical[, j] + rowSums(plogis(m_j))
        }
        rm(m, correction, p_two_stage)
      }

      mean_two_stage <- sum_two_stage / n_map_draws
      p_dynamical <- sum_dynamical / n_map_draws
      variance <- pmax(0, (sum_sq_two_stage - n_map_draws * mean_two_stage ^ 2) /
                         (n_map_draws - 1))
      results$two_stage_mortality[ok, ] <- mean_two_stage
      results$two_stage_sd_pp[ok, ] <- 100 * sqrt(variance)
      results$dynamical_mortality[ok, ] <- p_dynamical
      results$difference_pp[ok, ] <- 100 * (mean_two_stage - p_dynamical)
      results$correction_mean[ok, ] <- matrix(correction_mean[, 1], length(ok))
    }
  })
  n_outside_cells <- n_outside / length(map_years)
  if (n_outside_cells > 0) {
    warning(sprintf("%s: %.0f of %i mapped cells lie outside the mesh; ",
                    type, n_outside_cells,
                    sum(!is.na(cell_country))),
            "the correction there is 0 with SD 0")
  }

  for (q in quantities) {
    r <- rast(mask, nlyrs = length(map_years))
    full <- matrix(NA_real_, ncell(mask), length(map_years))
    full[grid$cells, ] <- results[[q]]
    values(r) <- full
    names(r) <- map_years
    writeRaster(r, raster_file(type, q), overwrite = TRUE, datatype = "FLT4S",
                gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  }
  write.csv(bind_cols(tibble(insecticide_type = type,
                             rho = train_k$rho[1]),
                      fit_summary(fit),
                      tibble(n_cells_outside_mesh = n_outside_cells,
                             time_fit_s = time_fit[["elapsed"]],
                             time_map_s = time_map[["elapsed"]],
                             peak_memory_gb = peak_memory_gb())),
            file.path(type_dir(type), "hyperparameters.csv"),
            row.names = FALSE)
  report("%s mapped in %.0f s (%.0f cells outside the mesh); peak memory %.1f GB",
         type, time_map[["elapsed"]], n_outside_cells, peak_memory_gb())
  quit(save = "no")
}


# figures --------------------------------------------------------------------------

hyperparameters <- bind_rows(lapply(types, function(type) {
  read.csv(file.path(type_dir(type), "hyperparameters.csv"))
}))
write.csv(hyperparameters, file.path(output_dir, "hyperparameters.csv"),
          row.names = FALSE)

# the look of R/fig_ir_maps.R: grey Africa background, thin grey borders, masked
# to the limits of Pf transmission and water bodies, one panel per year in two
# rows with the legend in the eighth slot
borders <- readRDS("data/clean/country_borders.RDS")
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")
africa_bg <- geom_sf(data = borders, linewidth = 0, fill = grey(0.75))
border_col <- grey(0.4)
country_borders <- geom_sf(data = borders, col = border_col, linewidth = 0.1,
                           fill = "transparent")
colourbar <- guide_colorbar(frame.colour = border_col, frame.linewidth = 0.1)

# the per-insecticide colours of R/fig_ir_maps.R
insecticides_plot <- c("Alpha-cypermethrin", "Deltamethrin",
                       "Lambda-cyhalothrin", "Permethrin", "Fenitrothion",
                       "Malathion", "Pirimiphos-methyl", "DDT", "Bendiocarb")
insecticides_col <- setNames(rev(scales::hue_pal()(length(insecticides_plot))),
                             insecticides_plot)

read_map <- function(type, quantity) {
  r <- rast(raster_file(type, quantity))
  names(r) <- map_years
  terra::mask(r, pf_water_mask)
}

# diverging scale centred at 0. Positive = more susceptible (higher mortality)
# under the two-stage model than under the dynamical model
diverging_scale <- function(name, limit) {
  scale_fill_gradient2(name = name, low = "#b2182b", mid = "white",
                       high = "#2166ac", midpoint = 0,
                       limits = c(-limit, limit), oob = scales::squish,
                       na.value = "transparent", guide = colourbar)
}

year_panels <- function(raster, fill_scale, title, subtitle, file) {
  years_list <- lapply(seq_len(nlyr(raster)), function(i) {
    ggplot() +
      africa_bg +
      geom_spatraster(data = raster[[i]]) +
      country_borders +
      fill_scale +
      facet_wrap(~lyr, nrow = 1, ncol = 1) +
      theme_ir_maps() +
      theme(plot.margin = unit(rep(0, 4), "cm"),
            legend.text.position = "left",
            legend.ticks = element_blank())
  })
  patchwork::wrap_plots(c(years_list, list(patchwork::guide_area()))) +
    patchwork::plot_layout(guides = "collect", nrow = 2) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste(strwrap(subtitle, width = 110), collapse = "\n"))
  ggsave(file, bg = "white", width = 13, height = 8, scale = 0.8, dpi = 300)
}

# common limits across types, so the maps can be compared between
# insecticides: a high quantile of the absolute value over all types, years
# and (masked) cells
pooled_quantile <- function(quantity, prob = 0.995) {
  values <- unlist(lapply(types, function(type) {
    v <- values(read_map(type, quantity), mat = FALSE)
    abs(v[!is.na(v)])
  }))
  quantile(values, prob, names = FALSE)
}
correction_limit <- ceiling(pooled_quantile("correction_mean") * 4) / 4
difference_limit <- ceiling(pooled_quantile("difference_pp") / 5) * 5
sd_limit <- ceiling(pooled_quantile("two_stage_sd_pp") / 5) * 5

model_note <- function(type) {
  T_k <- hyperparameters$T[hyperparameters$insecticide_type == type]
  sprintf(paste("target m + omega + xi (u and p are observation noise, not",
                "mapped); data to %i, later years are the AR(1) forecast",
                "of xi"), T_k)
}

for (type in insecticides_plot) {

  year_panels(
    read_map(type, "two_stage_mortality"),
    scale_fill_gradient(labels = scales::percent, name = "Susceptibility",
                        limits = c(0, 1), breaks = c(0, 0.5, 1),
                        high = insecticides_col[[type]], low = "white",
                        na.value = "transparent", guide = colourbar),
    title = sprintf("%s: two-stage model", type),
    subtitle = paste("Susceptibility of An. gambiae (s.l./s.s.) in WHO",
                     "bioassays; posterior mean of ilogit(m + omega + xi)",
                     sprintf("over %i paired draws;", n_map_draws),
                     model_note(type)),
    file = file.path(figure_dir, sprintf("%s_two_stage_ir_map.png", type)))

  year_panels(
    read_map(type, "two_stage_sd_pp"),
    # a coloured sequential palette, so that no SD reads as the grey of the
    # masked land behind it
    scale_fill_viridis_c(name = "SD<br>(% points)", limits = c(0, sd_limit),
                         option = "magma", direction = -1,
                         oob = scales::squish, na.value = "transparent",
                         guide = colourbar),
    title = sprintf("%s: two-stage model, posterior SD", type),
    subtitle = paste("Posterior SD of susceptibility, ilogit(m + omega + xi),",
                     "in percentage points, including the dynamical model's",
                     "uncertainty;", model_note(type)),
    file = file.path(figure_dir, sprintf("%s_two_stage_sd_map.png", type)))

  year_panels(
    read_map(type, "correction_mean"),
    diverging_scale("Correction<br>(logit)", correction_limit),
    title = sprintf("%s: second-stage correction", type),
    subtitle = paste("Posterior mean of omega + xi on the logit scale",
                     "(+ = more susceptible than the dynamical model);",
                     model_note(type)),
    file = file.path(figure_dir, sprintf("%s_correction_map.png", type)))

  year_panels(
    read_map(type, "difference_pp"),
    diverging_scale("Difference<br>(% points)", difference_limit),
    title = sprintf("%s: two-stage minus dynamical", type),
    subtitle = paste("Difference in posterior mean susceptibility",
                     "(including the shrinkage towards 50% from the",
                     "correction's variance), percentage points;",
                     model_note(type)),
    file = file.path(figure_dir, sprintf("%s_difference_map.png", type)))
}

# all types side by side, for one year, to compare the corrections' structure
compare_year <- 2020
all_types <- rast(lapply(insecticides_plot, function(type) {
  r <- read_map(type, "correction_mean")[[as.character(compare_year)]]
  names(r) <- type
  r
}))
ggplot() +
  africa_bg +
  geom_spatraster(data = all_types) +
  country_borders +
  facet_wrap(~lyr, ncol = 3) +
  diverging_scale("Correction<br>(logit)", correction_limit) +
  theme_ir_maps() +
  # a top margin, or the title is clipped
  theme(plot.margin = unit(c(0.3, 0, 0, 0), "cm"),
        legend.ticks = element_blank()) +
  labs(title = sprintf("Second-stage correction in %i", compare_year),
       subtitle = "Posterior mean of omega + xi, logit scale")
ggsave(file.path(figure_dir,
                 sprintf("correction_all_types_%i.png", compare_year)),
       bg = "white", width = 10, height = 10, scale = 0.8, dpi = 300)

# fitted hyperparameters per type
hyperparameters %>%
  select(insecticide_type,
         `omega range (km)` = range_omega, `omega SD` = sigma_omega,
         `eta range (km)` = range_eta, `eta SD` = sigma_eta, phi = phi,
         tau = tau, sigma_p = sigma_p) %>%
  pivot_longer(-insecticide_type) %>%
  mutate(name = factor(name, levels = unique(name)),
         insecticide_type = factor(insecticide_type,
                                   levels = rev(insecticides_plot))) %>%
  ggplot(aes(x = value, y = insecticide_type, colour = insecticide_type)) +
  geom_point(size = 2) +
  facet_wrap(~name, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = insecticides_col, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Second-stage hyperparameters, final model, fitted to all data") +
  theme_minimal()
ggsave(file.path(figure_dir, "map_hyperparameters.png"),
       bg = "white", width = 14, height = 3.5, dpi = 200)
report("figures written")
