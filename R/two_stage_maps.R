# Maps of the two-stage model (#21), the model behind the published maps and
# time series: the final correction model (R/two_stage_correction.R) fitted
# per insecticide type to ALL the bioassay data, on top of the full dynamical
# fit (R/fit_model.R), on the full prediction grid for every year from the
# baseline to 2030.
#
#   Rscript R/two_stage_maps.R prepare
#   Rscript R/two_stage_maps.R fit <type index>   (1-9, one process each)
#   Rscript R/two_stage_maps.R map <type index>   (1-9, one process each)
#   Rscript R/two_stage_maps.R map llin_effective (after the pyrethroids)
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
#            per type, the draws at its assays, the parameters the grid
#            recursion needs and the assays' key columns (to check that later
#            steps pair the draws with the same data) to
#            outputs/two_stage/maps/<type>/dynamical.rds;
#   fit      fit the final model for type k (fit_correction(), m_ref = the
#            posterior mean logit at the assays, per-type rho, t0 = 1995,
#            T = the type's last data year), check that its meshes cover every
#            cell of the prediction mask, and save it without its TMB object
#            to outputs/two_stage/maps/<type>/fit.rds, and its hyperparameters
#            to hyperparameters.csv there. The fit and dynamical.rds are what
#            the predictions of R/two_stage_predictions.R read, in the figure
#            scripts too;
#   map      on every mask cell and year, from n_map_draws paired draws
#            (dynamical draw d with a latent draw shifted for it, the cut
#            posterior), the posterior mean and SD of the predicted fraction
#            susceptible (bioassay mortality), ilogit(m + omega + xi), written
#            as R/predict.R writes the dynamical model's, one file per year:
#              outputs/two_stage/ir_maps/<type>/ir_<year>_susceptibility.tif
#              outputs/two_stage/ir_maps/<type>/ir_<year>_susceptibility_sd.tif
#            and, at the panel years of the figures, to outputs/two_stage/maps/
#            <type>/:
#              dynamical_mortality  posterior mean of ilogit(m), from the same
#                                   draws
#              correction_mean      posterior mean of omega + xi (logit)
#              local_mean, local_sd posterior mean and SD of ilogit(m + omega
#                                   + xi + u + p), the same draws plus fresh
#                                   u and p (local_terms), for n*_local only
#            Beyond T, xi is the AR(1) forecast. `map llin_effective` writes
#            llin_effective's mean and SD as the sums of the pyrethroids',
#            weighted by temporary/ingredient_weights.RDS. Summing the SDs
#            takes the pyrethroids' errors as perfectly correlated, the
#            largest SD their marginal SDs allow (n_star(), #40);
#   figures  the two-stage-specific maps in figures/two_stage/ (n*_pop and
#            n*_local (#40), correction, difference from the dynamical
#            model), in the layout of R/fig_ir_maps.R, which maps the
#            two-stage model's mortality.

arguments <- commandArgs(trailingOnly = TRUE)
step <- arguments[1]
stopifnot(length(arguments) == if (step %in% c("fit", "map")) 2 else 1)

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
ir_map_file <- function(output, year, quantity = "susceptibility") {
  ir_map_files(output, year, quantity)
}

# every year, as R/predict.R; the panel years of R/fig_ir_maps.R
end_year <- 2030
years_all <- baseline_year:end_year
map_years <- c(2000, 2005, 2010, 2015, 2020, 2025, 2030)

# 1000 paired draws (every other one of the 2000), in batches of 100: the sums
# behind the means and SD are accumulated batch by batch, so memory is set by
# chunk_size x map_batch_size x years, not by the number of draws. The Monte
# Carlo SE of a mean mortality, or of the difference map, is the posterior SD
# / sqrt(1000): ~0.3 percentage points where the SD is 10, ~1.3 where it is 40
# (far from data, late years). The draws are joint fields, so this error is
# spatially smooth
n_map_draws <- 1000
map_batch_size <- 100
chunk_size <- 10000

# forked workers, each mapping a group of the cells, or the environment
# variable IR_CUBE_MAP_WORKERS. Each held ~5 GB at 4 workers (its cells'
# covariates, and a batch of draws for a chunk)
n_map_workers <- as.integer(Sys.getenv("IR_CUBE_MAP_WORKERS", "1"))

# the types weighted into llin_effective
llin_weights <- unlist(readRDS("temporary/ingredient_weights.RDS"))

# The iid terms drawn fresh for n*_local (#40): u (pixel-year) and p (pixel)
# are either data error (sentinel sites, repeated habitats) or real local
# variation the model cannot predict, and the data cannot tell which. n*_pop
# leaves them out (the upper bound on the map's worth), n*_local draws them
# (the lower bound). Any survey-level term is data error either way, and is
# never drawn
local_terms <- c("u", "p")

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
    identical(fit_env$classes_index, classes_index)
  )
  # the fit's own design matrix, not one built from the default options
  # (R/validation_covariates.R); the check of the grid's covariates below
  # confirms that the current rasters reproduce it
  x_cell_years <- fit_env$x_cell_years
  cell_years_index <- fit_env$cell_years_index
  fold <- list(draws = fit_env$draws, options = fit_env$model_options,
               x_cells_init = fit_env$x_cells_init)
  rm(fit_env)
  parameters <- dynamical_parameter_draws(fold, classes_index, types)
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
                 logit_init = logit_init,
                 keys = df[df$type_id == k, c("cell", "year", "died",
                                              "mosquito_number")]),
            file.path(type_dir(types[k]), "dynamical.rds"))
  }
  report("saved; peak memory %.1f GB", peak_memory_gb())
  quit(save = "no")
}


# fit: one type -------------------------------------------------------------------

if (step == "fit") {

  source("R/two_stage_correction.R")
  k <- as.integer(arguments[2])
  stopifnot(!is.na(k), k >= 1, k <= length(types))
  type <- types[k]
  source("R/two_stage_predictions.R")
  dynamical <- readRDS(file.path(type_dir(type), "dynamical.rds"))
  check_two_stage_data(dynamical, df, type)
  rows_k <- which(df$type_id == k)
  train_k <- tibble(lon = df$longitude[rows_k], lat = df$latitude[rows_k],
                    year = df$year[rows_k], cell = df$cell[rows_k],
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

  # every mapped cell must be inside both meshes, or its correction would be 0
  grid <- grid_cells()
  xy <- terra::xyFromCell(mask, grid$cells[!is.na(grid$country)])
  coords <- project_km(xy[, 1], xy[, 2])
  for (mesh in list(fit$mesh, fit$mesh_xi)) {
    stopifnot(all(Matrix::rowSums(mesh_basis(mesh, coords)) > 0.5))
  }
  fit$obj <- NULL
  saveRDS(fit, file.path(type_dir(type), "fit.rds"))
  write.csv(bind_cols(tibble(insecticide_type = type,
                             rho = train_k$rho[1]),
                      fit_summary(fit),
                      tibble(time_fit_s = time_fit[["elapsed"]],
                             peak_memory_fit_gb = peak_memory_gb())),
            file.path(type_dir(type), "hyperparameters.csv"),
            row.names = FALSE)
  report("%s fitted in %.0f s and saved; peak memory %.1f GB", type,
         time_fit[["elapsed"]], peak_memory_gb())
  quit(save = "no")
}

# map: one type, or llin_effective from the pyrethroids' maps ---------------------

# one or more layers over the mapped cells, compressed as R/predict.R
write_map <- function(r, file) {
  writeRaster(r, file, overwrite = TRUE, datatype = "FLT4S",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
}

# llin_effective: the pyrethroids' rasters weighted and summed (see the header)
if (step == "map" && arguments[2] == "llin_effective") {
  weighted_sum <- function(file) {
    Reduce(`+`, Map(function(type, w) w * rast(file(type)),
                    names(llin_weights), llin_weights))
  }
  dir.create(dirname(ir_map_file("llin_effective", years_all[1])),
             showWarnings = FALSE, recursive = TRUE)
  for (year in years_all) {
    for (quantity in c("susceptibility", "susceptibility_sd")) {
      write_map(weighted_sum(function(type) {
                  ir_map_file(type, year, quantity)
                }),
                ir_map_file("llin_effective", year, quantity))
    }
  }
  report("llin_effective written")
  quit(save = "no")
}

if (step == "map") {

  source("R/two_stage_correction.R")
  source("R/two_stage_predictions.R")
  k <- as.integer(arguments[2])
  stopifnot(!is.na(k), k >= 1, k <= length(types))
  type <- types[k]

  time_map <- system.time({
    # the node fields of each batch of draws (omega and xi at every year,
    # without the full latent vector) are small, so they are drawn up front
    setup <- two_stage_setup(type, years = years_all, data = df,
                             n_draws = n_map_draws,
                             batch_size = map_batch_size,
                             baseline_year = baseline_year)
    # the posterior mean of omega + xi at the panel years
    fields_mean <- correction_node_draws(setup$fit, map_years, mean = TRUE)

    grid <- grid_cells()
    mapped <- which(!is.na(grid$country))
    n_mapped <- length(mapped)
    panel <- match(map_years, years_all)

    # the posterior mean and SD of mortality from sums over the draws (pmax()
    # keeps the dimensions of its first argument)
    moments <- function(sum, sum_sq) {
      mean <- sum / n_map_draws
      variance <- (sum_sq - n_map_draws * mean ^ 2) / (n_map_draws - 1)
      list(mean = mean, sd = sqrt(pmax(variance, 0)))
    }

    # for the cells `rows` of `cells` (two_stage_cells()), cells x years: the
    # posterior mean and SD of mortality at every year, and at the panel years
    # the panel quantities
    map_chunk <- function(cells, rows) {
      chunk <- two_stage_chunk(setup, cells, rows)
      zeros <- function(n_years) matrix(0, length(rows), n_years)
      sums <- list(p = zeros(length(years_all)),
                   p_sq = zeros(length(years_all)),
                   local = zeros(length(map_years)),
                   local_sq = zeros(length(map_years)),
                   dynamical = zeros(length(map_years)))
      # seeded by the chunk's first cell, so a rerun with the same
      # IR_CUBE_MAP_WORKERS reproduces the maps (forked workers otherwise
      # draw from fresh seeds). The chunks, and so the Monte Carlo noise,
      # change with the number of workers
      set.seed(string_seed(paste("local noise", type, chunk$cells[1])))
      for (b in seq_along(setup$batches)) {
        logit <- two_stage_logit_batch(setup, b, chunk)
        for (j in seq_along(years_all)) {
          p <- plogis(logit$two_stage[[j]])
          sums$p[, j] <- sums$p[, j] + colSums(p)
          sums$p_sq[, j] <- sums$p_sq[, j] + colSums(p ^ 2)
        }
        noise <- cell_noise(setup$fit, chunk, map_years,
                            length(setup$batches[[b]]), local_terms)
        for (i in seq_along(map_years)) {
          p <- plogis(logit$two_stage[[panel[i]]] + t(noise[[i]]))
          sums$local[, i] <- sums$local[, i] + colSums(p)
          sums$local_sq[, i] <- sums$local_sq[, i] + colSums(p ^ 2)
          sums$dynamical[, i] <- sums$dynamical[, i] +
            colSums(plogis(logit$dynamical[[panel[i]]]))
        }
        rm(logit, noise, p)
      }
      correction_mean <- project_cells(setup$fit, fields_mean, chunk,
                                       map_years)
      local <- moments(sums$local, sums$local_sq)
      list(rows = rows,
           values = c(moments(sums$p, sums$p_sq),
                      list(dynamical_mortality = sums$dynamical / n_map_draws,
                           correction_mean = sapply(correction_mean,
                                                    function(x) x[, 1]),
                           local_mean = local$mean, local_sd = local$sd)))
    }

    # the mapped cells in one contiguous group per forked worker, which shares
    # the setup and builds the covariates of its own cells (a few GB for the
    # whole grid), then maps them in chunks. Each worker saves its results to
    # a temporary file, read back one group at a time: returning them all at
    # once to the parent held every result twice and ran out of memory
    map_group <- function(group) {
      cells <- two_stage_cells(setup, grid$cells[mapped[group]],
                               grid$country[mapped[group]])
      invisible(gc())
      chunks <- split(seq_along(group),
                      ceiling(seq_along(group) / chunk_size))
      results <- lapply(seq_along(chunks), function(i) {
        chunk <- map_chunk(cells, chunks[[i]])
        chunk$rows <- group[chunks[[i]]]
        if (i %% 10 == 0) {
          report("cells %i-%i: %i of %i chunks", min(group), max(group), i,
                 length(chunks))
        }
        chunk
      })
      file <- tempfile(fileext = ".rds")
      saveRDS(results, file, compress = FALSE)
      file
    }
    groups <- split(seq_len(n_mapped),
                    ceiling(seq_len(n_mapped) /
                              ceiling(n_mapped / n_map_workers)))
    done <- parallel::mclapply(groups, map_group, mc.cores = n_map_workers,
                               mc.preschedule = FALSE)
    # a worker that fails returns a try-error, and one that is killed (e.g.
    # out of memory) returns NULL, rather than its file name
    failed <- !vapply(done, function(x) {
      is.character(x) && !inherits(x, "try-error")
    }, logical(1))
    if (any(failed)) {
      stop("mapping failed for ", sum(failed), " groups of cells: ",
           paste(unique(unlist(done[failed])), collapse = "; "))
    }
    files <- unlist(done)

    # assembled into mapped cells x years, freeing each chunk's results
    results <- NULL
    for (file in files) {
      done <- readRDS(file)
      unlink(file)
      for (chunk in done) {
        if (is.null(results)) {
          results <- lapply(chunk$values, function(v) {
            matrix(NA_real_, n_mapped, ncol(v))
          })
        }
        for (q in names(results)) {
          results[[q]][chunk$rows, ] <- chunk$values[[q]]
        }
      }
      rm(done)
      invisible(gc())
    }
  })

  # one raster per year (and quantity), as R/predict.R; one layer per panel
  # year for the panel quantities
  write_cells <- function(values, file, layer_names = NULL) {
    values <- as.matrix(values)
    r <- rast(mask, nlyrs = ncol(values))
    full <- matrix(NA_real_, ncell(mask), ncol(values))
    full[grid$cells[mapped], ] <- values
    values(r) <- full
    if (!is.null(layer_names)) names(r) <- layer_names
    write_map(r, file)
  }
  dir.create(dirname(ir_map_file(type, years_all[1])), showWarnings = FALSE,
             recursive = TRUE)
  for (j in seq_along(years_all)) {
    write_cells(results$mean[, j], ir_map_file(type, years_all[j]))
    write_cells(results$sd[, j],
                ir_map_file(type, years_all[j], "susceptibility_sd"))
  }
  for (q in setdiff(names(results), c("mean", "sd"))) {
    write_cells(results[[q]], raster_file(type, q), layer_names = map_years)
  }
  write.csv(tibble(insecticide_type = type, n_draws = n_map_draws,
                   time_map_s = time_map[["elapsed"]],
                   peak_memory_gb = peak_memory_gb()),
            file.path(output_dir, sprintf("map_summary_%s.csv", type)),
            row.names = FALSE)
  report("%s mapped in %.0f s; peak memory %.1f GB", type,
         time_map[["elapsed"]], peak_memory_gb())
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

# a quantity at the panel years, masked: n*_pop and n*_local (n_star(), from
# the posterior means and SDs of the type, or of the types llin_effective
# weights, without and with u and p) and the difference from the dynamical
# model (two-stage minus dynamical posterior mean, percentage points, from the
# same draws) from the yearly rasters, the others from the panel-year rasters
read_map <- function(type, quantity) {
  yearly <- function(stem, output = type) {
    rast(sapply(map_years, function(year) ir_map_file(output, year, stem)))
  }
  weights <- if (type == "llin_effective") llin_weights else setNames(1, type)
  n_star_map <- function(mean, sd) {
    n_star(lapply(names(weights), mean), lapply(names(weights), sd),
           rho = hyperparameters$rho[match(names(weights),
                                           hyperparameters$insecticide_type)],
           weights = weights)
  }
  r <- switch(quantity,
              n_star_pop = n_star_map(
                function(t) yearly("susceptibility", t),
                function(t) yearly("susceptibility_sd", t)),
              n_star_local = n_star_map(
                function(t) rast(raster_file(t, "local_mean")),
                function(t) rast(raster_file(t, "local_sd"))),
              difference_pp = 100 * (yearly("susceptibility") -
                                       rast(raster_file(type,
                                                        "dynamical_mortality"))),
              rast(raster_file(type, quantity)))
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

# panels in rows of four with the legend in the last slot, or in ncol columns
# with the legend at the side
year_panels <- function(raster, fill_scale, title, subtitle, file,
                        ncol = NULL) {
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
  if (is.null(ncol)) years_list <- c(years_list, list(patchwork::guide_area()))
  patchwork::wrap_plots(years_list) +
    patchwork::plot_layout(guides = "collect", nrow = 2, ncol = ncol) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste(strwrap(subtitle, width = 110), collapse = "\n"))
  ggsave(file, bg = "white", width = if (is.null(ncol)) 13 else 11, height = 8,
         scale = 0.8, dpi = 300)
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

model_note <- function(type) {
  T_k <- hyperparameters$T[hyperparameters$insecticide_type == type]
  sprintf(paste("target m + omega + xi, without u and p; data to %i, later",
                "years are the AR(1) forecast of xi"), T_k)
}

# the posterior uncertainty as n*, under both readings of u and p
# (local_terms). The posterior mean is mapped by R/fig_ir_maps.R
n_star_target <- c(
  pop = paste("ilogit(m + omega + xi), u and p taken as data error: the upper",
              "bound on the map's worth"),
  local = paste("ilogit(m + omega + xi + u + p), u and p taken as real local",
                "variation: the lower bound"))
for (type in c(insecticides_plot, "llin_effective")) {
  mix <- type == "llin_effective"
  rho <- hyperparameters$rho[hyperparameters$insecticide_type %in%
                               if (mix) names(llin_weights) else type]
  for (reading in names(n_star_target)) {
    year_panels(
      # not 2030, far beyond the data
      subset(read_map(type, paste0("n_star_", reading)), "2030",
             negate = TRUE),
      # binned, diverging at one bioassay on a log scale: orange where the
      # prediction is worth less than one bioassay (where data are needed),
      # purple where it is worth more. The 10-class PuOr without its two
      # palest, so that the bins either side of 1 differ clearly in hue
      binned_scale("fill", palette = function(x) {
                     RColorBrewer::brewer.pal(10, "PuOr")[-(5:6)][seq_along(x)]
                   },
                   name = "n*<br>(bioassays)", transform = "log10",
                   breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
                   labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"),
                   limits = c(1e-6, 1e6), oob = scales::squish,
                   na.value = "transparent",
                   guide = guide_coloursteps(frame.colour = border_col,
                                             frame.linewidth = 0.1)),
      title = sprintf("%s: two-stage model, uncertainty as n*_%s",
                      if (mix) "LLIN pyrethroids" else type, reading),
      subtitle = paste0(
        "n*_", reading, ": the number of independent bioassays of 100 ",
        "mosquitoes at the pixel-year",
        if (mix) ", split across the pyrethroids by weight,",
        " whose mean would estimate population mortality as precisely as the ",
        "posterior of ", n_star_target[[reading]], ". Bioassay noise: ",
        "within-pixel-year only, binomial with overdispersion rho ",
        paste(unique(sprintf("%.2f", range(rho))), collapse = "-"),
        ", from replicate bioassays.",
        if (mix) paste(" The pyrethroids' errors taken as perfectly",
                       "correlated, which gives the largest SD."),
        " Orange: worth less than one bioassay"),
      file = file.path(figure_dir,
                       sprintf("%s_n_star_%s_map.png", type, reading)),
      ncol = 3)
  }
}

for (type in insecticides_plot) {
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
