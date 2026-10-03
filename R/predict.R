# Predict the fraction susceptible (predicted bioassay mortality) on the full
# prediction grid, from the posterior draws of the fitted model
# (temporary/fitted_model.RData), for every insecticide type and for the
# effective susceptibility to the pyrethroids in single-AI LLINs.
#
#   Rscript R/predict.R
#
# Run with the greta environment of R/packages.R (the draws are a greta
# object), e.g.
#   R_LIBS=~/R/greta06-lib RETICULATE_PYTHON=.../greta06-env/bin/python \
#     OPENBLAS_NUM_THREADS=1 nice -n 10 Rscript R/predict.R
#
# The predictions are computed in plain R from the draws
# (R/dynamical_predictions.R, R/two_stage_map_functions.R), as
# R/two_stage_maps.R does. Countries and regions without data take their
# initial state from the hierarchical prior, drawn once per posterior draw
# with a fixed seed, so the maps are reproducible. Cells without a country in
# data/clean/country_raster.tif are NA. The maps are of predicted bioassay
# mortality, with the mortality floor if the fit has one.
#
# Writes, for each of the nine types and llin_effective, and each year from the
# baseline year to 2030:
#   outputs/ir_maps/<type>/ir_<year>_susceptibility.tif     posterior mean
#   outputs/ir_maps/<type>/ir_<year>_susceptibility_sd.tif  posterior SD
# llin_effective is the mortality of each draw weighted over the active
# ingredients by temporary/ingredient_weights.RDS.

source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")
# report()
source("R/two_stage_helpers.R")

# posterior draws to average over: an even subset of the 2000 thinned draws of
# paired_draw_index(). The Monte Carlo SE of a mean is the posterior SD /
# sqrt(n_draws)
n_draws <- 500

# the last year to predict to; covariates are carried forward from their last
# year
end_year <- 2030

# cells per chunk, and parallel workers (one output each at a time), or the
# environment variables IR_CUBE_PREDICT_CHUNK and IR_CUBE_PREDICT_WORKERS. The
# parent holds about 4 GB, most of it shared with the workers, and each worker
# about 2-2.5 GB more at chunks of 1000 cells (a cells x years mean and SD is
# 0.85 GB of that). Total memory (PSS) peaks at about 9 GB with 2 workers and
# 12.6 GB with 3; the whole run takes about 95 and 67 minutes
chunk_size <- as.integer(Sys.getenv("IR_CUBE_PREDICT_CHUNK", "1000"))
n_workers <- as.integer(Sys.getenv("IR_CUBE_PREDICT_WORKERS", "2"))

output_dir <- "outputs/ir_maps"

time_start <- Sys.time()


# posterior draws of the parameters ------------------------------------------

fit_env <- new.env()
load("temporary/fitted_model.RData", envir = fit_env)
baseline_year <- fit_env$baseline_year
types <- fit_env$types
classes_index <- fit_env$classes_index
countries <- fit_env$countries
regions <- fit_env$regions
df <- fit_env$df
fold <- list(draws = fit_env$draws, options = fit_env$model_options,
             x_cells_init = fit_env$x_cells_init)
options <- fold$options
rm(fit_env)
invisible(gc())

parameters <- dynamical_parameter_draws(fold, classes_index, types, df,
                                        options = options)

# the logit relative initial states for every country in the lookup. The
# prior draws for countries and regions without data are made for all 2000
# draws and then subset, so a draw has the same initial states whatever
# n_draws is
set.seed(1)
logit_init_all <- map_logit_init(parameters, countries, regions, df)
stopifnot(isTRUE(all.equal(logit_init_all[, countries, , drop = FALSE],
                           parameters$logit_init_relative,
                           check.attributes = FALSE)))

predict_draws <- round(seq(1, parameters$n_draws, length.out = n_draws))
parameters <- subset_draws(parameters, predict_draws)
logit_init <- logit_init_all[predict_draws, , , drop = FALSE]
rm(logit_init_all, fold)
invisible(gc())
report("%i draws of %i types", n_draws, length(types))


# the prediction grid --------------------------------------------------------

mask <- rast("data/clean/raster_mask.tif")
cells <- terra::cells(mask)
n_cells <- length(cells)
years_predict <- baseline_year:end_year
covariates <- map_covariates(cells, baseline_year, end_year,
                             options$selection_columns)

country_raster <- rast("data/clean/country_raster.tif")
cell_country <- as.character(terra::extract(country_raster,
                                            cells)$country_name)
cell_country_index <- match(cell_country, dimnames(logit_init)[[2]])
report("%i cells: %i without a country and %i in a country outside the lookup (both NA)",
       n_cells, sum(is.na(cell_country)),
       sum(!is.na(cell_country) & is.na(cell_country_index)))

chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))


# predictions ----------------------------------------------------------------

# each output as the weights of the types it combines
ingredient_weights <- unlist(readRDS("temporary/ingredient_weights.RDS"))
stopifnot(all(names(ingredient_weights) %in% types))
# llin_effective first, as it takes longest
outputs <- c(list(llin_effective = ingredient_weights),
             lapply(setNames(types, types), function(type) setNames(1, type)))

# posterior mean and SD of one output at every cell and year, written as one
# raster per year
predict_output <- function(output) {
  weights <- outputs[[output]]
  mean_out <- sd_out <- matrix(NA_real_, n_cells, length(years_predict))
  for (chunk in chunks) {
    ok <- chunk[!is.na(cell_country_index[chunk])]
    if (length(ok) == 0) next
    p <- NULL
    x_chunk <- map_x(covariates, ok, length(years_predict))
    for (type in names(weights)) {
      k <- match(type, types)
      dyn <- dynamical_logit_cells(
        parameters, k, matrix(logit_init[, cell_country_index[ok], k], n_draws),
        x_chunk, seq_along(years_predict),
        x_init = covariates$init[ok, , drop = FALSE])
      p_type <- lapply(dyn, function(x) weights[[type]] * plogis(x))
      p <- if (is.null(p)) p_type else Map(`+`, p, p_type)
      rm(dyn, p_type)
    }
    for (j in seq_along(years_predict)) {
      mean_out[ok, j] <- colMeans(p[[j]])
      sd_out[ok, j] <- col_sds(p[[j]])
    }
    rm(p, x_chunk)
  }

  # a fresh raster handle in this process
  template <- rast("data/clean/raster_mask.tif")
  dir.create(file.path(output_dir, output), showWarnings = FALSE,
             recursive = TRUE)
  write_layer <- function(values, file) {
    full <- rep(NA_real_, ncell(template))
    full[cells] <- values
    r <- rast(template)
    values(r) <- full
    writeRaster(r, file, overwrite = TRUE, datatype = "FLT4S",
                gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  }
  for (j in seq_along(years_predict)) {
    stem <- file.path(output_dir, output,
                      sprintf("ir_%i_susceptibility", years_predict[j]))
    write_layer(mean_out[, j], paste0(stem, ".tif"))
    write_layer(sd_out[, j], paste0(stem, "_sd.tif"))
  }
}

# forked workers, each writing its own outputs
done <- parallel::mclapply(names(outputs), function(output) {
  time <- system.time(predict_output(output))[["elapsed"]]
  report("%-18s %.0f s", output, time)
  output
}, mc.cores = n_workers, mc.preschedule = FALSE)
failed <- vapply(done, inherits, logical(1), "try-error")
if (any(failed)) {
  stop("prediction failed for ",
       paste(names(outputs)[failed], collapse = ", "), ": ",
       paste(unique(unlist(done[failed])), collapse = "; "))
}
report("done in %.0f min",
       as.numeric(difftime(Sys.time(), time_start, units = "mins")))
