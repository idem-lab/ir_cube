# Predictions of the two-stage model (#21) at any mask cells and years, from
# the per-type fits and paired dynamical draws that R/two_stage_maps.R saves
# (outputs/two_stage/maps/<type>/fit.rds and dynamical.rds). Used by the map
# step of R/two_stage_maps.R and by the figure scripts. Functions only; source
# from the repo root after R/packages.R. Sources the two-stage and dynamical
# prediction code it builds on.
#
# Draw d pairs dynamical draw d with a latent draw from N(mode, H^-1), shifted
# by the cut-posterior formula for it (correction_node_draws()). The target is
# m + omega + xi on the logit scale. u and p are observation-level noise: they
# enter as fresh draws, in predictions of new assays (noise = TRUE) and in the
# maps' n*_local (cell_noise()).

source("R/two_stage_helpers.R")
source("R/dynamical_predictions.R")
source("R/two_stage_correction.R")
source("R/two_stage_map_functions.R")

two_stage_maps_dir <- "outputs/two_stage/maps"

# The saved fit of `type` and the correction's node fields for n_draws of its
# 2000 paired dynamical draws (an even subset, by the scoring's rule), in
# batches of batch_size, at `years`; the fields are drawn from the stream
# seeded by `seed_key`. Two setups with the same n_draws and batch_size use the
# same dynamical draws, batch for batch, whatever the type, so their
# predictions can be combined draw by draw (e.g. over the LLIN pyrethroids).
# `data` is the modelled bioassays (df), which must be those the fit was paired
# with (check_two_stage_data())
two_stage_setup <- function(type, years, data, n_draws = 1000,
                            batch_size = 100,
                            seed_key = paste("two-stage draws", type),
                            baseline_year = 1995,
                            maps_dir = two_stage_maps_dir) {
  fit <- readRDS(file.path(maps_dir, type, "fit.rds"))
  dynamical <- readRDS(file.path(maps_dir, type, "dynamical.rds"))
  check_two_stage_data(dynamical, data, type)
  stopifnot(fit$t0 == baseline_year,
            max(abs(fit$m_ref - colMeans(dynamical$logit_train))) < 1e-12)
  draws <- thin_draws(matrix(seq_len(nrow(dynamical$logit_train))),
                      n_draws)[, 1]
  batches <- split(draws, ceiling(seq_along(draws) / batch_size))
  set.seed(string_seed(seed_key))
  fields <- lapply(batches, function(d) {
    f <- correction_node_draws(
      fit, years, length(d),
      m_draws_train = dynamical$logit_train[d, , drop = FALSE])
    f$theta <- NULL
    f
  })
  list(type = type,
       k = match(type, dynamical$parameters$types),
       fit = fit,
       years = years,
       baseline_year = baseline_year,
       batches = batches,
       fields = fields,
       n_draws = length(draws),
       parameters = lapply(batches, subset_draws,
                           parameters = dynamical$parameters),
       logit_init = dynamical$logit_init,
       design = dynamical$parameters$options$selection_columns)
}

# Stop unless the bioassays of `type` in `data` are those whose paired draws
# `dynamical` (dynamical.rds) holds, row for row: the draws, the fit and the
# u and p keys are matched to the data by row order
check_two_stage_data <- function(dynamical, data, type) {
  rows <- data[data$insecticide_type == type, names(dynamical$keys)]
  if (!isTRUE(all.equal(as.data.frame(rows), as.data.frame(dynamical$keys),
                        check.attributes = FALSE))) {
    stop("the ", type, " bioassays differ from those the two-stage fit was ",
         "paired with; rerun R/two_stage_maps.R from prepare")
  }
}

# The cells' covariates and coordinates: `cells` are mask cell numbers,
# `country` their country names (those of dimnames(setup$logit_init)[[2]]; for
# data cells, the fit's country of the cell, data_cell_country()). The
# covariates are kept in the compact form of map_covariates();
# two_stage_chunk() assembles them for a chunk of cells
two_stage_cells <- function(setup, cells, country,
                            mask = terra::rast("data/clean/raster_mask.tif")) {
  country_index <- match(country, dimnames(setup$logit_init)[[2]])
  stopifnot(length(country) == length(cells), !anyNA(country_index))
  xy <- terra::xyFromCell(mask, cells)
  list(cells = cells,
       country = country_index,
       covariates = map_covariates(cells, setup$baseline_year,
                                   max(setup$years), setup$design),
       coords = project_km(xy[, 1], xy[, 2]))
}

# The rows `rows` of `cells` (two_stage_cells()), with their covariates as the
# cells x years x n_covs array dynamical_logit_cells() takes, and the basis of
# the latent smooths at the cells (hsgp_basis(), from their coordinates in
# 1,000 km)
two_stage_chunk <- function(setup, cells, rows = seq_along(cells$cells)) {
  list(cells = cells$cells[rows],
       country = cells$country[rows],
       x = map_x(cells$covariates, rows,
                 max(setup$years) - setup$baseline_year + 1),
       x_init = cells$covariates$init[rows, , drop = FALSE],
       coords = cells$coords[rows, , drop = FALSE],
       basis = hsgp_basis(cells$coords[rows, , drop = FALSE] / 1000,
                          setup$parameters[[1]]$options$smooth))
}

# The correction omega + xi at the cells of `chunk` (two_stage_chunk()) for
# every year of the setup, as project_correction() computes it for the rows
# of a cells x years table, but projecting omega and each year's xi once per
# cell rather than once per cell-year: a list by year of cells x draws. With
# noise = TRUE, fresh draws of the iid terms are added (cell_noise()), as
# project_correction() does. The meshes cover the prediction mask (checked by
# the fit step of R/two_stage_maps.R)
project_cells <- function(fit, fields, chunk, years, noise = FALSE) {
  omega <- as.matrix(mesh_basis(fit$mesh, chunk$coords) %*% fields$omega)
  A_xi <- mesh_basis(fit$mesh_xi, chunk$coords)
  out <- lapply(years, function(year) {
    if (year <= fit$t0) return(omega)
    omega + as.matrix(A_xi %*% fields$xi[[as.character(year)]])
  })
  if (noise) out <- Map(`+`, out, cell_noise(fit, chunk, years, ncol(omega)))
  out
}

# Fresh draws of the iid terms named in `terms` (iid_noise()) at the cells of
# `chunk` and `years`, one per level of each key in the cells x years table: a
# list by year of cells x n_draws
cell_noise <- function(fit, chunk, years, n_draws,
                       terms = names(correction_iid_terms)) {
  n_cells <- length(chunk$cells)
  fresh <- iid_noise(fit, tibble(cell = rep(chunk$cells, length(years)),
                                 year = rep(years, each = n_cells)),
                     n_draws, terms)
  lapply(seq_along(years), function(j) {
    fresh[(j - 1) * n_cells + seq_len(n_cells), , drop = FALSE]
  })
}

# The logit of the dynamical model, m, and of the two-stage model, m + omega +
# xi (plus fresh u and p if noise = TRUE), for the draws of batch b at the
# cells of `chunk` (two_stage_chunk()): lists by year of draws x cells
two_stage_logit_batch <- function(setup, b, chunk, noise = FALSE) {
  draws <- setup$batches[[b]]
  year_index <- setup$years - setup$baseline_year + 1
  m <- dynamical_logit_cells(
    setup$parameters[[b]], setup$k,
    matrix(setup$logit_init[draws, chunk$country, setup$k], length(draws)),
    chunk$x, year_index, x_init = chunk$x_init,
    smooths = cell_smooths(setup$parameters[[b]], chunk$basis))
  correction <- project_cells(setup$fit, setup$fields[[b]], chunk,
                              setup$years, noise = noise)
  out <- list(dynamical = list(), two_stage = list())
  for (j in seq_along(setup$years)) {
    m_j <- m[[as.character(year_index[j])]]
    out$dynamical[[j]] <- m_j
    out$two_stage[[j]] <- m_j + t(correction[[j]])
  }
  out
}

# Draws of weighted averages of mortality over cells, for every year of the
# setup: `cells` are mask cell numbers, `country` their country names (as for
# two_stage_cells()), and `weights` a cells x groups matrix (columns summing to
# 1 for a weighted mean, or the identity for each cell alone). Returns draws x
# groups x years arrays of the two-stage model (ilogit(m + omega + xi), with
# fresh u and p if noise = TRUE) and of the dynamical model (ilogit(m)), with
# the groups and years as dimnames. Cells of zero weight in every group are
# skipped
two_stage_weighted_draws <- function(setup, cells, country, weights,
                                     noise = FALSE) {
  weights <- as.matrix(weights)
  stopifnot(nrow(weights) == length(cells))
  rows <- which(rowSums(weights != 0) > 0)
  w <- weights[rows, , drop = FALSE]
  dims <- c(setup$n_draws, ncol(w), length(setup$years))
  names <- list(NULL, colnames(weights), setup$years)
  out <- list(two_stage = array(NA_real_, dims, names),
              dynamical = array(NA_real_, dims, names))
  chunk <- two_stage_chunk(setup,
                           two_stage_cells(setup, cells[rows], country[rows]))
  start <- 0
  for (b in seq_along(setup$batches)) {
    index <- start + seq_along(setup$batches[[b]])
    logit <- two_stage_logit_batch(setup, b, chunk, noise = noise)
    for (j in seq_along(setup$years)) {
      out$two_stage[index, , j] <- plogis(logit$two_stage[[j]]) %*% w
      out$dynamical[index, , j] <- plogis(logit$dynamical[[j]]) %*% w
    }
    start <- max(index)
  }
  out
}

# The selection design (options$selection_columns) of the dynamical fit that
# the two-stage fits were built on
two_stage_design <- function(type, maps_dir = two_stage_maps_dir) {
  dynamical <- readRDS(file.path(maps_dir, type, "dynamical.rds"))
  dynamical$parameters$options$selection_columns
}

# The country of each data cell (cell_id) in the dynamical model: that of the
# cell's first record (dynamical_lookups()), by name
data_cell_country <- function(df) {
  lookups <- dynamical_lookups(df)
  lookups$levels$countries[lookups$cell_country_lookup]
}

# Posterior predictive draws of the mortality of new assays at every training
# assay (row of df), for n_draws of the paired draws, from both models, for
# the within-sample residual checks:
#   two_stage  ilogit(m + omega + xi + fresh u + fresh p), draws x rows: u is
#              shared by the assays of a pixel-year and p by those of a pixel,
#              within a draw (predict_correction())
#   dynamical  ilogit(m), draws x rows
# and the overdispersion each model scores new assays at, draws x rows: the
# two-stage model the external per-type rho (rho_lookup()), the dynamical
# model its own rho_types
assay_draws <- function(df, types, n_draws = 1000,
                        maps_dir = two_stage_maps_dir) {
  out <- list(two_stage = matrix(NA_real_, n_draws, nrow(df)),
              dynamical = matrix(NA_real_, n_draws, nrow(df)),
              rho_two_stage = matrix(rho_for_record(df, rho_lookup()),
                                     n_draws, nrow(df), byrow = TRUE),
              rho_dynamical = matrix(NA_real_, n_draws, nrow(df)))
  for (k in seq_along(types)) {
    rows <- which(df$type_id == k)
    fit <- readRDS(file.path(maps_dir, types[k], "fit.rds"))
    dynamical <- readRDS(file.path(maps_dir, types[k], "dynamical.rds"))
    check_two_stage_data(dynamical, df, types[k])
    stopifnot(identical(dynamical$parameters$types[k], types[k]))
    draws <- thin_draws(matrix(seq_len(nrow(dynamical$logit_train))),
                        n_draws)[, 1]
    stopifnot(length(draws) == n_draws)
    m <- dynamical$logit_train[draws, , drop = FALSE]
    train <- tibble(lon = df$longitude[rows], lat = df$latitude[rows],
                    year = df$year[rows], cell = df$cell[rows])
    set.seed(string_seed(paste("assay draws", types[k])))
    out$two_stage[, rows] <- plogis(predict_correction(
      fit, train, m_draws_train = m, m_draws_new = m, n_draws = n_draws))
    out$dynamical[, rows] <- plogis(m)
    out$rho_dynamical[, rows] <- dynamical$parameters$rho_types[draws, k]
    rm(fit, dynamical, m)
    invisible(gc())
  }
  out
}
