# Second-stage geostatistical correction to the dynamical model (#21).
#
# A Gaussian model for the dynamical model's residuals on the logit scale,
# conditioning on the dynamical model without ever updating it (a cut
# posterior):
#
#   lambda_i = m_i + omega(s_i) + xi(s_i, t_i) + u[pixel-year(i)] + p[pixel(i)]
#
# m is the dynamical model's posterior mean logit at the training assays
# (m_ref). omega is a Matern field correcting the initial conditions; xi
# accumulates an AR(1)-in-time, Matern-in-space field of annual anomalies eta
# from xi(., t0) = 0, so that forecasts of the correction plateau; u (per
# pixel-year) and p (per pixel) are iid.
#
# The prediction target (what the maps show) is m + omega + xi. u and p are
# observation-level noise, shared by the assays of a pixel-year or a pixel:
# they are fitted so that this noise stays out of omega and xi, and enter
# predictions of new assays only as fresh draws from N(0, tau^2) and
# N(0, sigma_p^2), never as their fitted values, even at pixel-years or
# pixels with training data.
#
# Fitting, per insecticide type (fit_correction()):
#   stage A  z, v = empirical logit and its variance, inflated by the per-type
#            replicate-based rho; hyperparameters by penalised maximum
#            marginal likelihood (TMB, tmb/two_stage_correction.cpp), latent
#            posterior exactly Gaussian given them (fit_stage_a());
#   stage B  PQL on the beta-binomial counts from the stage-A fit, with one
#            re-estimation of the hyperparameters (R/two_stage_pql.R).
# The result's latent posterior is N(mode, H^-1); predictions draw from it.
#
# API (source this file from the repo root; it sources R/two_stage_pql.R):
#
#   meshes <- build_correction_meshes(coords_km, prediction_mask_coords())
#                                                  list(omega, xi): the final
#                                                  model's meshes, covering the
#                                                  sites and the mask
#   fit <- fit_correction(train, t0, T, meshes)    stage A then stage B.
#     train: lon, lat (or x_km, y_km), year, cell, m (= m_ref), died,
#     mosquito_number, rho. T: last year xi is represented (forecast beyond).
#     Stops if stage A does not converge; check fit$opt$convergence for the
#     stage-B re-estimation. fit$hyper: sigma_omega, range_omega, sigma_eta,
#     range_eta, phi, persistence, tau (SD of u), sigma_p (and the kappas).
#     fit$stage_b: PQL diagnostics, incl. hyper_a and objective_a (stage A).
#   fit_summary(fit)                               one-row tibble for tables
#   fields <- correction_node_draws(fit, years, n_draws, m_draws_train = NULL,
#                                   mean = FALSE)
#     (a) joint latent draws (theta), with the cut-posterior shift per
#     dynamical draw when m_draws_train (n_draws x n_obs) is given, and the
#     node fields: omega (nodes x draws) and xi[[year]] (xi nodes x draws) for
#     every requested year, AR(1)-forecast beyond T. mean = TRUE gives the
#     posterior mean instead (the mode, mean forecast; one column).
#   project_correction(fit, fields, new, noise = FALSE)
#     (b) points x draws correction omega + xi at the rows of `new`
#     (lon/lat or x_km/y_km, year, and cell if noise = TRUE). noise = TRUE
#     adds fresh u per distinct pixel-year and fresh p per distinct pixel of
#     `new` (shared by its rows within a draw; caller's RNG), as for new
#     assays; maps use noise = FALSE. Warns for points outside the mesh
#     (their fields are 0).
#   predict_correction(fit, new, m_draws_train, m_draws_new, n_draws)
#     draws x points m_draw + omega + xi + fresh u + fresh p: the predictive
#     distribution of new assays' logit mortality, in batches (CV scoring).
#     map = TRUE also returns the same draws without u and p (the map).
#   Saving without refitting: fit$obj <- NULL; saveRDS(fit, file). Everything
#     above works on the reloaded fit (it keeps H, H_chol, meshes and levels);
#     only correction_adfun() needs tmb_data and par_list, which it keeps too.
#
# Speed depends heavily on the BLAS used by CHOLMOD's supernodal
# factorisation: run with OpenBLAS (reference BLAS is ~10x slower), and the
# template must be compiled with TMBad (CppAD is ~100x slower).

suppressMessages({
  library(TMB)
  library(Matrix)
  library(fmesher)
  library(sf)
  library(dplyr)
})
source("R/two_stage_pql.R")

correction_template <- "tmb/two_stage_correction.cpp"

# compile the template if its shared object is missing or stale, and load it
load_correction_template <- function(path = correction_template) {
  dll_path <- TMB::dynlib(sub("\\.cpp$", "", path))
  if (!file.exists(dll_path) || file.mtime(dll_path) < file.mtime(path)) {
    TMB::compile(path, flags = "-O2", framework = "TMBad")
  }
  dll_name <- basename(sub("\\.cpp$", "", path))
  if (!dll_name %in% names(getLoadedDLLs())) {
    dyn.load(dll_path)
  }
  dll_name
}

# empirical logit of bioassay mortality, and its approximate sampling variance
# inflated by the beta-binomial design effect 1 + (n - 1) rho
empirical_logit <- function(died, n, rho) {
  z <- log((died + 0.5) / (n - died + 0.5))
  v <- (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
  list(z = z, v = v)
}

# The iid noise terms, each as the key grouping assays into its levels, and
# the name of its SD in fit$hyper. Adding or removing a term here is all the
# model and the predictions need
correction_iid_terms <- list(
  u = list(key = function(d) paste(d$cell, d$year), sd = "tau"),
  p = list(key = function(d) as.character(d$cell), sd = "sigma_p")
)


# coordinates and meshes ----------------------------------------------------

# Albers equal-area conic for Africa, in km, so that the Matern range means
# the same distance everywhere and the range priors can be stated in km
africa_equal_area_crs <- paste(
  "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
)

project_km <- function(lon, lat, crs = africa_equal_area_crs) {
  xy <- sf::sf_project(from = "+proj=longlat +datum=WGS84 +no_defs", to = crs,
                       pts = cbind(lon, lat))
  colnames(xy) <- c("x_km", "y_km")
  xy
}

# projected coordinates of a data frame: x_km and y_km if present (simulated
# data), otherwise lon and lat projected
coords_km <- function(df) {
  if (all(c("x_km", "y_km") %in% names(df))) {
    cbind(x_km = df$x_km, y_km = df$y_km)
  } else {
    project_km(df$lon, df$lat)
  }
}

# Projected centres of the prediction mask's cells, on a coarse grid (every
# `fact`-th cell), for the outer boundary of the meshes. Callers pass them to
# build_correction_meshes()
prediction_mask_coords <- function(file = "data/clean/raster_mask.tif",
                                   fact = 10) {
  mask <- terra::aggregate(terra::rast(file), fact = fact, fun = "min",
                           na.rm = TRUE)
  xy <- terra::xyFromCell(mask, terra::cells(mask))
  project_km(xy[, 1], xy[, 2])
}

# One fmesher mesh with a node at every site, sites closer than `cutoff`
# merged, and triangles of at most max_edge_inner km between them inside a
# buffered non-convex hull of the sites. Resolution follows data density. If
# that inner mesh exceeds max_nodes, the cutoff grows (coarsening the densest
# clusters first) until it does not. Around it, a coarse outer region (edges of
# at most max_edge_outer km) out to outer_buffer km beyond the sites and beyond
# `outer_coords` (the prediction mask, prediction_mask_coords(), so that every
# map cell is inside the mesh), which keeps the SPDE's boundary effects away
# from the data and the maps
build_correction_mesh <- function(coords_km, outer_coords, max_edge_inner = 250,
                                  max_edge_outer = 800, cutoff = 30,
                                  inner_buffer = 200, outer_buffer = 1500,
                                  max_nodes = 2500, cutoff_growth = 1.2) {
  coords_km <- unique(as.matrix(coords_km))
  # a smooth (coarse resolution) hull avoids many short boundary segments
  inner <- suppressWarnings(
    fmesher::fm_nonconvex_hull(coords_km, convex = inner_buffer,
                               resolution = c(60, 60), format = "fm")
  )
  # the cutoff is chosen on the mesh without the mask
  repeat {
    mesh <- fmesher::fm_mesh_2d(loc = coords_km, boundary = list(inner),
                                max.edge = c(max_edge_inner, max_edge_outer),
                                cutoff = cutoff,
                                offset = c(-0.01, outer_buffer))
    if (mesh$n <= max_nodes) break
    cutoff <- cutoff * cutoff_growth
  }
  outer <- suppressWarnings(
    fmesher::fm_nonconvex_hull(rbind(coords_km, as.matrix(outer_coords)),
                               convex = outer_buffer,
                               resolution = c(60, 60), format = "fm")
  )
  mesh <- fmesher::fm_mesh_2d(loc = coords_km, boundary = list(inner, outer),
                              max.edge = c(max_edge_inner, max_edge_outer),
                              cutoff = cutoff)
  attr(mesh, "cutoff") <- cutoff
  mesh
}

# The meshes of the final model (omega5000_xi2500 in doc/two_stage_plan.md):
# omega on a fine mesh (15 km cutoff, inner edges of at most 150 km, at most
# 5000 nodes in the data region), xi on one of at most 2500. xi's latent
# dimension is multiplied by the number of years, which is why its mesh is
# coarser. Both cover `outer_coords`, the prediction mask
# (prediction_mask_coords())
build_correction_meshes <- function(coords_km, outer_coords) {
  list(omega = build_correction_mesh(coords_km, outer_coords, cutoff = 15,
                                     max_edge_inner = 150, max_nodes = 5000),
       xi = build_correction_mesh(coords_km, outer_coords, max_nodes = 2500))
}

# FEM matrices in the form R_inla::spde_t expects (c0 lumped, as in INLA)
correction_fem <- function(mesh) {
  fem <- fmesher::fm_fem(mesh, order = 2)
  as_dgc <- function(x) as(as(as(x, "dMatrix"), "generalMatrix"),
                           "CsparseMatrix")
  list(M0 = as_dgc(fem$c0), M1 = as_dgc(fem$g1), M2 = as_dgc(fem$g2))
}

# Matern precision at the nodes with marginal SD sigma, as in the template
matern_precision_r <- function(fem, kappa, sigma) {
  tau_spde <- 1 / (sigma * kappa * sqrt(4 * pi))
  tau_spde ^ 2 * (kappa ^ 4 * fem$M0 + 2 * kappa ^ 2 * fem$M1 + fem$M2)
}

# sparse projection from points to mesh nodes
mesh_basis <- function(mesh, coords) {
  A <- fmesher::fm_basis(mesh, loc = coords)
  as(as(A, "generalMatrix"), "CsparseMatrix")
}

# observations -> vec(x), x the xi nodes in years t0 + 1, ..., T: an
# observation in year t touches column t - t0, one in year t0 none
xi_design <- function(mesh_xi, coords, year, t0, T) {
  trip <- summary(mesh_basis(mesh_xi, coords))
  year_col <- year[trip$i] - t0
  keep <- year_col > 0 & year_col <= T - t0
  Matrix::sparseMatrix(i = trip$i[keep],
                       j = trip$j[keep] + (year_col[keep] - 1) * mesh_xi$n,
                       x = trip$x[keep],
                       dims = c(length(year), mesh_xi$n * (T - t0)))
}

# PC priors: P(range < 50 km) = 0.05 and P(sigma > 1) = 0.05 for both Matern
# fields, P(sd > 1) = 0.05 for each iid term, and persistence 1 / (1 - phi) ~
# lognormal(log 5, 0.62): a median of five years, 95% interval ~1.5-17
correction_priors <- function(range0 = 50, alpha_range = 0.05, sigma0 = 1,
                              alpha_sigma = 0.05, iid_sigma0 = 1,
                              alpha_iid = 0.05, persistence_median = 5,
                              persistence_sdlog = 0.62) {
  list(pc_omega = c(range0, alpha_range, sigma0, alpha_sigma),
       pc_eta = c(range0, alpha_range, sigma0, alpha_sigma),
       pc_iid = c(iid_sigma0, alpha_iid),
       persistence_prior = c(log(persistence_median), persistence_sdlog))
}


# fitting ----------------------------------------------------------------------

# the TMB object; fix_hyper = TRUE holds the hyperparameters at their values
# in `parameters` (used by the simulation check)
correction_adfun <- function(data, parameters, fix_hyper = FALSE,
                             silent = TRUE) {
  map <- list()
  if (fix_hyper) {
    hyper <- setdiff(names(parameters), c("w_omega", "x", "iid"))
    map <- lapply(parameters[hyper], function(x) factor(rep(NA, length(x))))
  }
  TMB::MakeADFun(data = data, parameters = parameters, map = map,
                 random = c("w_omega", "x", "iid"),
                 DLL = load_correction_template(), silent = silent)
}

# Stage A, or a refit of the hyperparameters to given (z, v). train as for
# fit_correction(); if it has z and v they are the response, otherwise the
# empirical logit is. meshes from build_correction_meshes()
fit_stage_a <- function(train, t0, T = max(train$year), meshes,
                        priors = correction_priors(),
                        start = list(),
                        control = list(eval.max = 1000, iter.max = 500),
                        silent = TRUE) {

  start_time <- Sys.time()
  stopifnot(T > t0, all(train$year >= t0 & train$year <= T))
  if (!all(c("z", "v") %in% names(train))) {
    response <- empirical_logit(train$died, train$mosquito_number, train$rho)
    train$z <- response$z
    train$v <- response$v
  }
  coords <- coords_km(train)
  n_obs <- nrow(train)

  # iid levels: each term's distinct keys, stacked into one vector
  iid_levels <- lapply(correction_iid_terms, function(term) {
    unique(term$key(train))
  })
  offset <- cumsum(c(0, lengths(iid_levels)))
  iid_index <- matrix(0L, n_obs, length(iid_levels))
  for (k in seq_along(iid_levels)) {
    iid_index[, k] <- match(correction_iid_terms[[k]]$key(train),
                            iid_levels[[k]]) + offset[k]
  }

  A_omega <- mesh_basis(meshes$omega, coords)
  A_xi <- xi_design(meshes$xi, coords, train$year, t0, T)
  data <- c(
    list(z = train$z, v = train$v, m = train$m, A_omega = A_omega,
         A_xi = A_xi, iid_index = iid_index - 1L,
         iid_term = rep(seq_along(iid_levels) - 1L, lengths(iid_levels)),
         spde = correction_fem(meshes$omega),
         spde_xi = correction_fem(meshes$xi)),
    priors
  )

  # ranges of a few hundred km, well inside the prior, and small SDs: the
  # correction should be modest if the dynamical model is good
  start <- modifyList(list(log_sigma_omega = log(0.5),
                           log_kappa_omega = log(sqrt(8) / 300),
                           log_sigma_eta = log(0.2),
                           log_kappa_eta = log(sqrt(8) / 300),
                           logit_phi = qlogis(0.8),
                           log_sigma_iid = rep(log(0.3), length(iid_levels))),
                      start)
  parameters <- c(list(w_omega = rep(0, meshes$omega$n),
                       x = matrix(0, meshes$xi$n, T - t0),
                       iid = rep(0, sum(lengths(iid_levels)))),
                  start[c("log_sigma_omega", "log_kappa_omega",
                          "log_sigma_eta", "log_kappa_eta", "logit_phi",
                          "log_sigma_iid")])

  obj <- correction_adfun(data, parameters, silent = silent)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
  optimise_time <- Sys.time()

  # latent mode and its Hessian at the optimum; evaluating fn at opt$par
  # leaves the inner problem solved at exactly these hyperparameters
  obj$fn(opt$par)
  par_full <- obj$env$last.par
  random <- obj$env$random
  block <- names(par_full)[random]
  blocks <- split(seq_along(random), factor(block, levels = unique(block)))
  blocks <- c(blocks[c("w_omega", "x")],
              lapply(seq_along(iid_levels), function(k) {
                blocks$iid[(offset[k] + 1):offset[k + 1]]
              }) |> setNames(names(iid_levels)))
  H <- obj$env$spHess(par_full, random = TRUE)
  H <- Matrix::forceSymmetric(as(H, "CsparseMatrix"), uplo = "L")

  # observations -> latent vector, in the latent order, for the
  # cut-posterior shift and PQL
  A_iid <- Matrix::sparseMatrix(i = rep(seq_len(n_obs), ncol(iid_index)),
                                j = as.vector(iid_index), x = 1,
                                dims = c(n_obs, sum(lengths(iid_levels))))

  par_list <- obj$env$parList(par = par_full)
  phi <- plogis(par_list$logit_phi)
  hyper <- list(sigma_omega = exp(par_list$log_sigma_omega),
                range_omega = sqrt(8) / exp(par_list$log_kappa_omega),
                kappa_omega = exp(par_list$log_kappa_omega),
                sigma_eta = exp(par_list$log_sigma_eta),
                range_eta = sqrt(8) / exp(par_list$log_kappa_eta),
                kappa_eta = exp(par_list$log_kappa_eta),
                phi = phi,
                persistence = 1 / (1 - phi))
  for (k in seq_along(iid_levels)) {
    hyper[[correction_iid_terms[[k]]$sd]] <- exp(par_list$log_sigma_iid[k])
  }

  structure(
    list(t0 = t0, T = T, n_years = T - t0,
         mesh = meshes$omega, mesh_xi = meshes$xi,
         fem_xi = data$spde_xi,
         obj = obj, opt = opt,
         max_gradient = max(abs(obj$gr(opt$par))),
         par_list = par_list, tmb_data = data, priors = priors,
         control = control, hyper = hyper,
         mode = par_full[random], blocks = blocks, H = H,
         H_chol = Matrix::Cholesky(H, perm = TRUE, LDL = FALSE, super = TRUE),
         A_latent = cbind(A_omega, A_xi, A_iid),
         precision_obs = 1 / data$v, m_ref = train$m,
         iid_levels = iid_levels, n_obs = n_obs,
         timings = c(optimise = as.numeric(optimise_time - start_time,
                                           units = "secs"),
                     total = as.numeric(Sys.time() - start_time,
                                        units = "secs"))),
    class = "correction_fit"
  )
}

# The final model: stage A, then stage B (PQL) from it. Arguments as for
# fit_stage_a(); pql_args go to fit_correction_pql()
fit_correction <- function(train, t0, T = max(train$year), meshes,
                           priors = correction_priors(),
                           control = list(eval.max = 1000, iter.max = 500),
                           pql_args = list()) {
  fit_a <- fit_stage_a(train[setdiff(names(train), c("z", "v"))], t0 = t0,
                       T = T, meshes = meshes, priors = priors,
                       control = control)
  if (fit_a$opt$convergence != 0) {
    stop("stage A did not converge: ", fit_a$opt$message)
  }
  do.call(fit_correction_pql, c(list(fit_a, train), pql_args))
}

# One row summarising a fit, from fit$hyper and fit$stage_b
fit_summary <- function(fit) {
  b <- fit$stage_b
  get_b <- function(name) if (is.null(b[[name]])) NA else b[[name]]
  hyper <- fit$hyper[c("sigma_omega", "range_omega", "sigma_eta", "range_eta",
                       "phi", "persistence",
                       vapply(correction_iid_terms, `[[`, "", "sd"))]
  tibble::as_tibble(c(
    list(n_train = fit$n_obs),
    lapply(fit$iid_levels, length) |>
      setNames(paste0("n_levels_", names(fit$iid_levels))),
    list(t0 = fit$t0, T = fit$T, mesh_nodes = fit$mesh$n,
         mesh_cutoff_km = attr(fit$mesh, "cutoff"),
         mesh_xi_nodes = fit$mesh_xi$n),
    hyper,
    list(objective = fit$opt$objective,
         convergence = fit$opt$convergence,
         nlminb_message = fit$opt$message,
         iterations = fit$opt$iterations,
         max_gradient = fit$max_gradient,
         stage_a_objective = get_b("objective_a"),
         pql_passes = get_b("passes_first"),
         pql_converged = get_b("converged_first"),
         pql_damped = get_b("damped_first"),
         pql_rms_move = get_b("rms_move"),
         pql_max_move = get_b("max_move"),
         pql_refit = get_b("refit"),
         pql_passes_second = get_b("passes_second"),
         pql_converged_second = get_b("converged_second"),
         pql_rms_move_second = get_b("rms_move_second"),
         pql_n_clamped = get_b("n_clamped"),
         time_stage_a_s = get_b("time_stage_a"),
         time_pql_s = get_b("time_pql"))
  ))
}


# prediction ------------------------------------------------------------------------

# Cut-posterior shift of the latent mode when the offset at the training
# observations changes from m_ref to m_new: the mode solves
# H theta = A' D (z - m), so it moves by H^-1 A' D (m_ref - m_new). m_new may
# have one column per dynamical draw; the factor of H is reused for all
correction_mode_shift <- function(fit, m_new) {
  rhs <- Matrix::crossprod(fit$A_latent,
                           fit$precision_obs * (fit$m_ref - as.matrix(m_new)))
  as.matrix(Matrix::solve(fit$H_chol, rhs, system = "A"))
}

# n draws from N(0, H^-1) from the permuted Cholesky factor P H P' = L L':
# P' L^-T e has covariance H^-1
sample_latent_deviation <- function(H_chol, n_latent, n) {
  e <- matrix(rnorm(n_latent * n), n_latent, n)
  as.matrix(Matrix::solve(H_chol, Matrix::solve(H_chol, e, system = "Lt"),
                          system = "Pt"))
}

# (a) Joint latent draws and the node fields per year; see the API above.
# Years at or before t0 have xi = 0, years in (t0, T] read the fitted x, and
# later years run eta forward from eta_T = x_T - x_{T-1} by the AR(1) with
# fresh Matern innovations (or their mean, 0), accumulating into xi.
correction_node_draws <- function(fit, years, n_draws, m_draws_train = NULL,
                                  mean = FALSE) {
  if (mean) {
    theta <- matrix(fit$mode, ncol = 1)
  } else {
    theta <- sample_latent_deviation(fit$H_chol, length(fit$mode), n_draws) +
      fit$mode
    if (!is.null(m_draws_train)) {
      stopifnot(nrow(m_draws_train) == n_draws,
                ncol(m_draws_train) == fit$n_obs)
      theta <- theta + correction_mode_shift(fit, t(m_draws_train))
    }
  }
  n <- ncol(theta)
  n_xi <- fit$mesh_xi$n
  x <- theta[fit$blocks$x, , drop = FALSE]
  x_year <- function(year) {
    if (year <= fit$t0) return(matrix(0, n_xi, n))
    x[(year - fit$t0 - 1) * n_xi + seq_len(n_xi), , drop = FALSE]
  }

  xi <- list()
  for (year in unique(years[years <= fit$T])) {
    xi[[as.character(year)]] <- x_year(year)
  }
  max_horizon <- max(c(0, years - fit$T))
  if (max_horizon > 0) {
    phi <- fit$hyper$phi
    if (!mean) {
      Q_eta <- matern_precision_r(fit$fem_xi, fit$hyper$kappa_eta,
                                  fit$hyper$sigma_eta)
      Q_eta_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta),
                                     perm = TRUE, LDL = FALSE, super = TRUE)
    }
    xi_t <- x_year(fit$T)
    eta <- xi_t - x_year(fit$T - 1)
    for (h in seq_len(max_horizon)) {
      eta <- phi * eta
      if (!mean) {
        eta <- eta + sqrt(1 - phi ^ 2) *
          sample_latent_deviation(Q_eta_chol, n_xi, n)
      }
      xi_t <- xi_t + eta
      if ((fit$T + h) %in% years) xi[[as.character(fit$T + h)]] <- xi_t
    }
  }
  list(theta = theta, omega = theta[fit$blocks$w_omega, , drop = FALSE],
       xi = xi)
}

# (b) The correction omega + xi at the rows of `new`, points x draws, plus
# fresh noise draws if noise = TRUE; see the API above
project_correction <- function(fit, fields, new, noise = FALSE) {
  coords <- coords_km(new)
  A_omega <- mesh_basis(fit$mesh, coords)
  outside <- Matrix::rowSums(A_omega) < 0.5
  if (any(outside)) {
    warning(sum(outside), " prediction points lie outside the mesh; ",
            "their spatial fields are set to zero")
  }
  out <- as.matrix(A_omega %*% fields$omega)
  A_xi <- mesh_basis(fit$mesh_xi, coords)
  for (year in unique(new$year[new$year > fit$t0])) {
    rows <- which(new$year == year)
    out[rows, ] <- out[rows, ] + as.matrix(
      A_xi[rows, , drop = FALSE] %*% fields$xi[[as.character(year)]])
  }
  if (noise) {
    out <- out + iid_noise(fit, new, ncol(out))
  }
  out
}

# Fresh draws of the iid terms (correction_iid_terms) at the rows of `new`
# (cell and year), rows x n_draws: one per level of each term's key, shared by
# the rows with that level within a draw (caller's RNG)
iid_noise <- function(fit, new, n_draws) {
  out <- 0
  for (term in correction_iid_terms) {
    key <- term$key(new)
    levels <- unique(key)
    fresh <- matrix(rnorm(length(levels) * n_draws, 0, fit$hyper[[term$sd]]),
                    length(levels))
    out <- out + fresh[match(key, levels), , drop = FALSE]
  }
  out
}

# Draws x points predictive draws of m + omega + xi + fresh u + fresh p at the
# rows of `new` (lon/lat or x_km/y_km, year, cell, and m, used when there are
# no dynamical draws): the logit mortality of new assays, before the
# beta-binomial assay noise. With m_draws_train (K x n_obs) and m_draws_new
# (K x nrow(new)), the paired dynamical draws, draw d uses dynamical draw
# ((d - 1) mod K) + 1: its cut-posterior shift and its m at the new points.
# Generated in batches to bound memory. With map = TRUE, returns
# list(draws, map): map holds the same draws without u and p, m + omega + xi,
# whose inverse logit is what the maps show; draws are unchanged by asking
# for it
predict_correction <- function(fit, new, m_draws_train = NULL,
                               m_draws_new = NULL, n_draws = 1000,
                               batch_size = 100, map = FALSE) {
  stopifnot(is.null(m_draws_train) == is.null(m_draws_new),
            is.null(m_draws_new) || ncol(m_draws_new) == nrow(new))
  draws <- matrix(NA_real_, n_draws, nrow(new))
  map_draws <- if (map) draws else NULL
  for (batch in split(seq_len(n_draws), ceiling(seq_len(n_draws) /
                                                 batch_size))) {
    if (is.null(m_draws_train)) {
      fields <- correction_node_draws(fit, new$year, length(batch))
      m_new <- new$m
    } else {
      k <- ((batch - 1) %% nrow(m_draws_train)) + 1
      fields <- correction_node_draws(fit, new$year, length(batch),
                                      m_draws_train[k, , drop = FALSE])
      m_new <- t(m_draws_new[k, , drop = FALSE])
    }
    # the noise is added to the correction before m, in the order
    # project_correction(noise = TRUE) adds it, so draws match it exactly
    correction <- project_correction(fit, fields, new)
    if (map) map_draws[batch, ] <- t(m_new + correction)
    correction <- correction + iid_noise(fit, new, length(batch))
    draws[batch, ] <- t(m_new + correction)
  }
  if (map) list(draws = draws, map = map_draws) else draws
}
