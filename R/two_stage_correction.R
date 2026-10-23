# Second-stage geostatistical correction to the dynamical model (#21).
#
# A latent Gaussian model for the bioassay counts, beta-binomial on the
# dynamical model's logit prediction plus spatial and noise terms,
# conditioning on the dynamical model without ever updating it (a cut
# posterior):
#
#   lambda_i = m_i + omega(s_i) + xi(s_i, t_i) + u[pixel-year(i)] + p[pixel(i)]
#              + sy[study-year(i)]
#
# m is the dynamical model's posterior mean logit at the training assays
# (m_ref). omega is a Matern field correcting the initial conditions; xi
# accumulates an AR(1)-in-time, Matern-in-space field of annual anomalies eta
# from xi(., t0) = 0, so that forecasts of the correction plateau; u (per
# pixel-year), p (per pixel) and sy (per study-year) are iid.
#
# The prediction target (what the maps show) is m + omega + xi. u, p and sy
# are observation-level noise, shared by the assays of a pixel-year, a pixel
# or a study in a year: they are fitted so that this noise stays out of omega
# and xi, and enter predictions of new assays only as fresh draws from
# N(0, tau^2), N(0, sigma_p^2) and N(0, sigma_sy^2), never as their fitted
# values, even at pixel-years, pixels or study-years with training data.
#
# Fitting, per insecticide type (fit_correction(), #34):
#   Laplace  the hyperparameters theta at the penalised maximum of TMB's
#            Laplace approximation to the marginal likelihood of the
#            beta-binomial counts (tmb/two_stage_correction.cpp), started from
#            the priors' medians;
#   CCD      the posterior of theta by INLA's central composite design
#            (ccd_points()); predictions mix the conditional latent posteriors
#            N(mode(theta_j), H(theta_j)^-1) over its points theta_j.
#
# API (source this file from the repo root):
#
#   meshes <- build_correction_meshes(coords_km, prediction_mask_coords())
#                                                  list(omega, xi): the final
#                                                  model's meshes, covering the
#                                                  sites and the mask
#   fit <- fit_correction(train, t0, T, meshes)    Laplace, then CCD.
#     train: lon, lat (or x_km, y_km), year, cell, study, m (= m_ref), died,
#     mosquito_number, rho. T: last year xi is represented (forecast beyond).
#     Stops if the optimiser does not converge. fit$hyper: sigma_omega,
#     range_omega, sigma_eta, range_eta, phi, persistence, tau (SD of u),
#     sigma_p, sigma_sy (and the kappas) at theta_hat; fit$ccd: the CCD points
#     (theta), their weights and the hyperparameters at each (hyper).
#   fit_summary(fit)                               one-row tibble for tables
#   fields <- correction_node_draws(fit, years, n_draws, m_draws_train = NULL)
#     (a) the node fields of joint latent draws mixed over the CCD points,
#     with the cut-posterior shift per dynamical draw when m_draws_train
#     (n_draws x n_obs) is given: omega (nodes x draws) and xi[[year]] (xi
#     nodes x draws) for every requested year, AR(1)-forecast beyond T, and
#     point, each draw's CCD point. The projection (b) is linear, so the
#     posterior mean of the correction is the projection of the draws' mean.
#   project_correction(fit, fields, new, noise = FALSE)
#     (b) points x draws correction omega + xi at the rows of `new`
#     (lon/lat or x_km/y_km, year, and cell and study if noise = TRUE).
#     noise = TRUE adds fresh u per distinct pixel-year, p per distinct pixel
#     and sy per distinct study-year of `new` (shared by its rows within a
#     draw; caller's RNG), as for new assays; maps use noise = FALSE. Warns for
#     points outside the mesh (their fields are 0). One set of fields can be
#     projected to several sets of points.
#   predict_correction(fit, new, m_draws_train, m_draws_new, n_draws)
#     draws x points m_draw + omega + xi + fresh noise, (a) then (b): the
#     predictive distribution of new assays' logit mortality, mixed over the
#     CCD points. noise = FALSE leaves out the noise: the map, m + omega + xi.
#   Saving: saveRDS(fit, file). Everything above works on the reloaded fit
#     (it keeps H_chol, meshes, levels, tmb_data and par_list); the TMB
#     object the mixture over theta needs is rebuilt once per session
#     (correction_obj()).
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
# africa_equal_area_crs, project_km() and prediction_mask_coords(), shared with
# the latent smooths of the dynamical model
source("R/projection.R")

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

# The iid noise terms, each as the key grouping assays into its levels, and
# the name of its SD in fit$hyper. Adding or removing a term here is all the
# model and the predictions need. sy is noise shared by the bioassays of a
# study in a year (the data's study label, #54), so that survey noise is not
# fitted by xi and mapped (#55). A point with no study (a new bioassay at a
# map cell) is its own study, so its draw is one per cell-year
correction_iid_terms <- list(
  u = list(key = function(d) paste(d$cell, d$year), sd = "tau"),
  p = list(key = function(d) as.character(d$cell), sd = "sigma_p"),
  sy = list(key = function(d) {
    paste(if (is.null(d$study)) paste("cell", d$cell) else d$study, d$year)
  }, sd = "sigma_sy")
)


# coordinates and meshes ----------------------------------------------------

# projected coordinates of a data frame: x_km and y_km if present (simulated
# data), otherwise lon and lat projected
coords_km <- function(df) {
  if (all(c("x_km", "y_km") %in% names(df))) {
    cbind(x_km = df$x_km, y_km = df$y_km)
  } else {
    project_km(df$lon, df$lat)
  }
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

# Priors (#34). The second stage only fills gaps the dynamical model leaves:
# finer-scale level error, and gradual changes of its trend.
#
# Ranges of omega and eta: normal on log range, with 95% of the mass in
# 50-1,500 km (median 274 km; 99.7% below 3,000 km), the form of R-INLA's
# non-PC SPDE priors (Lindgren & Rue 2015). The dynamical model's own fields
# have 1,500 km ranges, so continent-wide corrections are its job. (PC priors,
# Fuglstad et al. 2019, favour long ranges, against this intent.)
#
# AR(1) coefficient phi of the annual anomalies eta: normal on logit phi, with
# a 99% interval for the persistence L = 1 / (1 - phi) = 1 + exp(logit phi) of
# 3-10 years (median 5.2 years). Resistance does not jump up and down from
# year to year, and longer persistences are neither identified nor needed
# (omega holds static offsets).
#
# Half-normal priors on the SDs (Gelman 2006), so the cost of a large SD grows
# quadratically: P(sigma_omega > 1) = 0.05 (the level correction);
# P(sigma_eta > 0.05) = 0.05, so the trend correction xi has a prior SD of
# about 0.2 logit per decade at the median persistence, up to about 0.8 logit
# over 30 years at the 95% point of sigma_eta; and P(sd > 2) = 0.05 for the
# iid noise terms, permissive so that noise is never what forces signal into
# the structured fields.
correction_priors <- function(range_interval = c(50, 1500),
                              persistence_interval = c(3, 10),
                              sigma_omega = 1, sigma_eta = 0.05,
                              sigma_iid = 2, alpha_sigma = 0.05) {
  # a central interval, at `level`, of a normal prior
  normal_prior <- function(interval, level) {
    c(mean(interval), diff(interval) / (2 * qnorm((1 + level) / 2)))
  }
  half_normal_scale <- function(sigma0) sigma0 / qnorm(1 - alpha_sigma / 2)
  list(range_prior = normal_prior(log(range_interval), 0.95),
       phi_prior = normal_prior(log(persistence_interval - 1), 0.99),
       sd_scale = half_normal_scale(c(sigma_omega, sigma_eta, sigma_iid)))
}


# fitting ----------------------------------------------------------------------

# the hyperparameters, on the optimiser's scale, in the template's order
correction_hyper_names <- c("log_sigma_omega", "log_range_omega",
                            "log_sigma_eta", "log_range_eta", "logit_phi",
                            "log_sigma_iid")

# the TMB object, with the latent effects integrated out by TMB's Laplace
# approximation
correction_adfun <- function(data, parameters, silent = TRUE) {
  TMB::MakeADFun(data = data, parameters = parameters,
                 random = c("w_omega", "x", "iid"),
                 DLL = load_correction_template(), silent = silent)
}

# The TMB object of a fit's Laplace-approximate marginal posterior of theta,
# built once per session (it is not saved with the fit, and building it
# reloads the template), started at the fit's latent mode and theta
correction_obj <- local({
  cache <- new.env()
  function(fit) {
    if (is.null(cache[[fit$id]])) {
      cache[[fit$id]] <- correction_adfun(fit$tmb_data, fit$par_list)
    }
    cache[[fit$id]]
  }
})

# W, the negative second derivative in lambda of the beta-binomial log
# likelihood of `died` of n with mean ilogit(lambda) and intra-class
# correlation rho (a = c p, b = c (1 - p), c = (1 - rho) / rho)
beta_binomial_weight <- function(lambda, died, n, rho) {
  p <- plogis(lambda)
  s <- p * (1 - p)
  c <- (1 - rho) / rho
  a <- c * p
  b <- c * (1 - p)
  d1 <- c * (digamma(died + a) - digamma(n - died + b) - digamma(a) +
               digamma(b))
  d2 <- c ^ 2 * (trigamma(died + a) + trigamma(n - died + b) - trigamma(a) -
                   trigamma(b))
  -(d2 * s ^ 2 + d1 * s * (1 - 2 * p))
}

# H, the Hessian of the negative log joint density in the latent effects, from
# the TMB object obj at the full parameter vector par (latent effects and
# theta)
latent_hessian <- function(obj, par) {
  H <- obj$env$spHess(par, random = TRUE)
  Matrix::forceSymmetric(as(H, "CsparseMatrix"), uplo = "L")
}

# The conditional fit at theta (optimiser scale): the latent mode solved at
# theta, the factor of its Hessian H (updated from fit's symbolic
# factorisation when it has one: the sparsity pattern does not depend on
# theta), the hyperparameters, and W at the mode as precision_obs, for the
# linearised cut-posterior shift
conditional_fit <- function(fit, obj, theta) {
  obj$fn(theta)
  par_full <- obj$env$last.par
  H <- latent_hessian(obj, par_full)
  fit$par_list <- obj$env$parList(par = par_full)
  fit$hyper <- correction_hyper(fit$par_list)
  fit$mode <- par_full[obj$env$random]
  fit$H_chol <- if (is.null(fit$H_chol)) {
    Matrix::Cholesky(H, perm = TRUE, LDL = FALSE, super = TRUE)
  } else {
    Matrix::update(fit$H_chol, H)
  }
  data <- obj$env$data
  lambda <- fit$m_ref + as.vector(fit$A_latent %*% fit$mode)
  fit$precision_obs <- beta_binomial_weight(lambda, data$died,
                                            data$n_mosquito, data$rho)
  fit
}

# fit$hyper from the template's parameter list
correction_hyper <- function(par_list) {
  range_omega <- exp(par_list$log_range_omega)
  range_eta <- exp(par_list$log_range_eta)
  hyper <- list(sigma_omega = exp(par_list$log_sigma_omega),
                range_omega = range_omega,
                kappa_omega = sqrt(8) / range_omega,
                sigma_eta = exp(par_list$log_sigma_eta),
                range_eta = range_eta,
                kappa_eta = sqrt(8) / range_eta,
                phi = plogis(par_list$logit_phi),
                persistence = 1 + exp(par_list$logit_phi))
  for (k in seq_along(correction_iid_terms)) {
    hyper[[correction_iid_terms[[k]]$sd]] <- exp(par_list$log_sigma_iid[k])
  }
  hyper
}

# The model's TMB data and starting parameters for train (as fit_correction()),
# and the parts of a fit that do not depend on theta: the meshes, the iid
# levels, the latent blocks and the map from observations to the latent
# vector. The start is the latent field at 0 and theta at the priors' medians
# for the ranges and persistence, with small SDs (the correction should be
# modest if the dynamical model is good)
correction_setup <- function(train, t0, T, meshes, priors) {
  stopifnot(T > t0, all(train$year >= t0 & train$year <= T),
            !is.null(train$study))
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
    list(m = train$m, A_omega = A_omega, A_xi = A_xi,
         iid_index = iid_index - 1L,
         iid_term = rep(seq_along(iid_levels) - 1L, lengths(iid_levels)),
         spde = correction_fem(meshes$omega),
         spde_xi = correction_fem(meshes$xi),
         died = train$died, n_mosquito = train$mosquito_number,
         rho = train$rho),
    priors
  )
  parameters <- list(w_omega = rep(0, meshes$omega$n),
                     x = matrix(0, meshes$xi$n, T - t0),
                     iid = rep(0, sum(lengths(iid_levels))),
                     log_sigma_omega = log(0.5),
                     log_range_omega = priors$range_prior[1],
                     log_sigma_eta = log(0.03),
                     log_range_eta = priors$range_prior[1],
                     logit_phi = priors$phi_prior[1],
                     log_sigma_iid = rep(log(0.3), length(iid_levels)))

  # the latent vector's blocks, in the template's order
  n_x <- meshes$xi$n * (T - t0)
  blocks <- c(list(w_omega = seq_len(meshes$omega$n),
                   x = meshes$omega$n + seq_len(n_x)),
              lapply(seq_along(iid_levels), function(k) {
                meshes$omega$n + n_x + (offset[k] + 1):offset[k + 1]
              }) |> setNames(names(iid_levels)))
  # observations -> latent vector, for the cut-posterior shift
  A_iid <- Matrix::sparseMatrix(i = rep(seq_len(n_obs), ncol(iid_index)),
                                j = as.vector(iid_index), x = 1,
                                dims = c(n_obs, sum(lengths(iid_levels))))
  list(data = data, parameters = parameters,
       fit = structure(
         list(t0 = t0, T = T, mesh = meshes$omega, mesh_xi = meshes$xi,
              fem_xi = data$spde_xi, tmb_data = data,
              blocks = blocks, A_latent = cbind(A_omega, A_xi, A_iid),
              m_ref = train$m, iid_levels = iid_levels, n_obs = n_obs),
         class = "correction_fit"))
}

# The CCD design in d dimensions, as points x d in the standardised z scale:
# the centre, the 2d axial points at distance sqrt(d), and a two-level
# factorial of resolution at least V (full for d <= 4; the half fraction with
# the last factor the product of the others for d = 5-7; for d = 8 the quarter
# fraction with G = ABCD and H = ABEF), all scaled by f0
ccd_design <- function(d, f0 = 1.1) {
  stopifnot(d >= 2, d <= 8)
  base <- if (d <= 4) d else if (d <= 7) d - 1 else d - 2
  full <- as.matrix(expand.grid(rep(list(c(-1, 1)), base)))
  factorial <- if (d <= 4) {
    full
  } else if (d <= 7) {
    cbind(full, apply(full, 1, prod))
  } else {
    cbind(full, full[, 1] * full[, 2] * full[, 3] * full[, 4],
          full[, 1] * full[, 2] * full[, 5] * full[, 6])
  }
  axial <- rbind(diag(sqrt(d), d), diag(-sqrt(d), d))
  unname(rbind(0, axial, factorial) * f0)
}

# The normalised CCD weights of the points of ccd_design(d, f0), given the log
# posterior density at each relative to the mode: the density times the
# design weight, 1 at the centre and exp(d f0^2 / 2) / ((n_p - 1)(f0^2 - 1)) at
# the other n_p - 1 points, which all lie at radius f0 sqrt(d). If the
# posterior is Gaussian, the weighted points have its mean and covariance
ccd_weights <- function(log_density, d, f0 = 1.1) {
  n_points <- length(log_density)
  design_weight <- c(1, rep(exp(d * f0 ^ 2 / 2) /
                              ((n_points - 1) * (f0 ^ 2 - 1)), n_points - 1))
  weight <- design_weight * exp(log_density)
  weight / sum(weight)
}

# INLA's central composite design integration over theta (Rue, Martino &
# Chopin 2009, sec. 6.5), on the optimiser's scale, at the mode theta_hat of
# the marginal posterior (obj, as correction_obj()): theta = theta_hat +
# V S^1/2 z for the eigen-decomposition V S V' of the inverse Hessian there,
# each axis stretched separately on its two sides so that the log density
# falls by 1 at z = +-sqrt(2), as it would for a Gaussian (INLA's skewness
# correction), with ccd_weights(). Returns the points (theta) and their
# weights (summing to 1)
ccd_points <- function(obj, theta_hat, f0 = 1.1) {
  log_post <- function(theta) -obj$fn(theta)
  log_post_hat <- log_post(theta_hat)
  hessian <- optimHess(theta_hat, obj$fn, obj$gr)
  covariance <- solve((hessian + t(hessian)) / 2)
  d <- length(theta_hat)
  eigen_cov <- eigen(covariance, symmetric = TRUE)
  to_theta <- function(z) {
    as.vector(theta_hat + eigen_cov$vectors %*% (sqrt(eigen_cov$values) * z))
  }
  stretch <- matrix(1, d, 2, dimnames = list(NULL, c("negative", "positive")))
  for (i in seq_len(d)) {
    for (side in 1:2) {
      z <- rep(0, d)
      z[i] <- c(-1, 1)[side] * sqrt(2)
      drop <- log_post_hat - log_post(to_theta(z))
      stretch[i, side] <- if (is.finite(drop) && drop > 0) sqrt(1 / drop) else 1
    }
  }
  design <- ccd_design(d, f0)
  z <- design * ifelse(design < 0, rep(stretch[, 1], each = nrow(design)),
                       rep(stretch[, 2], each = nrow(design)))
  theta <- t(apply(z, 1, to_theta))
  colnames(theta) <- names(theta_hat)
  log_density <- apply(theta, 1, log_post) - log_post_hat
  # a failed inner solve would otherwise give NaN weights
  stopifnot(all(is.finite(log_density)))
  list(theta = theta, weight = ccd_weights(log_density, d, f0))
}

# The template's parameter list of fit with theta (optimiser scale) in place
# of its hyperparameters
theta_par_list <- function(fit, theta) {
  par_list <- fit$par_list
  index <- 0
  for (name in correction_hyper_names) {
    size <- length(par_list[[name]])
    par_list[[name]] <- unname(theta[index + seq_len(size)])
    index <- index + size
  }
  par_list
}

# The final model (#34). theta at the penalised maximum of TMB's Laplace
# approximation to the marginal likelihood of the beta-binomial counts,
# started from correction_setup(); then the posterior of theta by CCD
# (ccd_points()), with the hyperparameters at each point in fit$ccd$hyper.
# Stops if the optimiser does not converge, before the CCD's inner solves.
# The fit itself is the conditional fit at theta_hat: latent posterior
# N(mode, H^-1), with precision_obs W there for the linearised cut-posterior
# shift (conditional_fit()). Predictions mix over the CCD points
# (correction_node_draws(), predict_correction()). train, t0, T and meshes as
# in the API above
fit_correction <- function(train, t0, T = max(train$year), meshes,
                           priors = correction_priors()) {
  start_time <- Sys.time()
  elapsed <- function() as.numeric(Sys.time() - start_time, units = "secs")
  setup <- correction_setup(train, t0, T, meshes, priors)
  obj <- correction_adfun(setup$data, setup$parameters)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(eval.max = 1000, iter.max = 500))
  if (opt$convergence != 0) {
    stop("the Laplace fit did not converge: ", opt$message)
  }
  max_gradient <- max(abs(obj$gr(opt$par)))
  time_laplace <- elapsed()
  ccd <- ccd_points(obj, setNames(opt$par, names(obj$par)))
  fit <- conditional_fit(setup$fit, obj, opt$par)
  ccd$hyper <- bind_rows(lapply(seq_len(nrow(ccd$theta)), function(j) {
    as_tibble(correction_hyper(theta_par_list(fit, ccd$theta[j, ])))
  }))
  fit$opt <- opt
  fit$max_gradient <- max_gradient
  fit$ccd <- ccd
  fit$id <- sprintf("%s-%i", format(Sys.time(), "%Y%m%d%H%M%OS6"),
                    Sys.getpid())
  fit$timings <- c(laplace = time_laplace, ccd = elapsed() - time_laplace,
                   total = elapsed())
  fit
}

# One row summarising a fit: the hyperparameters at theta_hat and their CCD
# posterior means (suffix _ccd), and the optimisation's diagnostics
fit_summary <- function(fit) {
  names_hyper <- c("sigma_omega", "range_omega", "sigma_eta", "range_eta",
                   "phi", "persistence",
                   vapply(correction_iid_terms, `[[`, "", "sd"))
  ccd_mean <- colSums(as.matrix(fit$ccd$hyper[names_hyper]) * fit$ccd$weight)
  tibble::as_tibble(c(
    list(n_train = fit$n_obs),
    lapply(fit$iid_levels, length) |>
      setNames(paste0("n_levels_", names(fit$iid_levels))),
    list(t0 = fit$t0, T = fit$T, mesh_nodes = fit$mesh$n,
         mesh_cutoff_km = attr(fit$mesh, "cutoff"),
         mesh_xi_nodes = fit$mesh_xi$n),
    fit$hyper[names_hyper],
    as.list(setNames(ccd_mean, paste0(names_hyper, "_ccd"))),
    list(objective = fit$opt$objective,
         convergence = fit$opt$convergence,
         nlminb_message = fit$opt$message,
         iterations = fit$opt$iterations,
         max_gradient = fit$max_gradient,
         ccd_points = nrow(fit$ccd$theta),
         time_laplace_s = fit$timings[["laplace"]],
         time_ccd_s = fit$timings[["ccd"]])
  ))
}


# prediction ------------------------------------------------------------------------

# Cut-posterior shift of the latent mode when the offset at the training
# observations changes from m_ref to m_new, by one Newton step from the mode:
# to first order the score there changes by A' W (m_ref - m_new), with W the
# beta-binomial weight at the mode (precision_obs), so the mode moves by
# H^-1 A' W (m_ref - m_new). m_new may have one column per dynamical draw; the
# factor of H is reused for all
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

# The node fields omega and xi per year of n_draws joint latent draws at the
# fit's own theta. Years at or before t0 have xi = 0, years in (t0, T] read
# the fitted x, and later years run eta forward from eta_T = x_T - x_{T-1} by
# the AR(1) with fresh Matern innovations, accumulating into xi.
correction_node_draws_at <- function(fit, years, n_draws,
                                     m_draws_train = NULL) {
  latent <- sample_latent_deviation(fit$H_chol, length(fit$mode), n_draws) +
    fit$mode
  if (!is.null(m_draws_train)) {
    stopifnot(nrow(m_draws_train) == n_draws,
              ncol(m_draws_train) == fit$n_obs)
    latent <- latent + correction_mode_shift(fit, t(m_draws_train))
  }
  n_xi <- fit$mesh_xi$n
  x <- latent[fit$blocks$x, , drop = FALSE]
  x_year <- function(year) {
    if (year <= fit$t0) return(matrix(0, n_xi, n_draws))
    x[(year - fit$t0 - 1) * n_xi + seq_len(n_xi), , drop = FALSE]
  }

  xi <- list()
  for (year in unique(years[years <= fit$T])) {
    xi[[as.character(year)]] <- x_year(year)
  }
  max_horizon <- max(c(0, years - fit$T))
  if (max_horizon > 0) {
    phi <- fit$hyper$phi
    Q_eta <- matern_precision_r(fit$fem_xi, fit$hyper$kappa_eta,
                                fit$hyper$sigma_eta)
    Q_eta_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta),
                                   perm = TRUE, LDL = FALSE, super = TRUE)
    xi_t <- x_year(fit$T)
    eta <- xi_t - x_year(fit$T - 1)
    for (h in seq_len(max_horizon)) {
      eta <- phi * eta + sqrt(1 - phi ^ 2) *
        sample_latent_deviation(Q_eta_chol, n_xi, n_draws)
      xi_t <- xi_t + eta
      if ((fit$T + h) %in% years) xi[[as.character(fit$T + h)]] <- xi_t
    }
  }
  list(omega = latent[fit$blocks$w_omega, , drop = FALSE], xi = xi)
}

# The point of p(theta | y) each of n_draws draws uses: the CCD points
# allotted by systematic resampling of their weights, in random order, so that
# no point is paired with a contiguous run of the (MCMC-ordered) dynamical
# draws (caller's RNG)
ccd_allocation <- function(weight, n_draws) {
  u <- (seq_len(n_draws) - runif(1)) / n_draws
  sample(findInterval(u, cumsum(weight) / sum(weight)) + 1)
}

# The conditional fits at the CCD points `points` of fit: calls use(j, fit_j)
# for each. TMB starts each inner solve from the best mode so far
for_ccd_points <- function(fit, points, use) {
  obj <- correction_obj(fit)
  for (j in points) use(j, conditional_fit(fit, obj, fit$ccd$theta[j, ]))
}

# (a) The node fields per year of joint latent draws, mixed over theta; see
# the API above. Draw d uses CCD point j(d), allotted by the weights, and its
# conditional fit (correction_node_draws_at()); the fields carry `point`, the
# point of each draw
correction_node_draws <- function(fit, years, n_draws, m_draws_train = NULL) {
  out <- NULL
  point <- ccd_allocation(fit$ccd$weight, n_draws)
  for_ccd_points(fit, unique(point), function(j, fit_j) {
    draws <- which(point == j)
    fields <- correction_node_draws_at(
      fit_j, years, length(draws),
      if (is.null(m_draws_train)) NULL else
        m_draws_train[draws, , drop = FALSE])
    if (is.null(out)) {
      out <<- list(omega = matrix(0, nrow(fields$omega), n_draws),
                   xi = lapply(fields$xi, function(x) {
                     matrix(0, nrow(x), n_draws)
                   }))
    }
    out$omega[, draws] <<- fields$omega
    for (year in names(fields$xi)) {
      out$xi[[year]][, draws] <<- fields$xi[[year]]
    }
  })
  out$point <- point
  out
}

# (b) The correction omega + xi at the rows of `new`, points x draws, plus
# fresh noise draws if noise = TRUE (each draw with the SDs of its CCD point,
# fields$point, when the fields carry it); see the API above
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
    out <- out + iid_noise(fit, new, ncol(out), fields$point)
  }
  out
}

# Fresh draws of the iid terms (correction_iid_terms, or those named in
# `terms`) at the rows of `new` (cell, year and study), rows x n_draws: one per
# level of each term's key, shared by the rows with that level within a draw
# (caller's RNG). With `point`, the CCD point of each draw, each draw uses that
# point's SDs
iid_noise <- function(fit, new, n_draws, point = NULL,
                      terms = names(correction_iid_terms)) {
  out <- 0
  for (term in correction_iid_terms[terms]) {
    key <- term$key(new)
    levels <- unique(key)
    sd <- if (is.null(point)) fit$hyper[[term$sd]] else
      fit$ccd$hyper[[term$sd]][point]
    fresh <- matrix(rnorm(length(levels) * n_draws), length(levels)) *
      rep(sd, each = length(levels))
    out <- out + fresh[match(key, levels), , drop = FALSE]
  }
  out
}

# Draws x points predictive draws of m + omega + xi + fresh noise (u, p, sy)
# at the rows of `new` (lon/lat or x_km/y_km, year, cell, study, and m, used
# when there are no dynamical draws): the logit mortality of new assays,
# before the beta-binomial assay noise, mixed over theta
# (correction_node_draws(), project_correction()). With m_draws_train
# (K x n_obs) and m_draws_new (K x nrow(new)), the paired dynamical draws,
# draw d uses dynamical draw ((d - 1) mod K) + 1: its cut-posterior shift and
# its m at the new points. With noise = FALSE, the same without the noise:
# m + omega + xi, whose inverse logit is what the maps show
predict_correction <- function(fit, new, m_draws_train = NULL,
                               m_draws_new = NULL, n_draws = 1000,
                               noise = TRUE) {
  stopifnot(is.null(m_draws_train) == is.null(m_draws_new),
            is.null(m_draws_new) || ncol(m_draws_new) == nrow(new))
  m_train <- NULL
  m_new <- new$m
  if (!is.null(m_draws_train)) {
    k <- ((seq_len(n_draws) - 1) %% nrow(m_draws_train)) + 1
    m_train <- m_draws_train[k, , drop = FALSE]
    m_new <- t(m_draws_new[k, , drop = FALSE])
  }
  fields <- correction_node_draws(fit, unique(new$year), n_draws, m_train)
  t(m_new + project_correction(fit, fields, new, noise = noise))
}
