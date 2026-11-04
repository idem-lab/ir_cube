# Error of the one-step cut-posterior shift on the final fit (#69). The second
# stage is fitted once, at the mechanistic model's posterior mean logit at the
# training bioassays (m_ref). Each mechanistic draw's prediction m_d moves the
# latent mode by one Newton step from it (correction_mode_shift(),
# R/two_stage_correction.R), instead of a refit. This compares that step with
# an exact refit of the latent mode at m_d at the same hyperparameters (theta
# at its mode), as check 3 of R/check_two_stage.R does on simulated data.
#
#   Rscript R/check_mode_shift.R
#
# Run after R/two_stage_maps.R, which saves each type's fit and paired
# mechanistic draws. For each type, n_check of the 2,000 draws, evenly spaced;
# the error is reported on the latent correction at the training bioassays
# (A v, logit scale), against the size of the shift itself.

source("R/packages.R")
source("R/two_stage_predictions.R")

n_check <- 4
check_types <- c("Deltamethrin", "Permethrin", "DDT", "Bendiocarb",
                 "Pirimiphos-methyl")

check_type <- function(type) {
  fit <- readRDS(file.path(two_stage_maps_dir, type, "fit.rds"))
  dynamical <- readRDS(file.path(two_stage_maps_dir, type, "dynamical.rds"))
  draws <- round(seq(1, nrow(dynamical$logit_train), length.out = n_check))
  m_new <- dynamical$logit_train[draws, , drop = FALSE]
  shifted <- fit$mode + correction_mode_shift(fit, t(m_new))
  bind_rows(lapply(seq_along(draws), function(i) {
    data <- fit$tmb_data
    data$m <- m_new[i, ]
    obj <- correction_adfun(data, fit$par_list)
    invisible(obj$fn(fit$opt$par))
    refit <- obj$env$last.par[obj$env$random]
    at_data <- function(v) as.vector(fit$A_latent %*% v)
    shift <- at_data(refit - fit$mode)
    error <- at_data(shifted[, i] - refit)
    tibble(insecticide_type = type, draw = draws[i],
           rms_m_change = sqrt(mean((m_new[i, ] - fit$m_ref) ^ 2)),
           rms_shift = sqrt(mean(shift ^ 2)),
           rms_error = sqrt(mean(error ^ 2)),
           max_abs_error = max(abs(error)),
           error_over_shift = rms_error / rms_shift)
  }))
}

results <- bind_rows(lapply(check_types, function(type) {
  report("%s", type)
  check_type(type)
}))
cat("\nOne-step shift against an exact refit, on the correction at the",
    "training bioassays (logit scale): rms_m_change, the draw's change in m;",
    "rms_shift, the refit's change in the correction; rms_error and",
    "max_abs_error, the one-step shift's error\n")
print(results, n = Inf, width = Inf)
cat("\nover all draws: rms error", signif(sqrt(mean(results$rms_error ^ 2)), 3),
    "against rms shift", signif(sqrt(mean(results$rms_shift ^ 2)), 3),
    "; max error", signif(max(results$max_abs_error), 3), "\n")
