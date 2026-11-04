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
# mechanistic draws. For each type, n_check of the 2,000 draws, evenly spaced.
# The refits are the slow step: each type's are saved to
# temporary/mode_shift_check.rds as they finish, and reused on a rerun, so the
# summaries can change without refitting.
#
# Everything is on the logit scale at the training bioassays: the correction
# A v, and the prediction lambda = m + A v. The error is the one-step shift
# minus the exact one. It is compared with the size of the shift, and with the
# posterior SD of lambda at the fitted hyperparameters: over all the
# mechanistic draws (each with its one-step shift), plus the latent spread
# around the shifted mode, diag(A H^-1 A'). Its effect on calibration is split
# into
#   - bias: the error's mean over the checked draws, per bioassay. The one-step
#     shifts average exactly zero over all draws (they are linear in
#     m_ref - m_d); the exact shifts need not. A bias moves the intervals off
#     centre, which can only lower coverage;
#   - scale: the slope of the error on the exact shift, and the ratio of the
#     spread of lambda about its reference value under the one-step shift to
#     that under the exact one. A negative slope (ratio above 1) is an
#     undershoot: the correction cancels less of m's change than it should, so
#     lambda's draws spread more and intervals widen; a positive slope (ratio
#     below 1) narrows them.

source("R/packages.R")
source("R/two_stage_predictions.R")

n_check <- 10
check_types <- c("Deltamethrin", "Permethrin", "DDT", "Bendiocarb",
                 "Pirimiphos-methyl")
saved_file <- "temporary/mode_shift_check.rds"

# diag(A H^-1 A'): with P H P' = L L', a' H^-1 a = |L^-1 P a|^2; in chunks of
# bioassays to bound memory
latent_variance_at_data <- function(fit, chunk_size = 500) {
  At <- Matrix::t(fit$A_latent)
  chunks <- split(seq_len(ncol(At)), ceiling(seq_len(ncol(At)) / chunk_size))
  unlist(lapply(chunks, function(columns) {
    b <- Matrix::solve(fit$H_chol, as.matrix(At[, columns]), system = "P")
    Matrix::colSums(Matrix::solve(fit$H_chol, b, system = "L") ^ 2)
  }), use.names = FALSE)
}

# posterior SD of lambda = m + A v at the training bioassays, at the fitted
# hyperparameters: the variance over all mechanistic draws of m plus its
# one-step shift, plus the latent variance around the shifted mode
lambda_sd_at_data <- function(fit, m_draws, chunk_size = 250) {
  chunks <- split(seq_len(nrow(m_draws)),
                  ceiling(seq_len(nrow(m_draws)) / chunk_size))
  sums <- Reduce(`+`, lapply(chunks, function(rows) {
    m <- t(m_draws[rows, , drop = FALSE])
    lambda <- m + as.matrix(fit$A_latent %*%
                              correction_mode_shift(fit, m))
    cbind(rowSums(lambda), rowSums(lambda ^ 2))
  }))
  n <- nrow(m_draws)
  variance_draws <- (sums[, 2] - sums[, 1] ^ 2 / n) / (n - 1)
  sqrt(variance_draws + latent_variance_at_data(fit))
}

# the checked draws' one-step and exactly refitted latent modes, the
# correction's change under each at the bioassays, and lambda's SD there
check_type <- function(type) {
  fit <- readRDS(file.path(two_stage_maps_dir, type, "fit.rds"))
  dynamical <- readRDS(file.path(two_stage_maps_dir, type, "dynamical.rds"))
  draws <- round(seq(1, nrow(dynamical$logit_train), length.out = n_check))
  m_new <- dynamical$logit_train[draws, , drop = FALSE]
  one_step <- fit$mode + correction_mode_shift(fit, t(m_new))
  refit <- sapply(seq_along(draws), function(i) {
    data <- fit$tmb_data
    data$m <- m_new[i, ]
    obj <- correction_adfun(data, fit$par_list)
    invisible(obj$fn(fit$opt$par))
    obj$env$last.par[obj$env$random]
  })
  at_data <- function(v) as.matrix(fit$A_latent %*% (v - fit$mode))
  list(draws = draws,
       m_change = t(m_new) - fit$m_ref,
       one_step = one_step,
       refit = refit,
       shift_one_step = at_data(one_step),
       shift_exact = at_data(refit),
       sd_lambda = lambda_sd_at_data(fit, dynamical$logit_train))
}

saved <- if (file.exists(saved_file)) readRDS(saved_file) else list()
for (type in setdiff(check_types, names(saved))) {
  report("%s: %i refits", type, n_check)
  saved[[type]] <- check_type(type)
  saveRDS(saved, saved_file)
}

# per draw: the error's size, against the shift and against lambda's SD, and
# the scale of the one-step shift against the exact one
per_draw <- bind_rows(lapply(check_types, function(type) {
  s <- saved[[type]]
  bind_rows(lapply(seq_along(s$draws), function(i) {
    shift <- s$shift_exact[, i]
    error <- s$shift_one_step[, i] - shift
    lambda_change <- s$m_change[, i] + shift
    tibble(insecticide_type = type, draw = s$draws[i],
           rms_m_change = sqrt(mean(s$m_change[, i] ^ 2)),
           rms_shift = sqrt(mean(shift ^ 2)),
           rms_error = sqrt(mean(error ^ 2)),
           max_abs_error = max(abs(error)),
           rms_error_over_sd = sqrt(mean((error / s$sd_lambda) ^ 2)),
           max_error_over_sd = max(abs(error) / s$sd_lambda),
           slope = sum(error * shift) / sum(shift ^ 2),
           spread_ratio = sqrt(mean((lambda_change + error) ^ 2) /
                                 mean(lambda_change ^ 2)))
  }))
}))

# per type: lambda's SD, the error over it pooled over draws, and the bias:
# the error's mean over the checked draws per bioassay, over lambda's SD. With
# n_check draws its noise alone gives an rms of about the pooled rms error
# over sqrt(n_check), shown as noise_floor
per_type <- bind_rows(lapply(check_types, function(type) {
  s <- saved[[type]]
  error_over_sd <- (s$shift_one_step - s$shift_exact) / s$sd_lambda
  bias_over_sd <- rowMeans(error_over_sd)
  tibble(insecticide_type = type, bioassays = length(s$sd_lambda),
         median_sd_lambda = median(s$sd_lambda),
         rms_error_over_sd = sqrt(mean(error_over_sd ^ 2)),
         share_over_tenth_sd = mean(abs(error_over_sd) > 0.1),
         mean_bias_over_sd = mean(bias_over_sd),
         rms_bias_over_sd = sqrt(mean(bias_over_sd ^ 2)),
         noise_floor = sqrt(mean(error_over_sd ^ 2) / ncol(error_over_sd)),
         max_bias_over_sd = max(abs(bias_over_sd)),
         slope = sum((s$shift_one_step - s$shift_exact) * s$shift_exact) /
           sum(s$shift_exact ^ 2),
         spread_ratio = sqrt(mean((s$m_change + s$shift_one_step) ^ 2) /
                               mean((s$m_change + s$shift_exact) ^ 2)))
}))

cat("\nOne-step shift against an exact refit, per draw (logit scale, at the",
    "training bioassays): rms_m_change, the draw's change in m; rms_shift, the",
    "exact change in the correction; rms_error and max_abs_error, the one-step",
    "shift's error; rms_error_over_sd and max_error_over_sd, the error per",
    "bioassay over lambda's posterior SD there; slope, the error regressed on",
    "the exact shift (negative: undershoot); spread_ratio, the rms change in",
    "lambda under the one-step shift over that under the exact one\n")
print(per_draw, n = Inf, width = Inf)
cat("\nPer type, pooled over the draws: median_sd_lambda, lambda's median",
    "posterior SD; rms_error_over_sd; share_over_tenth_sd, the share of",
    "bioassay-draws with an error over a tenth of the SD; the bias (mean",
    "error over the draws per bioassay) over the SD: its mean, rms and",
    "maximum, with noise_floor the rms the draws' noise alone would give;",
    "slope and spread_ratio as above\n")
print(per_type, width = Inf)
