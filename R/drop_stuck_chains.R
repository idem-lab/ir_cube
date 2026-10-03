# Leave out MCMC chains that stopped moving (stuck_chains(), in
# R/dynamical_predictions.R) from a fitted fold or the full fit.
#
#   Rscript R/drop_stuck_chains.R <file> [<file> ...]
#
# A cross-validation fold (outputs/cv_draws/dynamical__*.rds): the rows of the
# stored predictions (p_draws, rho_type_draws, p_draws_before) that came from a
# stuck chain are dropped, and the chains recorded in `stuck_chains`, so that
# paired_draw_index() gives the draws the remaining rows correspond to. `draws`
# keeps every chain. The convergence diagnostics are recomputed without them.
#
# The full fit (temporary/fitted_model.RData): `draws` loses the stuck chains,
# in the raw free-state draws too, so calculate() and the plain-R predictions
# both leave them out; the original is kept as `draws_all_chains`, and the
# chains in `stuck_chains`.
#
# DROP_CHAINS="1,2" names the chains to drop instead, e.g. chains trapped in a
# minor mode of the posterior (lower log posterior), which are not stuck but
# disagree with the others.
#
# A file already processed, or with no stuck chain, is left as it is. Run with
# the greta library (R/packages.R), which the saved objects need.

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
})
source("R/dynamical_predictions.R")

chains_to_drop <- function(draws) {
  named <- Sys.getenv("DROP_CHAINS")
  if (nzchar(named)) {
    as.integer(strsplit(named, ",")[[1]])
  } else {
    stuck_chains(draws)
  }
}

for (file in commandArgs(trailingOnly = TRUE)) {

  if (grepl("\\.RData$", file)) {
    fit <- new.env()
    load(file, envir = fit)
    if (exists("stuck_chains", envir = fit, inherits = FALSE)) {
      cat(file, ": already processed, chains", fit$stuck_chains, "dropped\n")
      next
    }
    stuck <- chains_to_drop(fit$draws)
    if (length(stuck) == 0) {
      cat(file, ": no stuck chains\n")
      next
    }
    fit$draws_all_chains <- fit$draws
    fit$draws <- drop_chains(fit$draws, setdiff(seq_along(fit$draws), stuck))
    fit$stuck_chains <- stuck
    temporary_file <- paste0(file, ".tmp")
    save(list = ls(fit, all.names = TRUE), envir = fit,
         file = temporary_file)
    file.rename(temporary_file, file)
    cat(file, ": dropped chain(s)", stuck, "of", length(fit$draws_all_chains),
        "\n")
    next
  }

  fold <- readRDS(file)
  if (!is.null(fold$stuck_chains)) {
    cat(file, ": already processed, chains", fold$stuck_chains, "dropped\n")
    next
  }
  stuck <- chains_to_drop(fold$draws)
  if (length(stuck) == 0) {
    cat(file, ": no stuck chains\n")
    next
  }

  # the stored rows' draws, by the rule fit_fold() thinned them with
  stored_index <- paired_draw_index(fold, maximum = nrow(fold$p_draws))
  keep <- !draw_chain(fold$draws)[stored_index] %in% stuck
  fold$p_draws <- fold$p_draws[keep, , drop = FALSE]
  fold$rho_type_draws <- fold$rho_type_draws[keep, , drop = FALSE]
  if (!is.null(fold$p_draws_before)) {
    fold$p_draws_before <- fold$p_draws_before[keep, , drop = FALSE]
  }
  fold$stuck_chains <- stuck
  stopifnot(identical(stored_index[keep],
                      paired_draw_index(fold, maximum = length(keep))))

  kept_draws <- drop_chains(fold$draws, setdiff(seq_along(fold$draws), stuck))
  fold$convergence <- coda::gelman.diag(kept_draws, multivariate = FALSE,
                                        autoburnin = FALSE)$psrf
  fold$ess <- coda::effectiveSize(kept_draws)

  temporary_file <- paste0(file, ".tmp")
  saveRDS(fold, temporary_file)
  file.rename(temporary_file, file)
  cat(sprintf("%s: dropped chain(s) %s of %i; %i of %i stored draws kept\n",
              file, paste(stuck, collapse = ", "), length(fold$draws),
              sum(keep), length(keep)))
}
