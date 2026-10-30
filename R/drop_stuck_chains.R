# Leave out MCMC chains that stopped moving (stuck_chains(), in
# R/chain_checks.R) from the full fit. A fold's are dropped by
# R/run_one_fold.R before it is saved (drop_fold_chains()).
#
#   Rscript R/drop_stuck_chains.R <file> [<file> ...]
#
# The full fit (temporary/fitted_model.RData): `draws` loses the stuck chains,
# in the raw free-state draws too, so calculate() and the plain-R predictions
# both leave them out; the original is kept as `draws_all_chains`, and the
# chains in `stuck_chains`.
#
# A fit already processed, or with no stuck chain, is left as it is. Run with
# the greta library (R/packages.R), which the saved objects need.

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
})
source("R/dynamical_predictions.R")

for (file in commandArgs(trailingOnly = TRUE)) {
  fit <- new.env()
  load(file, envir = fit)
  # the chains dropped, if processed; the image also holds the function
  # stuck_chains() (R/chain_checks.R), which fit_model.R saves with it
  if (is.numeric(get0("stuck_chains", envir = fit, inherits = FALSE))) {
    cat(file, ": already processed, chains", fit$stuck_chains, "dropped\n")
    next
  }
  stuck <- stuck_chains(fit$draws)
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
}
