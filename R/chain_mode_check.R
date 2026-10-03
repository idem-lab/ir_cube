# Check a fitted fold or full fit for chains in different posterior modes, and
# for stuck chains, and say which chains to drop (R/drop_stuck_chains.R).
#
#   Rscript R/chain_mode_check.R <name> <file> <output dir> [threshold]
#
# <file> is a saved fold (outputs/cv_draws/dynamical__*.rds) or a full fit
# (temporary/fitted_model.RData). The constrained refit has two modes in the
# mortality floor, near 0.25-0.32 and near 0.001-0.002; where chains found
# both, the low-floor chains had a log posterior about 59-60 higher
# (doc/cv_run_plan.md, "Constrained refit"). So when the chains' mean floors
# differ by more than `threshold` (default 0.05), the chains whose mean floor
# is more than `threshold` above the lowest are marked to drop. Stuck chains
# (under half their draws distinct, as stuck_chains()) are marked too.
#
# Writes to <output dir>:
#   <name>.csv         one row: chains, per-chain mean floor, chains to drop,
#                      and rank Rhat, bulk and tail ESS (posterior) over every
#                      parameter, with all chains and with the kept chains
#   <name>_spread.csv  the parameters whose per-chain means disagree most,
#                      (max - min of the chain means) / mean within-chain sd
#   <name>.drop        the chains to drop, comma-separated (empty if none),
#                      for DROP_CHAINS
# A full fit already processed by drop_stuck_chains.R is checked on all its
# chains (draws_all_chains) and its recorded drop is reported.

suppressMessages({
  library(coda)
  library(posterior)
})

arguments <- commandArgs(trailingOnly = TRUE)
name <- arguments[1]
file <- arguments[2]
output_dir <- arguments[3]
threshold <- if (length(arguments) > 3) as.numeric(arguments[4]) else 0.05
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

already_dropped <- integer(0)
if (grepl("\\.RData$", file)) {
  fit <- new.env()
  load(file, envir = fit)
  draws <- if (exists("draws_all_chains", envir = fit, inherits = FALSE)) {
    already_dropped <- fit$stuck_chains
    fit$draws_all_chains
  } else {
    fit$draws
  }
  rm(fit)
} else {
  fold <- readRDS(file)
  draws <- fold$draws
  already_dropped <- as.integer(fold$stuck_chains)
  rm(fold)
}
invisible(gc())

chains <- lapply(draws, as.matrix)
rm(draws)
invisible(gc())
n_chains <- length(chains)
parameters <- colnames(chains[[1]])

floor_column <- grep("^mortality_floor", parameters, value = TRUE)[1]
floor_means <- vapply(chains, function(x) mean(x[, floor_column]), numeric(1))
floor_split <- diff(range(floor_means)) > threshold
high_floor <- if (floor_split) {
  which(floor_means - min(floor_means) > threshold)
} else {
  integer(0)
}
stuck <- which(vapply(chains, function(x) {
  nrow(unique(x)) / nrow(x) < 0.5
}, logical(1)))
drop <- sort(unique(c(high_floor, stuck)))
keep <- setdiff(seq_len(n_chains), drop)

# per-chain means against the within-chain sd
means <- vapply(chains, colMeans, numeric(length(parameters)))
sds <- vapply(chains, function(x) apply(x, 2, sd), numeric(length(parameters)))
spread <- (apply(means, 1, max) - apply(means, 1, min)) / rowMeans(sds)
spread[!is.finite(spread)] <- 0
worst <- order(-spread)[seq_len(min(10, length(spread)))]
spread_table <- data.frame(fit = name, parameter = parameters[worst],
                           spread = spread[worst],
                           matrix(means[worst, ], length(worst),
                                  dimnames = list(NULL, paste0("chain_",
                                                               seq_len(n_chains)))),
                           check.names = FALSE)
write.csv(spread_table, file.path(output_dir, paste0(name, "_spread.csv")),
          row.names = FALSE)

# rank-normalised diagnostics over the chains given
diagnostics <- function(which_chains) {
  if (length(which_chains) == 0) {
    return(c(rhat = NA, n_rhat_101 = NA, n_rhat_105 = NA, ess_bulk_min = NA,
             ess_bulk_median = NA, ess_tail_min = NA, worst = NA))
  }
  array <- simplify2array(chains[which_chains])
  if (length(dim(array)) == 2) {
    array <- array(array, c(dim(array), 1))
  }
  array <- aperm(array, c(1, 3, 2))
  dimnames(array) <- list(NULL, NULL, parameters)
  rhat_values <- apply(array, 3, rhat)
  ess_bulk_values <- apply(array, 3, ess_bulk)
  ess_tail_values <- apply(array, 3, ess_tail)
  rm(array)
  c(rhat = max(rhat_values, na.rm = TRUE),
    n_rhat_101 = sum(rhat_values > 1.01, na.rm = TRUE),
    n_rhat_105 = sum(rhat_values > 1.05, na.rm = TRUE),
    ess_bulk_min = min(ess_bulk_values, na.rm = TRUE),
    ess_bulk_median = median(ess_bulk_values, na.rm = TRUE),
    ess_tail_min = min(ess_tail_values, na.rm = TRUE),
    worst = parameters[which.max(rhat_values)])
}
all_chains <- diagnostics(seq_len(n_chains))
kept_chains <- if (length(drop) > 0) diagnostics(keep) else all_chains

row <- data.frame(
  fit = name, file = file, chains = n_chains,
  draws_per_chain = nrow(chains[[1]]), parameters = length(parameters),
  floor_per_chain = paste(sprintf("%.4f", floor_means), collapse = ", "),
  floor_split = floor_split,
  high_floor_chains = paste(high_floor, collapse = ","),
  stuck_chains = paste(stuck, collapse = ","),
  drop = paste(drop, collapse = ","),
  already_dropped = paste(already_dropped, collapse = ","),
  floor_kept = mean(floor_means[keep]),
  rhat_all = as.numeric(all_chains[["rhat"]]),
  worst_all = all_chains[["worst"]],
  rhat_kept = as.numeric(kept_chains[["rhat"]]),
  worst_kept = kept_chains[["worst"]],
  n_rhat_101_kept = as.numeric(kept_chains[["n_rhat_101"]]),
  n_rhat_105_kept = as.numeric(kept_chains[["n_rhat_105"]]),
  ess_bulk_min_kept = as.numeric(kept_chains[["ess_bulk_min"]]),
  ess_bulk_median_kept = as.numeric(kept_chains[["ess_bulk_median"]]),
  ess_tail_min_kept = as.numeric(kept_chains[["ess_tail_min"]]),
  top_spread = paste(sprintf("%s (%.1f)", parameters[worst[1:3]],
                             spread[worst[1:3]]), collapse = "; "))
write.csv(row, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
writeLines(paste(drop, collapse = ","), file.path(output_dir,
                                                  paste0(name, ".drop")))

options(width = 200)
print(t(row))
print(spread_table, digits = 3)
