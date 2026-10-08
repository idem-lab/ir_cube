# Initial values in each mortality-floor mode (#14, #37), for starting chains
# in both: the posterior means of a full fit's low-floor and high-floor chains,
# in the format fit_model.R caches (dynamical_inits_file), so that
# dynamical_inits() matches them to any model by name.
#
#   Rscript R/floor_mode_inits.R [fit file] [low chains] [high chains]
#
# Defaults: outputs/refit/fits/r2_main.RData (round 2 main fit), chains 1-3 in
# the low-floor mode and chain 4 in the high-floor mode, both from
# draws_all_chains (the draws before drop_stuck_chains.R dropped chain 4).
# Writes temporary/inits_floor_low.RDS and temporary/inits_floor_high.RDS;
# use them with
#   IR_CUBE_INITS=temporary/inits_floor_low.RDS,temporary/inits_floor_high.RDS
# which starts the first half of the chains in the low mode and the second
# half in the high mode (dynamical_chain_inits()). Plain R; about 1.5 GB.

arguments <- commandArgs(trailingOnly = TRUE)
file <- if (length(arguments) >= 1) arguments[1] else
  "outputs/refit/fits/r2_main.RData"
chains <- list(
  low = if (length(arguments) >= 2) {
    as.integer(strsplit(arguments[2], ",")[[1]])
  } else 1:3,
  high = if (length(arguments) >= 3) {
    as.integer(strsplit(arguments[3], ",")[[1]])
  } else 4L)

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
})
source("R/dynamical_predictions.R")

f <- new.env()
load(file, envir = f)
all_draws <- if (exists("draws_all_chains", envir = f, inherits = FALSE)) {
  f$draws_all_chains
} else {
  f$draws
}
dir.create("temporary", showWarnings = FALSE)

for (mode in names(chains)) {
  draws_matrix <- as.matrix(all_draws[chains[[mode]]])
  names <- unique(sub("\\[.*$", "", colnames(draws_matrix)))
  v <- lapply(setNames(nm = names), extract_parameter,
              draws_matrix = draws_matrix)
  # those of the non-centred model, if the fit centred some levels
  # (noncentred_draws()), the form dynamical_inits() reads
  v <- noncentred_draws(v, f$classes_index, f$model_options)
  # posterior means with the dimensions of the greta variables (vectors as
  # one-column matrices), as fit_model.R saves them
  means <- lapply(v, function(x) {
    m <- colMeans(matrix(x, nrow(x)))
    dims <- dim(x)[-1]
    if (length(dims) == 1) dims <- c(dims, 1)
    array(m, dims)
  })
  attr(means, "columns") <- colnames(f$x_cell_years)
  attr(means, "levels") <- dynamical_lookups(f$df)$levels
  out <- sprintf("temporary/inits_floor_%s.RDS", mode)
  saveRDS(means, out)
  cat(sprintf("%s: chains %s, mortality floor %.4f -> %s\n", mode,
              toString(chains[[mode]]), means$mortality_floor, out))
}
