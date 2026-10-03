# What separates the two mortality-floor modes (#14) of a full fit: predicted
# mortality at every assay and fitted country trajectories under each mode,
# with per-assay expected log likelihood and Pearson residuals.
#   Rscript R/floor_modes.R <fitted_model.RData> <low chains> <high chains> <label>
# e.g. Rscript R/floor_modes.R temporary/fitted_model.RData 1,2,3 4 round2
# Chains are given as comma-separated indices. Writes
# outputs/floor_modes_<label>.rds; figures are made by
# R/floor_modes_figures.R.
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
chains <- list(low = as.integer(strsplit(args[2], ",")[[1]]),
               high = as.integer(strsplit(args[3], ",")[[1]]))
label <- args[4]
n_per_chain <- 150

source("R/greta_setup.R"); start_greta(threads = 4)
suppressMessages({library(greta); library(dplyr); library(stringr)
  library(tidyr); library(Matrix)})
source("R/functions.R"); source("R/model_covariates.R")
source("R/dynamical_model.R"); source("R/dynamical_predictions.R")
f <- new.env(); load(file, envir = f)
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_iter <- nrow(as.matrix(f$draws[[1]]))
idx <- round(seq(1, n_iter, length.out = n_per_chain))
df <- f$df
n_times <- max(f$cell_years_index$year_id)
years <- 1994 + seq_len(n_times)
stopifnot(all(years[df$year_id] == df$year_start))

# net use (pyrethroid-weighted, as the model's selection column) at each
# cell-year, and its running sum from 1995
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_times)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
irs <- matrix(f$x_cell_years[x_row, "irs"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))
df$nets <- nets[cbind(df$cell_id, df$year_id)]
df$cum_nets <- cum_nets[cbind(df$cell_id, df$year_id)]
df$irs <- irs[cbind(df$cell_id, df$year_id)]

group_of_type <- ifelse(f$types == "DDT", "DDT",
                        ifelse(f$classes_index == 1, "Pyrethroids", "other"))
df$group <- group_of_type[df$type_id]

# trajectory rows: every year at each (cell, type) with pyrethroid or DDT
# assays, weighted by the assays there, averaged within country x group x year
ct <- df |> filter(group != "other") |>
  count(country_id, cell_id, type_id, group, name = "w")
traj_rows <- ct[rep(seq_len(nrow(ct)), each = n_times), ]
traj_rows$year_id <- rep(seq_len(n_times), nrow(ct))
traj_rows$key <- paste(traj_rows$country_id, traj_rows$group,
                       traj_rows$year_id)
keys <- unique(traj_rows$key)
W <- sparseMatrix(i = seq_len(nrow(traj_rows)),
                  j = match(traj_rows$key, keys), x = traj_rows$w)
W <- W %*% Diagonal(x = 1 / colSums(W))

out <- list()
for (m in names(chains)) {
  rows <- unlist(lapply(chains[[m]], function(k) (k - 1) * n_iter + idx))
  par <- dynamical_parameter_draws(fold, f$classes_index, f$types, df,
                                   draw_index = rows)
  p <- plogis(dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                          year_id),
                              df, f$x_cell_years, f$cell_years_index))
  rho <- par$rho_types
  ll <- matrix(NA_real_, nrow(p), ncol(p))
  for (i in seq_len(nrow(p))) {
    pr <- p[i, ]; rr <- rho[i, df$type_id]
    a <- pr * (1 / rr - 1); bb <- a * (1 - pr) / pr
    ll[i, ] <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                   beta = bb, log = TRUE)
  }
  pm <- colMeans(p); rm_ <- colMeans(rho)[df$type_id]
  n <- df$mosquito_number
  assay <- tibble(mode = m, row = seq_len(nrow(df)), pred = pm,
                  pred_lo = apply(p, 2, quantile, 0.05),
                  pred_hi = apply(p, 2, quantile, 0.95),
                  loglik = colMeans(ll),
                  pearson = (df$died - n * pm) /
                    sqrt(n * pm * (1 - pm) * (1 + (n - 1) * rm_)))
  cat(m, ": mean total loglik", sum(colMeans(ll)), "\n")
  rm(p, ll); gc()

  tp <- plogis(dynamical_logit(par, select(traj_rows, cell_id, country_id,
                                            type_id, year_id),
                               df, f$x_cell_years, f$cell_years_index))
  agg <- as.matrix(tp %*% W)
  rm(tp); gc()
  traj <- tibble(mode = m, key = keys, mean = colMeans(agg),
                 lo = apply(agg, 2, quantile, 0.05),
                 hi = apply(agg, 2, quantile, 0.95)) |>
    separate(key, c("country_id", "group", "year_id"), sep = " ",
             convert = TRUE)
  out[[m]] <- list(assay = assay, traj = traj,
                   floor = par$mortality_floor, rho = rho,
                   kappa = par$kappa_type,
                   effect_type = par$effect_type)
}
saveRDS(list(df = df, countries = f$countries, types = f$types,
             years = years, chains = chains, modes = out,
             net_country_year = tibble(
               country_id = rep(df$country_id[match(seq_len(nrow(nets)),
                                                    df$cell_id)], n_times),
               year_id = rep(seq_len(n_times), each = nrow(nets)),
               nets = c(nets), irs = c(irs))),
        sprintf("outputs/floor_modes_%s.rds", label))
cat("done\n")
