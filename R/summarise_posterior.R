# Posterior summary of the full fit's interpretable parameters, against their
# priors: the mortality floor, reversion per class, rho per type (with the
# external replicate-based estimate), the initial-state covariate coefficients
# per type, and the selection coefficients per type (beta_type; the effect on
# selection is exp(beta_type)).
#
#   Rscript R/summarise_posterior.R [fitted_model.RData] [output stem]
#
# Writes <stem>.csv (default outputs/posterior_summary.csv): per parameter the
# posterior mean, sd and 95% interval, the prior sd and 95% interval (from
# 4,000 prior simulations of the same greta arrays), the contraction
# 1 - posterior variance / prior variance, and rank-normalised Rhat and bulk
# and tail ESS over the ordered MCMC draws.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- if (length(arguments) >= 1) arguments[1] else
  "temporary/fitted_model.RData"
stem <- if (length(arguments) >= 2) arguments[2] else
  "outputs/posterior_summary"

source("R/greta_setup.R")
start_greta(threads = 2)
suppressMessages({
  library(dplyr)
  library(posterior)
})

fit <- new.env()
load(fit_file, envir = fit)

n_prior <- 4000

# the greta arrays to summarise, with dimnames for their elements
covariate_names <- colnames(fit$x_cell_years)
init_names <- fit$model_options$init_covariates
targets <- list(
  mortality_floor = list(fit$mortality_floor, NULL),
  reversion_rate = list(fit$reversion_rate, list(fit$classes)),
  rho_types = list(fit$rho_types, list(fit$types)),
  init_coef = list(fit$init_coef, list(init_names, fit$types)),
  beta_type = list(fit$beta_type, list(covariate_names, fit$types))
)
targets <- Filter(function(x) !is.null(x[[1]]), targets)

element_labels <- function(name, names) {
  if (is.null(names)) return(name)
  grid <- expand.grid(names, stringsAsFactors = FALSE)
  sprintf("%s[%s]", name, do.call(paste, c(grid, sep = ", ")))
}

rows <- list()
for (name in names(targets)) {
  target <- targets[[name]][[1]]
  labels <- element_labels(name, targets[[name]][[2]])

  # posterior: a deterministic function of each draw, keeping the chains
  post <- calculate(target, values = fit$draws)
  post_array <- as_draws_array(coda::as.mcmc.list(
    lapply(post, function(x) coda::mcmc(as.matrix(x)))))
  post_matrix <- as.matrix(do.call(rbind, lapply(post, as.matrix)))

  # prior: simulations of the same array from the model's priors
  prior <- calculate(target, nsim = n_prior)[[1]]
  prior_matrix <- matrix(prior, n_prior)

  stopifnot(ncol(post_matrix) == length(labels),
            ncol(prior_matrix) == length(labels))
  for (j in seq_along(labels)) {
    x <- post_matrix[, j]
    z <- prior_matrix[, j]
    v <- post_array[, , j]
    rows[[length(rows) + 1]] <- tibble(
      parameter = name,
      element = labels[j],
      mean = mean(x), sd = sd(x),
      lower = quantile(x, 0.025), upper = quantile(x, 0.975),
      prior_sd = sd(z),
      prior_lower = quantile(z, 0.025), prior_upper = quantile(z, 0.975),
      contraction = 1 - var(x) / var(z),
      rhat = rhat(v), ess_bulk = ess_bulk(v), ess_tail = ess_tail(v))
  }
}
summary_table <- bind_rows(rows)

# the external replicate-based rho per type
external <- read.csv("outputs/bioassay_rho_hierarchical.csv")
summary_table <- summary_table %>%
  mutate(type = ifelse(parameter == "rho_types",
                       sub("^rho_types\\[(.*)\\]$", "\\1", element),
                       NA_character_)) %>%
  left_join(select(external, type = insecticide_type,
                   external_rho = rho, external_lower = rho_lower,
                   external_upper = rho_upper),
            by = "type") %>%
  select(-type)

write.csv(summary_table, paste0(stem, ".csv"), row.names = FALSE)
options(width = 200)
print(as.data.frame(summary_table %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))
