# Estimate the bioassay overdispersion parameter rho from replicate bioassays.
#
# Bioassays carried out in the same 5km pixel, in the same year, with the same
# insecticide, are replicate measurements of a single population susceptible
# fraction, so the variation between them identifies the observation
# overdispersion without reference to any spatial model. That makes the
# resulting rho an external quantity: it can be compared with the rho estimated
# inside the dynamical model, and it sets the noise floor against which
# out-of-sample predictive performance is judged (see idem-lab/ir_cube#10).
#
# rho was previously estimated this way in fig_illustrate_bioassay_variability.R,
# but from the six most heavily sampled pixel-year-insecticide combinations
# only, which is enough to illustrate a figure but thin for anchoring a
# validation metric. This uses every replicated combination, and estimates rho
# separately for each insecticide class.
#
# The model is
#
#   died_gj ~ BetaBinomial(mosquito_number_gj, p_g, rho)
#   p_g     ~ Beta(a0, b0)
#
# with the group fractions p_g integrated out numerically rather than
# maximised over. Maximising jointly over one p_g per group would be subject to
# the incidental parameters problem: with most groups holding only two assays,
# the number of nuisance parameters grows with the sample size and the
# resulting rho is biased downwards.

source("R/validation_functions.R")

suppressMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

mask <- rast("data/clean/raster_mask.tif")
ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")

baseline_year <- 1995
final_data_year <- 2024

insecticides_keep <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin",
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")

# same subsetting as the model fitting and validation scripts
sample_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

df <- ir_africa %>%
  filter(insecticide_type %in% insecticides_keep) %>%
  group_by(insecticide_type) %>%
  filter(concentration == sample_mode(concentration)) %>%
  ungroup() %>%
  filter(
    year_start >= baseline_year,
    year_start <= final_data_year,
    !is.na(died),
    !is.na(mosquito_number),
    mosquito_number > 1
  ) %>%
  mutate(
    cell = cellFromXY(mask, as.matrix(select(., longitude, latitude)))
  ) %>%
  filter(!is.na(terra::extract(mask, cell)[, 1]))

# identify replicate groups: one pixel, one year, one insecticide
df <- df %>%
  group_by(cell, year_start, insecticide_type) %>%
  mutate(
    group_id = cur_group_id(),
    group_size = n()
  ) %>%
  ungroup()

replicated <- df %>%
  filter(group_size >= 2)

cat(sprintf("%i assays in %i replicated pixel-year-insecticide groups\n",
            nrow(replicated),
            n_distinct(replicated$group_id)))


# marginal likelihood ------------------------------------------------------

# Gauss-Legendre nodes and weights on (0, 1), for integrating over the group
# fraction. A fixed grid keeps the likelihood smooth in the parameters, which
# matters for the optimiser
gauss_legendre_01 <- function(n_nodes = 128) {
  # nodes of the Legendre polynomial by Newton iteration on the interval (-1, 1)
  i <- seq_len(n_nodes)
  x <- cos(pi * (i - 0.25) / (n_nodes + 0.5))
  for (iteration in 1:100) {
    p0 <- rep(1, n_nodes)
    p1 <- x
    for (j in 2:n_nodes) {
      p2 <- ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
      p0 <- p1
      p1 <- p2
    }
    derivative <- n_nodes * (x * p1 - p0) / (x ^ 2 - 1)
    step <- p1 / derivative
    x <- x - step
    if (max(abs(step)) < 1e-14) break
  }
  weights <- 2 / ((1 - x ^ 2) * derivative ^ 2)
  # map from (-1, 1) to (0, 1)
  list(node = (x + 1) / 2,
       weight = weights / 2)
}

quadrature <- gauss_legendre_01(128)

# negative log marginal likelihood, with parameters on unconstrained scales:
#   par[1] = logit(rho), par[2] = log(a0), par[3] = log(b0)
negative_log_likelihood <- function(par, died, mosquito_number, group_id,
                                    quadrature) {

  rho <- plogis(par[1])
  a0 <- exp(par[2])
  b0 <- exp(par[3])

  node <- quadrature$node
  n_nodes <- length(node)
  n_obs <- length(died)

  # log density of each observation at each quadrature node
  shape <- bb_shape(rep(node, each = n_obs),
                    rho)
  log_density <- lchoose(rep(mosquito_number, times = n_nodes),
                         rep(died, times = n_nodes)) +
    lbeta(rep(died, times = n_nodes) + shape$a,
          rep(mosquito_number - died, times = n_nodes) + shape$b) -
    lbeta(shape$a, shape$b)
  log_density <- matrix(log_density, nrow = n_obs, ncol = n_nodes)

  # sum the log densities within each group, at each node
  group_log_density <- rowsum(log_density, group_id)

  # prior density of the group fraction at each node, and the quadrature weight
  log_prior <- dbeta(node, a0, b0, log = TRUE) + log(quadrature$weight)

  # integrate over the group fraction, working on the log scale for stability
  integrand <- sweep(group_log_density, 2, log_prior, FUN = "+")
  maximum <- apply(integrand, 1, max)
  group_log_likelihood <- maximum +
    log(rowSums(exp(sweep(integrand, 1, maximum, FUN = "-"))))

  total <- sum(group_log_likelihood)
  if (!is.finite(total)) return(1e10)
  -total

}

# fit to one set of assays, returning rho with a standard error from the
# observed information
fit_rho <- function(data) {

  group_id <- match(data$group_id, unique(data$group_id))

  fit <- optim(
    par = c(qlogis(0.14), log(4), log(1.5)),
    fn = negative_log_likelihood,
    died = data$died,
    mosquito_number = data$mosquito_number,
    group_id = group_id,
    quadrature = quadrature,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 500)
  )

  rho <- plogis(fit$par[1])
  # delta method for the standard error of rho from that of logit(rho)
  standard_error_logit <- sqrt(diag(solve(fit$hessian)))[1]
  standard_error <- standard_error_logit * rho * (1 - rho)

  data.frame(
    rho = rho,
    rho_lower = plogis(fit$par[1] - 1.96 * standard_error_logit),
    rho_upper = plogis(fit$par[1] + 1.96 * standard_error_logit),
    standard_error = standard_error,
    n_assays = nrow(data),
    n_groups = length(unique(group_id)),
    convergence = fit$convergence
  )

}

# overall, and by insecticide class
rho_overall <- fit_rho(replicated) %>%
  mutate(insecticide_class = "all", .before = everything())

rho_by_class <- replicated %>%
  group_by(insecticide_class) %>%
  filter(n_distinct(group_id) >= 20) %>%
  group_split() %>%
  lapply(function(data) {
    fit_rho(data) %>%
      mutate(insecticide_class = data$insecticide_class[1],
             .before = everything())
  }) %>%
  bind_rows()

rho_estimates <- bind_rows(rho_overall, rho_by_class)

cat("\nestimated overdispersion:\n")
print(rho_estimates %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

write.csv(rho_estimates, "outputs/bioassay_rho.csv", row.names = FALSE)


# comparison with the previous estimate ------------------------------------

# the six most heavily sampled combinations, as used previously
top_six <- replicated %>%
  count(group_id, sort = TRUE) %>%
  slice(1:6) %>%
  pull(group_id)

rho_top_six <- replicated %>%
  filter(group_id %in% top_six) %>%
  fit_rho()

cat(sprintf("\nrho from the six most sampled groups only: %.3f (%.3f - %.3f), %i assays\n",
            rho_top_six$rho, rho_top_six$rho_lower, rho_top_six$rho_upper,
            rho_top_six$n_assays))
cat(sprintf("rho from all replicated groups:            %.3f (%.3f - %.3f), %i assays\n",
            rho_overall$rho, rho_overall$rho_lower, rho_overall$rho_upper,
            rho_overall$n_assays))


# are replicated pixel-years representative? -------------------------------

# The floor rests on these groups, but repeatedly sampled sites are plausibly
# sentinel sites, and may differ systematically from the rest of the data. This
# compares them on the quantities that matter for the floor
representativeness <- df %>%
  mutate(
    replicated = ifelse(group_size >= 2, "replicated", "single")
  ) %>%
  group_by(replicated) %>%
  summarise(
    assays = n(),
    pixel_years = n_distinct(group_id),
    mean_mortality = mean(died / mosquito_number),
    sd_mortality = sd(died / mosquito_number),
    median_mosquito_number = median(mosquito_number),
    median_year = median(year_start),
    countries = n_distinct(country_name),
    proportion_pyrethroid = mean(insecticide_class == "Pyrethroids"),
    .groups = "drop"
  )

cat("\nreplicated versus singly sampled pixel-years:\n")
print(representativeness %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

write.csv(representativeness,
          "outputs/bioassay_rho_representativeness.csv",
          row.names = FALSE)

# the same comparison by country, to show whether the replicated groups come
# from a narrow set of places
by_country <- df %>%
  mutate(replicated = group_size >= 2) %>%
  group_by(country_name) %>%
  summarise(
    assays = n(),
    proportion_replicated = mean(replicated),
    .groups = "drop"
  ) %>%
  arrange(desc(assays))

cat("\ntop countries by number of assays:\n")
print(head(by_country, 10))

representativeness_plot <- df %>%
  mutate(
    replicated = ifelse(group_size >= 2,
                        "replicated pixel-year",
                        "single assay")
  ) %>%
  ggplot(
    aes(x = died / mosquito_number,
        fill = replicated)
  ) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.6,
    position = "identity"
  ) +
  facet_wrap(~ insecticide_class, scales = "free_y") +
  scale_fill_manual(values = c("replicated pixel-year" = "#B2182B",
                               "single assay" = grey(0.5)),
                    name = "") +
  labs(
    x = "observed mortality",
    y = "density",
    title = "Are repeatedly sampled pixel-years representative?",
    subtitle = "the overdispersion estimate, and so the noise floor, rests on the replicated groups"
  ) +
  theme_minimal()

ggsave("figures/bioassay_rho_representativeness.png",
       representativeness_plot,
       bg = "white",
       width = 10,
       height = 6)
