# illustrate the inherent variability in biassay data at a set of locations with
# high sampling density

# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load the mask
mask <- rast("data/clean/raster_mask.tif")

insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  pull(insecticide)

insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")

df <- ir_africa %>%
  group_by(insecticide_type) %>%
  # subset to the most common concentration for each insecticide
  filter(
    concentration == sample_mode(concentration)
  ) %>%
  ungroup() %>%
  filter(
    # drop any from before when we have data on net coverage
    year_start >= baseline_year
  ) %>%
  mutate(
    # create an index to the simulation year (in 1-indexed integers)
    year_id = year_start - baseline_year + 1,
    # add on cell ids corresponding to these observations,
    cell = cellFromXY(mask,
                      as.matrix(select(., longitude, latitude)))
  ) %>%
  # drop a handful of datapoints missing covariates
  filter(
    !is.na(extract(mask, cell)[, 1])
  ) %>%
  # add an index to the vector of unique cells (now they have been subsetted)
  mutate(
    cell_id = match(cell, unique(cell))
  )

# subset df to the cell and year with the most bioassay results
df_most_sampled <- df %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  mutate(
    cell_year_type = sprintf("location %i\n(%i)\n%s",
                             cell_id,
                             year_start,
                             insecticide_type)
  ) %>%
  group_by(
    cell_year_type,
    cell_id,
    year_start,
    insecticide_type
  ) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(count)
  ) %>%
  slice(1:6) %>%
  # bind_cols(
  #   xyFromCell(mask, .$cell)
  # ) %>%
  # reverse_geocode(
  #   lat = y,
  #   long = x,
  #   method = 'osm',
  #   full_results = TRUE
  # ) %>%
  # mutate(
  #   precise_place = str_split_i(address, ",", 1),
  #   precise_place = case_when(
  #     grepl("Placodji", address) ~ "Placodji, Cotonou",
  #     grepl("Busia", address) ~ "Busia",
  #     .default = precise_place
  #   ),
  #   place = paste(precise_place, country, sep = ", ")
  # ) %>%
  select(
    cell_year_type,
    cell_id,
    year_start,
    insecticide_type,
    count,
    # place
  ) %>%
  left_join(
    mutate(df,
           index = row_number()),
    by = c("cell_id", "year_start", "insecticide_type")
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  )

# Estimate the bioassay overdispersion rho, hierarchically by insecticide type
# nested in class.
#
# This previously fitted a single rho to the six most sampled cell-year-type
# combinations, with a free fraction per combination. Two problems with using
# that for the validation noise floor. It used 81 assays where the whole
# replicated set has 9,420, and a single rho is badly misspecified: fitted
# separately, the nine insecticide types range from 0.093 (Lambda-cyhalothrin)
# to 0.256 (Alpha-cypermethrin), and heterogeneity *within* the pyrethroids
# (I2 = 97%) exceeds that between classes. Pooling to the class level forces
# Lambda-cyhalothrin and Alpha-cypermethrin to share 0.164, which overstates one
# floor by 76% and understates the other by 36%.
#
# The model, for group g (a cell x year x insecticide combination) with assays
# j, of insecticide type t in class c:
#
#   died_gj | p_g  ~ BetaBinomial(mosquito_number_gj, p_g, rho_t)
#   p_g            ~ Beta(a0_t, b0_t)
#   logit(rho_t)   ~ Normal(mu_c, sigma_within)
#   mu_c           ~ Normal(mu, sigma_between)
#
# The group fractions are sampled rather than maximised, so this does not suffer
# the incidental parameters problem that joint maximisation over one p_g per
# group would: most groups hold only two assays. Partial pooling across types
# shares strength toward the thinly sampled ones - Pirimiphos-methyl and
# Malathion have about 100 groups each, against 1,157 for Deltamethrin.
#
# Every replicated group is used, not just the six plotted below.
# restricted to the nine insecticides the model covers: `df` here carries all
# seventeen in the data, and the other eight have no class in `classes_index`
replicated_all <- df %>%
  filter(insecticide_type %in% types) %>%
  group_by(cell_id, year_start, insecticide_type) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(group_key = paste(cell_id, year_start, insecticide_type))

group_keys <- unique(replicated_all$group_key)
type_names <- sort(unique(replicated_all$insecticide_type))
class_of_type <- classes[classes_index][match(type_names, types)]
class_names <- sort(unique(class_of_type))

group_id <- match(replicated_all$group_key, group_keys)
type_of_group <- match(
  replicated_all$insecticide_type[match(group_keys, replicated_all$group_key)],
  type_names)
type_id <- match(replicated_all$insecticide_type, type_names)
class_id <- match(class_of_type, class_names)

n_groups <- length(group_keys)
n_types <- length(type_names)
n_classes <- length(class_names)

cat(sprintf("hierarchical rho: %i assays, %i groups, %i types, %i classes\n",
            nrow(replicated_all), n_groups, n_types, n_classes))

# Hyperparameters, non-centred: sampling the standardised deviations rather
# than the levels decouples them from the scale parameters, which a centred
# parameterisation makes hard to move when sigma is small.
mu <- normal(qlogis(0.15), 1)
sigma_between <- normal(0, 0.5, truncation = c(0, Inf))
sigma_within <- normal(0, 0.5, truncation = c(0, Inf))
z_class <- normal(0, 1, dim = n_classes)
z_type <- normal(0, 1, dim = n_types)
mu_class <- mu + sigma_between * z_class
logit_rho_type <- mu_class[class_id] + sigma_within * z_type
rho_type <- ilogit(logit_rho_type)

# The distribution of true fractions, per type. Parameterised by mean and
# concentration rather than by (a0, b0) directly, with a wide prior on the
# concentration: mortality against the organophosphates piles up at 1 - 90% of
# Pirimiphos-methyl assays are above 0.95 - so those types need a Beta skewed
# hard toward 1, which the maximum likelihood fit meets with a0 near 16. A
# lognormal(0, 1) prior on a0 spans only 0.14 to 7.1, and constraining it that
# way leaves rho to absorb the spread instead: it drove the first version of
# this fit to rho = 0.63 for those two types, against 0.13 and 0.19 by maximum
# likelihood.
beta_mean <- beta(1, 1, dim = n_types)
# sampled on the log scale, and kept as its own variable so initial values can
# be set on it: initials attach only to variable nodes, not to operations
log_beta_conc <- normal(log(5), 2, dim = n_types)
beta_conc <- exp(log_beta_conc)
a0 <- beta_mean * beta_conc
b0 <- (1 - beta_mean) * beta_conc
probs <- beta(a0[type_of_group], b0[type_of_group], dim = n_groups)

distribution(replicated_all$died) <- betabinomial_p_rho(
  N = replicated_all$mosquito_number,
  p = probs[group_id],
  rho = rho_type[type_id])

m <- model(rho_type, mu_class, sigma_within, sigma_between)

# initialised near the per-type maximum likelihood estimates, so the chains do
# not have to find the posterior from the prior in a 3,700-parameter model
rho_start <- c(`Alpha-cypermethrin` = 0.256, Bendiocarb = 0.131, DDT = 0.118,
               Deltamethrin = 0.155, Fenitrothion = 0.130,
               `Lambda-cyhalothrin` = 0.093, Malathion = 0.257,
               Permethrin = 0.176, `Pirimiphos-methyl` = 0.189)
mean_start <- c(`Alpha-cypermethrin` = 0.43, Bendiocarb = 0.86, DDT = 0.60,
                Deltamethrin = 0.68, Fenitrothion = 0.97,
                `Lambda-cyhalothrin` = 0.72, Malathion = 0.87,
                Permethrin = 0.59, `Pirimiphos-methyl` = 0.96)
conc_start <- c(`Alpha-cypermethrin` = 6.1, Bendiocarb = 3.1, DDT = 1.5,
                Deltamethrin = 3.6, Fenitrothion = 15.9,
                `Lambda-cyhalothrin` = 1.9, Malathion = 4.4,
                Permethrin = 2.6, `Pirimiphos-methyl` = 17.3)
start <- initials(
  mu = qlogis(0.15),
  sigma_between = 0.2,
  sigma_within = 0.3,
  z_type = as.numeric((qlogis(rho_start[type_names]) - qlogis(0.15)) / 0.3),
  beta_mean = as.numeric(mean_start[type_names]),
  log_beta_conc = as.numeric(log(conc_start[type_names])))

# seeded, so the overdispersion estimates the whole pipeline scores against can
# be reproduced from this script
set.seed(2026 - 9 - 30)
# Twice greta's default leapfrog steps (Lmin 5, Lmax 10). The binding
# constraint is mixing in the hierarchy variances, not warmup: at the default
# this fit reached worst Rhat 1.083 with minimum effective sample size 211 over
# 16,000 post-warmup draws. Longer trajectories cost proportionally more time
# per iteration but decorrelate those parameters.
rho_draws <- mcmc(m, n_samples = 4000, warmup = 5000, chains = 4,
                  sampler = hmc(Lmin = 10, Lmax = 20),
                  initial_values = replicate(4, start, simplify = FALSE))

worst_rhat <- max(coda::gelman.diag(rho_draws, multivariate = FALSE,
                                    autoburnin = FALSE)$psrf[, 1])
cat(sprintf("worst Rhat %.3f, minimum effective sample size %.0f\n",
            worst_rhat, min(coda::effectiveSize(rho_draws))))
if (worst_rhat > 1.05) {
  warning("the hierarchical rho fit has not converged (Rhat ", round(worst_rhat, 3),
          "); do not use these estimates")
}

rho_summary <- summary(rho_draws)
rho_type_mean <- rho_summary$statistics[paste0("rho_type[", seq_len(n_types), ",1]"),
                                        "Mean"]
rho_type_lower <- rho_summary$quantiles[paste0("rho_type[", seq_len(n_types), ",1]"),
                                        "2.5%"]
rho_type_upper <- rho_summary$quantiles[paste0("rho_type[", seq_len(n_types), ",1]"),
                                        "97.5%"]
names(rho_type_mean) <- type_names

# The joint draws, not only the summaries. The hierarchy correlates rho between
# types, so rebuilding them as independent logit-normals from each type's
# interval - which is what the noise-share interval used to do - understates the
# uncertainty in any quantity averaged over types (#12 review).
rho_type_draws <- do.call(rbind, lapply(rho_draws, function(chain) {
  chain[, paste0("rho_type[", seq_len(n_types), ",1]"), drop = FALSE]
}))
colnames(rho_type_draws) <- type_names
saveRDS(rho_type_draws, "outputs/bioassay_rho_type_draws.rds")

rho_hierarchical <- data.frame(
  insecticide_type = type_names,
  insecticide_class = class_of_type,
  rho = as.numeric(rho_type_mean),
  rho_lower = as.numeric(rho_type_lower),
  rho_upper = as.numeric(rho_type_upper),
  n_groups = as.numeric(table(factor(type_of_group,
                                     levels = seq_len(n_types)))),
  n_assays = as.numeric(table(factor(type_id, levels = seq_len(n_types)))),
  worst_rhat = worst_rhat)
write.csv(rho_hierarchical, "outputs/bioassay_rho_hierarchical.csv",
          row.names = FALSE)

cat("\nhierarchical rho by insecticide type:\n")
print(rho_hierarchical %>% mutate(across(where(is.numeric), ~ round(.x, 4))))

# the overall level, for the figure subtitle: the posterior mean of rho across
# types, weighted by the assays each contributes
rho_bayes <- weighted.mean(rho_hierarchical$rho, rho_hierarchical$n_assays)

# compute expected and observed statistics of these to elucidate the distribution
df_most_sampled_stats <- df_most_sampled %>%
  group_by(cell_year_type, insecticide_type) %>%
  summarise(
    # point estimate of suitability
    Susceptibility = sum(died) / sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    # 95% CIs under binomial (independence) sampling assumption
    binom_lower_100 = qbinom(0.025, 100, Susceptibility) / 100,
    binom_upper_100 = qbinom(0.975, 100, Susceptibility) / 100,
    # 95% CIs under betabinomial (non-independence) sampling assumption, using
    # the posterior mean estimated by our model
    # MLE is lower, maybe biased down
    # each example gets the rho of its own insecticide type, rather than one
    # value shared across all of them
    rho = rho_type_mean[insecticide_type],
    alpha = Susceptibility * (1 / rho - 1),
    beta = alpha * (1 - Susceptibility) / Susceptibility,
    betabinom_lower_100 = qbbinom(0.025, 100, alpha, beta) / 100,
    betabinom_upper_100 = qbbinom(0.975, 100, alpha, beta) / 100,
  )

colour_types <- scales::hue_pal(direction = -1)(9)
types_plot_id <- match(insecticides_plot_small, insecticides_plot)
colours_plot <- colour_types[types_plot_id]

set.seed(1)
df_most_sampled %>%
  # shuffle the x axes so the points don't overlap too much
  mutate(
    x_random = sample(seq(0.2, 1.8, length.out = n()))
  ) %>%
  ggplot(
    aes(
      xmin = 0,
      xmax = 2,
      group = cell_year_type
    )
  ) +
  geom_rect(
    aes(
      ymax = betabinom_upper_100,
      ymin = betabinom_lower_100,
    ),
    data = df_most_sampled_stats,
    fill = grey(0.9)
  ) +
  geom_rect(
    aes(
      ymax = binom_upper_100,
      ymin = binom_lower_100,
    ),
    data = df_most_sampled_stats,
    fill = grey(0.7)
  ) +
  geom_rect(
    aes(
      ymax = Susceptibility,
      ymin = Susceptibility
    ),
    data = df_most_sampled_stats,
    linewidth = 1,
    colour = grey(0.4)
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = mosquito_number,
      x = x_random,
      fill = insecticide_type
    ),
    shape = 21
  ) +
  facet_wrap(~cell_year_type,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    limits = c(0, 2)
  ) +
  scale_fill_manual(
    values = colours_plot
  ) +
  xlab("") +
  theme_minimal() +
  guides(
    fill = "none",
    size = guide_legend(title = "No. tested")
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggtitle(
    "Significant variability in replicated bioassay results",
    sprintf("95%s sampling intervals under binomial (mid grey) and betabinomial
(light grey) sampling, at N=100. Correlation rho is estimated hierarchically by
insecticide type (%s here; %s to %s across all nine types).",
            "%",
            paste(sprintf("%s %.2f", insecticides_plot_small,
                          rho_type_mean[insecticides_plot_small]),
                  collapse = ", "),
            round(min(rho_hierarchical$rho), 2),
            round(max(rho_hierarchical$rho), 2))
  )

ggsave(
  "figures/bioassay_variability.png",
  bg = "white",
  scale = 0.9,
  width = 8,
  height = 6
)

# Multiple discriminating concentration bioassay results against Pyrethroids at
# a single locations (5km grid cells) and years, as an example of the inherent
# variability in this type of data. Dark grey horizontal lines give the (sample
# size weighted) mean of the susceptibility estimates, medium grey bands give
# the 95% sampling interval of a binomial distribution, with this mean and 100
# samples, light grey bands give the 95% sampling interval of a beta-binomial
# distribution, with this mean, 100 samples, and correlation parameter estimate
# from each dataset independently.

# Using this estimate of the dispersion in susceptibility bioassays (and
# assuming each was sampled from a single location), plot how the statistical
# power to estimate the population susceptibility depends on the total number of
# mosquitoes assayed, and the number of unique collections they are collated
# from

tibble(
  n = seq(50, 500, by = 10)
) %>%
  rowwise() %>%
  mutate(
    moe_1 = moe_betabinomial(n, rho = rho_bayes),
    moe_3 = moe_betabinomial_cluster(n, rho = rho_bayes, clusters = 3),
    moe_5 = moe_betabinomial_cluster(n, rho = rho_bayes, clusters = 5),
    moe_10 = moe_betabinomial_cluster(n, rho = rho_bayes, clusters = 10),
    moe_50 = moe_betabinomial_cluster(n, rho = rho_bayes, clusters = 50),
    moe_independent = moe_binomial(n)
  ) %>%
  pivot_longer(
    cols = starts_with("moe"),
    names_to = "clusters",
    values_to = "moe",
    names_prefix = "moe_"
  ) %>%
  mutate(
    clusters = factor(
      clusters,
      levels = c(na.omit(unique(as.numeric(clusters))), "independent")
    )
  ) %>%
  ggplot(
    aes(
      x = n,
      y = moe,
      group = clusters,
      colour = clusters
    )
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 0.5)) +
  theme_minimal() +
  scale_x_continuous(breaks = c(50, 100, 250, 500)) +
  ylab("Margin of error") +
  xlab("Number of mosquitoes") +
  ggtitle(
    "Statistical power of cluster-stratified susceptibility bioassays",
  )

ggsave("figures/cluster_sampling_power.png",
       bg = "white",
       width = 6,
       height = 5)
