# see how variability in bioassays data looks for resistant colony mosquitoes
# against net samples (only brand new nets, to try to limit variability)

# greta first, so python starts before terra and sf are attached
source("R/greta_setup.R")
start_greta()
source("R/packages.R")
source("R/functions.R")

# load survival data for a resistance lab colony of resistant An. coluzzii
# against fresh dual-AI nets
ig2 <- "data/raw/durability" %>%
  list.files(
    recursive = TRUE,
    full.names = TRUE) %>%
  tibble(path = .) %>%
  rowwise() %>%
  mutate(
    file = basename(path),
    country = str_split_1(file, "_")[1],
    phase = str_split_1(file, "_")[2],
    type = str_split_1(file, "_")[3],
    type = str_remove(type, "^Bionet "),
    type = str_remove(type, ".csv$")
  ) %>%
  filter(
    type == "Interceptor G2"
  ) %>%
  mutate(
    data = list(read_csv(path))
  ) %>%
  unnest(
    data
  ) %>%
  filter(
    # only resistant vectors
    strain == "An. coluzzii VKPER",
    # only nets that have not been distributed
    phase == "baseline",
    is.na(loc)
    # strain == "An. gambiae Kisumu",
    # wash == 0
  ) %>%
  select(country,
         phase,
         type,
         strain,
         mosquito_number = n,
         died = n_dead24h) %>%
  mutate(
    phase_id = match(phase, unique(phase)),
    Susceptibility = died / mosquito_number
  )
  
rho <- variable(0, 1)
prob <- variable(0, 1, dim = max(ig2$phase_id))

distribution(ig2$died) <- betabinomial_p_rho(N = ig2$mosquito_number,
                                                p = prob[ig2$phase_id],
                                                rho = rho)
m <- model(rho)  
rho_draws <- mcmc(m)
o <- greta::opt(m, adjust = FALSE)
rho_ml <- o$par$rho
rho_bayes <- summary(rho_draws)$statistics["Mean"]

# compute expected and observed statistics of these to elucidate the distribution
ig2_stats <- ig2 %>%
  group_by(phase) %>%
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
    # rho = rho_ml,
    rho = rho_bayes,
    alpha = Susceptibility * (1 / rho - 1),
    beta = alpha * (1 - Susceptibility) / Susceptibility,
    betabinom_lower_100 = qbbinom(0.025, 100, alpha, beta) / 100,
    betabinom_upper_100 = qbbinom(0.975, 100, alpha, beta) / 100,
  )

set.seed(1)
ig2 %>%
  # shuffle the x axes so the points don't overlap too much
  mutate(
    x_random = sample(seq(0.2, 1.8, length.out = n()))
  ) %>%
  ggplot(
    aes(
      xmin = 0,
      xmax = 2,
      group = phase
    )
  ) +
  geom_rect(
    aes(
      ymax = betabinom_upper_100,
      ymin = betabinom_lower_100,
    ),
    data = ig2_stats,
    fill = grey(0.9)
  ) +
  geom_rect(
    aes(
      ymax = binom_upper_100,
      ymin = binom_lower_100,
    ),
    data = ig2_stats,
    fill = grey(0.7)
  ) +
  geom_rect(
    aes(
      ymax = Susceptibility,
      ymin = Susceptibility
    ),
    data = ig2_stats,
    linewidth = 1,
    colour = grey(0.4)
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = mosquito_number,
      x = x_random,
    ),
    fill = "green",
    shape = 21
  ) +
  facet_wrap(~phase,
             nrow = 1,
             strip.position = "bottom") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    limits = c(0, 2)
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
    "Significant variability in replicated bioassay results with colony mosquitoes",
    sprintf("95%s sampling intervals under binomial (mid grey) and betabinomial
(light grey; estimated correlation rho=%s) intervals assume sample size N=100.",
            "%",
            round(rho_bayes, 2))
  )

ggsave(
  "figures/colony_net_bioassay_variability.png",
  bg = "white",
  scale = 0.9,
  width = 3,
  height = 6
)
