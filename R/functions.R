# misc functions

# list of the countries in MTM data that are in malaria-endemic Africa
mtm_africa_countries <- function() {
  africa_countries <- c("Burkina Faso",
                        "Botswana",
                        "Mali",
                        "Kenya",
                        "Malawi",
                        "Cameroon",
                        "Niger",
                        "Benin",
                        "Côte d’Ivoire",
                        "Togo",
                        "Nigeria",
                        "United Republic of Tanzania",
                        "Ethiopia",
                        "Mozambique",
                        "Uganda",
                        "Ghana",
                        "Senegal",
                        "Democratic Republic of the Congo",
                        "Somalia",
                        "South Africa",
                        "Zambia",
                        "Eritrea",
                        "Sudan",
                        "Sierra Leone",
                        "Angola",
                        "Namibia",
                        "Djibouti",
                        "Rwanda",
                        "Chad",
                        "Liberia",
                        "Zimbabwe",
                        "Eswatini",
                        "Guinea",
                        "Equatorial Guinea",
                        "Congo",
                        "Central African Republic",
                        "Madagascar",
                        "Comoros",
                        "Sao Tome and Principe",
                        "Mauritania",
                        "Cabo Verde",
                        "Gabon",
                        "Burundi",
                        "Guinea-Bissau",
                        "Gambia")
}


# given the percent mortality, infer the sample size as the integer closest to
# 100 that implies an integer number that died. If one can't be found, within
# the limits, set to the default instead
infer_sample_size <- function(mortality_percent,
                              limits = c(50, 200),
                              default = 100,
                              tol = 1e-12) {
  
  # try integer values between the limits
  possible_integers <- seq(limits[1], limits[2])
  distance <- abs(possible_integers - default)
  
  # for all of these, see if the mortality rate gives an integer
  died <- possible_integers %*% t(mortality_percent / 100)
  is_integer <- abs(died - round(died)) < tol
  
  # find the value closest to 100 that gives an integer
  valid_distance <- sweep(ifelse(is_integer, 0, Inf),
                          1,
                          distance,
                          FUN = "+")
  
  # for any that don't have a match, put in the index for 100
  missing_index <- which(possible_integers == default)
  
  # function to find the smallest, without dropping values
  closest <- function(distances) {
    closest_index <- which.min(distances)
    if (length(closest_index) == 0) {
      closest_index <- missing_index
    } 
    closest_index
  }

  closest_index <- apply(valid_distance, 2, closest)
  
  # missing <- !apply(is_integer, 2, any)
  # closest_index[missing] <- missing_index
  
  # look up the corresponding integers
  possible_integers[closest_index]
  
}

# compute the mode of a sample
sample_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# given a vector of sample sizes N, predicted proportions p, and overdispersion
# parameter rho (on unit interval, where 0 means no overdispersioan and 1 means
# maximum), returna betatbinomial-distributed variable greta array
betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p

  # define betabinomial according to the greta interface
  beta_binomial(size = N, alpha = a, beta = b)
  
}

nearish <- function(x, y) {
  dplyr::near(x, y, tol = 1e-2)
}

# given a raster file, force it to be written to disk and not in memory (save to
# disk and reload). filename is the name fo the file to write it to, dots is any
# other arguments for writeraster
force_to_disk <- function(raster, filename = tempfile(fileext = ".tif"), ...) {
  terra::writeRaster(x = raster,
                     filename = filename,
                     overwrite = TRUE,
                     ...)
  terra::rast(filename)
}

# random RNG seed (copied from greta, used to fix downstream seeds, conditional
# on top-level randomness)
get_seed <- function ()  {
  sample.int(n = 2^30, size = 1)
}

# get rid of all the craps in plotting the IR maps
theme_ir_maps <- function() {
  theme_minimal() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.title = element_markdown(),
          legend.title = element_markdown())
}

# vectorised quantile function of the betabinomial distribution, by sampling 
qbbinom_once <- function(p, size, alpha = 1, beta = 1, nsims = 1e5) {
  sims <- extraDistr::rbbinom(nsims, size, alpha, beta)
  quantile(sims, p)
}
qbbinom <- Vectorize(qbbinom_once, c("p", "size", "alpha", "beta"))


# Compute the statistical power to estimate the population susceptibility
# (margin of error of the susceptibility estimate) under a betabinomial
# distribution, where: n = sample size, p = probability of a positive; rho =
# overdispersion; alpha = significance level (0.05 = MoE for 95% CI)
moe_betabinomial <- function(n = 100, p = 0.5, rho = 0.1, alpha = 0.05) {
  
  # alpha and beta of the betabinomial distribution (avoiding name clashes)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  
  # compute variance of the number dead
  num <- n * a * b * (a + b + n)
  denom <- (a + b) ^ 2 * (a + b + 1)
  var <- num / denom
  
  # to get the variance of the proportion, multiply the random variable for the
  # number dead by 1/n; compute the variance of this new random variable. From
  # normal distribution closure under affine transformation the variance of the
  # proportion is this variance divided by 1/n^2
  var_prop <- var / (n ^ 2)
  
  # get the standard deviation, and the standard error of the estimate
  sd_prop <- sqrt(var_prop)
  
  # compute Z score (for a two-tailed test at this significance level)
  z <- qnorm(1 - alpha / 2)
  
  # moe is z score times the standard deviation
  z * sd_prop
  
}

# Compute the statistical power to estimate the population susceptibility
# (margin of error of the susceptibility estimate) under a binomial
# distribution, where: n = sample size, p = probability of a positive; alpha =
# significance level (0.05 = MoE for 95% CI)
moe_binomial <- function(n = 100, p = 0.5, alpha = 0.05) {
  
  # compute variance of the number dead
  var <- n * p * (1 - p)
  
  # to get the variance of the proportion, multiply the random variable for the
  # number dead by 1/n; compute the variance of this new random variable. From
  # normal distribution closure under affine transformation the variance of the
  # proportion is this variance divided by 1/n^2
  var_prop <- var / (n ^ 2)
  
  # get the standard deviation, and the standard error of the estimate
  sd_prop <- sqrt(var_prop)
  
  # compute Z score (for a two-tailed test at this significance level)
  z <- qnorm(1 - alpha / 2)
  
  # moe is z score times the standard deviation
  z * sd_prop
  
}

# Compute the statistical power to estimate the population susceptibility
# (margin of error of the susceptibility estimate) under a sum of independent
# betabinomial samples, where: n = sample size, p = probability of a positive;
# rho = overdispersion; alpha = significance level (0.05 = MoE for 95% CI),
# clusters = number of independent samples (n is divided between these), and
# nsim = number of Monte Carlo simulations to perform to estimate the MoE
moe_betabinomial_cluster <- function(n = 100, p = 0.5, rho = 0.1, alpha = 0.05,
                                     clusters = 3, nsim = 1e4) {
  
  # alpha and beta of the betabinomial distribution (avoiding name clashes)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  
  # sample random cluster fractions from a beta distribution
  probs <- rbeta(n = nsim * clusters,
                 shape1 = a,
                 shape2 = b)
  
  # divide up the sample size between these clusters
  cluster_sample_sizes <- rep(n %/% clusters, clusters)
  which_get_extra <- seq_len(n %% clusters)
  cluster_sample_sizes[which_get_extra] <- cluster_sample_sizes[which_get_extra] + 1 
  
  # repeat these to match the probs (need 'each', so it aligns with matrix
  # dimension setting later)
  sample_sizes <- rep(cluster_sample_sizes,
                      each = nsim)
  
  # simulate successes in each cluster/sim
  sims <- rbinom(n = nsim * clusters,
                 size = sample_sizes,
                 prob = probs)
  sims <- matrix(sims,
                 nrow = nsim,
                 ncol = clusters)
  
  # sims of numbers of successes overall, and a Monte Carlo estimate of the
  # standard deviation of the estimates
  successes <- rowSums(sims)
  sd_prop <- sd(successes / n)
  
  # compute Z score (for a two-tailed test at this significance level)
  z <- qnorm(1 - alpha / 2)
  
  # moe is z score times the standard deviation
  z * sd_prop
  
}

# given a SpatRaster representing a space-time cube, pad with additional years
# back to an earlier baseline, repeating the earliest value in the cube. cube
# must have contiguous years, each layer correctly named using the naming
# format: name_year, e.g.: irs_2001 or nets_2010. This will return a SpatRaster
# with earlier years added, using the same naming convention, repeating the
# earliest year in the SpatRaster
pre_pad_cube <- function(cube, baseline_year = 1995) {
  # work out the naming system and first year
  first_layer_name <- names(cube)[1]
  earliest_year <- str_sub(first_layer_name, start = -4L)
  prefix <- str_remove(first_layer_name, earliest_year)
  n_years_pad <- as.numeric(earliest_year) - baseline_year
  
  # check the baseline year
  if (!(n_years_pad > 0)) {
    warning("baseline_year must be earlier than the first year in the raster",
            call. = FALSE)
    return(cube)
  }
  
  # make some padding
  pad_layer <- cube[[1]]
  padding <- replicate(n_years_pad, pad_layer, simplify = FALSE) %>%
    do.call(c, .)
  # name the padding years
  years_padding <- as.numeric(earliest_year) - rev(seq_len(n_years_pad))
  names(padding) <- paste0(prefix, years_padding)
  
  # prepend and return
  c(padding, cube)
  
}


# given a SpatRaster representing a space-time cube, pad with additional years
# forward to a later end year, repeating the last value in the cube. cube must
# have contiguous years, each layer correctly named using the naming format:
# name_year, e.g.: irs_2001 or nets_2010. This will return a SpatRaster with
# later years added, using the same naming convention, repeating the last year
# in the SpatRaster
post_pad_cube <- function(cube, end_year = 2030) {
  n_layers_orig <- terra::nlyr(cube)
  # work out the naming system and first year
  last_layer_name <- names(cube)[n_layers_orig]
  latest_year <- str_sub(last_layer_name, start = -4L)
  prefix <- str_remove(last_layer_name, latest_year)
  n_years_pad <- end_year - as.numeric(latest_year)
  
  # check the end year
  if (!(n_years_pad > 0)) {
    warning("end_year must be later than the last year in the raster",
            call. = FALSE)
    return(cube)
  }
  
  # make some padding
  pad_layer <- cube[[n_layers_orig]]
  padding <- replicate(n_years_pad, pad_layer, simplify = FALSE) %>%
    do.call(c, .)
  # name the padding years
  years_padding <- as.numeric(latest_year) + seq_len(n_years_pad)
  names(padding) <- paste0(prefix, years_padding)
  
  # append and return
  c(cube, padding)
  
}

# read in the UN statistical division's geoscheme and turn into a lookup for
# African region.

# why the ever-loving fuck is this semi-colon delimited?
country_region_lookup <- function() {
  read_delim("data/raw/UNSD — Methodology.csv",
             delim = ";",
             col_select = any_of(c("Country or Area",
                                   "Region Name",
                                   "Sub-region Name",
                                   "Intermediate Region Name")),
             col_types = "c") %>%
    filter(`Region Name` == "Africa") %>%
    select(-"Region Name") %>%
    mutate(
      # combine two levels of regions
      region = case_when(
        `Sub-region Name` == "Northern Africa" ~ `Sub-region Name`,
        `Sub-region Name` == "Sub-Saharan Africa" ~ `Intermediate Region Name`
      ),
      country_name = printable_country_name(`Country or Area`)
    ) %>%
    select(
      country_name,
      region
    )
}

# our printable versions of some UNSD country names
printable_country_name <- function(country) {
  case_when(
    country == "United Republic of Tanzania" ~ "Tanzania",
    country == "Democratic Republic of the Congo" ~ "DR Congo",
    country == "Central African Republic" ~ "CAR",
    country == "Sao Tome and Principe" ~ "Sao Tome & Principe",
    .default = country
  )
}


# the nine modelled insecticide types in plotting order (by class, then name),
# and the colours the main figures use for them
insecticides_plot_order <- c("Alpha-cypermethrin",
                             "Deltamethrin",
                             "Lambda-cyhalothrin",
                             "Permethrin",
                             "Fenitrothion",
                             "Malathion",
                             "Pirimiphos-methyl",
                             "DDT",
                             "Bendiocarb")
insecticide_colours <- function() {
  setNames(rev(scales::hue_pal()(length(insecticides_plot_order))),
           insecticides_plot_order)
}
