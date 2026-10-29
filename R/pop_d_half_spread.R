# The density d_half of the population transform at which the population
# covariate is most spread across the modelled bioassays, so that the data
# carry the most information about population (#37).
#
#   Rscript R/pop_d_half_spread.R
#
# At each modelled bioassay (modelled_bioassays(), R/bioassay_subset.R) the
# encounter transform 1 - exp(-d log(2) / d_half) of the population density d
# at its pixel and year (pop_density_matrix(), R/model_covariates.R), and its
# sd over the assays, for d_half on a grid. Prints the d_half where the sd is
# largest, and the sd at d_half 50, 200 and 270.

suppressMessages({
  library(dplyr)
  library(stringr)
  library(terra)
})
source("R/functions.R")
source("R/model_covariates.R")
source("R/bioassay_subset.R")
baseline_year <- 1995
modelled <- modelled_bioassays(baseline_year)
d <- pop_density_matrix(modelled$unique_cells, baseline_year,
                        max(modelled$df$year))
density <- d[cbind(modelled$df$cell_id, modelled$df$year_id)]

d_half <- sort(unique(c(round(exp(seq(log(5), log(1000), length.out = 60))),
                        50, 200, 270)))
spread <- vapply(d_half, function(h) sd(-expm1(-density * log(2) / h)),
                 numeric(1))
cat(sprintf("sd largest at d_half %g: %.3f\n", d_half[which.max(spread)],
            max(spread)))
print(data.frame(d_half, sd = spread)[d_half %in% c(50, 200, 270), ],
      row.names = FALSE, digits = 3)
