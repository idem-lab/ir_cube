# The exposure w of a conventionally treated ITN (cITN) relative to an LLIN,
# for the pyrethroid-only net use covariate (#26): the ratio of the median
# durations over which the two net types meet the WHO efficacy criteria.
#
# Source: Tungu P, Magesa S, Maxwell C, Malima R, Masue D, Sudi W, Myamba J,
# Pigeon O, Rowland M (2016). Evaluation of alpha-cypermethrin long-lasting
# insecticidal nets (Interceptor LN) against conventionally treated nets in
# households in north-eastern Tanzania. Parasites & Vectors 9:189.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4831182/
# Interceptor LN and conventionally treated nets (CTN, alpha-cypermethrin),
# three years of household use.
#
# The model uses w = 0.25 (main) and w = 0.47 (sensitivity) as fixed values
# (selection_design() in R/model_covariates.R; R/prep_net_use_pyrethroid.R).
#
#   Rscript R/net_type_weight.R
# writes outputs/net_type_weight.csv.

# nets meeting the WHO efficacy criteria (knockdown or mortality) by months of
# use; all nets are assumed to pass when new
efficacy <- data.frame(
  type = rep(c("LLIN", "CTN"), each = 4),
  months = rep(c(0, 12, 24, 36), 2),
  pass = c(30, 29, 27, 26,
           30, 16, 4, 0),
  n = c(30, 30, 30, 30,
        30, 27, 30, 30)
)
efficacy$years <- efficacy$months / 12

# alpha-cypermethrin content, mg/m2, by months of use
content <- data.frame(
  type = rep(c("LLIN", "CTN"), each = 3),
  months = rep(c(0, 12, 36), 2),
  mg_m2 = c(204, 117, 42,
            32, 9.6, 1.3)
)
content$years <- content$months / 12

# median effective duration, years: where a binomial GLM of passing on years
# of use crosses 0.5, -a / b
median_duration <- function(data) {
  fit <- glm(cbind(pass, n - pass) ~ years, family = binomial, data = data)
  unname(-coef(fit)[1] / coef(fit)[2])
}

durations <- function(data) {
  sapply(c(LLIN = "LLIN", CTN = "CTN"), function(type) {
    median_duration(data[data$type == type, ])
  })
}

# main: with the assumed passes at 0 months; the LLIN median extrapolates
# beyond 36 months
with_zero <- durations(efficacy)
# without them
without_zero <- durations(efficacy[efficacy$months > 0, ])

# sensitivity: the ratio of exponential half-lives of insecticide content. The
# main fit is a log-link quasi-Poisson GLM of content on years; the others fit
# the log of content by least squares, use the 0 and 36 month values only, or
# the 0 and 12 month values only
half_lives <- function(data) {
  slope <- c(
    glm_quasipoisson = unname(coef(glm(mg_m2 ~ years,
                                       family = quasipoisson(link = "log"),
                                       data = data))[2]),
    lm_log = unname(coef(lm(log(mg_m2) ~ years, data = data))[2]),
    months_0_36 = diff(log(data$mg_m2[data$months %in% c(0, 36)])) / 3,
    months_0_12 = diff(log(data$mg_m2[data$months %in% c(0, 12)])) / 1
  )
  log(2) / -slope
}
half_life <- sapply(c(LLIN = "LLIN", CTN = "CTN"), function(type) {
  half_lives(content[content$type == type, ])
})

results <- rbind(
  data.frame(method = "efficacy_glm_with_t0", llin_years = with_zero[["LLIN"]],
             citn_years = with_zero[["CTN"]]),
  data.frame(method = "efficacy_glm_without_t0",
             llin_years = without_zero[["LLIN"]],
             citn_years = without_zero[["CTN"]]),
  data.frame(method = paste0("content_half_life_", rownames(half_life)),
             llin_years = half_life[, "LLIN"],
             citn_years = half_life[, "CTN"])
)
results$w <- results$citn_years / results$llin_years
rownames(results) <- NULL

print(results, digits = 3)
cat(sprintf("\nmain: w = %.2f (median durations: cITN %.2f, LLIN %.2f years)\n",
            results$w[1], results$citn_years[1], results$llin_years[1]))
cat(sprintf("sensitivity: w = %.2f (content half-lives, range %.2f-%.2f)\n",
            results$w[3], min(results$w[-(1:2)]), max(results$w[-(1:2)])))

dir.create("outputs", showWarnings = FALSE)
write.csv(results, "outputs/net_type_weight.csv", row.names = FALSE)
