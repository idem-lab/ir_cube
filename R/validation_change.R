# The in-sample change figure (#36, fig_variance_explained.R): skill against
# no change, for the change in mortality between a pair of single bioassays at
# one pixel, predicted by the full-data fit's maps.
#
# The target. At each forecasting origin (cuts 2014 and 2018), every pair of
# one assay in the five-year window before the cut and one in the five-year
# window after it (the fold definitions, R/validation_folds.R), at the same
# pixel and of the same insecticide, at least min_gap years apart. Its change,
# and every prediction of it, is divided by the years between its two assays:
# a rate per year. Shorter gaps mostly measure assay noise. A pair of a 2013
# and a 2018 assay enters at both origins, and is scored once. Each
# pixel-insecticide counts once (its pairs share its weight).
#
# Scored as 100 (1 - sum w (P - A)^2 / sum w A^2), the skill against
# predicting no change (R/validation_skill.R). The ceiling is the expected
# score of the true change in population mortality given the pair's assay
# noise; the local bars the expected score of the mean of 1, 2 or 3
# independent pairs at the same pixel.
#
# Reads outputs/full_fit_maps_at_assays.csv. Writes outputs/cv_change_pairs.csv
# (one row per target pair) and outputs/cv_change_skill.csv (the bars).

source("R/validation_skill.R")
source("R/validation_folds.R")

set.seed(2026 - 10 - 6)

cuts <- c("2014", "2018")
min_gap <- 3

rho_source <- rho_lookup()

record_columns <- function(data) {
  data %>%
    transmute(cell, insecticide_type, year, died, mosquito_number)
}

# the assays of one origin's two windows
windows <- lapply(setNames(cuts, cuts), function(cut) {
  fold <- temporal_forecasting_folds[[cut]]
  list(before = record_columns(fold$before),
       holdout = record_columns(fold$test))
})

# every pair of a before and a holdout assay at one pixel, of one insecticide
pairs <- bind_rows(lapply(cuts, function(cut) {
  inner_join(windows[[cut]]$before, windows[[cut]]$holdout,
             by = c("cell", "insecticide_type"), suffix = c("_b", "_a"),
             relationship = "many-to-many") %>%
    mutate(cut = as.integer(cut), .before = 1)
})) %>%
  mutate(gap = year_a - year_b,
         rate = (died_a / mosquito_number_a - died_b / mosquito_number_b) /
           gap)

# a pair entering at both origins is scored once, as its first origin's copy.
# Identical assays (#31) make the values ambiguous, so copies are matched by
# the pair's assays and the rank among identical pairs at the origin
pair_columns <- c("cell", "insecticide_type", "year_b", "died_b",
                  "mosquito_number_b", "year_a", "died_a",
                  "mosquito_number_a")
pairs <- pairs %>%
  group_by(cut, across(all_of(pair_columns))) %>%
  mutate(copy = row_number()) %>%
  group_by(across(all_of(c(pair_columns, "copy")))) %>%
  mutate(in_sample = row_number() == 1) %>%
  ungroup() %>%
  select(-copy)

n_all_gaps <- sum(pairs$in_sample)
pairs <- pairs %>% filter(gap >= min_gap)

# the data floor at rho-hat; an assay of one mosquito carries none
rho_hat <- setNames(rho_source$table$rho, rho_source$table$key)
pairs$floor <- pair_floor(pairs, rho_hat)
pairs <- pairs %>% filter(!is.na(floor))

# predicted rates: the full-data fit's maps. In one run the full fit and the
# folds are made from the same data, so every pair has a full-fit map
full_fit_file <- "outputs/full_fit_maps_at_assays.csv"
full_fit <- read.csv(full_fit_file)
at <- function(year) {
  tibble(cell = pairs$cell, insecticide_type = pairs$insecticide_type,
         year = year) %>%
    left_join(full_fit, by = c("cell", "year", "insecticide_type"),
              relationship = "many-to-one")
}
full_b <- at(pairs$year_b)
full_a <- at(pairs$year_a)
if (anyNA(full_b) || anyNA(full_a)) {
  stop("pairs have no full-fit map at one of their assays' pixel-years: ",
       full_fit_file, " is from different data than the folds")
}
pairs <- pairs %>%
  mutate(full_dynamical = (full_a$dynamical - full_b$dynamical) / gap,
         full_two_stage = (full_a$two_stage - full_b$two_stage) / gap)

# each pixel-insecticide counts once (w_in_sample, 0 for the second copy of a
# pair)
pairs <- pairs %>%
  group_by(cell, insecticide_type) %>%
  mutate(w_in_sample = in_sample / sum(in_sample)) %>%
  ungroup()

pair_counts <- function(p, n_all) {
  tibble(n_pairs = nrow(p),
         pixel_insecticides = n_distinct(paste(p$cell, p$insecticide_type)),
         pixels = n_distinct(p$cell),
         median_gap = median(p$gap),
         share_of_all_pairs = nrow(p) / n_all)
}
counts <- pair_counts(filter(pairs, in_sample), n_all_gaps) %>%
  mutate(panel = "in_sample", .before = 1)
with(counts, cat(sprintf(paste(
  "%i pairs of single assays %i+ years apart (%.0f%% of all pairs) at",
  "%i pixel-insecticides in %i pixels; median gap %g years\n"),
  n_pairs, min_gap, 100 * share_of_all_pairs, pixel_insecticides,
  pixels, median_gap)))

write.csv(pairs %>% select(cut, cell, insecticide_type, year_b,
                           year_a, gap, died_b, mosquito_number_b,
                           died_a, mosquito_number_a,
                           rate, floor, in_sample, w_in_sample,
                           starts_with("full_")),
          "outputs/cv_change_pairs.csv", row.names = FALSE)


# the bars --------------------------------------------------------------------------

# in sample, also the mean of 2 or 3 local pairs (the pooled figure)
in_sample <- change_statistic(c("full_dynamical", "full_two_stage"),
                              weight = "w_in_sample", pools = 1:3)
change <- bootstrap_summary(pairs, function(pairs, rho) {
  d <- in_sample(pairs, rho)
  setNames(d, sub("^full_", "", names(d)))
}, rho_replicates()) %>%
  mutate(panel = "in_sample", source = "full-data fit (outputs/two_stage/maps)",
         .before = 1) %>%
  left_join(counts, by = "panel")

write.csv(change, "outputs/cv_change_skill.csv", row.names = FALSE)

cat("\nskill against no change (%), pairs of single bioassays:\n")
print(as.data.frame(change %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(panel, bar, value) %>%
  pivot_wider(names_from = bar, values_from = value)), row.names = FALSE)

