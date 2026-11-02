# The change panels of the main cross-validation figure (#36, D and E;
# fig_variance_explained.R): skill against no change, for the change in
# mortality between a pair of single bioassays at one pixel.
#
# The target. At each forecasting origin (cuts 2014 and 2018), every pair of
# one assay in the five-year window before the cut (in the fold's training
# data) and one in the five-year holdout after it, at the same pixel and of the
# same insecticide, at least min_gap years apart. Its change, and every
# prediction of it, is divided by the years between its two assays: a rate per
# year. Shorter gaps mostly measure assay noise. Each pixel-insecticide counts
# once (its pairs share its weight), whichever origins it appears at.
#
# Scored as 100 (1 - sum w (P - A)^2 / sum w A^2), the skill against
# predicting no change (R/validation_skill.R). The ceiling is the expected
# score of the true change in population mortality given the pair's assay
# noise; the local bar the expected score of an independent second pair at the
# same pixel.
#
#   D  in sample  the full-data fit's maps at both assays' pixel-years
#                 (outputs/full_fit_maps_at_assays.csv), which saw both
#   E  forecast   each fold's maps: at the before assay, fitted to it; at the
#                 holdout assay, forecast
#
# A pair of a 2013 and a 2018 assay enters at both origins: 2013 is in both
# before windows and 2018 in both holdouts. In E its two copies are different
# forecasts; in D they are one prediction, so D scores each such pair once.
# Otherwise D and E share their pairs, so the step from D to E changes only
# whether the model saw the local data. E has no bioassay-based bar: no
# bioassay exists yet in a forecast year.
#
# Inputs: the dynamical forecasting folds (outputs/cv_draws, for the
# before-window predictions p_draws_before), outputs/cv_scores.csv (the holdout
# maps) and outputs/full_fit_maps_at_assays.csv. The two-stage model's map at
# the before window is read from its fold (map_before, in the order of the
# dynamical fold's before_df).
#
# Replaces the change between pooled window means scored here before (#28).
#
# Writes outputs/cv_change_pairs.csv (one row per target pair) and
# outputs/cv_change_skill.csv (the bars).

source("R/validation_skill.R")

set.seed(2026 - 10 - 6)

cuts <- c("2014", "2018")
min_gap <- 3

rho_source <- rho_lookup()
scores <- require_map_columns(read.csv("outputs/cv_scores.csv",
                                   colClasses = c(fold = "character")))

record_columns <- function(data) {
  data %>%
    transmute(cell, insecticide_type, year, died, mosquito_number)
}

# the assays of one origin's two windows, with each model's map at them
window_assays <- function(cut) {

  fold_file <- function(model) {
    file.path(draws_dir, sprintf("%s__temporal_forecasting__%s.rds", model,
                                 cut))
  }
  dynamical <- readRDS(fold_file("dynamical"))
  if (is.null(dynamical$p_draws_before)) {
    stop(basename(fold_file("dynamical")), " carries no before-window ",
         "predictions; it predates the change score and has to be refitted")
  }
  stopifnot(ncol(dynamical$p_draws_before) == nrow(dynamical$before_df))
  before <- record_columns(dynamical$before_df) %>%
    mutate(dynamical = colMeans(dynamical$p_draws_before))

  # the holdout: the maps as scored (#32), on the fold's held-out records
  holdout <- record_columns(dynamical$test_df)
  for (model in c("dynamical", "two_stage")) {
    rows <- scores %>%
      filter(experiment == paste0("temporal_forecasting_", cut),
             .data$model == !!model)
    stopifnot(identical(record_key(rows), record_key(dynamical$test_df)))
    holdout[[model]] <- rows$map
  }
  rm(dynamical)
  invisible(gc())

  two_stage <- readRDS(fold_file("two_stage"))
  stopifnot(length(two_stage$map_before) == nrow(before))
  before$two_stage <- two_stage$map_before

  list(before = before, holdout = holdout)
}

windows <- lapply(setNames(cuts, cuts), window_assays)

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

# D's copy of each pair: the first origin's. Identical assays (#31) make the
# values ambiguous, so copies are matched by the pair's assays and the rank
# among identical pairs at the origin
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

n_all_gaps <- c(D = sum(pairs$in_sample), E = nrow(pairs))
pairs <- pairs %>% filter(gap >= min_gap)

# the data floor at rho-hat; an assay of one mosquito carries none
rho_hat <- setNames(rho_source$table$rho, rho_source$table$key)
pairs$floor <- pair_floor(pairs, rho_hat)
pairs <- pairs %>% filter(!is.na(floor))

# predicted rates: the fold's maps (E), and the full-data fit's (D). D and E
# have to share their pairs: in one run the full fit and the folds are made
# from the same data, so every pair has a full-fit map
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
  mutate(fold_dynamical = (dynamical_a - dynamical_b) / gap,
         fold_two_stage = (two_stage_a - two_stage_b) / gap,
         full_dynamical = (full_a$dynamical - full_b$dynamical) / gap,
         full_two_stage = (full_a$two_stage - full_b$two_stage) / gap)

# each pixel-insecticide counts once, in E (w) and in D (w_in_sample, 0 for
# the second copy of a pair)
pairs <- pairs %>%
  group_by(cell, insecticide_type) %>%
  mutate(w = 1 / n(),
         w_in_sample = in_sample / sum(in_sample)) %>%
  ungroup()

pair_counts <- function(p, n_all) {
  tibble(n_pairs = nrow(p),
         pixel_insecticides = n_distinct(paste(p$cell, p$insecticide_type)),
         pixels = n_distinct(p$cell),
         median_gap = median(p$gap),
         share_of_all_pairs = nrow(p) / n_all)
}
counts <- bind_rows(
  D = pair_counts(filter(pairs, in_sample), n_all_gaps[["D"]]),
  E = pair_counts(pairs, n_all_gaps[["E"]]), .id = "panel")
for (panel in c("D", "E")) {
  with(counts[counts$panel == panel, ], cat(sprintf(paste(
    "%s: %i pairs of single assays %i+ years apart (%.0f%% of all pairs) at",
    "%i pixel-insecticides in %i pixels; median gap %g years\n"),
    panel, n_pairs, min_gap, 100 * share_of_all_pairs, pixel_insecticides,
    pixels, median_gap)))
}

write.csv(pairs %>% select(cut, cell, insecticide_type, year_b,
                           year_a, gap, died_b, mosquito_number_b,
                           died_a, mosquito_number_a,
                           rate, floor, w, in_sample, w_in_sample,
                           starts_with("fold_"), starts_with("full_")),
          "outputs/cv_change_pairs.csv", row.names = FALSE)


# the bars --------------------------------------------------------------------------

# D and E on the same resamples of pixels, each with its own weights
in_sample <- change_statistic(c("full_dynamical", "full_two_stage"),
                              weight = "w_in_sample")
forecast <- change_statistic(c("fold_dynamical", "fold_two_stage"))
both_panels <- function(pairs, rho) {
  d <- in_sample(pairs, rho)
  e <- forecast(pairs, rho)[c("ceiling", "fold_dynamical", "fold_two_stage")]
  c(setNames(d, paste0("D|", sub("^full_", "", names(d)))),
    setNames(e, paste0("E|", sub("^fold_", "", names(e)))))
}
change <- bootstrap_summary(pairs, both_panels, rho_replicates()) %>%
  separate(bar, c("panel", "bar"), sep = "\\|") %>%
  mutate(source = ifelse(panel == "D", "full-data fit (outputs/two_stage/maps)",
                         "temporal_forecasting 2014, 2018"),
         .after = panel) %>%
  left_join(counts, by = "panel")

write.csv(change, "outputs/cv_change_skill.csv", row.names = FALSE)

cat("\nskill against no change (%), pairs of single bioassays:\n")
print(as.data.frame(change %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(panel, bar, value) %>%
  pivot_wider(names_from = bar, values_from = value)), row.names = FALSE)

