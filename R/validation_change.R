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
# D and E share their pairs, so the step from D to E changes only whether the
# model saw the local data. E has no bioassay-based bar: no bioassay exists
# yet in a forecast year.
#
# Inputs: the dynamical forecasting folds (outputs/cv_draws, for the
# before-window predictions p_draws_before), outputs/cv_scores.csv (the holdout
# maps) and outputs/full_fit_maps_at_assays.csv. The two-stage model's map at
# the before window is read from its fold (map_before, aligned with before_df)
# where the fold has it; folds fitted before that was added do not, and its E
# bar is then a placeholder.
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
    transmute(cell, insecticide_type, year_start, died, mosquito_number)
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
  before$two_stage <- if (is.null(two_stage$map_before)) NA_real_ else {
    stopifnot(identical(record_key(two_stage$before_df), record_key(before)))
    two_stage$map_before
  }

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
  mutate(gap = year_start_a - year_start_b,
         rate = (died_a / mosquito_number_a - died_b / mosquito_number_b) /
           gap)

n_all_gaps <- nrow(pairs)
pairs <- pairs %>% filter(gap >= min_gap)

# the data floor at rho-hat; an assay of one mosquito carries none
rho_hat <- setNames(rho_source$table$rho, rho_source$table$key)
pairs$floor <- pair_floor(pairs, rho_hat)
pairs <- pairs %>% filter(!is.na(floor))

# predicted rates: the fold's maps (E), and the full-data fit's (D)
pairs <- pairs %>%
  mutate(fold_dynamical = (dynamical_a - dynamical_b) / gap,
         fold_two_stage = (two_stage_a - two_stage_b) / gap)

full_fit_file <- "outputs/full_fit_maps_at_assays.csv"
if (file.exists(full_fit_file)) {
  full_fit <- read.csv(full_fit_file)
  at <- function(year) {
    tibble(cell = pairs$cell, insecticide_type = pairs$insecticide_type,
           year_start = year) %>%
      left_join(full_fit, by = c("cell", "year_start", "insecticide_type"),
                relationship = "many-to-one")
  }
  full_b <- at(pairs$year_start_b)
  full_a <- at(pairs$year_start_a)
  pairs <- pairs %>%
    mutate(full_dynamical = (full_a$dynamical - full_b$dynamical) / gap,
           full_two_stage = (full_a$two_stage - full_b$two_stage) / gap)
  # D and E have to share their pairs. In one run the full fit and the folds
  # are made from the same data, and nothing is dropped here
  missing <- is.na(pairs$full_dynamical) | is.na(pairs$full_two_stage)
  if (any(missing)) {
    warning(sum(missing), " of ", nrow(pairs), " pairs have no full-fit map ",
            "at one of their assays' pixel-years and are dropped from D and E ",
            "alike: ", full_fit_file, " is from different data than the folds")
    pairs <- pairs[!missing, ]
  }
} else {
  warning("no ", full_fit_file, "; panel D's model bars are placeholders")
  pairs <- pairs %>% mutate(full_dynamical = NA_real_,
                            full_two_stage = NA_real_)
}

# each pixel-insecticide counts once
pairs <- pairs %>%
  group_by(cell, insecticide_type) %>%
  mutate(w = 1 / n()) %>%
  ungroup()

counts <- with(pairs, tibble(
  n_pairs = length(gap),
  pixel_insecticides = n_distinct(paste(cell, insecticide_type)),
  pixels = n_distinct(cell),
  median_gap = median(gap),
  share_of_all_pairs = length(gap) / n_all_gaps))
cat(sprintf(paste("%i pairs of single assays %i+ years apart (%.0f%% of all",
                  "pairs) at %i pixel-insecticides in %i pixels; median gap",
                  "%g years\n"),
            counts$n_pairs, min_gap, 100 * counts$share_of_all_pairs,
            counts$pixel_insecticides, counts$pixels, counts$median_gap))

write.csv(pairs %>% select(cut, cell, insecticide_type, year_start_b,
                           year_start_a, gap, died_b, mosquito_number_b,
                           died_a, mosquito_number_a,
                           rate, floor, w, starts_with("fold_"),
                           starts_with("full_")),
          "outputs/cv_change_pairs.csv", row.names = FALSE)


# the bars --------------------------------------------------------------------------

predictions <- c("full_dynamical", "full_two_stage", "fold_dynamical",
                 "fold_two_stage")
summary <- bootstrap_summary(pairs, change_statistic(predictions),
                             rho_replicates())

two_stage_exact <- !anyNA(pairs$fold_two_stage) &&
  all(scores$map_exact[grepl("^temporal_forecasting", scores$experiment) &
                         scores$model == "two_stage"])

change <- bind_rows(
  summary %>%
    filter(bar %in% c("ceiling", "local", "full_dynamical",
                      "full_two_stage")) %>%
    mutate(panel = "D", source = "full-data fit (outputs/two_stage/maps)"),
  summary %>%
    filter(bar %in% c("ceiling", "fold_dynamical", "fold_two_stage")) %>%
    mutate(panel = "E", source = "temporal_forecasting 2014, 2018")
) %>%
  mutate(map_exact = ifelse(bar == "fold_two_stage", two_stage_exact, NA),
         bar = sub("^(full|fold)_", "", bar),
         # a bar with no estimate is a placeholder, and drawn as one
         source = ifelse(is.na(estimate), "placeholder", source)) %>%
  relocate(panel, source, .before = 1) %>%
  bind_cols(counts)

write.csv(change, "outputs/cv_change_skill.csv", row.names = FALSE)

cat("\nskill against no change (%), pairs of single bioassays:\n")
print(as.data.frame(change %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(panel, bar, value) %>%
  pivot_wider(names_from = bar, values_from = value)), row.names = FALSE)


# the local bar against actual second pairs ---------------------------------------
