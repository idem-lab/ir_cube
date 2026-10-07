# The full-data fit's maps at the pixel-year of every assay it was fitted to:
# the in-sample predictions of the main cross-validation figure (#36, panels A
# and D, fig_variance_explained.R).
#
#   Rscript R/full_fit_maps_at_assays.R [bioassay data file]
#
# The data file defaults to data/clean/all_gambiae_complex_data.RDS. It must be
# the one the full fit was made from: the saved draws are in the order of the
# assays it gives, and the script stops if they are not (the two-stage fit
# stores its assays' pixel-years in order, which are checked).
#
# Per insecticide type, from what R/two_stage_maps.R saved in
# outputs/two_stage/maps/<type>/: dynamical.rds (the full dynamical fit's 2000
# posterior logit draws at the type's assays, in the order of df) and fit.rds
# (the two-stage correction fitted to all the type's assays). For each distinct
# pixel-year with an assay:
#
#   dynamical  posterior mean of ilogit(m), from all 2000 draws: the
#              dynamical_mortality map at that pixel-year
#   two_stage  posterior mean of ilogit(m + omega + xi), at the pixel centre,
#              from n_draws paired draws (dynamical draw d with a latent draw
#              shifted for it, the cut posterior), as two_stage_maps.R maps
#              it. No u or p: the two_stage_mortality map at that pixel-year
#
# The map rasters hold only the panel years (2000, 2005, ...), so the values
# are recomputed here rather than read off them. The Monte Carlo SE of a
# two-stage value is its posterior SD / sqrt(n_draws), under 0.5 percentage
# points where the SD is 10.
#
# Peak memory ~2 GB (one type's fit at a time); ~15 min. Writes
# outputs/full_fit_maps_at_assays.csv.

arguments <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(arguments) > 0) arguments[1] else
  "data/clean/all_gambiae_complex_data.RDS"

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/two_stage_helpers.R")
source("R/two_stage_correction.R")

# the modelled assays, as R/validation_folds.R and R/fit_model.R build them
mask <- rast("data/clean/raster_mask.tif")
df <- subset_modelled_bioassays(readRDS(data_file), mask,
                                baseline_year = 1995, final_data_year = 2024)
types <- unique(df$insecticide_type)

maps_dir <- "outputs/two_stage/maps"
n_draws <- 400
batch_size <- 100

maps_for_type <- function(type) {

  rows <- which(df$insecticide_type == type)
  dynamical <- readRDS(file.path(maps_dir, type, "dynamical.rds"))
  fit <- readRDS(file.path(maps_dir, type, "fit.rds"))

  # the saved draws are in the order of the df they were computed from; that
  # df has to be this one, row for row, or the draws land on the wrong assays.
  # The fit's u levels are its assays' pixel-years in order of appearance
  stopifnot(
    ncol(dynamical$logit_train) == length(rows),
    fit$n_obs == length(rows),
    max(abs(fit$m_ref - colMeans(dynamical$logit_train))) < 1e-12,
    identical(fit$iid_levels$u,
              unique(paste(df$cell[rows], df$year_start[rows])))
  )

  # one row per pixel-year: m depends on the pixel-year alone, so any of its
  # assays' columns will do
  first <- rows[!duplicated(paste(df$cell[rows], df$year_start[rows]))]
  column <- match(first, rows)
  xy <- terra::xyFromCell(mask, df$cell[first])
  new <- tibble(lon = xy[, 1], lat = xy[, 2], year = df$year_start[first])

  draws <- thin_draws(matrix(seq_len(nrow(dynamical$logit_train))),
                      n_draws)[, 1]
  batches <- split(draws, ceiling(seq_along(draws) / batch_size))
  set.seed(string_seed(paste("full fit maps at assays", type)))
  total <- numeric(length(first))
  for (batch in batches) {
    m <- dynamical$logit_train[batch, column, drop = FALSE]
    fields <- correction_node_draws(
      fit, unique(new$year), length(batch),
      m_draws_train = dynamical$logit_train[batch, , drop = FALSE])
    total <- total + rowSums(plogis(t(m) + project_correction(fit, fields,
                                                              new)))
  }

  report("%-18s %5i assays, %5i pixel-years; peak memory %.1f GB", type,
         length(rows), length(first), peak_memory_gb())
  tibble(cell = df$cell[first], year_start = df$year_start[first],
         insecticide_type = type,
         dynamical = colMeans(plogis(dynamical$logit_train[, column,
                                                          drop = FALSE])),
         two_stage = total / length(draws))
}

maps <- bind_rows(lapply(types, maps_for_type))
write.csv(maps, "outputs/full_fit_maps_at_assays.csv", row.names = FALSE)
report("wrote outputs/full_fit_maps_at_assays.csv: %i pixel-years",
       nrow(maps))
