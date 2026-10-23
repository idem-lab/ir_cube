# The full-data fit's maps at the pixel-year of every assay it was fitted to:
# the in-sample predictions of the main cross-validation figure (#36, panels A
# and D, fig_variance_explained.R), read off the published yearly rasters:
#
#   dynamical  outputs/ir_maps (R/predict.R): posterior mean of ilogit(m)
#   two_stage  outputs/two_stage/ir_maps (R/two_stage_maps.R): posterior mean
#              of ilogit(m + omega + xi), without u, p and sy
#
#   Rscript R/full_fit_maps_at_assays.R
#
# Run after both, on the same full fit. Writes
# outputs/full_fit_maps_at_assays.csv.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")

# the modelled assays, as R/fit_model.R builds them
invisible(list2env(modelled_bioassays(1995, 2024), environment()))

maps <- df %>%
  distinct(insecticide_type, year_start, cell) %>%
  group_by(insecticide_type, year_start) %>%
  group_modify(function(keys, group) {
    at <- function(dir) {
      file <- file.path(dir, group$insecticide_type,
                        sprintf("ir_%i_susceptibility.tif", group$year_start))
      terra::extract(rast(file), keys$cell)[, 1]
    }
    mutate(keys,
           dynamical = at("outputs/ir_maps"),
           two_stage = at("outputs/two_stage/ir_maps"))
  }) %>%
  ungroup() %>%
  select(cell, year_start, insecticide_type, dynamical, two_stage)
stopifnot(!anyNA(maps$dynamical), !anyNA(maps$two_stage))

write.csv(maps, "outputs/full_fit_maps_at_assays.csv", row.names = FALSE)
cat(sprintf("wrote outputs/full_fit_maps_at_assays.csv: %i pixel-years\n",
            nrow(maps)))
