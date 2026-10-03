# Start greta's python, for scripts that use greta.
#
# Scripts that build, sample or calculate() through a greta model call
# start_greta() before sourcing R/packages.R. Two ordering constraints
# (doc/cv_run_plan.md, section 2):
# - TensorFlow will not change its thread count after initialisation, so
#   `threads` is set before python comes up.
# - python must initialise before terra and sf are attached: they load the
#   system XML libraries, against which the conda environment's pyexpat is then
#   resolved, and tensorflow_probability fails to import.
# Plain-R scripts do not call it, and so never start python.
#
# Definitions only; sourcing this file starts nothing.

start_greta <- function(threads = NULL) {
  suppressMessages(library(greta))
  if (!is.null(threads)) {
    tensorflow_module <- reticulate::import("tensorflow")
    tensorflow_module$config$threading$set_intra_op_parallelism_threads(
      as.integer(threads))
    tensorflow_module$config$threading$set_inter_op_parallelism_threads(
      as.integer(threads))
  }
  check_greta_fill()
}

# greta must fill subassignments into greta arrays column-major, as R does:
# greta >= 43f9c52 (0.6.0.9000, greta-dev/greta#844, fixed in #847). This
# initialises python if it is not up yet. Run once per session.
check_greta_fill <- function() {
  if (isTRUE(getOption("ir_cube.greta_fill_checked"))) {
    return(invisible(TRUE))
  }
  g <- greta::zeros(4, 2)
  g[c(1, 3), ] <- matrix(1:4, 2)
  stopifnot(
    "greta fills subassignments row-major; install greta >= 43f9c52 (greta-dev/greta#847)" =
      identical(unname(greta::calculate(g)[[1]]), rbind(c(1, 3), 0, c(2, 4), 0))
  )
  options(ir_cube.greta_fill_checked = TRUE)
  invisible(TRUE)
}
