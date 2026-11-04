# Cached figure data (#70). Some figures draw from the two-stage posterior,
# which takes tens of minutes, then plot in seconds. Their scripts compute the
# plotted summaries once, through cached_figure_data(), and keep them in
# temporary/, so that changes to the plot alone rerun in seconds.
#
# The cache holds a fingerprint of everything the summaries were computed
# from:
#   - the code of the compute function (parsed, so its comments don't count)
#     and its arguments: settings and data go in as arguments, never as
#     globals the function reads;
#   - the R files the script sources, and those they source in turn
#     (r_sources()), whole: any edit to them counts;
#   - the contents of the input files, and of every file in the input
#     directories.
# A cache whose fingerprint matches is read; otherwise the script reports what
# changed, recomputes and saves. Package versions are not part of it.

# the files `files` source by a literal source("R/...") call, recursively,
# with `files`
r_sources <- function(files) {
  found <- character()
  while (length(files) > 0) {
    found <- union(found, files)
    calls <- unlist(lapply(files, function(file) {
      regmatches(readLines(file),
                 gregexpr("(?<=source\\(\")R/[^\"]+\\.R(?=\"\\))",
                          readLines(file), perl = TRUE))
    }))
    files <- setdiff(unique(calls), found)
  }
  sort(found)
}

figure_fingerprint <- function(compute, args, sources, inputs) {
  is_dir <- dir.exists(inputs)
  files <- sort(c(r_sources(sources), inputs[!is_dir],
                  list.files(inputs[is_dir], recursive = TRUE,
                             full.names = TRUE)))
  c(code = rlang::hash(deparse(compute, control = c("keepNA", "keepInteger",
                                                    "niceNames"))),
    args = rlang::hash(args),
    setNames(rlang::hash_file(files), files))
}

# compute's result, from `file` if its fingerprint matches, otherwise computed
# with `args` and saved there. `sources` are the R files the script sources;
# `inputs` the files and directories compute reads
cached_figure_data <- function(file, compute, args = list(),
                               sources = character(), inputs = character()) {
  stopifnot(all(file.exists(c(sources, inputs))))
  fingerprint <- figure_fingerprint(compute, args, sources, inputs)
  if (file.exists(file)) {
    cache <- readRDS(file)
    if (identical(cache$fingerprint, fingerprint)) {
      message("figure data from ", file)
      return(cache$data)
    }
    keys <- union(names(fingerprint), names(cache$fingerprint))
    changed <- keys[fingerprint[keys] %in% NA |
                      cache$fingerprint[keys] %in% NA |
                      fingerprint[keys] != cache$fingerprint[keys]]
    message(file, " is stale (changed: ",
            paste(stats::na.omit(changed), collapse = ", "),
            "); recomputing")
  } else {
    message("no ", file, "; computing")
  }
  data <- do.call(compute, args)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(fingerprint = fingerprint, data = data), file)
  data
}

# the inputs of figure data drawn from the two-stage posterior of `types`
# (R/two_stage_predictions.R): their saved fits and paired dynamical draws,
# the external overdispersion, and the prepared data (covariates, masks)
two_stage_figure_inputs <- function(types, maps_dir = "outputs/two_stage/maps") {
  c(file.path(maps_dir, types, "fit.rds"),
    file.path(maps_dir, types, "dynamical.rds"),
    "outputs/bioassay_rho_hierarchical.csv",
    "data/clean")
}
