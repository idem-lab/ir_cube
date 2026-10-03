# Covariates of the dynamical model that are built the same way for the fits,
# the folds and the maps. Functions only; needs terra. Source from the repo
# root.


# initial-state covariates (#19) -----------------------------------------------

# Static predictors of the logit initial fraction susceptible, which let the
# initial state vary within a country. Each is standardised over every cell of
# the mask, so the fits, folds and maps share them:
#   population    population in 2000, the earliest layer of the population
#                 cube, transformed as design$init_pop says: "log", named
#                 log_pop_2000 (the fits before the refit); otherwise pop_2000,
#                 with no trend (init_pop_layer())
#   all_crops     the "all crops" yield layer of crop_group_scaled.tif, as
#                 log(x + 1e-4). The layer is very skewed (56% zeros, mean
#                 0.003, sd 0.012), so standardised untransformed it reaches
#                 47 sd at the data cells; 1e-4 is about the 10th percentile
#                 of its non-zero values, and on this scale the data cells
#                 span -0.7 to 4.1 sd
init_covariate_names <- function(design = selection_design()) {
  design <- complete_selection_design(design)
  c(if (design$init_pop == "log") "log_pop_2000" else "pop_2000", "all_crops")
}

# The unstandardised population layer of the initial-state covariates
init_pop_layer <- function(design) {
  pop <- rast("data/clean/pop_cube.tif")[["pop_2000"]]
  if (design$init_pop == "log") {
    return(log(pop))
  }
  if (design$init_pop == "raw") {
    return(pop)
  }
  d <- pop_zero_filled(pop) /
    terra::cellSize(rast("data/clean/raster_mask.tif"), unit = "km")
  switch(design$init_pop,
         encounter = 1 - exp(-d * log(2) / design$pop_d_half),
         saturating = d / (d + design$pop_d_half))
}

init_covariate_layers <- function(design = selection_design()) {
  design <- complete_selection_design(design)
  mask <- rast("data/clean/raster_mask.tif")
  pop <- init_pop_layer(design)
  crops <- log(rast("data/clean/crop_group_scaled.tif")[["all crops"]] + 1e-4)
  layers <- terra::mask(c(pop, crops), mask)
  names(layers) <- init_covariate_names(design)
  moments <- terra::global(layers, c("mean", "sd"), na.rm = TRUE)
  (layers - moments$mean) / moments$sd
}

# The standardised initial-state covariates at mask cells `cells`, as a
# cells x covariates matrix with named columns
init_covariate_matrix <- function(cells, design = selection_design(),
                                  layers = init_covariate_layers(design)) {
  x <- as.matrix(terra::extract(layers, cells))
  colnames(x) <- names(layers)
  x
}


# selection covariates (#23) ---------------------------------------------------

# How the selection design matrix is built. Net use and IRS enter as they are,
# with any hinge columns. Population and the crop layers are proxies for
# insecticide use outside vector control; with a trend they enter only as
# products g(t) s(x, t), so that there is no selection from them where g is 0:
# g_dom for population (domestic and public-health use), g_ag for the crops
# (agricultural use).
#   pop          transform of population s(x, t), each on 0-1, with d the
#                population density in people per km2 (pop_density_matrix()):
#                "encounter": 1 - exp(-d log(2) / pop_d_half), the chance of
#                encountering at least one source of exposure when sources
#                are Poisson in number with mean proportional to d;
#                "saturating": d / (d + pop_d_half);
#                "raw": pop_scaled_cube.tif, population per cell min-max
#                scaled; "log": pop_log_scaled_cube.tif, log population
#                min-max scaled (R/prep_rasters.R)
#   pop_d_half   the density at which "encounter" and "saturating" are 0.5
#   init_pop     the transform of the initial-state population covariate
#                (#19), one of those of pop
#   hinges       named list of knots, e.g. list(nets = c(0.18, 0.35)): hinge
#                columns min(x, k) for any of nets, irs, pop (after its
#                transform), knots on the 0-1 scale, each with a positive
#                coefficient, so the effect can saturate but not decrease
#   trend_pop    g_dom(t), multiplying population and its hinges, and
#   trend_crops  g_ag(t), multiplying the crops, each one of: "none";
#                "linear_0_1", linear in the year, 0 in trend_years[1] and 1
#                in trend_years[2]; or a regions x years matrix of g, rownames
#                the regions of country_region_lookup() and colnames every
#                year the design is built for (cell_regions())
#   trend_years  the years where a linear trend is 0 and 1; the first should
#                be the baseline year. After the second it continues on the
#                same line (1.17 in 2030)
#   nets         the net use column: "pyrethroid", use times the share of
#                nets that are pyrethroid-only, conventional ITNs weighted by
#                net_w (#26; R/prep_net_use_pyrethroid.R); "all", all net use
#   net_w        the exposure of a conventional ITN relative to an LLIN
#                (R/net_type_weight.R)
#   net_use_source  the net use layers: "run06", MITN run 06, the run of the
#                net crop by type; "legacy", those of net_use_cube.tif
# Net use, IRS and their hinges have no trend. The settings of the fits before
# the refit are selection_design_untrended().
selection_design <- function(pop = c("encounter", "saturating", "raw",
                                     "log"),
                             pop_d_half = 50,
                             init_pop = pop,
                             hinges = list(),
                             trend_pop = "linear_0_1",
                             trend_crops = "linear_0_1",
                             trend_years = c(1995, 2025),
                             nets = c("pyrethroid", "all"),
                             net_w = 0.25,
                             net_use_source = c("run06", "legacy")) {
  pop <- match.arg(pop)
  nets <- match.arg(nets)
  net_use_source <- match.arg(net_use_source)
  init_pop <- match.arg(init_pop, c("encounter", "saturating", "raw", "log"))
  stopifnot(is.list(hinges),
            all(names(hinges) %in% c("nets", "irs", "pop")),
            all(vapply(hinges, function(k) {
              is.numeric(k) && all(k > 0 & k < 1)
            }, logical(1))),
            is.numeric(pop_d_half), length(pop_d_half) == 1,
            pop_d_half > 0,
            is.numeric(net_w), length(net_w) == 1, net_w >= 0, net_w <= 1,
            is.numeric(trend_years), length(trend_years) == 2,
            trend_years[2] > trend_years[1])
  for (trend in list(trend_pop, trend_crops)) {
    if (is.matrix(trend)) {
      stopifnot(is.numeric(trend), !is.null(rownames(trend)),
                !anyNA(suppressWarnings(as.integer(colnames(trend)))))
    } else {
      stopifnot(length(trend) == 1, trend %in% c("none", "linear_0_1"))
    }
  }
  list(pop = pop, pop_d_half = pop_d_half, init_pop = init_pop,
       hinges = hinges, trend_pop = trend_pop,
       trend_crops = trend_crops, trend_years = trend_years,
       nets = nets, net_w = net_w,
       net_use_source = net_use_source)
}

# The design of the fits before the refit
selection_design_untrended <- function() {
  selection_design(pop = "raw", init_pop = "log", trend_pop = "none",
                   trend_crops = "none", nets = "all",
                   net_use_source = "legacy")
}

# A saved design, completed and checked as selection_design() would build it
# (it ignores settings since removed)
complete_selection_design <- function(design) {
  do.call(selection_design,
          design[intersect(names(design), names(formals(selection_design)))])
}

# A linear trend g(t) from g_start in trend_years[1] to 1 in trend_years[2],
# continuing on the same line, as a regions x years matrix for trend_pop or
# trend_crops: the sensitivity fits with g(1995) = 0.37 (doc/cv_run_plan.md)
linear_trend_matrix <- function(g_start, years = 1995:2030,
                                trend_years = c(1995, 2025)) {
  regions <- unique(na.omit(country_region_lookup()$region))
  g <- g_start + (1 - g_start) * (years - trend_years[1]) / diff(trend_years)
  matrix(g, length(regions), length(years), byrow = TRUE,
         dimnames = list(regions, years))
}

# whether a trend (design$trend_pop or design$trend_crops) is on
has_trend <- function(trend) {
  is.matrix(trend) || trend != "none"
}

# the columns of the selection design matrix, in order: nets, irs, pop, the
# hinge columns, then the crops; with their trend, the pop and population
# hinge columns are named "<name>:g_dom" (e.g. "pop_enc:g_dom",
# "pop_enc_min_0.5:g_dom") and the crops "<name>:g_ag" ("all crops:g_ag")
selection_column_names <- function(design = selection_design()) {
  design <- complete_selection_design(design)
  pop_name <- c(raw = "pop", log = "log_pop",
                encounter = "pop_enc", saturating = "pop_sat")[[design$pop]]
  trended <- function(name, trend, suffix) {
    if (has_trend(trend)) paste0(name, ":", suffix) else name
  }
  hinge_names <- unlist(lapply(names(design$hinges), function(name) {
    if (name == "pop") {
      trended(paste0(pop_name, "_min_", design$hinges[[name]]),
              design$trend_pop, "g_dom")
    } else {
      paste0(name, "_min_", design$hinges[[name]])
    }
  }))
  c("nets", "irs", trended(pop_name, design$trend_pop, "g_dom"), hinge_names,
    trended(selection_crop_names, design$trend_crops, "g_ag"))
}

# the crop columns: the crop group totals, then cotton, vegetables and rice
selection_crop_names <- c("all crops", "cereal crops", "root crops",
                          "pulse crops", "oil crops", "fibre crops",
                          "other crops", "cotton", "vegetables", "rice")

# A cube as a cells x years matrix at mask cells `cells`, years
# baseline_year..end_year, padded back to baseline_year and forward to
# end_year by repeating its first and last layers
read_padded_cube <- function(cube, cells, baseline_year, end_year) {
  if (is.character(cube)) {
    cube <- rast(cube)
  }
  cube <- suppressWarnings(pre_pad_cube(cube, baseline_year))
  cube <- suppressWarnings(post_pad_cube(cube, end_year))
  years <- as.numeric(str_sub(names(cube), start = -4L))
  cube <- cube[[years >= baseline_year & years <= end_year]]
  stopifnot(identical(as.numeric(str_sub(names(cube), start = -4L)),
                      as.numeric(baseline_year:end_year)))
  as.matrix(terra::extract(cube, cells))
}

# Population density, people per km2, at mask cells `cells` for years
# baseline_year..end_year (padded as read_padded_cube()), cells x years.
# prep_rasters.R filled the pixels WorldPop has as empty (0 or NA on land)
# with 1 person per cell, relevelled by a factor close to 1 in each year;
# these are set back to 0. The fill is each layer's most common value near 1
# (51,408 pixels in every year, against at most 3 for any other value there).
pop_density_matrix <- function(cells, baseline_year, end_year) {
  pop <- pop_zero_filled(rast("data/clean/pop_cube.tif"))
  area <- terra::cellSize(rast("data/clean/raster_mask.tif"), unit = "km")
  area_cells <- terra::extract(area, cells)[, 1]
  read_padded_cube(pop, cells, baseline_year, end_year) / area_cells
}

# Layers of pop_cube.tif with the fill for empty pixels set back to 0
pop_zero_filled <- function(pop) {
  fill <- vapply(seq_len(terra::nlyr(pop)), function(i) {
    v <- terra::values(pop[[i]], mat = FALSE)
    v <- v[!is.na(v) & v > 0.99 & v < 1.01]
    u <- unique(v)
    n <- tabulate(match(v, u))
    stopifnot(max(n) > 10000)
    u[which.max(n)]
  }, numeric(1))
  out <- pop * (pop != fill)
  names(out) <- names(pop)
  out
}

# The population column s(x, t) of the design, cells x years
selection_pop_matrix <- function(cells, baseline_year, end_year, design) {
  switch(
    design$pop,
    raw = read_padded_cube("data/clean/pop_scaled_cube.tif", cells,
                           baseline_year, end_year),
    log = {
      file <- "data/clean/pop_log_scaled_cube.tif"
      if (!file.exists(file)) {
        stop(file, " not found; run R/prep_rasters.R")
      }
      read_padded_cube(file, cells, baseline_year, end_year)
    },
    encounter = {
      d <- pop_density_matrix(cells, baseline_year, end_year)
      -expm1(-d * log(2) / design$pop_d_half)
    },
    saturating = {
      d <- pop_density_matrix(cells, baseline_year, end_year)
      d / (d + design$pop_d_half)
    }
  )
}

# The region of each mask cell in `cells`, from the country raster and the
# UNSD lookup as predict.R assigns them; NA outside any country
cell_regions <- function(cells) {
  country <- as.character(terra::extract(
    rast("data/clean/country_raster.tif"), cells)$country_name)
  lookup <- country_region_lookup()
  lookup$region[match(country, lookup$country_name)]
}

# A trend g(t) (design$trend_pop or design$trend_crops) at mask cells `cells`
# for years baseline_year..end_year, as a cells x years matrix (see
# selection_design()). With a supplied regional trend, cells outside any
# region are NA.
selection_trend_matrix <- function(cells, baseline_year, end_year, trend,
                                   design) {
  years <- baseline_year:end_year
  if (is.matrix(trend)) {
    missing_years <- setdiff(years, as.integer(colnames(trend)))
    if (length(missing_years) > 0) {
      stop("the supplied trend has no values for ",
           toString(missing_years))
    }
    region <- cell_regions(cells)
    missing_regions <- setdiff(na.omit(unique(region)), rownames(trend))
    if (length(missing_regions) > 0) {
      stop("the supplied trend has no row for ", toString(missing_regions))
    }
    return(trend[match(region, rownames(trend)),
                 match(years, as.integer(colnames(trend))),
                 drop = FALSE])
  }
  stopifnot(trend == "linear_0_1")
  g <- (years - design$trend_years[1]) / diff(design$trend_years)
  g <- pmax(g, 0)
  matrix(g, length(cells), length(years), byrow = TRUE)
}

# The time-varying selection covariates at mask cells `cells`, years
# baseline_year..end_year, as a cells x years x columns array: nets, irs and
# pop from the cubes (padded as read_padded_cube()), then the hinge columns,
# with trend_pop, pop and its hinges multiplied by g_dom, and with
# trend_crops, the crop products g_ag x crop after them (selection_static()
# then has none).
selection_time_varying <- function(cells, baseline_year, end_year,
                                   design = selection_design()) {
  design <- complete_selection_design(design)
  columns <- selection_column_names(design)
  n_time_varying <- 3 + length(unlist(design$hinges)) +
    if (has_trend(design$trend_crops)) length(selection_crop_names) else 0
  out <- array(NA_real_,
               c(length(cells), end_year - baseline_year + 1, n_time_varying),
               dimnames = list(NULL, baseline_year:end_year,
                               columns[seq_len(n_time_varying)]))
  out[, , 1] <- read_padded_cube(net_use_file(design), cells, baseline_year,
                                 end_year)
  out[, , 2] <- read_padded_cube("data/clean/irs_coverage_scaled_cube.tif",
                                 cells, baseline_year, end_year)
  out[, , 3] <- selection_pop_matrix(cells, baseline_year, end_year, design)
  j <- 3
  for (name in names(design$hinges)) {
    for (k in design$hinges[[name]]) {
      j <- j + 1
      out[, , j] <- pmin(out[, , match(name, c("nets", "irs", "pop"))], k)
    }
  }
  if (has_trend(design$trend_pop)) {
    g <- selection_trend_matrix(cells, baseline_year, end_year,
                                design$trend_pop, design)
    pop_columns <- c(3, 3 + which(rep(names(design$hinges),
                                      lengths(design$hinges)) == "pop"))
    for (p in pop_columns) {
      out[, , p] <- out[, , p] * g
    }
  }
  if (has_trend(design$trend_crops)) {
    g <- selection_trend_matrix(cells, baseline_year, end_year,
                                design$trend_crops, design)
    flat <- selection_flat(cells)
    for (i in seq_len(ncol(flat))) {
      j <- j + 1
      out[, , j] <- flat[, i] * g
    }
  }
  stopifnot(j == n_time_varying)
  out
}

# The net use cube of the design (design$nets, design$net_w,
# design$net_use_source)
net_use_file <- function(design) {
  if (design$nets == "all" && design$net_use_source == "legacy") {
    return("data/clean/net_use_cube.tif")
  }
  file <- if (design$nets == "all") {
    "data/clean/net_use_run06_cube.tif"
  } else {
    sprintf("data/clean/net_use_pyrethroid%s_cube_w%.2f.tif",
            if (design$net_use_source == "legacy") "_legacy" else "",
            design$net_w)
  }
  if (!file.exists(file)) {
    stop(file, " not found; run R/prep_net_use_pyrethroid.R")
  }
  file
}

# The crop layers at mask cells `cells`, cells x selection_crop_names
selection_flat <- function(cells) {
  crops_all <- rast("data/clean/crop_scaled.tif")
  flat <- as.matrix(terra::extract(
    c(rast("data/clean/crop_group_scaled.tif"),
      crops_all[[c("cotton", "vegetables", "rice")]]), cells))
  stopifnot(identical(colnames(flat), selection_crop_names))
  flat
}

# The static columns of the design at mask cells `cells`: the crop layers
# without trend_crops, and none (a cells x 0 matrix) with it, when the crop
# products are among the time-varying columns
selection_static <- function(cells, design = selection_design()) {
  design <- complete_selection_design(design)
  if (has_trend(design$trend_crops)) {
    return(matrix(numeric(0), length(cells), 0))
  }
  selection_flat(cells)
}

# The selection design matrix at mask cells `cells` (cell_id = position in
# `cells`) for years baseline_year..end_year: x_cell_years, one row per
# (cell_id, year_id), cell-major, and its cell_years_index, as
# build_dynamical_model() takes them.
selection_design_matrix <- function(cells, baseline_year, end_year,
                                    design = selection_design()) {
  time_varying <- selection_time_varying(cells, baseline_year, end_year,
                                         design)
  flat <- selection_static(cells, design)
  n_years <- dim(time_varying)[2]
  # cells x years x columns to (years x cells) x columns, year fastest
  long <- matrix(aperm(time_varying, c(2, 1, 3)),
                 ncol = dim(time_varying)[3])
  x_cell_years <- cbind(long, flat[rep(seq_along(cells), each = n_years), ,
                                   drop = FALSE])
  colnames(x_cell_years) <- selection_column_names(design)
  list(x_cell_years = x_cell_years,
       cell_years_index = tibble(
         cell_id = rep(seq_along(cells), each = n_years),
         year_id = rep(seq_len(n_years), length(cells))))
}

# The covariates of `design` on their own scales at mask cells `cells`, for
# figures and summaries: a tibble of cell_id, year_id and one column per
# covariate (nets, irs, pop, the crops), with population min-max scaled and no
# trends or hinges
covariate_extract <- function(cells, baseline_year, end_year,
                              design = selection_design()) {
  design$pop <- "raw"
  design$init_pop <- NULL
  design$trend_pop <- "none"
  design$trend_crops <- "none"
  design$hinges <- list()
  selection <- selection_design_matrix(cells, baseline_year, end_year,
                                       design)
  bind_cols(selection$cell_years_index,
            as_tibble(selection$x_cell_years))
}
