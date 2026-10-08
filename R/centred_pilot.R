# Pilot of the centred parameterisation (#48): fits to the full data compared
# on effective samples per gradient evaluation and per hour.
#
#   Rscript R/centred_pilot.R fit <label> [threads]
#     one fit, with the model options in IR_CUBE_MODEL_OPTIONS and the sampler
#     settings in IR_CUBE_MCMC_SETTINGS (as R/fit_model.R; by default
#     dynamical_mcmc_settings(warmup = 1000, n_samples = 1000)), from the
#     cached initial values (dynamical_inits_files()), recording the number
#     of leapfrog steps, step size and time of every burst. Writes
#     temporary/centred_pilot/<label>.rds.
#   Rscript R/centred_pilot.R compare <fit> <fit> ...
#     for each fit, a label of the above or a fit_model.R image
#     (.../temporary/fitted_model.RData): the adapted step size, the time per
#     gradient, and the bulk and tail ESS and R-hat of the quantities every
#     parameterisation shares (the variables common to all fits, and
#     beta_class, beta_type and rho_types), with the bulk ESS per 1,000
#     gradient evaluations and per hour. For an image, the number of leapfrog
#     steps is taken as the mean of its range, and the hours are those of the
#     whole run (sampling_time), if recorded, or given as <fit>=<hours>. Plain
#     R, with greta attached to read an image.
#
# e.g. a local pilot, about 6 GB and 8 threads each:
#   IR_CUBE_MCMC_SETTINGS='dynamical_mcmc_settings(warmup = 1000, n_samples = 1000)' \
#     Rscript R/centred_pilot.R fit noncentred 8 &
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(centred = centred_options_data_informed())' \
#   IR_CUBE_MCMC_SETTINGS='dynamical_mcmc_settings(warmup = 1000, n_samples = 1000, Lmin = 15, Lmax = 30)' \
#     Rscript R/centred_pilot.R fit centred 8 &
#   wait; Rscript R/centred_pilot.R compare noncentred centred
#
# Run the fits with the greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
mode <- arguments[1]
pilot_dir <- "temporary/centred_pilot"


# a fit ----------------------------------------------------------------------

# windowed_hmc() with a record of each burst: its iterations, number of
# leapfrog steps, step size, mean acceptance and elapsed seconds, in the
# sampler's burst_log (one data frame per burst)
logged_windowed_hmc <- function(settings) {
  sampler <- windowed_hmc(Lmin = settings$Lmin, Lmax = settings$Lmax,
                          accept_target = settings$accept_target)
  # R6 evaluates `inherit` when an object is created, so the parent class is
  # bound to its own name here
  windowed_class <- sampler$class
  sampler$class <- R6::R6Class(
    "logged_windowed_hmc_sampler",
    inherit = windowed_class,
    public = list(
      burst_l = NA_integer_,
      burst_log = list(),
      sampler_parameter_values = function() {
        values <- super$sampler_parameter_values()
        self$burst_l <- values$hmc_l
        values
      },
      run_burst = function(n_samples, thin = 1L) {
        start <- proc.time()[["elapsed"]]
        super$run_burst(n_samples, thin)
        self$burst_log[[length(self$burst_log) + 1]] <- data.frame(
          iterations = n_samples,
          leapfrog_steps = self$burst_l,
          epsilon = self$parameters$epsilon[1],
          accept = self$mean_accept_stat,
          seconds = proc.time()[["elapsed"]] - start)
      }
    )
  )
  sampler
}

fit_pilot <- function(label, threads) {
  source("R/greta_setup.R")
  start_greta(threads = threads)
  source("R/dynamical_model.R")
  suppressMessages({
    sink("/dev/null")
    source("R/validation_folds.R")
    source("R/validation_covariates.R")
    sink()
  })
  source("R/dynamical_predictions.R")
  settings <- eval(str2lang(Sys.getenv(
    "IR_CUBE_MCMC_SETTINGS",
    "dynamical_mcmc_settings(warmup = 1000, n_samples = 1000)")))
  cat("options:", Sys.getenv("IR_CUBE_MODEL_OPTIONS", "defaults"), "\n")
  print(settings[c("n_chains", "warmup", "n_samples", "Lmin", "Lmax")])

  built <- build_dynamical_model(train_df = df,
                                 df = df,
                                 x_cell_years = x_cell_years,
                                 cell_years_index = cell_years_index,
                                 classes_index = classes_index,
                                 types = types,
                                 options = model_options,
                                 x_cells_init = x_cells_init)
  inits <- dynamical_chain_inits(dynamical_inits_files(), built$variables,
                                 levels = built$lookups$levels,
                                 columns = colnames(x_cell_years),
                                 n_chains = settings$n_chains,
                                 options = built$options,
                                 classes_index = classes_index)

  # as run_dynamical_mcmc(), with the logged sampler; mcmc() matches the
  # initial values to the variables by name in the calling frame
  run <- function(m, variables) {
    list2env(variables, environment())
    mcmc(m,
         chains = settings$n_chains,
         initial_values = inits,
         warmup = settings$warmup,
         sampler = logged_windowed_hmc(settings),
         n_samples = settings$n_samples,
         pb_update = settings$pb_update)
  }
  elapsed <- system.time(draws <- run(built$model, built$variables))
  sampler <- attr(draws, "model_info")$samplers[[1]]
  bursts <- do.call(rbind, sampler$burst_log)
  in_warmup <- cumsum(bursts$iterations) <= settings$warmup
  bursts$phase <- ifelse(in_warmup, "warmup", "sampling")

  values <- as.matrix(draws)
  dir.create(pilot_dir, showWarnings = FALSE, recursive = TRUE)
  file <- file.path(pilot_dir, paste0(label, ".rds"))
  saveRDS(list(label = label,
               options = built$options,
               settings = settings,
               threads = threads,
               n_chains = length(draws),
               values = values,
               shared = shared_quantities(values, classes_index, types,
                                          built$options),
               bursts = bursts,
               epsilon = sampler$parameters$epsilon[1],
               diag_sd = sampler$parameters$diag_sd,
               metric_log = sampler$metric_log,
               elapsed = elapsed[["elapsed"]]),
          file)
  cat(sprintf("%s: %.0f s; written to %s\n", label, elapsed[["elapsed"]],
              file))
}


# the comparison -------------------------------------------------------------

# The quantities every parameterisation shares, from the draws `values` (draws
# x columns, chains stacked) of a model with `options`: beta_class, beta_type
# and rho_types, as draws x elements named as greta names them
shared_quantities <- function(values, classes_index, types, options) {
  names <- unique(sub("\\[.*$", "", colnames(values)))
  v <- lapply(setNames(nm = names), extract_parameter, draws_matrix = values)
  terms <- dynamical_terms_draws(v, classes_index, types,
                                 terms = c("beta_class", "beta_type",
                                           "rho_types"),
                                 options = options)
  shared <- lapply(names(terms), function(name) {
    x <- matrix(terms[[name]], nrow(values))
    index <- arrayInd(seq_len(ncol(x)), dim(terms[[name]])[-1])
    colnames(x) <- sprintf("%s[%s]", name,
                           apply(index, 1, paste, collapse = ","))
    x
  })
  do.call(cbind, shared)
}

# What the comparison needs from one fit: a pilot fit's label, or a
# fit_model.R image
pilot_inputs <- function(fit) {
  if (!grepl("\\.RData$", fit)) {
    pilot <- readRDS(file.path(pilot_dir, paste0(fit, ".rds")))
    sampling <- pilot$bursts[pilot$bursts$phase == "sampling", ]
    return(list(label = pilot$label,
                settings = pilot$settings,
                n_chains = pilot$n_chains,
                values = pilot$values,
                shared = pilot$shared,
                epsilon = pilot$epsilon,
                diag_sd = pilot$diag_sd,
                gradients = sum(sampling$iterations *
                                  sampling$leapfrog_steps),
                total_gradients = sum(pilot$bursts$iterations *
                                        pilot$bursts$leapfrog_steps),
                n_samples = sum(sampling$iterations),
                accept = weighted.mean(sampling$accept, sampling$iterations),
                sampling_hours = sum(sampling$seconds) / 3600,
                total_hours = sum(pilot$bursts$seconds) / 3600))
  }
  image <- new.env()
  load(fit, envir = image)
  settings <- image$settings
  sampler <- attr(image$draws, "model_info")$samplers[[1]]
  accept <- sampler$accept_history
  sampling_rows <- seq(nrow(accept) - settings$n_samples + 1, nrow(accept))
  values <- as.matrix(image$draws)
  mean_l <- (settings$Lmin + settings$Lmax) / 2
  list(label = basename(sub("/temporary/fitted_model\\.RData$", "", fit)),
       settings = settings,
       n_chains = length(image$draws),
       values = values,
       shared = shared_quantities(values, image$classes_index, image$types,
                                  image$model_options),
       epsilon = sampler$parameters$epsilon[1],
       diag_sd = sampler$parameters$diag_sd,
       gradients = settings$n_samples * mean_l,
       total_gradients = (settings$warmup + settings$n_samples) * mean_l,
       n_samples = settings$n_samples,
       accept = mean(accept[sampling_rows, ]),
       sampling_hours = NA,
       total_hours = if (!is.null(image$sampling_time)) {
         image$sampling_time[["elapsed"]] / 3600
       } else {
         NA
       })
}

compare_pilots <- function(fits) {
  suppressMessages({
    library(greta)
    library(dplyr)
    library(stringr)
  })
  source("R/dynamical_predictions.R")
  # a fit given as <fit>=<hours> has its run time from there
  inputs <- lapply(strsplit(fits, "="), function(fit) {
    input <- pilot_inputs(fit[1])
    if (length(fit) == 2) {
      input$total_hours <- as.numeric(fit[2])
    }
    input
  })
  # the variables every fit has, with the same meaning in each
  variable_names <- function(input) {
    unique(sub("\\[.*$", "", colnames(input$values)))
  }
  common <- Reduce(intersect, lapply(inputs, variable_names))
  centring_variables <- c("beta_class_raw", "beta_type_raw", "rho_type_raw",
                          "beta_class_centred", "beta_type_centred",
                          "logit_rho_type")
  common <- setdiff(common, centring_variables)

  rows <- lapply(inputs, function(input) {
    columns <- sub("\\[.*$", "", colnames(input$values)) %in% common
    quantities <- cbind(input$values[, columns, drop = FALSE], input$shared)
    n_iterations <- nrow(quantities) / input$n_chains
    draws <- posterior::as_draws_array(
      array(quantities, c(n_iterations, input$n_chains, ncol(quantities)),
            dimnames = list(NULL, NULL, colnames(quantities))))
    summary <- posterior::summarise_draws(draws, "rhat", "ess_bulk",
                                          "ess_tail")
    # gradient evaluations per chain while sampling: L per iteration
    gradients <- input$gradients
    mean_l <- gradients / input$n_samples
    # greta's step for parameter i is epsilon diag_sd_i / sum(diag_sd), so
    # this is the step in units of its posterior sd
    step_per_sd <- input$epsilon / sum(input$diag_sd)
    data.frame(
      fit = input$label,
      L = sprintf("%d-%d", input$settings$Lmin, input$settings$Lmax),
      mean_L = mean_l,
      epsilon = input$epsilon,
      step_per_sd = step_per_sd,
      trajectory_per_sd = step_per_sd * mean_l,
      accept = input$accept,
      ms_per_gradient = 3.6e6 * input$total_hours / input$total_gradients,
      sampling_hours = input$sampling_hours,
      total_hours = input$total_hours,
      quantities = nrow(summary),
      ess_bulk_min = min(summary$ess_bulk),
      ess_bulk_median = median(summary$ess_bulk),
      ess_tail_min = min(summary$ess_tail),
      rhat_max = max(summary$rhat),
      rhat_over_1.01 = sum(summary$rhat > 1.01),
      ess_min_per_1000_gradients = 1000 * min(summary$ess_bulk) / gradients,
      ess_median_per_1000_gradients = 1000 * median(summary$ess_bulk) /
        gradients,
      ess_min_per_hour = min(summary$ess_bulk) / input$total_hours,
      ess_median_per_hour = median(summary$ess_bulk) / input$total_hours,
      worst = summary$variable[which.min(summary$ess_bulk)])
  })
  out <- do.call(rbind, rows)
  print(t(format(out, digits = 3)), quote = FALSE)
  invisible(out)
}


if (identical(mode, "fit")) {
  fit_pilot(arguments[2],
            threads = if (length(arguments) >= 3) as.integer(arguments[3]))
} else if (identical(mode, "compare")) {
  compare_pilots(arguments[-1])
} else {
  stop("usage: Rscript R/centred_pilot.R fit <label> [threads] | ",
       "compare <fit> ...")
}
