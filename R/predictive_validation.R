# Run out-of-sample predictive validation experiments

# load packages, functions, the cross-validation fold definitions and the null
# model definitions
source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/null_models.R")


# binomial deviance
betabinom_dev <- function(died, mosquito_number, predicted, rho = 0.14) {
  
  # reparameterise from prediction and overdispersion ot the beta parameters
  a <- predicted * (1 / rho - 1)
  b <- a * (1 - predicted) / predicted
  
  log_probs <- extraDistr::dbbinom(x = died,
                                   size = mosquito_number,
                                   alpha = a,
                                   beta = b,
                                   log = TRUE)
  -2 * sum(log_probs)
}

binom_dev <- function(died, mosquito_number, predicted) {
  log_probs <- dbinom(x = died,
                      size = mosquito_number,
                      prob = predicted,
                      log = TRUE)
  -2 * sum(log_probs)
}

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted) ^ 2))
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}


# define the intercept only null model

# sanity check that glm preds and manually calculating mean are numerically similar
spatial_interpolation$training %>% 
  group_by(insecticide_type) %>% 
  mutate(prop_pred = emplog_prop(died,mosquito_number)) %>% 
  summarise(prop_pred = mean(prop_pred),
            mean_pred = mean(died/mosquito_number)) %>% 
  mutate(glm_pred = predict(
    glm(
      cbind(died, (mosquito_number-died)) ~ insecticide_type, 
      data = spatial_interpolation$training, 
      family = stats::binomial), 
    newdata = as.tibble(insecticide_type),
    type = "response"))

# predict type-specific mean predictions
spatial_interpolation_intercept_preds <- predict(
  glm(
    cbind(died, (mosquito_number-died)) ~ insecticide_type, 
    data = spatial_interpolation$training, 
    family = stats::binomial), 
  newdata = spatial_interpolation$test,
  type = "response")

spatial_interpolation_intercept_error <- spatial_interpolation$test %>% 
  mutate(observed = prop(died, mosquito_number),
         predicted = spatial_interpolation_intercept_preds) %>% 
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  group_by(year_start) %>% 
  summarise(
    pred_error_intercept = mean(pred_error),
    bias_intercept = mean(predicted - observed)) 

spatial_interpolation_intercept_error_overall <- spatial_interpolation$test %>% 
  mutate(observed = prop(died, mosquito_number),
         predicted = spatial_interpolation_intercept_preds) %>% 
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  summarise(
    pred_error_intercept = mean(pred_error),
    bias_intercept = mean(predicted - observed)) 

temporal_forecasting_intercept_preds <- predict(
  glm(
    cbind(died, (mosquito_number-died)) ~ insecticide_type, 
    data = temporal_forecasting$training, 
    family = stats::binomial), 
  newdata = temporal_forecasting$test,
  type = "response")

temporal_forecasting_intercept_error <- temporal_forecasting$test %>% 
  mutate(observed = prop(died, mosquito_number),
         predicted = temporal_forecasting_intercept_preds) %>% 
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  group_by(year_start) %>% 
  summarise(pred_error_intercept = mean(pred_error),
            bias_intercept = mean(predicted - observed)) 

get_spatial_extrapolation_intercept_pred_error <- function(test,training) {
  pred <- predict(
    glm(
      cbind(died, (mosquito_number-died)) ~ insecticide_type, 
      data = training, 
      family = stats::binomial), 
    newdata = test,
    type = "response")
  
  test %>% 
    mutate(predicted = pred,
           observed = prop(died, mosquito_number)) %>% 
    group_by(row_number()) %>% 
    mutate(pred_error = betabinom_dev(died = died,
                                  mosquito_number = mosquito_number,
                                  predicted = predicted)) %>% 
    ungroup() %>% 
    summarise(pred_error_intercept = mean(pred_error),
              #full_deviance = sum(pred_error),
              bias_intercept = mean(predicted - observed)) 
}

spatial_extrapolation_intercept_error <- mapply(
  get_spatial_extrapolation_intercept_pred_error,
  spatial_extrapolation$test,
  spatial_extrapolation$training,
  SIMPLIFY = FALSE
) %>%
  do.call(
    bind_rows, .
  ) %>% 
  mutate(country_name = countries_to_validate)



# Do grid search to find the optimal number of nearest neighbours for spatial
# interpolation and temporal forecasting

# for a fairer test, use an 'internal' test set from the training data for
# selecting the optimal NN number, such that the approach is not cheating wrt to
# optimising performance on test fold
test_internal_train <- spatial_interpolation$training %>%
  mutate(for_optim = FALSE)

set.seed(111)
test_internal_train[sample(nrow(test_internal_train),100),'for_optim'] <- TRUE


spatial_interpolation_optimal_nn <- test_internal_train %>%
  filter(for_optim) %>% 
  mutate(
    mortality = prop(died, mosquito_number)
  ) %>%
  predict_null_optimal_nn(
    training_data = test_internal_train %>%
      filter(!for_optim)
  ) %>%
  filter(nn_is_optimal) %>%
  slice(1) %>%
  pull(n_neighbours)
  
# build full prediction onto the test fold
spatial_interpolation_preds <- spatial_interpolation$test %>%
  mutate(
    mortality = prop(died, mosquito_number)
  ) %>%
  predict_null_optimal_nn(
    training_data = spatial_interpolation$training,
    n_nearest_neighbour_range = c(spatial_interpolation_optimal_nn,
                                  spatial_interpolation_optimal_nn) 
  )


# for temporal forecasting, use up to 3 years prior to enable prediction to the
# third year into the future
test_internal_train <- temporal_forecasting$training %>%
  mutate(for_optim = FALSE)

set.seed(111)
test_internal_train[sample(nrow(test_internal_train),100),'for_optim'] <- TRUE


temporal_forecasting_optimal_nn <- test_internal_train %>%
  filter(for_optim) %>% 
  mutate(
    mortality = prop(died, mosquito_number)
  ) %>%
  predict_null_optimal_nn(
    training_data = test_internal_train %>%
      filter(!for_optim)
  ) %>%
  filter(nn_is_optimal) %>%
  slice(1) %>%
  pull(n_neighbours)

temporal_forecasting_preds <- temporal_forecasting$test %>%
  mutate(
    mortality = prop(died, mosquito_number)
  ) %>%
  predict_null_optimal_nn(
    training_data = temporal_forecasting$training,
    n_years_prior = 3,
    n_nearest_neighbour_range = c(temporal_forecasting_optimal_nn,
                                  temporal_forecasting_optimal_nn)
  )

# for spatial extrapolation, need to run it for each country and compute the average rmses to
# RMSEs to identify the optimal number of neighbours

each_country_optimal_nn_internal <- function(training) {
  
  test_internal_train <- training %>%
    mutate(for_optim = FALSE)
  
  set.seed(111)
  test_internal_train[sample(nrow(test_internal_train),100),'for_optim'] <- TRUE
  
  test_internal_train %>%
    filter(for_optim) %>%
    mutate(
      mortality = prop(died, mosquito_number)
    ) %>%
    predict_null_optimal_nn(
      training_data = test_internal_train %>%
        filter(!for_optim)
    )
}

spatial_extrapolation_internal_preds <- lapply(
  FUN = each_country_optimal_nn_internal,
  spatial_extrapolation$training
) %>%
  do.call(
    bind_rows, .
  )

# get and plot the overall rmses
spatial_extrapolation_internal_pred_errors <- spatial_extrapolation_internal_preds %>%
  group_by(
    n_neighbours
  ) %>%
  summarise(
    pred_error = betabinom_dev(died = died,
                               mosquito_number = mosquito_number,
                               predicted = predicted)
  )

plot(pred_error ~ n_neighbours,
     data = spatial_extrapolation_internal_pred_errors,
     type = "b")

# find the optimum
spatial_extrapolation_optimal_nn <- spatial_extrapolation_internal_pred_errors %>%
  filter(pred_error == min(pred_error)) %>%
  slice(1) %>%
  pull(n_neighbours)

# make full predictions onto test fold
each_country_optimal_nn <- function(test, training) {
  test %>%
    mutate(
      mortality = prop(died, mosquito_number)
    ) %>%
    predict_null_optimal_nn(
      training_data = training,
      n_nearest_neighbour_range = rep(spatial_extrapolation_optimal_nn,2)
    )
}

spatial_extrapolation_preds <- mapply(
  FUN = each_country_optimal_nn,
  spatial_extrapolation$test,
  spatial_extrapolation$training,
  SIMPLIFY = FALSE
) %>%
  do.call(
    bind_rows, .
  )

# get and plot the overall rmses
spatial_extrapolation_pred_errors <- spatial_extrapolation_preds %>%
  group_by(
    n_neighbours
  ) %>%
  summarise(
    pred_error = betabinom_dev(died = died,
                               mosquito_number = mosquito_number,
                               predicted = predicted)
  )

# overwrite the optimal flag
spatial_extrapolation_preds <- spatial_extrapolation_preds %>%
  mutate(
    nn_is_optimal = n_neighbours == spatial_extrapolation_optimal_nn
  )

# subset each of these validation experiments to the optimal numbers of
# neighbours
optimal_nn_preds <- bind_rows(
  spatial_interpolation = spatial_interpolation_preds,
  spatial_extrapolation = spatial_extrapolation_preds,
  temporal_forecasting = temporal_forecasting_preds,
  .id = "experiment"
) %>%
  filter(
    nn_is_optimal
  )

# record the optima
optimal_nn <- optimal_nn_preds %>%
  group_by(
    experiment
  ) %>%
  summarise(
    n_neighbours = n_neighbours[1]
  )
optimal_nn

# save optimal nn numbers
write_csv(optimal_nn,"outputs/optimal_nn.csv")

# compute overall bias and RMSE for each of these and tabulate
optimal_nn_preds %>%
  group_by(
    experiment
  ) %>%
  summarise(
    pred_error = betabinom_dev(died = died,
                               mosquito_number = mosquito_number,
                               predicted = predicted),
    bias = mean(predicted - observed),
    .groups = "drop"
  )

# only for spatial interpolation
optimal_nn_preds %>%
  filter(
    experiment == "spatial_interpolation"
  ) %>%
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  summarise(
    pred_error_nn = mean(pred_error),
    bias_nn = mean(predicted - observed),
    .groups = "drop"
  ) %>% 
  mutate(experiment = "spatial_interpolation") -> interp_result

# save null model CV results as csvs
write_csv(interp_result,"outputs/interp_result.csv")

# split these by country (spatial extrapolation) and years ahead (temporal
# forecasting)
optimal_nn_preds %>%
  filter(
    experiment == "spatial_extrapolation"
  ) %>%
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  group_by(
    country_name
  ) %>%
  summarise(
    #full_deviance = sum(pred_error),
    pred_error = mean(pred_error),
    bias = mean(predicted - observed),
    .groups = "drop"
  ) %>% 
  mutate(experiment = "spatial_extrapolation") -> extrap_result

# save null model CV results as csvs
write_csv(extrap_result,"outputs/extrap_result.csv")

optimal_nn_preds %>%
  filter(
    experiment == "temporal_forecasting"
  ) %>%
  group_by(row_number()) %>% 
  mutate(pred_error = betabinom_dev(died = died,
                                mosquito_number = mosquito_number,
                                predicted = predicted)) %>% 
  ungroup() %>% 
  group_by(
    year_start
  ) %>%
  summarise(
    pred_error = mean(pred_error),
    bias = mean(predicted - observed),
    .groups = "drop"
  ) %>% 
  mutate(experiment = "temporal_forecasting") -> temp_forecast_result
# save null model CV results as csvs
write_csv(temp_forecast_result,"outputs/temp_forecast_result.csv")
# These are substantially more negative (overpredicting susceptibility) 3y into
# the future.


# # Compute spatial extrapolation and interpolation ability based on Penny's
# # published maps, subsetted to (2010 to 2017) and only in the regions covered
# # https://doi.org/10.6084/m9.figshare.9912623
# 
# hancock_preds <- terra::rast("data/clean/hancock_2020_predictions.tif")
# 
# # extract predictions form Hancock layers for a given a fixed insecticide type
# # and year, and vectors of coordinates
# pred_hancock <- function(insecticide_type,
#                          year,
#                          latitude,
#                          longitude) {
#   
#   # form the coordinates matrix
#   coords <- cbind(longitude, latitude)
#   
#   # extract all the values, for all layers and reshape these to get the row
#   # number, insecticide_type, and year
#   values <- hancock_preds %>%
#     terra::extract(coords) %>%
#     as_tibble() %>%
#     mutate(
#       row = row_number()
#     ) %>%
#     pivot_longer(
#       cols = !any_of("row"),
#       names_to = c("insecticide_type", "year"),
#       names_pattern = "(.*)_(.*)",
#       values_to = "prediction"
#     ) %>%
#     mutate(
#       insecticide_type = case_when(
#         insecticide_type == "Alphacypermethrin" ~ "Alpha-cypermethrin",
#         insecticide_type == "Lambdacyhalothrin" ~ "Lambda-cyhalothrin-cypermethrin",
#         .default = insecticide_type
#       ),
#       year = as.numeric(year)
#     )
# 
#   # list the data we want, and get it from the list of values
#   targets <- tibble(
#     row = seq_along(longitude),
#     insecticide_type = insecticide_type,
#     year = year
#   ) %>%
#   left_join(
#     values,
#     by = c("row", "insecticide_type", "year")
#   ) %>%
#     pull(prediction)
# 
# }
# 
# # for the full set of optimal nn predictions, extract the prediction from
# # Hancock and append. Then subset to only those records with Hancock predictions
# # and compute null and Hancock validation metrics
# 
# 
# hancock_test_set <- optimal_nn_preds %>%
#   # pre-emptively subset to Hancock prediction years and insecticides
#   filter(
#     year_start %in% 2006:2017,
#     insecticide_type %in% c("Alpha-cypermethrin",
#                             "DDT",
#                             "Deltamethrin",
#                             "Lambda-cyhalothrin",
#                             "Permethrin") 
#   ) %>%
#   # do extraction in batches of insecticide types and year, so we can extract
#   # from one raster layer at a time
#   mutate(
#     predicted_hancock = pred_hancock(insecticide_type = insecticide_type,
#                                      year = year_start,
#                                      latitude = latitude,
#                                      longitude = longitude)
#   ) %>%
#   # drop any test data for which hancock made no prediction
#   filter(
#     !is.na(predicted_hancock)
#   ) %>%
#   mutate(
#     # clamp greater-than-1 predictions to 1
#     predicted_hancock = pmin(predicted_hancock, max(predicted_hancock[predicted_hancock < 1]))
#   )
# 
# # summarise these predictions and the nn predictions
# plot(hancock_test_set$predicted_hancock ~ hancock_test_set$predicted)
# range(hancock_test_set$predicted_hancock)
# 
# # split these by country (spatial extrapolation) and years ahead (temporal
# # forecasting)
# hancock_test_set %>%
#   filter(
#     experiment == "spatial_interpolation"
#   ) %>%
#   mutate(predicted_intercept = predict(
#     glm(
#       cbind(died, (mosquito_number-died)) ~ insecticide_type, 
#       family = stats::binomial), 
#     type = "response")) %>% 
#   group_by(
#     year_start
#   ) %>%
#   summarise(
#     pred_error_nn = betabinom_dev(died = died,
#                                   mosquito_number = mosquito_number,
#                                   predicted = predicted),
#     pred_error_hancock = betabinom_dev(died = died,
#                                   mosquito_number = mosquito_number,
#                                   predicted = predicted_hancock),
#     bias_nn = mean(predicted - observed),
#     bias_hancock = mean(predicted_hancock - observed),
#     pred_error_intercept = betabinom_dev(died = died,
#                                          mosquito_number = mosquito_number,
#                                          predicted = predicted_intercept),
#     bias_intercept = mean(predicted_intercept - observed),
#     .groups = "drop"
#   ) %>% 
#   mutate(experiment = "hancock_spatial_interpolation") -> hancock_interp_result
# # At spatial interpolation, Hancock et al. is consistently better than the
# # nearest neighbour heuristic on both prediction error, and
# # similar in terms of bias (though generally overestimating susceptibility)
# 
# # save null model CV results as csvs
# write_csv(hancock_interp_result,"outputs/hancock_interp_result.csv")
# 
# hancock_test_set %>%
#   filter(
#     experiment == "spatial_extrapolation"
#   ) %>%
#   mutate(predicted_intercept = predict(
#     glm(
#       cbind(died, (mosquito_number-died)) ~ insecticide_type, 
#       family = stats::binomial), 
#     type = "response")) %>% 
#   group_by(
#     country_name
#   ) %>%
#   summarise(
#     pred_error_nn = betabinom_dev(died = died,
#                                   mosquito_number = mosquito_number,
#                                   predicted = predicted),
#     pred_error_hancock = betabinom_dev(died = died,
#                                        mosquito_number = mosquito_number,
#                                        predicted = predicted_hancock),
#     bias_nn = mean(predicted - observed),
#     bias_hancock = mean(predicted_hancock - observed),
#     pred_error_intercept = betabinom_dev(died = died,
#                                          mosquito_number = mosquito_number,
#                                          predicted = predicted_intercept),
#     bias_intercept = mean(predicted_intercept - observed),
#     .groups = "drop"
#   ) %>% 
#   mutate(experiment = "hancock_spatial_extrapolation") -> hancock_extrap_result
# # At spatial extrapolation Hancock et al is consistently better than the nearest
# # neighbour heuristic (41 nearest neighbour) on prediction error and better on
# # bias for most countries; however Hancock et al has consistent bias of
# # overestimating susceptibility, whereas null model doesn't seem to be biased in
# # one direction on average
# 
# # save null model CV results as csvs
# write_csv(hancock_extrap_result,"outputs/hancock_extrap_result.csv")
