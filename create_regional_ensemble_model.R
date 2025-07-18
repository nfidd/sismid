# Create Regional Ensemble Model
# Strategy: Fit individual ARIMA models per region, then ensemble
# This allows each region to have its own dynamics while benefiting from ensemble averaging

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("lubridate")
library("hubUtils")

set.seed(406)
data(flu_data_hhs)

print("=== Creating Regional Ensemble Model ===")
print("Strategy: Individual ARIMA per region with sqrt transformation")
print("Innovation: Region-specific parameters with ensemble robustness")

# Test on validation data first
validation_start <- as.Date("2016-10-22")
validation_end <- as.Date("2017-10-14")

validation_dates <- flu_data_hhs |>
  filter(origin_date >= validation_start & origin_date <= validation_end) |>
  pull(origin_date) |>
  unique() |>
  sort()

print(paste("Testing on", length(validation_dates), "validation dates"))

# Function to generate regional ensemble forecasts
generate_regional_ensemble_forecast <- function(forecast_date, flu_data) {
  # Prepare data
  train_data <- flu_data |>
    filter(origin_date <= forecast_date) |>
    mutate(
      # Apply sqrt transformation
      wili_sqrt = sqrt(wili)
    ) |>
    as_tsibble(index = origin_date, key = location)
  
  # Check if we have enough data
  min_obs <- train_data |>
    as_tibble() |>
    group_by(location) |>
    summarise(n = n()) |>
    pull(n) |>
    min()
  
  if (min_obs < 20) {
    return(NULL)
  }
  
  # Fit models - use automatic ARIMA selection per region
  models <- train_data |>
    model(
      # Three different models to ensemble
      arima_auto = ARIMA(wili_sqrt),
      arima_seasonal = ARIMA(wili_sqrt ~ pdq(2,1,0) + PDQ(1,0,0)),
      ets = ETS(wili_sqrt)
    )
  
  # Generate forecasts from each model
  forecasts <- models |>
    forecast(h = 4)
  
  # Extract samples for ensemble
  samples_list <- list()
  
  for (model_name in c("arima_auto", "arima_seasonal", "ets")) {
    model_samples <- models |>
      select(!!model_name) |>
      generate(h = 4, times = 33) |>  # 33 samples per model = 99 total
      mutate(
        # Transform back from sqrt
        .sim = .sim^2,
        .sim = pmax(.sim, 0),  # Ensure non-negative
        model = model_name
      ) |>
      as_tibble()
    
    samples_list[[model_name]] <- model_samples
  }
  
  # Combine all samples
  all_samples <- bind_rows(samples_list) |>
    group_by(location, origin_date) |>
    mutate(
      horizon = row_number(),
      target_end_date = origin_date
    ) |>
    ungroup() |>
    mutate(
      origin_date = forecast_date,
      value = .sim,
      target = "ili perc"
    )
  
  # Calculate quantiles from combined samples
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- all_samples |>
    group_by(location, origin_date, horizon, target_end_date) |>
    reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
    mutate(
      output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
      target = "ili perc",
      output_type = "quantile"
    ) |>
    select(-quantile)
  
  return(forecast_quantiles)
}

# Test on validation set
validation_results <- list()
successful_dates <- 0

for (i in seq_along(validation_dates)) {
  date <- validation_dates[i]
  if (i %% 10 == 1) {
    print(paste("Processing validation date", i, "of", length(validation_dates)))
  }
  
  tryCatch({
    forecast <- generate_regional_ensemble_forecast(date, flu_data_hhs)
    if (!is.null(forecast)) {
      validation_results[[length(validation_results) + 1]] <- forecast
      successful_dates <- successful_dates + 1
    }
  }, error = function(e) {
    # Silently skip errors during validation
  })
}

print(paste("\nSuccessfully generated forecasts for", successful_dates, "dates"))

if (length(validation_results) > 0) {
  validation_forecasts <- bind_rows(validation_results)
  
  # Calculate WIS for validation
  print("\n=== Calculating Validation WIS ===")
  
  # Get actual values
  actual_values <- flu_data_hhs |>
    filter(origin_date %in% unique(validation_forecasts$target_end_date)) |>
    select(origin_date, location, wili) |>
    rename(target_end_date = origin_date, observation = wili)
  
  # Join forecasts with actuals
  validation_eval <- validation_forecasts |>
    left_join(actual_values, by = c("target_end_date", "location")) |>
    filter(!is.na(observation))
  
  # Calculate WIS
  wis_by_horizon <- validation_eval |>
    pivot_wider(names_from = output_type_id, values_from = value, names_prefix = "q") |>
    mutate(
      # Interval scores
      is_50 = pmax(0, pmax(q0.25 - observation, observation - q0.75)) * 2/0.5,
      is_80 = pmax(0, pmax(q0.1 - observation, observation - q0.9)) * 2/0.8,
      is_95 = pmax(0, pmax(q0.025 - observation, observation - q0.975)) * 2/0.95,
      # WIS
      wis = (is_50 + is_80 + is_95) / 3 + 0.5 * abs(q0.5 - observation)
    ) |>
    group_by(horizon) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      n = n()
    )
  
  print(wis_by_horizon)
  
  overall_wis <- mean(wis_by_horizon$mean_wis)
  print(paste("\nOverall Validation WIS:", round(overall_wis, 3)))
  
  # Compare to our VAR(2) sqrt baseline of 0.208
  improvement <- (0.208 - overall_wis) / 0.208 * 100
  print(paste("Improvement over VAR(2) sqrt:", round(improvement, 1), "%"))
  
  if (overall_wis < 0.25) {
    print("\n*** This model shows promise! Proceeding with full implementation ***")
  }
  
} else {
  print("No validation forecasts generated")
}