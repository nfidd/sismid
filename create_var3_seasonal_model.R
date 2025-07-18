# Create VAR(3) with seasonal indicators model
# Building on our learnings: simple VAR models work best, sqrt transformation is optimal
# New innovation: Include seasonal indicators directly in the VAR model

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("lubridate")
library("hubUtils")

set.seed(406)
data(flu_data_hhs)

# Define sqrt transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

print("=== Creating VAR(3) with Seasonal Indicators Model ===")
print("Innovation: Direct seasonal modeling within VAR framework")

# Test on validation data first
validation_start <- as.Date("2016-10-22")
validation_end <- as.Date("2017-10-14")

# Get validation dates
validation_dates <- flu_data_hhs |>
  filter(origin_date >= validation_start & origin_date <= validation_end) |>
  pull(origin_date) |>
  unique() |>
  sort()

print(paste("Testing on", length(validation_dates), "validation dates"))

# Function to add seasonal features
add_seasonal_features <- function(data) {
  data |>
    mutate(
      week_num = week(origin_date),
      # Create seasonal indicator: 1 during flu season, 0 otherwise
      flu_season = as.numeric(week_num >= 40 | week_num <= 20),
      # Create smooth seasonal transition using cosine
      seasonal_smooth = cos(2 * pi * week_num / 52)
    )
}

# Function to generate forecasts
generate_var3_seasonal_forecast <- function(forecast_date, flu_data) {
  # Prepare data with seasonal features
  train_data <- flu_data |>
    filter(origin_date <= forecast_date) |>
    add_seasonal_features() |>
    select(origin_date, location, wili, flu_season, seasonal_smooth) |>
    pivot_wider(
      names_from = location, 
      values_from = wili,
      names_prefix = "loc_"
    ) |>
    as_tsibble(index = origin_date)
  
  # Check if we have enough data
  if (nrow(train_data) < 20) {
    return(NULL)
  }
  
  # Fit VAR(3) with sqrt transformation and seasonal indicators
  model_fit <- train_data |>
    model(
      var3_seasonal = VAR(vars(
        `loc_HHS Region 1` = my_sqrt(`loc_HHS Region 1`),
        `loc_HHS Region 2` = my_sqrt(`loc_HHS Region 2`),
        `loc_HHS Region 3` = my_sqrt(`loc_HHS Region 3`),
        `loc_HHS Region 4` = my_sqrt(`loc_HHS Region 4`),
        `loc_HHS Region 5` = my_sqrt(`loc_HHS Region 5`),
        `loc_HHS Region 6` = my_sqrt(`loc_HHS Region 6`),
        `loc_HHS Region 7` = my_sqrt(`loc_HHS Region 7`),
        `loc_HHS Region 8` = my_sqrt(`loc_HHS Region 8`),
        `loc_HHS Region 9` = my_sqrt(`loc_HHS Region 9`),
        `loc_HHS Region 10` = my_sqrt(`loc_HHS Region 10`),
        `loc_US National` = my_sqrt(`loc_US National`)
      ) ~ AR(3) + flu_season + seasonal_smooth)
    )
  
  # Create future seasonal features
  future_dates <- seq(from = forecast_date + 7, by = 7, length.out = 4)
  future_seasonal <- data.frame(
    origin_date = future_dates,
    flu_season = as.numeric(week(future_dates) >= 40 | week(future_dates) <= 20),
    seasonal_smooth = cos(2 * pi * week(future_dates) / 52)
  ) |>
    as_tsibble(index = origin_date)
  
  # Generate forecasts with seasonal features
  forecast_samples <- model_fit |>
    generate(
      new_data = future_seasonal,
      times = 100, 
      bootstrap = TRUE
    ) |>
    mutate(horizon = row_number()) |>
    as_tibble()
  
  # Convert to long format and prepare output
  location_cols <- paste0("loc_", c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                                    "HHS Region 4", "HHS Region 5", "HHS Region 6",
                                    "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                                    "HHS Region 10", "US National"))
  
  forecast_long <- forecast_samples |>
    select(.rep, origin_date, horizon, all_of(location_cols)) |>
    pivot_longer(
      cols = all_of(location_cols),
      names_to = "location",
      values_to = "value",
      names_prefix = "loc_"
    ) |>
    mutate(
      target_end_date = origin_date,
      origin_date = forecast_date,
      target = "ili perc"
    )
  
  # Calculate quantiles
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- forecast_long |>
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
for (i in seq_along(validation_dates)) {
  date <- validation_dates[i]
  if (i %% 10 == 1) {
    print(paste("Processing validation date", i, "of", length(validation_dates)))
  }
  
  tryCatch({
    forecast <- generate_var3_seasonal_forecast(date, flu_data_hhs)
    if (!is.null(forecast)) {
      validation_results[[length(validation_results) + 1]] <- forecast
    }
  }, error = function(e) {
    print(paste("Error for date", date, ":", e$message))
  })
}

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
  
  # Simple WIS calculation
  wis_by_horizon <- validation_eval |>
    pivot_wider(names_from = output_type_id, values_from = value, names_prefix = "q") |>
    mutate(
      # Interval scores for each interval
      is_50 = pmax(0, pmax(q0.25 - observation, observation - q0.75)) * 2/0.5,
      is_80 = pmax(0, pmax(q0.1 - observation, observation - q0.9)) * 2/0.8,
      is_95 = pmax(0, pmax(q0.025 - observation, observation - q0.975)) * 2/0.95,
      # WIS components
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
  
} else {
  print("No validation forecasts generated")
}