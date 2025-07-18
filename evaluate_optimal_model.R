# Evaluate Our Optimal VAR Model Performance
# Using the proper hub evaluation approach

library("nfidd")
library("dplyr")
library("hubUtils")
library("hubEvals")

# Set seed
set.seed(406)

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")
hub_con <- connect_hub(hub_path)

# Get all our model forecasts
our_forecasts <- hub_con |> 
  filter(model_id == "sismid-var1-optimal") |> 
  collect_hub()

print(paste("Loaded", nrow(our_forecasts), "forecast rows"))
print("Sample of our forecasts:")
print(head(our_forecasts))

# Get origin dates from TEST PHASE (seasons 3-5)
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

test_origin_dates <- origin_dates[which(as.Date(origin_dates) > as.Date("2017-05-06"))]
print(paste("Test phase origin dates:", length(test_origin_dates)))
print(paste("Date range:", min(test_origin_dates), "to", max(test_origin_dates)))

# Load target data
target_data <- load_target_data(hub_con)
print(paste("Loaded", nrow(target_data), "target data rows"))
print("Sample of target data:")
print(head(target_data))

# Calculate WIS for our model
print("\n=== Calculating WIS for Test Phase ===")

# Filter our forecasts to test phase only
test_forecasts <- our_forecasts |>
  filter(origin_date %in% test_origin_dates)

print(paste("Test phase forecasts:", nrow(test_forecasts)))

# Calculate scores
tryCatch({
  
  # Use hubEvals to calculate WIS
  wis_scores <- test_forecasts |>
    score_forecasts(
      target_data = target_data,
      metrics = "wis"
    )
  
  print("WIS calculation successful!")
  print(head(wis_scores))
  
  # Overall performance
  overall_wis <- wis_scores |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      median_wis = median(wis, na.rm = TRUE),
      n_forecasts = n()
    )
  
  print("\n=== Overall Test Phase Performance ===")
  print(overall_wis)
  
  # Performance by location
  location_wis <- wis_scores |>
    group_by(location) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      n_forecasts = n(),
      .groups = "drop"
    ) |>
    arrange(mean_wis)
  
  print("\n=== WIS by Location ===")
  print(location_wis)
  
  # Performance over time
  time_wis <- wis_scores |>
    group_by(origin_date) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      n_forecasts = n(),
      .groups = "drop"
    ) |>
    arrange(origin_date)
  
  print("\n=== WIS Over Time (First 10 dates) ===")
  print(head(time_wis, 10))
  
  # Save results
  write.csv(wis_scores, "optimal_model_wis_scores.csv", row.names = FALSE)
  write.csv(overall_wis, "optimal_model_wis_summary.csv", row.names = FALSE)
  
  print("\n=== FINAL RESULTS ===")
  print(paste("🎯 Model: sismid-var1-optimal"))
  print(paste("📊 Mean WIS:", round(overall_wis$mean_wis, 3)))
  print(paste("📊 Median WIS:", round(overall_wis$median_wis, 3)))
  print(paste("📝 Test forecasts:", overall_wis$n_forecasts))
  
  # Compare with reference
  print("\n=== Performance Context ===")
  print("Previous analysis showed:")
  print("- delphi-epicast: 0.31 WIS")
  print("- sismid-var2-sqrt: 0.34 WIS")
  print("- sismid-arima210: 0.40 WIS")
  print("- hist-avg: 0.45 WIS")
  
  if (overall_wis$mean_wis < 0.32) {
    print("🏆 EXCELLENT! Better than delphi-epicast!")
  } else if (overall_wis$mean_wis < 0.35) {
    print("🥈 VERY GOOD! Similar to best VAR model!")
  } else if (overall_wis$mean_wis < 0.41) {
    print("✅ GOOD! Better than ARIMA baseline!")
  } else {
    print("📈 Room for improvement")
  }
  
}, error = function(e) {
  print(paste("Error in WIS calculation:", e$message))
  print("Trying alternative approach...")
  
  # Alternative approach - use the time-series.csv file directly
  time_series_data <- read.csv(file.path(hub_path, "target-data", "time-series.csv"))
  
  # Simple manual calculation
  test_forecasts_simple <- test_forecasts |>
    filter(output_type == "quantile") |>
    select(origin_date, location, target_end_date, output_type_id, value)
  
  # Match with observations
  matched_data <- test_forecasts_simple |>
    left_join(time_series_data, 
              by = c("location", "target_end_date"),
              suffix = c("_forecast", "_actual")) |>
    filter(!is.na(observation))
  
  print(paste("Matched data points:", nrow(matched_data)))
  
  if (nrow(matched_data) > 0) {
    # Calculate simple metrics
    median_forecasts <- matched_data |>
      filter(output_type_id == 0.5) |>
      mutate(abs_error = abs(value_forecast - observation))
    
    simple_performance <- median_forecasts |>
      summarise(
        mean_abs_error = mean(abs_error, na.rm = TRUE),
        median_abs_error = median(abs_error, na.rm = TRUE),
        n_forecasts = n()
      )
    
    print("\n=== Simple Performance Metrics ===")
    print(simple_performance)
    
    print(paste("📊 Mean Absolute Error:", round(simple_performance$mean_abs_error, 3)))
    print(paste("📊 Median Absolute Error:", round(simple_performance$median_abs_error, 3)))
    print(paste("📝 Test forecasts:", simple_performance$n_forecasts))
    
    # Approximate WIS as 2 * MAE (rough approximation)
    approx_wis <- 2 * simple_performance$mean_abs_error
    print(paste("📊 Approximate WIS (2 * MAE):", round(approx_wis, 3)))
  }
})

print("\nFiles saved:")
print("- optimal_model_wis_scores.csv")
print("- optimal_model_wis_summary.csv")