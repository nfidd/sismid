# Evaluate WIS Performance of Our Submitted Model (Simple Version)
# Check how our optimal VAR(1) + season(52) model performs

library("nfidd")
library("dplyr")
library("tidyr")
library("hubUtils")

# Hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Load target data
target_data <- read.csv(file.path(hub_path, "target-data", "target-data-raw.csv"))

# Load our model forecasts
model_folder <- file.path(hub_path, "model-output", "sismid-var1-optimal")
forecast_files <- list.files(model_folder, pattern = "\\.csv$", full.names = TRUE)

print(paste("Found", length(forecast_files), "forecast files"))

# Load all forecasts
all_forecasts <- list()
for (file in forecast_files) {
  df <- read.csv(file)
  all_forecasts[[length(all_forecasts) + 1]] <- df
}

forecasts_combined <- bind_rows(all_forecasts)

print(paste("Total forecast rows:", nrow(forecasts_combined)))

# Convert target data to correct format
target_data_formatted <- target_data |>
  mutate(
    origin_date = as.Date(origin_date),
    target_end_date = as.Date(target_end_date),
    value = wili,
    target = "ili perc"
  ) |>
  select(origin_date, location, target_end_date, target, value)

# Convert forecasts to correct format
forecasts_formatted <- forecasts_combined |>
  mutate(
    origin_date = as.Date(origin_date),
    target_end_date = as.Date(target_end_date),
    model_id = "sismid-var1-optimal"
  ) |>
  filter(output_type == "quantile") |>
  select(origin_date, location, target_end_date, target, output_type, 
         output_type_id, value, model_id)

print("Sample of forecasts:")
print(head(forecasts_formatted))

print("Sample of target data:")
print(head(target_data_formatted))

# Calculate WIS manually
print("\n=== Calculating WIS ===")

# Check date ranges
print("Forecast date range:")
print(paste("Min:", min(forecasts_formatted$origin_date)))
print(paste("Max:", max(forecasts_formatted$origin_date)))

print("Target date range:")
print(paste("Min:", min(target_data_formatted$origin_date)))
print(paste("Max:", max(target_data_formatted$origin_date)))

# Align forecasts with target data
# For each forecast, find the matching target data
wis_data <- forecasts_formatted |>
  left_join(target_data_formatted, 
            by = c("location", "target_end_date", "target"),
            suffix = c("_forecast", "_target")) |>
  filter(!is.na(value_target))

print(paste("Matched forecast-target pairs:", nrow(wis_data)))

# Calculate WIS manually
if (nrow(wis_data) > 0) {
  
  # Convert to wide format for quantiles
  wis_wide <- wis_data |>
    select(location, target_end_date, output_type_id, value_forecast, value_target) |>
    pivot_wider(
      names_from = output_type_id,
      values_from = value_forecast,
      names_prefix = "q"
    )
  
  print("WIS calculation data structure:")
  print(head(wis_wide))
  print("Available quantiles:")
  print(names(wis_wide)[grepl("^q", names(wis_wide))])
  
  # Manual WIS calculation
  # WIS = (1/M) * sum of interval scores + absolute error to median
  
  # Define quantile levels
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  # Calculate interval scores for each prediction interval
  interval_scores <- wis_wide |>
    mutate(
      # Point forecast (median)
      median_forecast = q0.5,
      abs_error_median = abs(value_target - median_forecast),
      
      # 50% prediction interval (25th to 75th percentile)
      interval_50_score = ifelse(
        !is.na(q0.25) & !is.na(q0.75),
        pmax(q0.25 - value_target, 0) * 2 * 0.25 + 
        pmax(value_target - q0.75, 0) * 2 * 0.25 + 
        (q0.75 - q0.25),
        NA
      ),
      
      # 80% prediction interval (10th to 90th percentile)
      interval_80_score = ifelse(
        !is.na(q0.1) & !is.na(q0.9),
        pmax(q0.1 - value_target, 0) * 2 * 0.1 + 
        pmax(value_target - q0.9, 0) * 2 * 0.1 + 
        (q0.9 - q0.1),
        NA
      ),
      
      # 90% prediction interval (5th to 95th percentile)
      interval_90_score = ifelse(
        !is.na(q0.05) & !is.na(q0.95),
        pmax(q0.05 - value_target, 0) * 2 * 0.05 + 
        pmax(value_target - q0.95, 0) * 2 * 0.05 + 
        (q0.95 - q0.05),
        NA
      ),
      
      # 95% prediction interval (2.5th to 97.5th percentile)
      interval_95_score = ifelse(
        !is.na(q0.025) & !is.na(q0.975),
        pmax(q0.025 - value_target, 0) * 2 * 0.025 + 
        pmax(value_target - q0.975, 0) * 2 * 0.025 + 
        (q0.975 - q0.025),
        NA
      ),
      
      # 98% prediction interval (1st to 99th percentile)
      interval_98_score = ifelse(
        !is.na(q0.01) & !is.na(q0.99),
        pmax(q0.01 - value_target, 0) * 2 * 0.01 + 
        pmax(value_target - q0.99, 0) * 2 * 0.01 + 
        (q0.99 - q0.01),
        NA
      ),
      
      # Coverage indicators
      coverage_50 = ifelse(!is.na(q0.25) & !is.na(q0.75), 
                           value_target >= q0.25 & value_target <= q0.75, NA),
      coverage_80 = ifelse(!is.na(q0.1) & !is.na(q0.9), 
                           value_target >= q0.1 & value_target <= q0.9, NA),
      coverage_90 = ifelse(!is.na(q0.05) & !is.na(q0.95), 
                           value_target >= q0.05 & value_target <= q0.95, NA),
      coverage_95 = ifelse(!is.na(q0.025) & !is.na(q0.975), 
                           value_target >= q0.025 & value_target <= q0.975, NA)
    )
  
  # Calculate WIS as weighted average of interval scores
  wis_results <- interval_scores |>
    mutate(
      # Approximate WIS using available intervals
      wis = (abs_error_median + 
             ifelse(!is.na(interval_50_score), interval_50_score * 0.5, 0) +
             ifelse(!is.na(interval_80_score), interval_80_score * 0.8, 0) +
             ifelse(!is.na(interval_90_score), interval_90_score * 0.9, 0) +
             ifelse(!is.na(interval_95_score), interval_95_score * 0.95, 0) +
             ifelse(!is.na(interval_98_score), interval_98_score * 0.98, 0)) / 5
    )
  
  # Overall WIS summary
  overall_wis <- wis_results |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      median_wis = median(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error_median, na.rm = TRUE),
      coverage_50_pct = mean(coverage_50, na.rm = TRUE) * 100,
      coverage_80_pct = mean(coverage_80, na.rm = TRUE) * 100,
      coverage_90_pct = mean(coverage_90, na.rm = TRUE) * 100,
      coverage_95_pct = mean(coverage_95, na.rm = TRUE) * 100,
      n_forecasts = n()
    )
  
  print("\n=== Overall WIS Performance ===")
  print(overall_wis)
  
  # WIS by location
  location_wis <- wis_results |>
    group_by(location) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error_median, na.rm = TRUE),
      coverage_80_pct = mean(coverage_80, na.rm = TRUE) * 100,
      n_forecasts = n(),
      .groups = "drop"
    ) |>
    arrange(mean_wis)
  
  print("\n=== WIS by Location ===")
  print(location_wis)
  
  # Save results
  write.csv(wis_results, "optimal_var_wis_results.csv", row.names = FALSE)
  write.csv(overall_wis, "optimal_var_wis_summary.csv", row.names = FALSE)
  
  print("\n=== FINAL RESULTS ===")
  print(paste("🎯 Model: sismid-var1-optimal (VAR(1) + season(52))"))
  print(paste("📊 Overall WIS:", round(overall_wis$mean_wis, 3)))
  print(paste("📈 Mean Absolute Error:", round(overall_wis$mean_abs_error, 3)))
  print(paste("🎯 50% Coverage:", round(overall_wis$coverage_50_pct, 1), "%"))
  print(paste("🎯 80% Coverage:", round(overall_wis$coverage_80_pct, 1), "%"))
  print(paste("🎯 90% Coverage:", round(overall_wis$coverage_90_pct, 1), "%"))
  print(paste("🎯 95% Coverage:", round(overall_wis$coverage_95_pct, 1), "%"))
  print(paste("📝 Total forecasts evaluated:", overall_wis$n_forecasts))
  
  # Compare with baseline expectation
  print("\n=== Performance Context ===")
  print("Expected WIS ranges:")
  print("- Excellent: < 0.3")
  print("- Good: 0.3 - 0.4")
  print("- Acceptable: 0.4 - 0.5")
  print("- Poor: > 0.5")
  
  if (overall_wis$mean_wis < 0.3) {
    print("🏆 EXCELLENT PERFORMANCE!")
  } else if (overall_wis$mean_wis < 0.4) {
    print("🥈 GOOD PERFORMANCE!")
  } else if (overall_wis$mean_wis < 0.5) {
    print("✅ ACCEPTABLE PERFORMANCE")
  } else {
    print("⚠️ NEEDS IMPROVEMENT")
  }
  
} else {
  print("❌ No matching forecast-target pairs found!")
}

print("\nFiles saved:")
print("- optimal_var_wis_results.csv")
print("- optimal_var_wis_summary.csv")