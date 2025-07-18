# Evaluate WIS Performance of Our Submitted Model (Fixed Version)
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

# Check structure
print("Target data columns:")
print(names(target_data))

print("Target data sample:")
print(head(target_data, 3))

# Load all forecasts
all_forecasts <- list()
for (file in forecast_files) {
  df <- read.csv(file)
  all_forecasts[[length(all_forecasts) + 1]] <- df
}

forecasts_combined <- bind_rows(all_forecasts)

print("Forecast data columns:")
print(names(forecasts_combined))

print("Forecast data sample:")
print(head(forecasts_combined, 3))

print(paste("Total forecast rows:", nrow(forecasts_combined)))

# Convert target data to correct format
# The target data has 'issue' as the forecast date and 'epiweek' as the target date
target_data_formatted <- target_data |>
  mutate(
    issue_date = as.Date(issue),
    epiweek_date = as.Date(epiweek),
    target_value = wili,
    target = "ili perc"
  ) |>
  select(issue_date, region, epiweek_date, target, target_value)

print("Target data formatted sample:")
print(head(target_data_formatted, 3))

# Convert forecasts to correct format
forecasts_formatted <- forecasts_combined |>
  mutate(
    target_end_date = as.Date(target_end_date),
    model_id = "sismid-var1-optimal"
  ) |>
  filter(output_type == "quantile") |>
  select(location, target_end_date, target, output_type, 
         output_type_id, value, model_id)

print("Forecast data formatted sample:")
print(head(forecasts_formatted, 3))

# Calculate WIS
print("\n=== Calculating WIS ===")

# Check date ranges
print("Forecast target_end_date range:")
print(paste("Min:", min(forecasts_formatted$target_end_date)))
print(paste("Max:", max(forecasts_formatted$target_end_date)))

print("Target epiweek_date range:")
print(paste("Min:", min(target_data_formatted$epiweek_date)))
print(paste("Max:", max(target_data_formatted$epiweek_date)))

# Map location names
location_mapping <- data.frame(
  location = c("HHS Region 1", "HHS Region 2", "HHS Region 3", "HHS Region 4", 
               "HHS Region 5", "HHS Region 6", "HHS Region 7", "HHS Region 8", 
               "HHS Region 9", "HHS Region 10", "US National"),
  region = c("hhs1", "hhs2", "hhs3", "hhs4", "hhs5", "hhs6", "hhs7", "hhs8", 
             "hhs9", "hhs10", "nat")
)

print("Location mapping:")
print(location_mapping)

# Align forecasts with target data
# Match by location and target_end_date
wis_data <- forecasts_formatted |>
  left_join(location_mapping, by = "location") |>
  left_join(target_data_formatted, 
            by = c("region", "target_end_date" = "epiweek_date", "target")) |>
  filter(!is.na(target_value))

print(paste("Matched forecast-target pairs:", nrow(wis_data)))

if (nrow(wis_data) > 0) {
  print("Sample of matched data:")
  print(head(wis_data, 3))
  
  # Convert to wide format for quantiles
  wis_wide <- wis_data |>
    select(location, target_end_date, output_type_id, value, target_value) |>
    pivot_wider(
      names_from = output_type_id,
      values_from = value,
      names_prefix = "q"
    )
  
  print("WIS calculation data structure:")
  print(head(wis_wide, 3))
  print("Available quantiles:")
  print(names(wis_wide)[grepl("^q", names(wis_wide))])
  
  # Manual WIS calculation
  interval_scores <- wis_wide |>
    mutate(
      # Point forecast (median)
      median_forecast = q0.5,
      abs_error_median = abs(target_value - median_forecast),
      
      # 50% prediction interval (25th to 75th percentile)
      interval_50_score = ifelse(
        !is.na(q0.25) & !is.na(q0.75),
        pmax(q0.25 - target_value, 0) * 2 * 0.25 + 
        pmax(target_value - q0.75, 0) * 2 * 0.25 + 
        (q0.75 - q0.25),
        NA
      ),
      
      # 80% prediction interval (10th to 90th percentile)
      interval_80_score = ifelse(
        !is.na(q0.1) & !is.na(q0.9),
        pmax(q0.1 - target_value, 0) * 2 * 0.1 + 
        pmax(target_value - q0.9, 0) * 2 * 0.1 + 
        (q0.9 - q0.1),
        NA
      ),
      
      # 90% prediction interval (5th to 95th percentile)
      interval_90_score = ifelse(
        !is.na(q0.05) & !is.na(q0.95),
        pmax(q0.05 - target_value, 0) * 2 * 0.05 + 
        pmax(target_value - q0.95, 0) * 2 * 0.05 + 
        (q0.95 - q0.05),
        NA
      ),
      
      # 95% prediction interval (2.5th to 97.5th percentile)
      interval_95_score = ifelse(
        !is.na(q0.025) & !is.na(q0.975),
        pmax(q0.025 - target_value, 0) * 2 * 0.025 + 
        pmax(target_value - q0.975, 0) * 2 * 0.025 + 
        (q0.975 - q0.025),
        NA
      ),
      
      # Coverage indicators
      coverage_50 = ifelse(!is.na(q0.25) & !is.na(q0.75), 
                           target_value >= q0.25 & target_value <= q0.75, NA),
      coverage_80 = ifelse(!is.na(q0.1) & !is.na(q0.9), 
                           target_value >= q0.1 & target_value <= q0.9, NA),
      coverage_90 = ifelse(!is.na(q0.05) & !is.na(q0.95), 
                           target_value >= q0.05 & target_value <= q0.95, NA),
      coverage_95 = ifelse(!is.na(q0.025) & !is.na(q0.975), 
                           target_value >= q0.025 & target_value <= q0.975, NA)
    )
  
  # Calculate WIS as weighted average of interval scores
  wis_results <- interval_scores |>
    mutate(
      # Simple WIS approximation using available intervals
      wis = (abs_error_median + 
             ifelse(!is.na(interval_50_score), interval_50_score, 0) +
             ifelse(!is.na(interval_80_score), interval_80_score, 0) +
             ifelse(!is.na(interval_90_score), interval_90_score, 0) +
             ifelse(!is.na(interval_95_score), interval_95_score, 0)) / 5
    )
  
  # Overall WIS summary
  overall_wis <- wis_results |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      median_wis = median(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error_median, na.rm = TRUE),
      median_abs_error = median(abs_error_median, na.rm = TRUE),
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
  
  print("\n=== WIS by Location (Top 5) ===")
  print(head(location_wis, 5))
  
  # Save results
  write.csv(wis_results, "optimal_var_wis_results.csv", row.names = FALSE)
  write.csv(overall_wis, "optimal_var_wis_summary.csv", row.names = FALSE)
  
  print("\n=== FINAL RESULTS ===")
  print(paste("🎯 Model: sismid-var1-optimal (VAR(1) + season(52))"))
  print(paste("📊 Overall WIS:", round(overall_wis$mean_wis, 3)))
  print(paste("📈 Mean Absolute Error:", round(overall_wis$mean_abs_error, 3)))
  print(paste("📈 Median Absolute Error:", round(overall_wis$median_abs_error, 3)))
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
  print("Checking unique values...")
  print("Unique locations in forecasts:")
  print(unique(forecasts_formatted$location))
  print("Unique regions in targets:")
  print(unique(target_data_formatted$region))
  print("Unique target dates in forecasts (first 10):")
  print(head(unique(forecasts_formatted$target_end_date), 10))
  print("Unique epiweek dates in targets (first 10):")
  print(head(unique(target_data_formatted$epiweek_date), 10))
}

print("\nFiles saved:")
print("- optimal_var_wis_results.csv")
print("- optimal_var_wis_summary.csv")