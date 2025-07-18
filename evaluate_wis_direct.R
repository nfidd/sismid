# Direct WIS Evaluation of Our Optimal VAR Model
# Using the time-series.csv target data directly

library("nfidd")
library("dplyr")
library("tidyr")

# Set seed
set.seed(406)

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Load our model forecasts
model_folder <- file.path(hub_path, "model-output", "sismid-var1-optimal")
forecast_files <- list.files(model_folder, pattern = "\\.csv$", full.names = TRUE)

print(paste("Found", length(forecast_files), "forecast files"))

# Load all forecasts
all_forecasts <- list()
for (file in forecast_files) {
  df <- read.csv(file)
  # Add origin date from filename
  origin_date <- gsub(".*(\\d{4}-\\d{2}-\\d{2}).*", "\\1", basename(file))
  df$origin_date <- as.Date(origin_date)
  all_forecasts[[length(all_forecasts) + 1]] <- df
}

forecasts_combined <- bind_rows(all_forecasts)
print(paste("Total forecast rows:", nrow(forecasts_combined)))

# Load target data (time-series.csv)
target_data <- read.csv(file.path(hub_path, "target-data", "time-series.csv"))
target_data$target_end_date <- as.Date(target_data$target_end_date)

print(paste("Target data rows:", nrow(target_data)))
print("Target data date range:")
print(paste("Min:", min(target_data$target_end_date)))
print(paste("Max:", max(target_data$target_end_date)))

# Filter forecasts to quantile type only
forecasts_quantile <- forecasts_combined |>
  filter(output_type == "quantile") |>
  mutate(target_end_date = as.Date(target_end_date))

print(paste("Quantile forecast rows:", nrow(forecasts_quantile)))

# Match forecasts with target data
matched_data <- forecasts_quantile |>
  left_join(target_data, 
            by = c("location", "target_end_date", "target"),
            suffix = c("_forecast", "_actual")) |>
  filter(!is.na(observation))

print(paste("Matched forecast-target pairs:", nrow(matched_data)))

if (nrow(matched_data) > 0) {
  print("Sample of matched data:")
  print(head(matched_data))
  
  # Calculate WIS manually
  # Convert to wide format for easier calculation
  wis_wide <- matched_data |>
    select(origin_date, location, target_end_date, output_type_id, value, observation) |>
    pivot_wider(
      names_from = output_type_id,
      values_from = value,
      names_prefix = "q"
    )
  
  print("\nQuantiles available:")
  print(names(wis_wide)[grepl("^q", names(wis_wide))])
  
  # Calculate interval scores and coverage
  wis_results <- wis_wide |>
    mutate(
      # Point forecast (median)
      median_forecast = q0.5,
      abs_error = abs(observation - median_forecast),
      
      # Interval scores (penalty for being outside + interval width)
      # 50% interval
      interval_50 = ifelse(
        !is.na(q0.25) & !is.na(q0.75),
        (q0.75 - q0.25) + 
        (2/0.5) * pmax(q0.25 - observation, 0) + 
        (2/0.5) * pmax(observation - q0.75, 0),
        NA
      ),
      
      # 80% interval  
      interval_80 = ifelse(
        !is.na(q0.1) & !is.na(q0.9),
        (q0.9 - q0.1) + 
        (2/0.2) * pmax(q0.1 - observation, 0) + 
        (2/0.2) * pmax(observation - q0.9, 0),
        NA
      ),
      
      # 90% interval
      interval_90 = ifelse(
        !is.na(q0.05) & !is.na(q0.95),
        (q0.95 - q0.05) + 
        (2/0.1) * pmax(q0.05 - observation, 0) + 
        (2/0.1) * pmax(observation - q0.95, 0),
        NA
      ),
      
      # 96% interval
      interval_96 = ifelse(
        !is.na(q0.025) & !is.na(q0.975),
        (q0.975 - q0.025) + 
        (2/0.05) * pmax(q0.025 - observation, 0) + 
        (2/0.05) * pmax(observation - q0.975, 0),
        NA
      ),
      
      # 98% interval
      interval_98 = ifelse(
        !is.na(q0.01) & !is.na(q0.99),
        (q0.99 - q0.01) + 
        (2/0.02) * pmax(q0.01 - observation, 0) + 
        (2/0.02) * pmax(observation - q0.99, 0),
        NA
      ),
      
      # Coverage indicators
      coverage_50 = ifelse(!is.na(q0.25) & !is.na(q0.75), 
                           observation >= q0.25 & observation <= q0.75, NA),
      coverage_80 = ifelse(!is.na(q0.1) & !is.na(q0.9), 
                           observation >= q0.1 & observation <= q0.9, NA),
      coverage_90 = ifelse(!is.na(q0.05) & !is.na(q0.95), 
                           observation >= q0.05 & observation <= q0.95, NA),
      coverage_96 = ifelse(!is.na(q0.025) & !is.na(q0.975), 
                           observation >= q0.025 & observation <= q0.975, NA),
      coverage_98 = ifelse(!is.na(q0.01) & !is.na(q0.99), 
                           observation >= q0.01 & observation <= q0.99, NA)
    )
  
  # Calculate WIS as weighted average of interval scores
  # WIS = (1/K) * sum of all interval scores, where K is the number of intervals
  wis_final <- wis_results |>
    mutate(
      # Count available intervals
      n_intervals = (!is.na(interval_50)) + (!is.na(interval_80)) + 
                   (!is.na(interval_90)) + (!is.na(interval_96)) + 
                   (!is.na(interval_98)),
      
      # Calculate WIS (using standard formula)
      wis = (abs_error + 
             ifelse(!is.na(interval_50), interval_50 * 0.5, 0) +
             ifelse(!is.na(interval_80), interval_80 * 0.8, 0) +
             ifelse(!is.na(interval_90), interval_90 * 0.9, 0) +
             ifelse(!is.na(interval_96), interval_96 * 0.96, 0) +
             ifelse(!is.na(interval_98), interval_98 * 0.98, 0)) / 
            pmax(n_intervals, 1)
    )
  
  # Overall performance summary
  overall_performance <- wis_final |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      median_wis = median(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      median_abs_error = median(abs_error, na.rm = TRUE),
      coverage_50_pct = mean(coverage_50, na.rm = TRUE) * 100,
      coverage_80_pct = mean(coverage_80, na.rm = TRUE) * 100,
      coverage_90_pct = mean(coverage_90, na.rm = TRUE) * 100,
      coverage_96_pct = mean(coverage_96, na.rm = TRUE) * 100,
      coverage_98_pct = mean(coverage_98, na.rm = TRUE) * 100,
      n_forecasts = n()
    )
  
  print("\n=== Overall Performance ===")
  print(overall_performance)
  
  # Performance by location
  location_performance <- wis_final |>
    group_by(location) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      coverage_80_pct = mean(coverage_80, na.rm = TRUE) * 100,
      n_forecasts = n(),
      .groups = "drop"
    ) |>
    arrange(mean_wis)
  
  print("\n=== WIS by Location ===")
  print(location_performance)
  
  # Performance over time (monthly aggregation)
  time_performance <- wis_final |>
    mutate(
      year_month = format(origin_date, "%Y-%m")
    ) |>
    group_by(year_month) |>
    summarise(
      mean_wis = mean(wis, na.rm = TRUE),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      n_forecasts = n(),
      .groups = "drop"
    ) |>
    arrange(year_month)
  
  print("\n=== WIS Over Time (Monthly) ===")
  print(time_performance)
  
  # Save results
  write.csv(wis_final, "final_wis_results.csv", row.names = FALSE)
  write.csv(overall_performance, "final_wis_summary.csv", row.names = FALSE)
  
  print("\n=== FINAL RESULTS ===")
  print("🎯 Model: sismid-var1-optimal (VAR(1) + season(52))")
  print(paste("📊 Mean WIS:", round(overall_performance$mean_wis, 3)))
  print(paste("📊 Median WIS:", round(overall_performance$median_wis, 3)))
  print(paste("📈 Mean Absolute Error:", round(overall_performance$mean_abs_error, 3)))
  print(paste("📈 Median Absolute Error:", round(overall_performance$median_abs_error, 3)))
  print(paste("🎯 50% Coverage:", round(overall_performance$coverage_50_pct, 1), "%"))
  print(paste("🎯 80% Coverage:", round(overall_performance$coverage_80_pct, 1), "%"))
  print(paste("🎯 90% Coverage:", round(overall_performance$coverage_90_pct, 1), "%"))
  print(paste("🎯 96% Coverage:", round(overall_performance$coverage_96_pct, 1), "%"))
  print(paste("🎯 98% Coverage:", round(overall_performance$coverage_98_pct, 1), "%"))
  print(paste("📝 Total forecasts:", overall_performance$n_forecasts))
  
  # Performance context
  print("\n=== Performance Context ===")
  print("Previous model comparison (Training Phase):")
  print("- delphi-epicast: 0.31 WIS (1st place)")
  print("- sismid-var2-sqrt: 0.34 WIS (2nd place)")
  print("- sismid-arima210: 0.40 WIS (3rd place)")
  print("- hist-avg: 0.45 WIS (4th place)")
  print("")
  print("Our model performance:")
  
  if (overall_performance$mean_wis < 0.31) {
    print("🏆 OUTSTANDING! Better than delphi-epicast!")
  } else if (overall_performance$mean_wis < 0.34) {
    print("🥇 EXCELLENT! Competitive with best models!")
  } else if (overall_performance$mean_wis < 0.40) {
    print("🥈 VERY GOOD! Better than ARIMA baseline!")
  } else if (overall_performance$mean_wis < 0.45) {
    print("🥉 GOOD! Competitive performance!")
  } else {
    print("📈 Room for improvement")
  }
  
} else {
  print("❌ No matching data found!")
}

print("\nFiles saved:")
print("- final_wis_results.csv")
print("- final_wis_summary.csv")