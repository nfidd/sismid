# Evaluate WIS Performance of Our Submitted Model
# Check how our optimal VAR(1) + season(52) model performs

library("nfidd")
library("dplyr")
library("tidyr")
library("hubUtils")
library("covidHubUtils")

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

# Calculate WIS
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

# Calculate WIS using covidHubUtils
if (nrow(wis_data) > 0) {
  
  # Prepare data for WIS calculation
  wis_input <- wis_data |>
    select(
      model = model_id,
      location,
      target_end_date,
      target,
      type = output_type,
      quantile = output_type_id,
      value = value_forecast,
      truth = value_target
    ) |>
    mutate(
      type = "quantile",
      quantile = as.numeric(quantile)
    )
  
  print("WIS input data structure:")
  print(head(wis_input))
  
  # Calculate WIS
  tryCatch({
    wis_results <- covidHubUtils::score_forecasts(
      wis_input,
      metrics = "wis"
    )
    
    print("\n=== WIS Results ===")
    print(wis_results)
    
    # Summary statistics
    overall_wis <- wis_results |>
      summarise(
        mean_wis = mean(wis, na.rm = TRUE),
        median_wis = median(wis, na.rm = TRUE),
        sd_wis = sd(wis, na.rm = TRUE),
        n_forecasts = n()
      )
    
    print("\n=== Overall WIS Performance ===")
    print(overall_wis)
    
    # WIS by location
    location_wis <- wis_results |>
      group_by(location) |>
      summarise(
        mean_wis = mean(wis, na.rm = TRUE),
        median_wis = median(wis, na.rm = TRUE),
        n_forecasts = n(),
        .groups = "drop"
      ) |>
      arrange(mean_wis)
    
    print("\n=== WIS by Location ===")
    print(location_wis)
    
    # WIS over time
    time_wis <- wis_results |>
      group_by(target_end_date) |>
      summarise(
        mean_wis = mean(wis, na.rm = TRUE),
        n_forecasts = n(),
        .groups = "drop"
      ) |>
      arrange(target_end_date)
    
    print("\n=== WIS Over Time (First 10 dates) ===")
    print(head(time_wis, 10))
    
    # Save results
    write.csv(wis_results, "optimal_var_wis_results.csv", row.names = FALSE)
    write.csv(overall_wis, "optimal_var_wis_summary.csv", row.names = FALSE)
    
    print("\n=== Final Summary ===")
    print(paste("Model: sismid-var1-optimal (VAR(1) + season(52))"))
    print(paste("Overall WIS:", round(overall_wis$mean_wis, 3)))
    print(paste("Total forecasts evaluated:", overall_wis$n_forecasts))
    
  }, error = function(e) {
    print(paste("Error calculating WIS:", e$message))
    
    # Try alternative approach - simple WIS calculation
    print("Attempting manual WIS calculation...")
    
    # Convert to wide format for quantiles
    wis_wide <- wis_data |>
      select(location, target_end_date, output_type_id, value_forecast, value_target) |>
      pivot_wider(
        names_from = output_type_id,
        values_from = value_forecast,
        names_prefix = "q"
      )
    
    print("Manual WIS calculation structure:")
    print(head(wis_wide))
    
    # Calculate simple coverage and interval scores
    if ("q0.5" %in% names(wis_wide)) {
      simple_scores <- wis_wide |>
        mutate(
          abs_error_median = abs(value_target - q0.5),
          coverage_50 = ifelse(!is.na(q0.25) & !is.na(q0.75), 
                               value_target >= q0.25 & value_target <= q0.75, NA),
          coverage_80 = ifelse(!is.na(q0.1) & !is.na(q0.9), 
                               value_target >= q0.1 & value_target <= q0.9, NA),
          coverage_95 = ifelse(!is.na(q0.025) & !is.na(q0.975), 
                               value_target >= q0.025 & value_target <= q0.975, NA)
        )
      
      summary_scores <- simple_scores |>
        summarise(
          mean_abs_error = mean(abs_error_median, na.rm = TRUE),
          coverage_50_pct = mean(coverage_50, na.rm = TRUE),
          coverage_80_pct = mean(coverage_80, na.rm = TRUE),
          coverage_95_pct = mean(coverage_95, na.rm = TRUE),
          n_forecasts = n()
        )
      
      print("\n=== Simple Forecast Scores ===")
      print(summary_scores)
      
      print(paste("Approximate WIS (mean absolute error):", round(summary_scores$mean_abs_error, 3)))
      print(paste("50% coverage:", round(summary_scores$coverage_50_pct * 100, 1), "%"))
      print(paste("80% coverage:", round(summary_scores$coverage_80_pct * 100, 1), "%"))
      print(paste("95% coverage:", round(summary_scores$coverage_95_pct * 100, 1), "%"))
    }
  })
  
} else {
  print("No matching forecast-target pairs found!")
}

print("\nFiles saved:")
print("- optimal_var_wis_results.csv")
print("- optimal_var_wis_summary.csv")