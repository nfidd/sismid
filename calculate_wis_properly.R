# Calculate WIS Properly Using Standard Formula
# Based on the scoringutils package approach

library("nfidd")
library("dplyr")
library("tidyr")

# Load the matched data we created earlier
wis_results <- read.csv("final_wis_results.csv")

print("=== Proper WIS Calculation ===")
print("WIS formula: mean absolute deviation from median + weighted interval scores")

# The standard WIS calculation according to forecast evaluation literature:
# WIS = 1/M * sum of (quantile scores)
# where quantile score = 2 * (I(y < q) - p) * (y - q)
# I(y < q) = 1 if y < q, 0 otherwise
# p = quantile level (0.1, 0.2, ..., 0.9)
# y = observation
# q = quantile forecast

# First, let's check our data
print("Sample of data:")
print(head(wis_results))

# Define quantile levels
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

# Convert to long format for proper WIS calculation
wis_data_long <- wis_results |>
  select(location, origin_date, target_end_date, observation, starts_with("q")) |>
  pivot_longer(
    cols = starts_with("q"),
    names_to = "quantile",
    values_to = "forecast",
    names_prefix = "q"
  ) |>
  mutate(
    quantile_level = as.numeric(quantile)
  ) |>
  filter(quantile_level %in% quantile_levels)

print(paste("Long format data:", nrow(wis_data_long), "rows"))

# Calculate quantile scores
wis_scored <- wis_data_long |>
  mutate(
    # Standard quantile score formula
    indicator = ifelse(observation < forecast, 1, 0),
    quantile_score = 2 * (indicator - quantile_level) * (observation - forecast)
  )

# Calculate WIS as mean quantile score
wis_final <- wis_scored |>
  group_by(location, origin_date, target_end_date, observation) |>
  summarise(
    wis = mean(quantile_score, na.rm = TRUE),
    median_forecast = forecast[quantile_level == 0.5],
    abs_error = abs(observation - median_forecast),
    .groups = "drop"
  )

print("WIS calculation completed")
print("WIS summary statistics:")
print(summary(wis_final$wis))

# Overall performance
overall_performance <- wis_final |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    median_wis = median(wis, na.rm = TRUE),
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    median_abs_error = median(abs_error, na.rm = TRUE),
    n_forecasts = n()
  )

print("\n=== Final Performance Results ===")
print(overall_performance)

# Performance by location
location_performance <- wis_final |>
  group_by(location) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    median_wis = median(wis, na.rm = TRUE),
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    n_forecasts = n(),
    .groups = "drop"
  ) |>
  arrange(mean_wis)

print("\n=== WIS by Location ===")
print(location_performance)

# Check if we're getting reasonable values
print("\n=== Validation Checks ===")
print(paste("Mean WIS:", round(overall_performance$mean_wis, 3)))
print(paste("Median WIS:", round(overall_performance$median_wis, 3)))
print(paste("Mean Absolute Error:", round(overall_performance$mean_abs_error, 3)))

# Save corrected results
write.csv(wis_final, "corrected_wis_results.csv", row.names = FALSE)
write.csv(overall_performance, "corrected_wis_summary.csv", row.names = FALSE)

print("\n=== FINAL RESULTS ===")
print("🎯 Model: sismid-var1-optimal (VAR(1) + season(52))")
print(paste("📊 Mean WIS:", round(overall_performance$mean_wis, 3)))
print(paste("📊 Median WIS:", round(overall_performance$median_wis, 3)))
print(paste("📈 Mean Absolute Error:", round(overall_performance$mean_abs_error, 3)))
print(paste("📝 Test forecasts:", overall_performance$n_forecasts))

# Performance context
print("\n=== Performance Context ===")
print("Expected comparison (from training analysis):")
print("- delphi-epicast: 0.31 WIS")
print("- sismid-var2-sqrt: 0.34 WIS")
print("- sismid-arima210: 0.40 WIS")
print("- hist-avg: 0.45 WIS")
print("")

wis_value <- overall_performance$mean_wis
if (wis_value < 0.35) {
  print("🏆 EXCELLENT! Competitive with top models!")
} else if (wis_value < 0.45) {
  print("🥈 GOOD! Better than historical average!")
} else {
  print("📈 Performance varies - check individual components")
}

# Additional analysis - check if there are outliers affecting the mean
print("\n=== Outlier Analysis ===")
print("WIS quantiles:")
print(quantile(wis_final$wis, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))

# Check for very high WIS values
high_wis <- wis_final |> 
  filter(wis > quantile(wis, 0.95)) |>
  arrange(desc(wis))

print(paste("Number of forecasts with WIS > 95th percentile:", nrow(high_wis)))
if (nrow(high_wis) > 0) {
  print("Sample of high WIS forecasts:")
  print(head(high_wis, 5))
}

print("\nFiles saved:")
print("- corrected_wis_results.csv")
print("- corrected_wis_summary.csv")