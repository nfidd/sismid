# Debug WIS Calculation - Check what's causing high WIS values
# Let's examine the results more carefully

library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")

# Load the results we just calculated
wis_results <- read.csv("final_wis_results.csv")

print("=== WIS Results Summary ===")
print(summary(wis_results$wis))
print("Quantiles:")
print(quantile(wis_results$wis, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))

# Check for outliers
print("\n=== High WIS Values ===")
high_wis <- wis_results |> 
  filter(wis > 10) |>
  select(location, origin_date, target_end_date, wis, abs_error, observation, median_forecast)

print(paste("Number of forecasts with WIS > 10:", nrow(high_wis)))
print(head(high_wis, 10))

# Check the distribution by location
print("\n=== WIS Distribution by Location ===")
location_stats <- wis_results |>
  group_by(location) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    median_wis = median(wis, na.rm = TRUE),
    q95_wis = quantile(wis, 0.95, na.rm = TRUE),
    max_wis = max(wis, na.rm = TRUE),
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    n_forecasts = n()
  ) |>
  arrange(desc(mean_wis))

print(location_stats)

# Let's check the problem region (HHS Region 6)
print("\n=== HHS Region 6 Analysis ===")
region6 <- wis_results |> 
  filter(location == "HHS Region 6") |>
  arrange(desc(wis))

print("Top 10 worst forecasts for HHS Region 6:")
print(head(region6, 10))

# Check if there are any zero or negative target values
print("\n=== Target Value Analysis ===")
print("Target value summary:")
print(summary(wis_results$observation))

print("Number of zero targets:")
print(sum(wis_results$observation == 0, na.rm = TRUE))

print("Number of very small targets (< 0.1):")
print(sum(wis_results$observation < 0.1, na.rm = TRUE))

# Check if intervals are reasonable
print("\n=== Interval Analysis ===")
interval_analysis <- wis_results |>
  mutate(
    interval_50_width = q0.75 - q0.25,
    interval_80_width = q0.9 - q0.1,
    interval_90_width = q0.95 - q0.05
  )

print("Interval width summary:")
print(summary(interval_analysis$interval_50_width))
print(summary(interval_analysis$interval_80_width))
print(summary(interval_analysis$interval_90_width))

# Check for negative forecasts
print("\n=== Forecast Value Analysis ===")
print("Median forecast summary:")
print(summary(wis_results$median_forecast))

print("Number of negative median forecasts:")
print(sum(wis_results$median_forecast < 0, na.rm = TRUE))

# Let's recalculate WIS using the standard formula more carefully
print("\n=== Recalculating WIS with Proper Formula ===")

# Standard WIS formula: WIS = 1/2 * sum over all quantiles of interval scores
# For each quantile level α, interval score = width + 2/α * overprediction + 2/(1-α) * underprediction

# First, let's try a simpler approach - calculate WIS as average of quantile scores
quantile_levels <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
                    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

# Convert to long format for quantile score calculation
wis_long <- wis_results |>
  select(location, origin_date, target_end_date, observation, starts_with("q")) |>
  pivot_longer(cols = starts_with("q"), names_to = "quantile", values_to = "forecast") |>
  mutate(
    quantile_level = as.numeric(gsub("q", "", quantile))
  ) |>
  filter(quantile_level %in% quantile_levels)

# Calculate quantile scores
wis_long_scored <- wis_long |>
  mutate(
    # Quantile score = 2 * (indicator - quantile_level) * (observation - forecast)
    # where indicator = 1 if observation < forecast, 0 otherwise
    indicator = ifelse(observation < forecast, 1, 0),
    quantile_score = 2 * (indicator - quantile_level) * (observation - forecast)
  )

# Calculate WIS as mean of quantile scores
wis_corrected <- wis_long_scored |>
  group_by(location, origin_date, target_end_date, observation) |>
  summarise(
    wis_corrected = mean(quantile_score, na.rm = TRUE),
    median_forecast = forecast[quantile_level == 0.5],
    abs_error = abs(observation - median_forecast),
    .groups = "drop"
  )

print("Corrected WIS summary:")
print(summary(wis_corrected$wis_corrected))

corrected_overall <- wis_corrected |>
  summarise(
    mean_wis = mean(wis_corrected, na.rm = TRUE),
    median_wis = median(wis_corrected, na.rm = TRUE),
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    n_forecasts = n()
  )

print("\n=== Corrected Overall Performance ===")
print(corrected_overall)

# Compare with previous models
print("\n=== Performance Comparison ===")
print("Previous models (training phase):")
print("- delphi-epicast: 0.31 WIS")
print("- sismid-var2-sqrt: 0.34 WIS")
print("- sismid-arima210: 0.40 WIS")
print("- hist-avg: 0.45 WIS")
print("")
print(paste("Our model (corrected):", round(corrected_overall$mean_wis, 3), "WIS"))

if (corrected_overall$mean_wis < 0.35) {
  print("🏆 EXCELLENT! Competitive with best models!")
} else if (corrected_overall$mean_wis < 0.45) {
  print("🥈 GOOD! Better than historical average!")
} else {
  print("📈 Performance needs analysis")
}