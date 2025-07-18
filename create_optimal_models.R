# Create and compare the optimal models based on our analysis

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")

set.seed(406)
data(flu_data_hhs)

# Test date
test_date <- as.Date("2018-12-01")

# Prepare data
flu_data_wide <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

flu_data_long <- flu_data_hhs |>
  filter(origin_date <= test_date)

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

fourth_root_transform <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root_transform, inv_fourth_root)

print("=== Comparing Best Transformations ===")

# 1. VAR(2) with sqrt transformation (best VAR from our analysis)
print("1. VAR(2) with sqrt transformation")
var2_sqrt <- flu_data_wide |>
  model(
    var2_sqrt = VAR(vars(
      `HHS Region 1` = my_sqrt(`HHS Region 1`),
      `HHS Region 2` = my_sqrt(`HHS Region 2`),
      `HHS Region 3` = my_sqrt(`HHS Region 3`),
      `HHS Region 4` = my_sqrt(`HHS Region 4`),
      `HHS Region 5` = my_sqrt(`HHS Region 5`),
      `HHS Region 6` = my_sqrt(`HHS Region 6`),
      `HHS Region 7` = my_sqrt(`HHS Region 7`),
      `HHS Region 8` = my_sqrt(`HHS Region 8`),
      `HHS Region 9` = my_sqrt(`HHS Region 9`),
      `HHS Region 10` = my_sqrt(`HHS Region 10`),
      `US National` = my_sqrt(`US National`)
    ) ~ AR(2))
  )
print(paste("VAR(2) sqrt AICc:", round(glance(var2_sqrt)$AICc, 0)))

# 2. VAR(2) with fourth-root transformation (playground style)
print("2. VAR(2) with fourth-root transformation")
var2_fourth <- flu_data_wide |>
  model(
    var2_fourth = VAR(vars(
      `HHS Region 1` = my_fourth_root(`HHS Region 1`),
      `HHS Region 2` = my_fourth_root(`HHS Region 2`),
      `HHS Region 3` = my_fourth_root(`HHS Region 3`),
      `HHS Region 4` = my_fourth_root(`HHS Region 4`),
      `HHS Region 5` = my_fourth_root(`HHS Region 5`),
      `HHS Region 6` = my_fourth_root(`HHS Region 6`),
      `HHS Region 7` = my_fourth_root(`HHS Region 7`),
      `HHS Region 8` = my_fourth_root(`HHS Region 8`),
      `HHS Region 9` = my_fourth_root(`HHS Region 9`),
      `HHS Region 10` = my_fourth_root(`HHS Region 10`),
      `US National` = my_fourth_root(`US National`)
    ) ~ AR(2))
  )
print(paste("VAR(2) fourth-root AICc:", round(glance(var2_fourth)$AICc, 0)))

# 3. ARIMA(2,1,0) with fourth-root (playground approach)
print("3. ARIMA(2,1,0) with fourth-root transformation")
arima_fourth <- flu_data_long |>
  model(
    arima_fourth = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
  )
arima_fourth_metrics <- arima_fourth |>
  glance() |>
  group_by(.model) |>
  summarise(mean_aicc = mean(AICc), .groups = "drop")
print(paste("ARIMA(2,1,0) fourth-root AICc:", round(arima_fourth_metrics$mean_aicc, 0)))

# Model comparison
print("\n=== AICc Comparison ===")
comparison <- data.frame(
  Model = c("VAR(2) sqrt", "VAR(2) fourth-root", "ARIMA(2,1,0) fourth-root"),
  AICc = c(glance(var2_sqrt)$AICc, glance(var2_fourth)$AICc, arima_fourth_metrics$mean_aicc)
) |>
  arrange(AICc)

print(comparison)

best_model <- comparison$Model[1]
print(paste("\nBest model by AICc:", best_model))

# Generate forecasts for the best model
print("\n=== Generating Forecasts for Best Model ===")

if (best_model == "VAR(2) sqrt") {
  print("Generating VAR(2) sqrt forecasts...")
  
  # Generate forecasts
  forecast_samples <- var2_sqrt |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  # Convert to hub format
  location_cols <- c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                     "HHS Region 4", "HHS Region 5", "HHS Region 6",
                     "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                     "HHS Region 10", "US National")
  
  forecast_long <- forecast_samples |>
    pivot_longer(
      cols = all_of(location_cols),
      names_to = "location",
      values_to = "value"
    ) |>
    rename(target_end_date = origin_date) |>
    mutate(
      model_id = "sismid-var2-sqrt",
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
  
} else if (best_model == "VAR(2) fourth-root") {
  print("Generating VAR(2) fourth-root forecasts...")
  
  forecast_samples <- var2_fourth |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  location_cols <- c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                     "HHS Region 4", "HHS Region 5", "HHS Region 6",
                     "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                     "HHS Region 10", "US National")
  
  forecast_long <- forecast_samples |>
    pivot_longer(
      cols = all_of(location_cols),
      names_to = "location",
      values_to = "value"
    ) |>
    rename(target_end_date = origin_date) |>
    mutate(
      model_id = "sismid-var2-fourth",
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
  
} else {
  print("Generating ARIMA(2,1,0) fourth-root forecasts...")
  
  forecast_samples <- arima_fourth |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, location, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  forecast_long <- forecast_samples |>
    rename(
      target_end_date = origin_date,
      value = .sim
    ) |>
    mutate(
      model_id = "sismid-arima-fourth",
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
}

# Calculate quantiles
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

forecast_quantiles <- forecast_long |>
  group_by(model_id, location, origin_date, horizon, target_end_date) |>
  reframe(
    quantile_level = quantile_levels,
    value = quantile(value, quantile_levels, na.rm = TRUE)
  ) |>
  mutate(
    output_type = "quantile",
    output_type_id = quantile_level
  ) |>
  select(-quantile_level)

print(paste("Generated", nrow(forecast_quantiles), "forecast rows"))

# Save the best model forecasts
filename <- paste0("best_model_forecasts_", gsub("[^A-Za-z0-9]", "_", best_model), ".csv")
write.csv(forecast_quantiles, filename, row.names = FALSE)

print(paste("Saved forecasts to:", filename))

# Summary
print("\n=== Summary ===")
print(paste("Best performing model:", best_model))
print(paste("Best AICc:", round(comparison$AICc[1], 0)))
print("")
print("This represents Nick Reich's preferred approach:")
print("- Simple, interpretable models")
print("- Appropriate transformations")
print("- No over-complicated seasonality")
print("- Focus on what works in practice")

print("\n=== Next Steps ===")
print("1. Evaluate WIS performance of this model")
print("2. Compare against our current 0.466 WIS")
print("3. If better, generate full test phase forecasts")
print("4. Submit as our final model")