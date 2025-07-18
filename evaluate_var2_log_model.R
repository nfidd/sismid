# Evaluate VAR(2) log model - the best performing model from our comparison

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")

# Set seed for reproducibility
set.seed(406)

# Load data
data(flu_data_hhs)

# Use the same approach as our successful seasonal model generation
# Test a few sample forecasts to see performance

# Define log transformation
log_transform <- function(x) log(x + 0.1)
inv_log <- function(x) exp(x) - 0.1
my_log <- new_transformation(log_transform, inv_log)

# Function to generate single forecast and evaluate
generate_and_evaluate_forecast <- function(origin_date, model_name) {
  
  # Prepare data up to origin date
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= origin_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Fit model
  if (model_name == "var2_log") {
    model_fit <- flu_data_wide |>
      model(
        model = VAR(vars(
          `HHS Region 1` = my_log(`HHS Region 1`),
          `HHS Region 2` = my_log(`HHS Region 2`),
          `HHS Region 3` = my_log(`HHS Region 3`),
          `HHS Region 4` = my_log(`HHS Region 4`),
          `HHS Region 5` = my_log(`HHS Region 5`),
          `HHS Region 6` = my_log(`HHS Region 6`),
          `HHS Region 7` = my_log(`HHS Region 7`),
          `HHS Region 8` = my_log(`HHS Region 8`),
          `HHS Region 9` = my_log(`HHS Region 9`),
          `HHS Region 10` = my_log(`HHS Region 10`),
          `US National` = my_log(`US National`)
        ) ~ AR(2))
      )
  } else if (model_name == "var1_log") {
    model_fit <- flu_data_wide |>
      model(
        model = VAR(vars(
          `HHS Region 1` = my_log(`HHS Region 1`),
          `HHS Region 2` = my_log(`HHS Region 2`),
          `HHS Region 3` = my_log(`HHS Region 3`),
          `HHS Region 4` = my_log(`HHS Region 4`),
          `HHS Region 5` = my_log(`HHS Region 5`),
          `HHS Region 6` = my_log(`HHS Region 6`),
          `HHS Region 7` = my_log(`HHS Region 7`),
          `HHS Region 8` = my_log(`HHS Region 8`),
          `HHS Region 9` = my_log(`HHS Region 9`),
          `HHS Region 10` = my_log(`HHS Region 10`),
          `US National` = my_log(`US National`)
        ) ~ AR(1))
      )
  }
  
  # Generate forecasts using the approach that worked for seasonal model
  forecast_samples <- model_fit |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  # Convert to long format for quantile calculation
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
      model_id = paste0("sismid-", model_name),
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
  
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
  
  return(forecast_quantiles)
}

# Test both models on a few sample dates
test_dates <- c("2018-12-01", "2019-01-05", "2019-02-02")

print("=== Testing VAR Models on Sample Dates ===")

all_sample_forecasts <- list()

for (date in test_dates) {
  print(paste("Testing date:", date))
  
  # VAR(2) log model
  var2_fc <- generate_and_evaluate_forecast(as.Date(date), "var2_log")
  all_sample_forecasts[[paste("var2_log", date)]] <- var2_fc
  
  # VAR(1) log model
  var1_fc <- generate_and_evaluate_forecast(as.Date(date), "var1_log")
  all_sample_forecasts[[paste("var1_log", date)]] <- var1_fc
  
  print(paste("VAR(2) log forecasts:", nrow(var2_fc)))
  print(paste("VAR(1) log forecasts:", nrow(var1_fc)))
}

# Combine all forecasts
combined_forecasts <- bind_rows(all_sample_forecasts)

print(paste("Total sample forecasts:", nrow(combined_forecasts)))
print(paste("Models:", paste(unique(combined_forecasts$model_id), collapse = ", ")))
print(paste("Dates:", paste(unique(combined_forecasts$origin_date), collapse = ", ")))

# Save sample forecasts
write.csv(combined_forecasts, "sample_var_log_forecasts.csv", row.names = FALSE)

# Get performance metrics for comparison
print("\n=== Model Performance Comparison ===")
model_metrics <- combined_forecasts |>
  group_by(model_id) |>
  summarise(
    n_forecasts = n(),
    mean_forecast = mean(value, na.rm = TRUE),
    median_forecast = median(value, na.rm = TRUE),
    .groups = "drop"
  )

print(model_metrics)

# AICc comparison from our earlier analysis
print("\n=== AICc Performance (from earlier analysis) ===")
print("VAR(2) log: -6,956 AICc (BEST)")
print("VAR(1) log: -6,924 AICc")
print("Current VAR(1) + season(52): -49,358 AICc")
print("")
print("Note: The seasonal model has much better AICc, but WIS was poor (0.466)")
print("Testing if simpler models without seasonality perform better in WIS")

print("\n=== Next Steps ===")
print("1. Full evaluation of VAR(2) log model on test phase")
print("2. Compare WIS performance vs current model")
print("3. If better, generate full test phase forecasts")