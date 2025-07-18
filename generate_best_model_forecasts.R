# Generate forecasts for best performing models and evaluate WIS

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("hubData")
library("hubEvals")

# Set seed for reproducibility
set.seed(406)

# Load data
data(flu_data_hhs)

# Hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Get one test phase origin date to evaluate
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
log_transform <- function(x) log(x + 0.1)
inv_log <- function(x) exp(x) - 0.1
my_log <- new_transformation(log_transform, inv_log)

sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

fourth_root_transform <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root_transform, inv_fourth_root)

print("=== Testing Best Models with Forecast Generation ===")

# 1. VAR(2) with log transformation (best overall)
print("1. VAR(2) with log transformation")
var2_log_fit <- flu_data_wide |>
  model(
    var2_log = VAR(vars(
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

# 2. VAR(1) with log transformation
print("2. VAR(1) with log transformation")  
var1_log_fit <- flu_data_wide |>
  model(
    var1_log = VAR(vars(
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

# 3. ARIMA with fourth root transformation
print("3. ARIMA with fourth root transformation")
arima_fourth_fit <- flu_data_long |>
  model(
    arima_fourth = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
  )

# Generate forecasts for each model
print("\n=== Generating Forecasts ===")

# VAR(2) log forecasts
var2_log_fc <- forecast(var2_log_fit, h = 4, bootstrap = TRUE)
print(paste("VAR(2) log forecasts generated:", nrow(var2_log_fc)))

# VAR(1) log forecasts  
var1_log_fc <- forecast(var1_log_fit, h = 4, bootstrap = TRUE)
print(paste("VAR(1) log forecasts generated:", nrow(var1_log_fc)))

# ARIMA forecasts
arima_fourth_fc <- forecast(arima_fourth_fit, h = 4)
print(paste("ARIMA fourth root forecasts generated:", nrow(arima_fourth_fc)))

# Function to convert forecasts to hub format
convert_to_hub_format <- function(forecasts, model_name) {
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  # Location mapping
  loc_df <- data.frame(
    region = c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
               "hhs6", "hhs7", "hhs8", "hhs9", "hhs10"),
    location = c("US National", paste("HHS Region", 1:10))
  )
  
  if ("wili" %in% names(forecasts)) {
    # Univariate ARIMA format
    forecasts_hub <- forecasts |>
      as_tibble() |>
      rename(
        target_end_date = origin_date,
        value = wili
      ) |>
      mutate(
        model_id = model_name,
        target = "ili perc",
        horizon = row_number(),
        .by = c(location, target_end_date)
      ) |>
      mutate(origin_date = target_end_date - horizon * 7) |>
      left_join(loc_df, by = "location") |>
      select(-region)
  } else {
    # VAR format - need to pivot longer for locations
    forecasts_hub <- forecasts |>
      as_tibble() |>
      pivot_longer(
        cols = c(`HHS Region 1`, `HHS Region 2`, `HHS Region 3`, 
                 `HHS Region 4`, `HHS Region 5`, `HHS Region 6`,
                 `HHS Region 7`, `HHS Region 8`, `HHS Region 9`, 
                 `HHS Region 10`, `US National`),
        names_to = "location",
        values_to = "value"
      ) |>
      rename(target_end_date = origin_date) |>
      mutate(
        model_id = model_name,
        target = "ili perc",
        horizon = row_number(),
        .by = c(location, target_end_date)
      ) |>
      mutate(origin_date = target_end_date - horizon * 7)
  }
  
  # Convert to quantile format
  forecasts_quantiles <- forecasts_hub |>
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
  
  return(forecasts_quantiles)
}

# Convert forecasts to hub format
print("\n=== Converting Forecasts to Hub Format ===")

var2_log_hub <- convert_to_hub_format(var2_log_fc, "sismid-var2-log")
var1_log_hub <- convert_to_hub_format(var1_log_fc, "sismid-var1-log")
arima_fourth_hub <- convert_to_hub_format(arima_fourth_fc, "sismid-arima-fourth")

print(paste("VAR(2) log hub format:", nrow(var2_log_hub)))
print(paste("VAR(1) log hub format:", nrow(var1_log_hub)))
print(paste("ARIMA fourth hub format:", nrow(arima_fourth_hub)))

# Combine all forecasts
all_forecasts <- bind_rows(var2_log_hub, var1_log_hub, arima_fourth_hub)

# Get oracle output for evaluation
oracle_output <- connect_target_oracle_output(hub_path) |>
  collect()

# Evaluate each model
print("\n=== Model Evaluation ===")

if (nrow(all_forecasts) > 0) {
  model_scores <- score_model_out(all_forecasts, oracle_output)
  
  print("Model performance comparison:")
  print(model_scores |> 
        select(model_id, wis, interval_coverage_50, interval_coverage_90) |>
        arrange(wis))
  
  # Save results
  write.csv(model_scores, "best_models_wis_comparison.csv", row.names = FALSE)
  
  # Find best model
  best_model <- model_scores |> 
    slice_min(wis, n = 1)
  
  print(paste("\nBest model:", best_model$model_id, "with WIS:", round(best_model$wis, 3)))
  
  # Compare to current performance
  print("\nComparison to current model:")
  print(paste("Current VAR(1) + season(52):", 0.466, "WIS"))
  print(paste("New best model:", round(best_model$wis, 3), "WIS"))
  
  improvement <- (0.466 - best_model$wis) / 0.466 * 100
  print(paste("Improvement:", round(improvement, 1), "%"))
  
} else {
  print("No forecasts generated for evaluation")
}

print("\n=== Next Steps ===")
print("1. If significant improvement found, generate full test phase forecasts")
print("2. Submit best model to hub")
print("3. Compare against all existing models")