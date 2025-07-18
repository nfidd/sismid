# Test Multiple Model Configurations for Better WIS Performance

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("ggplot2")

# Set seed for reproducibility
set.seed(406)

# Load data
data(flu_data_hhs)

# Test on a representative forecast date
test_date <- as.Date("2018-12-01")

# Prepare data for different model types
flu_data_wide <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

flu_data_long <- flu_data_hhs |>
  filter(origin_date <= test_date)

print("=== Testing Multiple Model Configurations ===")
print(paste("Training data up to:", test_date))
print(paste("Rows in wide data:", nrow(flu_data_wide)))

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

# Test VAR models with different configurations
print("\n=== Testing VAR Model Configurations ===")

# 1. VAR(1) with sqrt transformation
print("1. VAR(1) with sqrt transformation")
var1_sqrt <- flu_data_wide |>
  model(
    var1_sqrt = VAR(vars(
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
    ) ~ AR(1))
  )
var1_sqrt_metrics <- glance(var1_sqrt)
print(var1_sqrt_metrics)

# 2. VAR(2) with sqrt transformation
print("\n2. VAR(2) with sqrt transformation")
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
var2_sqrt_metrics <- glance(var2_sqrt)
print(var2_sqrt_metrics)

# 3. VAR(1) with log transformation
print("\n3. VAR(1) with log transformation")
var1_log <- flu_data_wide |>
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
var1_log_metrics <- glance(var1_log)
print(var1_log_metrics)

# 4. VAR(2) with log transformation
print("\n4. VAR(2) with log transformation")
var2_log <- flu_data_wide |>
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
var2_log_metrics <- glance(var2_log)
print(var2_log_metrics)

# Test univariate ARIMA models
print("\n=== Testing Univariate ARIMA Models ===")

# 5. ARIMA with fourth root transformation
print("5. ARIMA with fourth root transformation")
arima_fourth <- flu_data_long |>
  model(
    arima_fourth = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
  )
arima_fourth_metrics <- arima_fourth |>
  glance() |>
  group_by(.model) |>
  summarise(
    mean_aicc = mean(AICc),
    mean_aic = mean(AIC),
    mean_bic = mean(BIC)
  )
print(arima_fourth_metrics)

# 6. ARIMA with log transformation
print("\n6. ARIMA with log transformation")
arima_log <- flu_data_long |>
  model(
    arima_log = ARIMA(my_log(wili) ~ pdq(2,1,0))
  )
arima_log_metrics <- arima_log |>
  glance() |>
  group_by(.model) |>
  summarise(
    mean_aicc = mean(AICc),
    mean_aic = mean(AIC),
    mean_bic = mean(BIC)
  )
print(arima_log_metrics)

# 7. ARIMA with automatic selection
print("\n7. ARIMA with automatic selection")
arima_auto <- flu_data_long |>
  model(
    arima_auto = ARIMA(my_log(wili))
  )
arima_auto_metrics <- arima_auto |>
  glance() |>
  group_by(.model) |>
  summarise(
    mean_aicc = mean(AICc),
    mean_aic = mean(AIC),
    mean_bic = mean(BIC)
  )
print(arima_auto_metrics)

# Compile results
print("\n=== MODEL COMPARISON SUMMARY ===")
model_comparison <- data.frame(
  Model = c("VAR(1) sqrt", "VAR(2) sqrt", "VAR(1) log", "VAR(2) log", 
            "ARIMA fourth", "ARIMA log", "ARIMA auto"),
  AICc = c(var1_sqrt_metrics$AICc, var2_sqrt_metrics$AICc, var1_log_metrics$AICc,
           var2_log_metrics$AICc, 
           arima_fourth_metrics$mean_aicc, arima_log_metrics$mean_aicc,
           arima_auto_metrics$mean_aicc)
) |>
  arrange(AICc)

print(model_comparison)

# Save results
write.csv(model_comparison, "model_configuration_comparison.csv", row.names = FALSE)

# Identify best models
best_var <- model_comparison |> filter(grepl("VAR", Model)) |> slice(1)
best_arima <- model_comparison |> filter(grepl("ARIMA", Model)) |> slice(1)

print("\n=== BEST MODELS BY TYPE ===")
print(paste("Best VAR model:", best_var$Model, "with AICc:", round(best_var$AICc, 2)))
print(paste("Best ARIMA model:", best_arima$Model, "with AICc:", round(best_arima$AICc, 2)))
print(paste("Overall best model:", model_comparison$Model[1], "with AICc:", round(model_comparison$AICc[1], 2)))

print("\n=== RECOMMENDATIONS ===")
print("1. Test the top 2-3 models with actual forecast generation")
print("2. Evaluate WIS performance on test data")
print("3. Consider ensemble approaches combining best models")