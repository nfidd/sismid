# Debug VAR Model

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")

# Load data
data(flu_data_hhs)

# Define transformations
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# Try a single forecast
origin_date <- as.Date("2015-10-24")

# Filter data up to origin date
train_data <- flu_data_hhs |>
  filter(origin_date <= !!origin_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

print("Train data structure:")
print(train_data)
print(paste("Number of rows:", nrow(train_data)))

# Check for missing values
print("\nMissing values per column:")
print(colSums(is.na(train_data)))

# Try fitting VAR model without transformation first
print("\n=== Fitting simple VAR model ===")
tryCatch({
  simple_var <- train_data |>
    model(
      var1 = VAR(vars(
        `HHS Region 1`, `HHS Region 2`, `HHS Region 3`, 
        `HHS Region 4`, `HHS Region 5`, `HHS Region 6`,
        `HHS Region 7`, `HHS Region 8`, `HHS Region 9`, 
        `HHS Region 10`, `US National`
      ) ~ AR(1))
    )
  print("Simple VAR(1) model fitted successfully!")
  
  # Try forecasting
  forecast_result <- forecast(simple_var, h = 4, bootstrap = TRUE)
  print("Forecast generated successfully!")
  
}, error = function(e) {
  print(paste("Error:", e$message))
  print("Full error:")
  print(e)
})

# Try with transformation
print("\n=== Fitting VAR with transformation ===")
tryCatch({
  var_with_transform <- train_data |>
    model(
      var2 = VAR(vars(
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
  print("VAR(2) with transformation fitted successfully!")
  
  # Try generating samples
  generated_samples <- generate(var_with_transform, h = 4, times = 10, bootstrap = TRUE)
  print("Samples generated successfully!")
  print(generated_samples)
  
}, error = function(e) {
  print(paste("Error:", e$message))
  print("Full error:")
  print(e)
})

# Alternative approach: Try VARIMA
print("\n=== Trying VARIMA approach ===")
tryCatch({
  # Prepare data in different format
  train_data_long <- flu_data_hhs |>
    filter(origin_date <= !!origin_date)
  
  # Try individual ARIMA models for comparison
  arima_models <- train_data_long |>
    model(
      arima = ARIMA(my_fourth_root(wili) ~ pdq(2,1,0))
    )
  
  print("Individual ARIMA models fitted successfully!")
  
  # Generate forecasts
  arima_forecasts <- generate(arima_models, h = 4, times = 10, bootstrap = TRUE)
  print("ARIMA forecasts generated successfully!")
  
}, error = function(e) {
  print(paste("Error:", e$message))
})