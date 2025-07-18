# Debug VAR2sqrt forecast generation
library("dplyr")
library("readr")
library("lubridate")
library("tsibble")
library("fable")
library("tidyr")
library("here")

hub_path <- here::here("sismid-ili-forecasting-sandbox")
flu_data <- read_csv(file.path(hub_path, "target-data", "time-series.csv"))

# Test with a specific date
test_date <- as.Date("2016-11-05")

# Filter data up to origin date
train_data <- flu_data |>
  filter(target_end_date < test_date) |>
  mutate(
    week_ending = as.Date(target_end_date),
    ili = sqrt(observation)
  ) |>
  as_tsibble(index = week_ending, key = location) |>
  filter(!is.na(ili))

print("Train data summary:")
summary_data <- train_data |> 
  as_tibble() |>
  group_by(location) |> 
  summarise(n = n(), first = min(week_ending), last = max(week_ending))
print(summary_data)

# Check data structure
print("\nData structure:")
print(head(train_data))

# Try fitting model and generating forecasts
tryCatch({
  model_fit <- train_data |>
    model(var_model = VAR(ili ~ AR(2)))
  print("Model fitted successfully!")
  print(model_fit)
  
  # Generate forecasts
  forecasts <- model_fit |>
    forecast(h = 4) |>
    mutate(
      # Transform back from sqrt
      .mean = .mean^2,
      .distribution = distributional::dist_transformed(
        .distribution, 
        transform = list(forward = function(x) x^2, inverse = sqrt),
        inverse_dist = .distribution
      )
    )
  print("\nForecasts generated successfully!")
  print(forecasts)
  
  # Extract quantiles
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  
  forecast_df <- forecasts |>
    hilo(quantiles * 100) |>
    unpack_hilo(names(quantiles)) |>
    select(location, week_ending, .mean, all_of(names(quantiles)))
  
  print("\nForecast dataframe created!")
  print(head(forecast_df))
  
}, error = function(e) {
  print(paste("Error:", e$message))
  print(traceback())
})