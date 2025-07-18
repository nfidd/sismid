# Create Simple ARIMA Ensemble Model
# Strategy: Multiple ARIMA models with different specs, then average forecasts

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("readr")

set.seed(406)
data(flu_data_hhs)

print("=== Creating Simple ARIMA Ensemble Model ===")
print("Strategy: 5 ARIMA models, average their forecast quantiles")

# Get all forecast dates
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter to valid dates (after sufficient training data)
min_data_date <- as.Date("2016-03-12")  # Use same as VAR models
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("Total dates to process:", length(valid_dates)))

# Simple transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Generate ensemble forecast for a single date
generate_ensemble_forecast <- function(forecast_date) {
  # Prepare training data
  train_data <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week, key = location)
  
  # Check minimum observations
  min_obs <- train_data |>
    as_tibble() |>
    group_by(location) |>
    summarise(n = n()) |>
    pull(n) |>
    min()
  
  if (min_obs < 15) {
    return(NULL)
  }
  
  # Fit 5 different ARIMA models
  model_fit <- train_data |>
    model(
      arima1 = ARIMA(wili ~ pdq(2,1,0)),           # AR(2)
      arima2 = ARIMA(wili ~ pdq(1,1,1)),           # ARIMA(1,1,1)
      arima3 = ARIMA(wili),                        # Auto ARIMA
      arima4 = ARIMA(my_sqrt(wili) ~ pdq(2,1,0)),  # Sqrt transform
      arima5 = ARIMA(wili ~ pdq(0,1,2))            # MA(2)
    )
  
  # Generate forecasts
  fc_list <- list()
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  for (model_name in c("arima1", "arima2", "arima3", "arima4", "arima5")) {
    tryCatch({
      # Extract single model
      single_model <- model_fit |>
        select(location, !!model_name)
      
      # Generate forecast
      fc <- single_model |> forecast(h = 4)
      
      # Extract quantiles
      fc_quantiles <- fc |>
        as_tibble() |>
        mutate(
          target_end_date = as.Date(week),
          origin_date = forecast_date,
          target = "wk inc flu hosp"
        ) |>
        select(origin_date, target, target_end_date, location, wili)
      
      # Get quantiles from distribution
      model_forecasts <- list()
      for (i in 1:nrow(fc_quantiles)) {
        row <- fc_quantiles[i, ]
        dist <- row$wili[[1]]
        
        if (!is.null(dist)) {
          q_values <- quantile(dist, probs = quantile_levels)
          
          for (j in seq_along(quantile_levels)) {
            model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
              origin_date = row$origin_date,
              target = row$target,
              target_end_date = row$target_end_date,
              location = row$location,
              output_type = "quantile",
              output_type_id = quantile_levels[j],
              value = max(q_values[j], 0),
              model = model_name
            )
          }
        }
      }
      
      if (length(model_forecasts) > 0) {
        fc_list[[model_name]] <- do.call(rbind, model_forecasts)
      }
      
    }, error = function(e) {
      # Skip failed models
    })
  }
  
  if (length(fc_list) == 0) {
    return(NULL)
  }
  
  # Combine and average forecasts
  all_forecasts <- do.call(rbind, fc_list)
  
  # Average quantiles across models
  ensemble_forecast <- all_forecasts |>
    group_by(origin_date, target, target_end_date, location, output_type, output_type_id) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    mutate(value = pmax(value, 0))  # Ensure non-negative
  
  return(ensemble_forecast)
}

# Generate forecasts for all dates
print("\n=== Generating Ensemble Forecasts ===")
all_forecasts <- list()
successful <- 0

for (i in seq_along(valid_dates)) {
  date <- valid_dates[i]
  
  if (i %% 10 == 0) {
    print(paste("Progress:", i, "of", length(valid_dates)))
  }
  
  forecast <- generate_ensemble_forecast(as.Date(date))
  if (!is.null(forecast)) {
    all_forecasts[[length(all_forecasts) + 1]] <- forecast
    successful <- successful + 1
  }
}

print(paste("\nSuccessfully generated forecasts for", successful, "dates"))

# Save results
if (length(all_forecasts) > 0) {
  combined_forecasts <- do.call(rbind, all_forecasts)
  
  # Create directories
  model_folder <- file.path(hub_path, "model-output", "sismid-arima5ensemble")
  if (!dir.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Create metadata
  metadata_path <- file.path(hub_path, "model-metadata", "sismid-arima5ensemble.yml")
  metadata_text <- c(
    'team_abbr: "sismid"',
    'model_abbr: "arima5ensemble"',
    'designated_model: false',
    'model_details: "Ensemble of 5 ARIMA models: AR(2), ARIMA(1,1,1), Auto ARIMA, sqrt-AR(2), MA(2). Quantiles averaged across models."'
  )
  writeLines(metadata_text, metadata_path)
  
  # Save forecast files
  print("\nSaving forecast files...")
  forecast_groups <- combined_forecasts |>
    group_by(origin_date) |>
    group_split()
  
  saved <- 0
  for (group in forecast_groups) {
    origin_date <- group$origin_date[1]
    filename <- paste0(origin_date, "-sismid-arima5ensemble.csv")
    write_csv(group, file.path(model_folder, filename))
    saved <- saved + 1
    
    if (saved %% 20 == 0) {
      print(paste("Saved", saved, "files"))
    }
  }
  
  print(paste("\nTotal files saved:", saved))
  print("\n=== ARIMA Ensemble Model Ready ===")
  print("Model: sismid-arima5ensemble")
  print("5 ARIMA models: AR(2), ARIMA(1,1,1), Auto, sqrt-AR(2), MA(2)")
  print("Quantiles averaged across all successful models")
} else {
  print("No forecasts generated!")
}