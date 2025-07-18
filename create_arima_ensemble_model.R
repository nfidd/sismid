# Create ARIMA Ensemble Model - 10 different ARIMA structures
# Each model captures different aspects of the data

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("readr")

set.seed(406)
data(flu_data_hhs)

print("=== Creating ARIMA Ensemble Model ===")
print("Strategy: 10 different ARIMA models with various structures")
print("Each model votes equally in the ensemble")

# Get all dates
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter to valid dates
min_data_date <- as.Date("2015-10-24") + 20 * 7
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("Total dates to process:", length(valid_dates)))

# Define 10 different ARIMA model specifications
model_specs <- list(
  arima1 = ~ ARIMA(wili ~ pdq(2,1,0)),                    # Simple AR
  arima2 = ~ ARIMA(wili ~ pdq(1,1,1)),                    # Classic ARIMA
  arima3 = ~ ARIMA(wili ~ pdq(2,1,2)),                    # Higher order
  arima4 = ~ ARIMA(wili ~ pdq(3,1,0)),                    # AR(3)
  arima5 = ~ ARIMA(wili ~ pdq(0,1,2)),                    # Pure MA
  arima6 = ~ ARIMA(wili ~ pdq(1,0,1) + PDQ(1,0,0)),      # Seasonal
  arima7 = ~ ARIMA(wili ~ pdq(2,1,1) + PDQ(1,0,0)),      # Seasonal ARIMA
  arima8 = ~ ARIMA(wili),                                 # Auto ARIMA
  arima9 = ~ ARIMA(sqrt(wili) ~ pdq(2,1,0)),             # Sqrt transform
  arima10 = ~ ARIMA(log(wili + 0.1) ~ pdq(1,1,1))        # Log transform
)

# Generate forecasts function
generate_arima_ensemble_forecast <- function(forecast_date) {
  # Prepare data
  train_data <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    as_tsibble(index = origin_date, key = location)
  
  # Check for sufficient data
  min_obs <- train_data |>
    as_tibble() |>
    group_by(location) |>
    summarise(n = n()) |>
    pull(n) |>
    min()
  
  if (min_obs < 20) {
    return(NULL)
  }
  
  # Fit all 10 models
  all_forecasts <- list()
  
  for (i in 1:10) {
    model_name <- names(model_specs)[i]
    model_spec <- model_specs[[i]]
    
    tryCatch({
      # Fit model
      model_fit <- train_data |>
        model(!!model_name := !!model_spec)
      
      # Generate 10 samples per model (total 100 across ensemble)
      model_samples <- model_fit |>
        generate(h = 4, times = 10, bootstrap = TRUE) |>
        mutate(
          # Back-transform if needed
          .sim = case_when(
            model_name == "arima9" ~ .sim^2,
            model_name == "arima10" ~ exp(.sim) - 0.1,
            TRUE ~ .sim
          ),
          .sim = pmax(.sim, 0),  # Ensure non-negative
          model = model_name
        ) |>
        as_tibble()
      
      all_forecasts[[model_name]] <- model_samples
      
    }, error = function(e) {
      # Skip models that fail
    })
  }
  
  if (length(all_forecasts) == 0) {
    return(NULL)
  }
  
  # Combine all model forecasts
  combined_samples <- bind_rows(all_forecasts) |>
    mutate(
      value = .sim,
      target = "wk inc flu hosp",
      target_end_date = as.Date(origin_date),
      origin_date = forecast_date
    )
  
  # Calculate quantiles from combined ensemble
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- combined_samples |>
    group_by(location, target_end_date) |>
    summarise(
      quantile_values = list(quantile(value, quantile_levels, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    unnest_longer(quantile_values, indices_to = "quantile_idx") |>
    mutate(
      output_type_id = quantile_levels[quantile_idx],
      value = pmax(quantile_values, 0),
      origin_date = forecast_date,
      target = "wk inc flu hosp",
      output_type = "quantile"
    ) |>
    select(origin_date, target, target_end_date, location, output_type, output_type_id, value)
  
  return(forecast_quantiles)
}

# Generate all forecasts
print("\n=== Generating Ensemble Forecasts ===")
all_forecasts <- list()
successful <- 0

for (i in seq_along(valid_dates)) {
  date <- valid_dates[i]
  
  if (i %% 10 == 0) {
    print(paste("Progress:", i, "of", length(valid_dates)))
  }
  
  forecast <- generate_arima_ensemble_forecast(as.Date(date))
  if (!is.null(forecast)) {
    all_forecasts[[length(all_forecasts) + 1]] <- forecast
    successful <- successful + 1
  }
}

print(paste("\nSuccessfully generated forecasts for", successful, "dates"))

if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  # Create output directory
  model_folder <- file.path(hub_path, "model-output", "sismid-arima10claude")
  if (!dir.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Create metadata
  metadata_path <- file.path(hub_path, "model-metadata", "sismid-arima10claude.yml")
  metadata_text <- c(
    'team_abbr: "sismid"',
    'model_abbr: "arima10claude"',
    'designated_model: false',
    'model_details: "Ensemble of 10 ARIMA models with different structures developed with Claude. Includes AR, MA, ARIMA, seasonal variants, and transformations."'
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
    filename <- paste0(origin_date, "-sismid-arima10claude.csv")
    write_csv(group, file.path(model_folder, filename))
    saved <- saved + 1
    
    if (saved %% 20 == 0) {
      print(paste("Saved", saved, "files"))
    }
  }
  
  print(paste("\nTotal files saved:", saved))
  print("\n=== Model Ready for Submission ===")
  print("Model: sismid-arima10claude")
  print("10 ARIMA models in ensemble:")
  print("- ARIMA(2,1,0), ARIMA(1,1,1), ARIMA(2,1,2)")
  print("- ARIMA(3,1,0), ARIMA(0,1,2)")
  print("- Seasonal: ARIMA(1,0,1)(1,0,0)[52], ARIMA(2,1,1)(1,0,0)[52]")
  print("- Auto ARIMA, sqrt-ARIMA(2,1,0), log-ARIMA(1,1,1)")
}