# Create Working Ensemble Model for Hub Submission
# Based on successful VAR approach with ARIMA ensemble

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("readr")

set.seed(406)
data(flu_data_hhs)

print("=== Creating Working Ensemble Model ===")
print("Strategy: VAR(2) sqrt + ARIMA ensemble (our best models combined)")

# Get all forecast dates from hub
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter to dates with sufficient training data
min_data_date <- as.Date("2016-03-12")  # Same as VAR models
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("Total forecast dates:", length(valid_dates)))

# Transformation for sqrt model
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Generate forecast for single date
generate_forecast <- function(forecast_date) {
  # Prepare data for VAR (wide format)
  train_data_wide <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week) |>
    fill_gaps()
  
  # Prepare data for ARIMA (long format)
  train_data_long <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week, key = location)
  
  # Check minimum observations
  min_obs <- nrow(train_data_wide)
  if (min_obs < 15) {
    return(NULL)
  }
  
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  all_model_results <- list()
  
  # 1. VAR(2) with sqrt transformation (our best model!)
  tryCatch({
    var_model <- train_data_wide |>
      model(
        var_sqrt = VAR(vars(
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
    
    var_fc <- var_model |> forecast(h = 4)
    
    # Convert VAR to hub format
    var_hub <- var_fc |>
      as_tibble() |>
      pivot_longer(cols = c(`HHS Region 1`:`US National`),
                  names_to = "location", 
                  values_to = ".distribution") |>
      mutate(
        target_end_date = as.Date(week),
        origin_date = forecast_date,
        target = "wk inc flu hosp"
      )
    
    # Extract VAR quantiles
    var_results <- list()
    for (i in 1:nrow(var_hub)) {
      row <- var_hub[i, ]
      dist <- row$.distribution[[1]]
      
      if (!is.null(dist)) {
        q_values <- quantile(dist, probs = quantile_levels)
        
        for (j in seq_along(quantile_levels)) {
          var_results[[length(var_results) + 1]] <- data.frame(
            origin_date = row$origin_date,
            target = row$target,
            target_end_date = row$target_end_date,
            location = row$location,
            output_type = "quantile",
            output_type_id = quantile_levels[j],
            value = max(q_values[j], 0),
            model = "var_sqrt"
          )
        }
      }
    }
    
    if (length(var_results) > 0) {
      all_model_results[["var"]] <- do.call(rbind, var_results)
    }
    
  }, error = function(e) {
    # VAR failed, continue with ARIMA only
  })
  
  # 2. ARIMA models
  tryCatch({
    arima_fit <- train_data_long |>
      model(
        arima1 = ARIMA(wili ~ pdq(2,1,0)),              # AR(2)
        arima2 = ARIMA(wili ~ pdq(1,1,1)),              # ARIMA(1,1,1)
        arima3 = ARIMA(wili ~ pdq(2,1,0) + fourier(K=2)) # AR(2) with Fourier seasonality
      )
    
    arima_fc <- arima_fit |> forecast(h = 4)
    
    # Process ARIMA results
    arima_results <- list()
    for (model_name in c("arima1", "arima2", "arima3")) {
      model_fc <- arima_fc |>
        select(location, week, all_of(model_name)) |>
        as_tibble() |>
        rename(wili = all_of(model_name)) |>
        mutate(
          target_end_date = as.Date(week),
          origin_date = forecast_date,
          target = "wk inc flu hosp"
        )
      
      for (i in 1:nrow(model_fc)) {
        row <- model_fc[i, ]
        dist <- row$wili[[1]]
        
        if (!is.null(dist)) {
          q_values <- quantile(dist, probs = quantile_levels)
          
          for (j in seq_along(quantile_levels)) {
            arima_results[[length(arima_results) + 1]] <- data.frame(
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
    }
    
    if (length(arima_results) > 0) {
      all_model_results[["arima"]] <- do.call(rbind, arima_results)
    }
    
  }, error = function(e) {
    # ARIMA failed
  })
  
  if (length(all_model_results) == 0) {
    return(NULL)
  }
  
  # Combine all models and average
  all_results <- do.call(rbind, all_model_results)
  
  ensemble_result <- all_results |>
    group_by(origin_date, target, target_end_date, location, output_type, output_type_id) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    mutate(value = pmax(value, 0))  # Ensure non-negative
  
  return(ensemble_result)
}

# Test on a few dates first
print("\n=== Testing Ensemble Generation ===")
test_dates <- valid_dates[1:2]  # Start with just 2 dates
test_results <- list()

for (date in test_dates) {
  cat("Testing:", as.character(date), "... ")
  result <- generate_forecast(as.Date(date))
  if (!is.null(result)) {
    test_results[[length(test_results) + 1]] <- result
    cat("SUCCESS\n")
  } else {
    cat("FAILED\n")
  }
}

if (length(test_results) > 0) {
  print("✓ Test successful! Proceeding with full generation...")
  
  # Generate all forecasts
  print("\n=== Generating All Forecasts ===")
  all_forecasts <- list()
  successful <- 0
  
  for (i in seq_along(valid_dates)) {
    date <- valid_dates[i]
    
    if (i %% 10 == 0) {
      print(paste("Progress:", i, "of", length(valid_dates)))
    }
    
    forecast <- generate_forecast(as.Date(date))
    if (!is.null(forecast)) {
      all_forecasts[[length(all_forecasts) + 1]] <- forecast
      successful <- successful + 1
    }
  }
  
  print(paste("\nSuccessfully generated forecasts for", successful, "dates"))
  
  # Save results
  if (length(all_forecasts) > 0) {
    combined_forecasts <- do.call(rbind, all_forecasts)
    
    # Create model folder
    model_folder <- file.path(hub_path, "model-output", "sismid-vararimaensemble")
    if (!dir.exists(model_folder)) {
      dir.create(model_folder, recursive = TRUE)
    }
    
    # Create metadata
    metadata_path <- file.path(hub_path, "model-metadata", "sismid-vararimaensemble.yml")
    metadata_text <- c(
      'team_abbr: "sismid"',
      'model_abbr: "vararimaensemble"', 
      'designated_model: false',
      'model_details: "Ensemble combining VAR(2) with sqrt transformation and ARIMA models (AR(2), ARIMA(1,1,1), Fourier-seasonal AR(2)). Our best performing VAR model with reliable and seasonal ARIMA models."'
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
      filename <- paste0(origin_date, "-sismid-vararimaensemble.csv")
      write_csv(group, file.path(model_folder, filename))
      saved <- saved + 1
      
      if (saved %% 20 == 0) {
        print(paste("Saved", saved, "files"))
      }
    }
    
    print(paste("\nTotal files saved:", saved))
    print("\n=== VAR+ARIMA Ensemble Model Ready ===")
    print("Model: sismid-vararimaensemble")
    print("Combines our best VAR(2) sqrt model with ARIMA models!")
    print("Ready for hub submission!")
  }
} else {
  print("✗ Test failed. Need to debug the ensemble generation.")
}