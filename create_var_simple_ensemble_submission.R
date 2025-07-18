# Create VAR + Simple Models Ensemble Submission
# Combines VAR(2) sqrt with individual ARIMA/ETS models

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("lubridate")
library("vars")

set.seed(406)

# Load data
data(flu_data_hhs)

# Define transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Test dates for submission
test_dates <- seq(as.Date("2017-11-04"), as.Date("2019-05-18"), by = "week")

print("=== Creating VAR + Simple Models Ensemble ===")
print(paste("Number of forecast dates:", length(test_dates)))

# Store all forecasts
all_forecasts <- list()

# Process each origin date
for (i in seq_along(test_dates)) {
  origin_date <- test_dates[i]
  
  tryCatch({
    # Prepare data for individual models (long format)
    train_long <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      mutate(week = yearweek(origin_date)) |>
      as_tsibble(index = week, key = location)
    
    # Prepare data for VAR (wide format)
    train_wide <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      dplyr::select(origin_date, location, wili) |>
      pivot_wider(names_from = location, values_from = wili) |>
      arrange(origin_date)
    
    # Apply transformation for VAR
    numeric_cols <- setdiff(names(train_wide), "origin_date")
    train_wide_transformed <- train_wide
    train_wide_transformed[numeric_cols] <- lapply(train_wide[numeric_cols], sqrt)
    
    # Fit VAR(2) using vars package
    ts_data <- ts(train_wide_transformed[numeric_cols], frequency = 52)
    var_model <- VAR(ts_data, p = 2, type = "const")
    var_forecast <- predict(var_model, n.ahead = 4)
    
    # Fit individual models using fable
    individual_models <- train_long |>
      model(
        arima = ARIMA(my_sqrt(wili)),
        ets = ETS(my_sqrt(wili))
      )
    
    # Generate individual forecasts
    individual_fc <- individual_models |> forecast(h = 4)
    
    # Process forecasts for each location
    locations <- unique(train_long$location)
    quantile_levels <- seq(0.01, 0.99, by = 0.01)
    
    for (loc in locations) {
      # Get VAR forecasts for this location
      var_pred <- var_forecast$fcst[[loc]][, "fcst"]
      var_se <- var_forecast$fcst[[loc]][, "SE"]
      
      # Get individual model forecasts
      arima_fc <- individual_fc |> 
        filter(location == loc, .model == "arima")
      ets_fc <- individual_fc |> 
        filter(location == loc, .model == "ets")
      
      # Process each horizon
      for (h in 1:4) {
        target_date <- origin_date + h * 7
        
        # VAR quantiles (using normal approximation)
        var_mean <- inv_sqrt(var_pred[h])
        var_sd <- var_se[h] * 2  # Adjust SD for back-transformation
        var_quantiles <- qnorm(quantile_levels, mean = var_mean, sd = var_sd)
        
        # Individual model quantiles
        if (nrow(arima_fc) >= h) {
          arima_dist <- arima_fc$wili[[h]]
          arima_quantiles <- inv_sqrt(quantile(arima_dist, probs = quantile_levels))
        } else {
          arima_quantiles <- var_quantiles  # Fallback
        }
        
        if (nrow(ets_fc) >= h) {
          ets_dist <- ets_fc$wili[[h]]
          ets_quantiles <- inv_sqrt(quantile(ets_dist, probs = quantile_levels))
        } else {
          ets_quantiles <- var_quantiles  # Fallback
        }
        
        # Ensemble: weighted average (VAR 50%, ARIMA 30%, ETS 20%)
        ensemble_quantiles <- 0.5 * var_quantiles + 
                             0.3 * arima_quantiles + 
                             0.2 * ets_quantiles
        
        # Ensure non-negative
        ensemble_quantiles <- pmax(ensemble_quantiles, 0)
        
        # Store forecasts in hub format
        for (k in seq_along(quantile_levels)) {
          all_forecasts[[length(all_forecasts) + 1]] <- data.frame(
            model_id = "sismid-ensemble",
            location = loc,
            origin_date = origin_date,
            horizon = h,
            target_end_date = target_date,
            value = ensemble_quantiles[k],
            output_type_id = quantile_levels[k],
            target = "ili perc",
            output_type = "quantile",
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (i %% 10 == 0) {
      print(paste("Completed", i, "of", length(test_dates), "forecasts"))
    }
    
  }, error = function(e) {
    # Fallback to simple ARIMA if ensemble fails
    tryCatch({
      train_long <- flu_data_hhs |>
        filter(origin_date <= !!origin_date) |>
        mutate(week = yearweek(origin_date)) |>
        as_tsibble(index = week, key = location)
      
      simple_model <- train_long |>
        model(arima = ARIMA(wili ~ pdq(2,1,1)))
      
      simple_fc <- simple_model |> forecast(h = 4)
      
      for (loc in unique(train_long$location)) {
        loc_fc <- simple_fc |> filter(location == loc)
        
        for (h in 1:4) {
          target_date <- origin_date + h * 7
          
          if (nrow(loc_fc) >= h) {
            dist <- loc_fc$wili[[h]]
            q_values <- quantile(dist, probs = quantile_levels)
            q_values <- pmax(q_values, 0)
            
            for (k in seq_along(quantile_levels)) {
              all_forecasts[[length(all_forecasts) + 1]] <- data.frame(
                model_id = "sismid-ensemble",
                location = loc,
                origin_date = origin_date,
                horizon = h,
                target_end_date = target_date,
                value = q_values[k],
                output_type_id = quantile_levels[k],
                target = "ili perc",
                output_type = "quantile",
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }, error = function(e2) {
      print(paste("Fallback also failed at", origin_date))
    })
  })
}

# Combine and save
if (length(all_forecasts) > 0) {
  submission_df <- do.call(rbind, all_forecasts)
  
  print("\nSubmission summary:")
  print(paste("Total rows:", nrow(submission_df)))
  print(paste("Unique origin dates:", n_distinct(submission_df$origin_date)))
  print(paste("Unique locations:", n_distinct(submission_df$location)))
  
  # Show sample
  print("\nSample forecasts:")
  sample_df <- submission_df |>
    filter(origin_date == submission_df$origin_date[1],
           location == "US National",
           output_type_id %in% c(0.1, 0.5, 0.9)) |>
    arrange(horizon, output_type_id)
  print(sample_df)
  
  # Save submission
  filename <- "ensemble_submission.csv"
  write.csv(submission_df, filename, row.names = FALSE)
  
  print(paste("\nSubmission saved to:", filename))
  print("\nEnsemble composition:")
  print("- 50% VAR(2) sqrt (captures cross-regional patterns)")
  print("- 30% ARIMA (automatic selection per location)")
  print("- 20% ETS (exponential smoothing)")
  print("\nThis ensemble provides robustness through model diversity!")
  
} else {
  print("No forecasts generated!")
}

print("\n=== Ensemble Submission Complete ===")