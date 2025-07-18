# Test ensemble and fundamentally different model architectures

library("nfidd")
library("dplyr")
library("tidyr")
library("forecast")
library("randomForest")
library("hubEvals")

# Load data
data(flu_data_hhs)

print("=== Testing Ensemble and Advanced Models ===")

# Validation setup
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")

# Helper function to create lagged features
create_features <- function(data, lags = 8) {
  features <- data
  
  # Add lagged values
  for (i in 1:lags) {
    features[[paste0("lag_", i)]] <- lag(data$wili, i)
  }
  
  # Add seasonal features
  features$week_of_year <- lubridate::week(data$origin_date)
  features$month <- lubridate::month(data$origin_date)
  features$cos_week <- cos(2 * pi * features$week_of_year / 52)
  features$sin_week <- sin(2 * pi * features$week_of_year / 52)
  
  # Add trend
  features$time_trend <- as.numeric(data$origin_date - min(data$origin_date))
  
  return(features)
}

# Test different model architectures
models_to_test <- list(
  "ensemble_arima_rf" = "Ensemble: ARIMA + Random Forest",
  "sarima_auto" = "Auto SARIMA per location",
  "random_forest" = "Random Forest with engineered features",
  "hierarchical_arima" = "Hierarchical ARIMA with national trend"
)

results <- list()

for (model_name in names(models_to_test)) {
  print(paste("\n--- Testing", models_to_test[[model_name]], "---"))
  
  model_forecasts <- list()
  
  for (i in seq_along(validation_dates[1:10])) {  # Test on subset first
    origin_date <- validation_dates[i]
    
    # Get training data
    train_data <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      select(origin_date, location, wili) |>
      arrange(origin_date)
    
    # Loop through locations
    for (loc in unique(train_data$location)) {
      loc_data <- train_data |> filter(location == loc)
      
      tryCatch({
        
        if (model_name == "ensemble_arima_rf") {
          # Ensemble: ARIMA + Random Forest
          
          # ARIMA component
          ts_data <- ts(loc_data$wili, frequency = 52)
          arima_model <- auto.arima(ts_data, seasonal = TRUE)
          arima_forecast <- forecast(arima_model, h = 4)
          
          # Random Forest component
          rf_data <- create_features(loc_data) |>
            na.omit() |>
            select(-origin_date, -location)
          
          if (nrow(rf_data) > 20) {
            rf_model <- randomForest(wili ~ ., data = rf_data, ntree = 100)
            
            # Create future features for RF prediction
            last_obs <- tail(rf_data, 1)
            rf_preds <- numeric(4)
            for (h in 1:4) {
              # Use last observation as base for features
              future_features <- last_obs
              future_features$week_of_year <- (future_features$week_of_year + h - 1) %% 52 + 1
              future_features$month <- ((future_features$month + h - 1) %% 12) + 1
              future_features$cos_week <- cos(2 * pi * future_features$week_of_year / 52)
              future_features$sin_week <- sin(2 * pi * future_features$week_of_year / 52)
              future_features$time_trend <- future_features$time_trend + h * 7
              
              rf_preds[h] <- predict(rf_model, future_features)
            }
            
            # Ensemble: average ARIMA and RF
            ensemble_preds <- (arima_forecast$mean + rf_preds) / 2
          } else {
            # Fall back to ARIMA only
            ensemble_preds <- arima_forecast$mean
          }
          
          # Use ARIMA prediction intervals as base
          pred_intervals <- cbind(
            arima_forecast$lower[, "95%"],
            arima_forecast$lower[, "80%"],
            ensemble_preds,
            arima_forecast$upper[, "80%"],
            arima_forecast$upper[, "95%"]
          )
          
        } else if (model_name == "sarima_auto") {
          # Auto SARIMA
          ts_data <- ts(loc_data$wili, frequency = 52)
          sarima_model <- auto.arima(ts_data, seasonal = TRUE, 
                                     stepwise = FALSE, approximation = FALSE)
          sarima_forecast <- forecast(sarima_model, h = 4)
          
          pred_intervals <- cbind(
            sarima_forecast$lower[, "95%"],
            sarima_forecast$lower[, "80%"],
            sarima_forecast$mean,
            sarima_forecast$upper[, "80%"],
            sarima_forecast$upper[, "95%"]
          )
          
        } else if (model_name == "random_forest") {
          # Random Forest only
          rf_data <- create_features(loc_data) |>
            na.omit() |>
            select(-origin_date, -location)
          
          if (nrow(rf_data) > 20) {
            rf_model <- randomForest(wili ~ ., data = rf_data, ntree = 200)
            
            rf_preds <- numeric(4)
            for (h in 1:4) {
              last_obs <- tail(rf_data, 1)
              last_obs$week_of_year <- (last_obs$week_of_year + h - 1) %% 52 + 1
              last_obs$month <- ((last_obs$month + h - 1) %% 12) + 1
              last_obs$cos_week <- cos(2 * pi * last_obs$week_of_year / 52)
              last_obs$sin_week <- sin(2 * pi * last_obs$week_of_year / 52)
              last_obs$time_trend <- last_obs$time_trend + h * 7
              
              rf_preds[h] <- predict(rf_model, last_obs)
            }
            
            # Approximate prediction intervals using recent volatility
            recent_residuals <- tail(rf_data$wili - predict(rf_model, rf_data), 20)
            pred_sd <- sd(recent_residuals, na.rm = TRUE)
            
            pred_intervals <- cbind(
              rf_preds - 1.96 * pred_sd,
              rf_preds - 1.28 * pred_sd,
              rf_preds,
              rf_preds + 1.28 * pred_sd,
              rf_preds + 1.96 * pred_sd
            )
          } else {
            next  # Skip if not enough data
          }
          
        } else if (model_name == "hierarchical_arima") {
          # Hierarchical: use national trend + location-specific ARIMA
          
          # Get national data for trend
          national_data <- train_data |>
            filter(location == "US National") |>
            pull(wili)
          
          if (loc == "US National") {
            # For national, use regular ARIMA
            ts_data <- ts(loc_data$wili, frequency = 52)
            arima_model <- auto.arima(ts_data, seasonal = TRUE)
            arima_forecast <- forecast(arima_model, h = 4)
            
            pred_intervals <- cbind(
              arima_forecast$lower[, "95%"],
              arima_forecast$lower[, "80%"],
              arima_forecast$mean,
              arima_forecast$upper[, "80%"],
              arima_forecast$upper[, "95%"]
            )
          } else {
            # For regions, model as national + local deviation
            if (length(national_data) == nrow(loc_data)) {
              local_deviation <- loc_data$wili - national_data
              
              # Model deviations
              ts_dev <- ts(local_deviation, frequency = 52)
              dev_model <- auto.arima(ts_dev, seasonal = TRUE)
              dev_forecast <- forecast(dev_model, h = 4)
              
              # Model national trend
              ts_national <- ts(national_data, frequency = 52)
              nat_model <- auto.arima(ts_national, seasonal = TRUE)
              nat_forecast <- forecast(nat_model, h = 4)
              
              # Combine
              combined_mean <- nat_forecast$mean + dev_forecast$mean
              combined_sd <- sqrt(nat_forecast$mean^2 + dev_forecast$mean^2) * 0.1
              
              pred_intervals <- cbind(
                combined_mean - 1.96 * combined_sd,
                combined_mean - 1.28 * combined_sd,
                combined_mean,
                combined_mean + 1.28 * combined_sd,
                combined_mean + 1.96 * combined_sd
              )
            } else {
              next  # Skip if data doesn't align
            }
          }
        }
        
        # Ensure non-negative predictions
        pred_intervals <- pmax(pred_intervals, 0)
        
        # Convert to quantile format
        quantiles <- c(0.025, 0.1, 0.5, 0.9, 0.975)
        
        for (h in 1:4) {
          target_date <- origin_date + h * 7
          
          for (q in seq_along(quantiles)) {
            model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
              origin_date = origin_date,
              target = "wk inc flu hosp",
              target_end_date = target_date,
              location = loc,
              output_type = "quantile",
              output_type_id = quantiles[q],
              value = pred_intervals[h, q]
            )
          }
        }
        
      }, error = function(e) {
        print(paste("Error for", model_name, "location", loc, ":", e$message))
      })
    }
    
    if (i %% 2 == 0) {
      print(paste("Completed", i, "of 10 validation dates"))
    }
  }
  
  # Evaluate if we have forecasts
  if (length(model_forecasts) > 0) {
    forecasts_df <- do.call(rbind, model_forecasts)
    
    # Calculate WIS
    hub_path <- here::here("sismid-ili-forecasting-sandbox")
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(validation_dates[1:10]) + 7,
             target_end_date <= max(validation_dates[1:10]) + 28) |>
      collect()
    
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    results[[model_name]] <- list(
      wis = mean_wis,
      forecasts = forecasts_df,
      description = models_to_test[[model_name]]
    )
    
    print(paste("Mean WIS:", round(mean_wis, 4)))
  }
}

# Compare results
print("\n=== Advanced Model Comparison ===")
print("Current VAR(2) sqrt baseline: 0.208 WIS")
print("")

for (model_name in names(results)) {
  wis <- results[[model_name]]$wis
  improvement <- ((0.208 - wis) / 0.208) * 100
  
  print(paste(
    results[[model_name]]$description, 
    "- WIS:", round(wis, 4),
    "- Change:", ifelse(improvement > 0, "+", ""), round(improvement, 1), "%"
  ))
}

# Identify best model
if (length(results) > 0) {
  wis_scores <- sapply(results, function(x) x$wis)
  best_model <- names(wis_scores)[which.min(wis_scores)]
  
  if (wis_scores[best_model] < 0.208) {
    print(paste("\nBEST MODEL:", results[[best_model]]$description))
    print("This model beats our current VAR(2) sqrt!")
    print("Ready to create full forecasts for second submission")
  } else {
    print("\nNo advanced model beats our VAR(2) sqrt baseline")
    print("Consider other approaches or stick with current submission")
  }
}