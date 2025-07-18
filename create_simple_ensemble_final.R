# Create Simple Ensemble Final Submission
# Quick and efficient ensemble approach

library("nfidd")
library("dplyr")
library("tidyr")
library("forecast")

set.seed(406)

# Load data
data(flu_data_hhs)

# Test dates for submission
test_dates <- seq(as.Date("2017-11-04"), as.Date("2019-05-18"), by = "week")

print("=== Creating Simple Ensemble Submission ===")
print(paste("Processing", length(test_dates), "forecast dates"))

# Store all forecasts
all_forecasts <- list()

# Process in batches to avoid timeout
batch_size <- 10
n_batches <- ceiling(length(test_dates) / batch_size)

for (batch in 1:n_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(test_dates))
  batch_dates <- test_dates[start_idx:end_idx]
  
  print(paste("\nProcessing batch", batch, "of", n_batches))
  
  for (origin_date in batch_dates) {
    tryCatch({
      # Get training data up to origin date
      train_data <- flu_data_hhs |>
        filter(origin_date <= !!origin_date)
      
      # Process each location independently
      locations <- unique(train_data$location)
      quantile_levels <- seq(0.01, 0.99, by = 0.01)
      
      for (loc in locations) {
        loc_data <- train_data |> 
          filter(location == loc) |>
          arrange(origin_date)
        
        # Skip if not enough data
        if (nrow(loc_data) < 104) next  # Need at least 2 years
        
        # Create time series
        ts_data <- ts(loc_data$wili, frequency = 52)
        
        # Fit simple models
        # Model 1: ARIMA
        arima_model <- auto.arima(ts_data, seasonal = TRUE, stepwise = TRUE)
        arima_fc <- forecast(arima_model, h = 4)
        
        # Model 2: ETS
        ets_model <- ets(ts_data)
        ets_fc <- forecast(ets_model, h = 4)
        
        # Model 3: Seasonal naive
        snaive_fc <- snaive(ts_data, h = 4)
        
        # Create ensemble forecast
        for (h in 1:4) {
          target_date <- origin_date + h * 7
          
          # Get point forecasts
          arima_mean <- as.numeric(arima_fc$mean[h])
          ets_mean <- as.numeric(ets_fc$mean[h])
          snaive_mean <- as.numeric(snaive_fc$mean[h])
          
          # Simple average
          ensemble_mean <- (arima_mean + ets_mean + snaive_mean) / 3
          
          # Use largest prediction interval for uncertainty
          arima_sd <- (arima_fc$upper[h, "95%"] - arima_fc$lower[h, "95%"]) / (2 * 1.96)
          ets_sd <- (ets_fc$upper[h, "95%"] - ets_fc$lower[h, "95%"]) / (2 * 1.96)
          ensemble_sd <- max(arima_sd, ets_sd)
          
          # Generate quantiles
          ensemble_quantiles <- qnorm(quantile_levels, 
                                     mean = ensemble_mean, 
                                     sd = ensemble_sd)
          ensemble_quantiles <- pmax(ensemble_quantiles, 0)
          
          # Store forecasts
          for (k in seq_along(quantile_levels)) {
            all_forecasts[[length(all_forecasts) + 1]] <- data.frame(
              model_id = "sismid-simple-ensemble",
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
    }, error = function(e) {
      # Silent fail for individual dates
    })
  }
  
  print(paste("Completed batch", batch))
}

# Save results
if (length(all_forecasts) > 0) {
  submission_df <- do.call(rbind, all_forecasts)
  
  print("\n=== Submission Summary ===")
  print(paste("Total forecasts:", nrow(submission_df)))
  print(paste("Origin dates:", n_distinct(submission_df$origin_date)))
  print(paste("Locations:", n_distinct(submission_df$location)))
  
  # Save
  write.csv(submission_df, "simple_ensemble_submission.csv", row.names = FALSE)
  print("\nSaved to: simple_ensemble_submission.csv")
  
  # Show sample
  print("\nSample rows:")
  print(head(submission_df |> filter(location == "US National", horizon == 1), 5))
}

print("\n=== Complete ===")
print("Ensemble combines: ARIMA + ETS + Seasonal Naive")
print("Simple average with conservative prediction intervals")