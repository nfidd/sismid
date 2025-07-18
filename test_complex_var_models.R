# Test more complex VAR models for potential second submission

library("nfidd")
library("dplyr")
library("vars")
library("hubEvals")

# Load and prepare data
data(flu_data_hhs)

# Use same validation setup as before
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- seq(as.Date("2017-11-04"), as.Date("2019-05-18"), by = "week")

print("=== Testing Complex VAR Models ===")
print(paste("Testing", length(validation_dates), "validation dates"))

# Test configurations
models_to_test <- list(
  "VAR(3)_sqrt" = list(lags = 3, transformation = "sqrt"),
  "VAR(4)_sqrt" = list(lags = 4, transformation = "sqrt"),
  "VAR(3)_log" = list(lags = 3, transformation = "log"),
  "VAR(4)_log" = list(lags = 4, transformation = "log"),
  "VAR(3)_fourth" = list(lags = 3, transformation = "fourth_root"),
  "VAR(2)_log" = list(lags = 2, transformation = "log")
)

results <- list()

for (model_name in names(models_to_test)) {
  config <- models_to_test[[model_name]]
  
  print(paste("\n--- Testing", model_name, "---"))
  
  # Apply transformation
  if (config$transformation == "sqrt") {
    transform_fn <- sqrt
    inverse_fn <- function(x) x^2
  } else if (config$transformation == "log") {
    transform_fn <- function(x) log(x + 0.01)  # Add small constant for zeros
    inverse_fn <- function(x) exp(x) - 0.01
  } else if (config$transformation == "fourth_root") {
    transform_fn <- function(x) x^0.25
    inverse_fn <- function(x) x^4
  }
  
  model_forecasts <- list()
  
  for (i in seq_along(validation_dates)) {
    origin_date <- validation_dates[i]
    
    # Get training data up to origin date
    train_data <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      select(origin_date, location, wili) |>
      pivot_wider(names_from = location, values_from = wili) |>
      arrange(origin_date)
    
    # Apply transformation
    numeric_cols <- setdiff(names(train_data), "origin_date")
    train_data[numeric_cols] <- lapply(train_data[numeric_cols], transform_fn)
    
    # Create time series
    ts_data <- ts(train_data[numeric_cols], frequency = 52)
    
    # Fit VAR model
    tryCatch({
      var_model <- VAR(ts_data, p = config$lags, type = "const")
      
      # Generate forecast
      forecast_obj <- predict(var_model, n.ahead = 4)
      
      # Extract forecasts and transform back
      for (loc in numeric_cols) {
        pred_values <- forecast_obj$fcst[[loc]][, "fcst"]
        pred_values <- pmax(inverse_fn(pred_values), 0)  # Ensure non-negative
        
        # Create forecast for each horizon
        for (h in 1:4) {
          target_date <- origin_date + h * 7
          
          # Bootstrap prediction intervals
          n_boot <- 1000
          boot_preds <- replicate(n_boot, {
            resid_sample <- sample(residuals(var_model)[, loc], 
                                   size = config$lags, replace = TRUE)
            pred_val <- pred_values[h] + rnorm(1, 0, sd(resid_sample))
            max(inverse_fn(pred_val), 0)
          })
          
          # Calculate quantiles
          quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
          pred_quantiles <- quantile(boot_preds, quantiles)
          
          # Store forecast
          for (q in seq_along(quantiles)) {
            model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
              origin_date = origin_date,
              target = "wk inc flu hosp",
              target_end_date = target_date,
              location = loc,
              output_type = "quantile",
              output_type_id = quantiles[q],
              value = pred_quantiles[q]
            )
          }
        }
      }
      
      if (i %% 5 == 0) {
        print(paste("Completed", i, "of", length(validation_dates), "forecasts"))
      }
      
    }, error = function(e) {
      print(paste("Error for", model_name, "at", origin_date, ":", e$message))
    })
  }
  
  # Combine forecasts
  if (length(model_forecasts) > 0) {
    forecasts_df <- do.call(rbind, model_forecasts)
    
    # Calculate WIS using hubEvals
    hub_path <- here::here("sismid-ili-forecasting-sandbox")
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(validation_dates) + 7,
             target_end_date <= max(validation_dates) + 28) |>
      collect()
    
    # Score the model
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    results[[model_name]] <- list(
      wis = mean_wis,
      forecasts = forecasts_df
    )
    
    print(paste("Mean WIS for", model_name, ":", round(mean_wis, 4)))
  }
}

# Compare results
print("\n=== Model Comparison ===")
wis_scores <- sapply(results, function(x) x$wis)
wis_scores <- sort(wis_scores)

for (i in seq_along(wis_scores)) {
  model_name <- names(wis_scores)[i]
  wis <- wis_scores[i]
  print(paste(i, ".", model_name, "- WIS:", round(wis, 4)))
}

# Save best model if it beats our current best (0.208)
best_model <- names(wis_scores)[1]
best_wis <- wis_scores[1]

print(paste("\nBest model:", best_model, "with WIS:", round(best_wis, 4)))
print(paste("Current best VAR(2) sqrt: 0.208"))

if (best_wis < 0.208) {
  print("NEW BEST MODEL FOUND!")
  print("Saving forecasts for second submission...")
  
  # Save the best model's forecasts
  write.csv(results[[best_model]]$forecasts, 
            paste0("complex_model_forecasts_", gsub("[^A-Za-z0-9]", "_", best_model), ".csv"),
            row.names = FALSE)
  
  improvement <- ((0.208 - best_wis) / 0.208) * 100
  print(paste("Improvement:", round(improvement, 1), "%"))
} else {
  print("No improvement found over current VAR(2) sqrt model")
}