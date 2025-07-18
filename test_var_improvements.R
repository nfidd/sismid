# Test VAR Model Improvements
# Focus on variations of VAR that might beat baseline

library("nfidd")
library("dplyr")
library("tidyr")
library("vars")
library("forecast")
library("hubEvals")

set.seed(406)

# Load data
data(flu_data_hhs)

# Use same validation setup
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")

# Test on subset first
test_dates <- validation_dates[1:5]

print("=== Testing VAR Model Improvements ===")

# Function to fit VAR and generate hub forecasts
generate_var_forecasts <- function(train_data, config, origin_date) {
  # Apply transformation
  transform_fn <- config$transform
  inverse_fn <- config$inverse
  
  # Prepare data
  numeric_cols <- setdiff(names(train_data), "origin_date")
  train_data[numeric_cols] <- lapply(train_data[numeric_cols], transform_fn)
  
  # Create time series
  ts_data <- ts(train_data[numeric_cols], frequency = 52)
  
  # Add external regressors if specified
  exogen <- NULL
  if (!is.null(config$seasonality)) {
    if (config$seasonality == "fourier") {
      # Create Fourier terms
      exogen <- fourier(ts_data, K = config$K)
    } else if (config$seasonality == "dummy") {
      # Create seasonal dummies
      exogen <- seasonaldummy(ts_data)
    }
  }
  
  # Fit VAR model
  if (!is.null(exogen)) {
    var_model <- VAR(ts_data, p = config$lags, type = config$type, 
                     exogen = exogen)
  } else {
    var_model <- VAR(ts_data, p = config$lags, type = config$type)
  }
  
  # Generate forecasts
  if (!is.null(exogen)) {
    # Need future values of exogenous variables
    future_exogen <- if (config$seasonality == "fourier") {
      fourier(ts_data, K = config$K, h = 4)
    } else {
      # For seasonal dummies, need to calculate future periods
      future_periods <- (nrow(ts_data) + 1:4 - 1) %% 52 + 1
      future_dummy <- matrix(0, nrow = 4, ncol = 51)
      for (i in 1:4) {
        if (future_periods[i] < 52) {
          future_dummy[i, future_periods[i]] <- 1
        }
      }
      future_dummy
    }
    forecast_obj <- predict(var_model, n.ahead = 4, 
                           dumvar = future_exogen)
  } else {
    forecast_obj <- predict(var_model, n.ahead = 4)
  }
  
  # Extract forecasts and create hub format
  model_forecasts <- list()
  
  for (loc in numeric_cols) {
    pred_values <- forecast_obj$fcst[[loc]][, "fcst"]
    pred_values <- pmax(inverse_fn(pred_values), 0)
    
    # Get prediction errors
    se_values <- forecast_obj$fcst[[loc]][, "CI"]
    
    # Generate quantiles
    quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    
    for (h in 1:4) {
      target_date <- origin_date + h * 7
      
      # Use normal approximation for quantiles
      mean_val <- pred_values[h]
      # Estimate SD from prediction interval
      sd_val <- (se_values[h] / 1.96) * sqrt(h)  # Adjust for horizon
      
      for (q in quantiles) {
        q_val <- qnorm(q, mean = mean_val, sd = sd_val)
        q_val <- max(q_val, 0)
        
        model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
          origin_date = origin_date,
          target = "wk inc flu hosp",
          target_end_date = target_date,
          location = loc,
          output_type = "quantile",
          output_type_id = q,
          value = q_val
        )
      }
    }
  }
  
  return(model_forecasts)
}

# Test configurations
configs <- list(
  "VAR2_sqrt_baseline" = list(
    lags = 2, 
    transform = sqrt, 
    inverse = function(x) x^2,
    type = "const",
    seasonality = NULL
  ),
  "VAR3_sqrt" = list(
    lags = 3, 
    transform = sqrt, 
    inverse = function(x) x^2,
    type = "const",
    seasonality = NULL
  ),
  "VAR2_sqrt_fourier2" = list(
    lags = 2, 
    transform = sqrt, 
    inverse = function(x) x^2,
    type = "const",
    seasonality = "fourier",
    K = 2
  ),
  "VAR2_sqrt_trend" = list(
    lags = 2, 
    transform = sqrt, 
    inverse = function(x) x^2,
    type = "trend",
    seasonality = NULL
  ),
  "VAR2_log" = list(
    lags = 2, 
    transform = function(x) log(x + 0.01), 
    inverse = function(x) exp(x) - 0.01,
    type = "const",
    seasonality = NULL
  )
)

# Test each configuration
results <- list()

for (config_name in names(configs)) {
  print(paste("\n--- Testing", config_name, "---"))
  config <- configs[[config_name]]
  
  all_forecasts <- list()
  errors <- 0
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Get training data
      train_data <- flu_data_hhs |>
        dplyr::filter(origin_date <= !!origin_date) |>
        dplyr::select(origin_date, location, wili) |>
        pivot_wider(names_from = location, values_from = wili) |>
        arrange(origin_date)
      
      # Generate forecasts
      forecasts <- generate_var_forecasts(train_data, config, origin_date)
      all_forecasts <- c(all_forecasts, forecasts)
      
    }, error = function(e) {
      errors <<- errors + 1
      if (errors <= 2) {
        print(paste("Error:", e$message))
      }
    })
  }
  
  # Calculate WIS
  if (length(all_forecasts) > 0) {
    forecasts_df <- do.call(rbind, all_forecasts)
    
    hub_path <- here::here("sismid-ili-forecasting-sandbox")
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(test_dates) + 7,
             target_end_date <= max(test_dates) + 28) |>
      collect()
    
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    results[[config_name]] <- mean_wis
    
    print(paste("Mean WIS:", round(mean_wis, 4)))
  }
}

# Summary
print("\n=== VAR Model Comparison ===")
print(paste(rep("=", 50), collapse = ""))

# Sort by WIS
sorted_results <- sort(unlist(results))

for (i in seq_along(sorted_results)) {
  name <- names(sorted_results)[i]
  wis <- sorted_results[i]
  
  if (name == "VAR2_sqrt_baseline") {
    print(sprintf("%d. %-25s WIS: %.4f (BASELINE)", i, name, wis))
  } else {
    improvement <- (results[["VAR2_sqrt_baseline"]] - wis) / 
                  results[["VAR2_sqrt_baseline"]] * 100
    print(sprintf("%d. %-25s WIS: %.4f (%+.1f%%)", 
                  i, name, wis, improvement))
  }
}

# Check against 0.208 benchmark
print("\nComparison to reported baseline (0.208):")
best_wis <- min(sorted_results)
if (best_wis < 0.208) {
  print(paste("FOUND IMPROVEMENT! Best model:", names(sorted_results)[1]))
  print(paste("WIS:", round(best_wis, 4)))
  print("Run full validation to confirm.")
} else {
  print("No improvement found over baseline.")
}