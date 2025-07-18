# Test VARIMA models with Moving Average components
# Systematic approach to test different MA configurations

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("hubEvals")
library("feasts")

# Set seed for reproducibility
set.seed(406)

# Load data
data(flu_data_hhs)

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Use same validation setup as other tests
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")

# Test on first 10 dates for initial evaluation
test_dates <- validation_dates[1:10]

print("=== Testing VARIMA Models with MA Components ===")
print(paste("Testing", length(test_dates), "validation dates"))

# Prepare data function - need wide format for VAR
prepare_var_data <- function(data, origin_date) {
  data |>
    filter(origin_date <= !!origin_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week) |>
    fill_gaps()
}

# Store results
model_results <- list()

# Test different VARIMA configurations
varima_configs <- list(
  "VARIMA(2,0,1)" = list(p = 2, d = 0, q = 1),
  "VARIMA(2,0,2)" = list(p = 2, d = 0, q = 2),
  "VARIMA(1,0,1)" = list(p = 1, d = 0, q = 1),
  "VARIMA(3,0,1)" = list(p = 3, d = 0, q = 1),
  "VARIMA(2,1,1)" = list(p = 2, d = 1, q = 1),
  "VARIMA(2,0,1)_seasonal" = list(p = 2, d = 0, q = 1, seasonal = TRUE)
)

for (model_name in names(varima_configs)) {
  config <- varima_configs[[model_name]]
  print(paste("\n--- Testing", model_name, "---"))
  
  model_forecasts <- list()
  errors <- 0
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Prepare training data
      train_data <- prepare_var_data(flu_data_hhs, origin_date)
      
      # Fit VARIMA model
      if (model_name == "VARIMA(2,0,1)_seasonal") {
        # VAR with MA component and Fourier seasonality
        model_fit <- train_data |>
          model(
            varima = VAR(vars(
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
            ) ~ AR(p = config$p) + MA(q = config$q) + fourier(K = 2))
          )
      } else if (config$d == 0) {
        # Standard VARIMA without differencing
        model_fit <- train_data |>
          model(
            varima = VAR(vars(
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
            ) ~ AR(p = config$p) + MA(q = config$q))
          )
      } else {
        # VARIMA with differencing - use ARIMA per location
        model_fit <- train_data |>
          pivot_longer(cols = -week, names_to = "location", values_to = "wili") |>
          as_tsibble(index = week, key = location) |>
          model(
            varima = ARIMA(my_sqrt(wili) ~ pdq(config$p, config$d, config$q))
          )
      }
      
      # Generate forecasts
      fc <- model_fit |>
        forecast(h = 4)
      
      # Extract quantiles
      fc_quantiles <- fc |>
        as_tibble() |>
        mutate(
          origin_date = origin_date,
          target = "wk inc flu hosp",
          target_end_date = as.Date(week),
          output_type = "quantile"
        )
      
      # Generate quantiles from distribution
      quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
      
      for (j in 1:nrow(fc_quantiles)) {
        row <- fc_quantiles[j, ]
        
        # Extract quantiles from distribution
        if (".distribution" %in% names(row)) {
          dist <- row$.distribution[[1]]
          q_values <- quantile(dist, probs = quantile_levels)
          
          for (k in seq_along(quantile_levels)) {
            model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
              origin_date = origin_date,
              target = "wk inc flu hosp",
              target_end_date = row$target_end_date,
              location = row$location,
              output_type = "quantile",
              output_type_id = quantile_levels[k],
              value = max(q_values[k], 0)  # Ensure non-negative
            )
          }
        }
      }
      
      if (i %% 3 == 0) {
        print(paste("Completed", i, "of", length(test_dates), "forecasts"))
      }
      
    }, error = function(e) {
      errors <- errors + 1
      if (errors <= 3) {  # Only print first 3 errors
        print(paste("Error for", model_name, "at", origin_date, ":", e$message))
      }
    })
  }
  
  # Calculate WIS if we have forecasts
  if (length(model_forecasts) > 0) {
    forecasts_df <- do.call(rbind, model_forecasts)
    
    # Get oracle data
    hub_path <- here::here("sismid-ili-forecasting-sandbox")
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(test_dates) + 7,
             target_end_date <= max(test_dates) + 28) |>
      collect()
    
    # Score the model
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    median_wis <- median(scores$wis, na.rm = TRUE)
    
    model_results[[model_name]] <- list(
      mean_wis = mean_wis,
      median_wis = median_wis,
      n_forecasts = nrow(forecasts_df),
      n_errors = errors
    )
    
    print(paste("Mean WIS:", round(mean_wis, 4), 
                "| Median WIS:", round(median_wis, 4)))
  } else {
    print(paste("No successful forecasts for", model_name))
  }
}

# Summary comparison
print("\n=== VARIMA Model Comparison ===")
print("Baseline VAR(2) sqrt: WIS = 0.208")
print(paste(rep("-", 50), collapse = ""))

# Create comparison dataframe
comparison_df <- do.call(rbind, lapply(names(model_results), function(x) {
  data.frame(
    Model = x,
    Mean_WIS = model_results[[x]]$mean_wis,
    Median_WIS = model_results[[x]]$median_wis,
    N_Forecasts = model_results[[x]]$n_forecasts,
    Improvement_pct = round((0.208 - model_results[[x]]$mean_wis) / 0.208 * 100, 1)
  )
}))

comparison_df <- comparison_df[order(comparison_df$Mean_WIS), ]

print(comparison_df)

# Save results
write.csv(comparison_df, "varima_ma_comparison.csv", row.names = FALSE)

# Identify best model
best_model <- comparison_df$Model[1]
best_wis <- comparison_df$Mean_WIS[1]

if (best_wis < 0.208) {
  print(paste("\nBEST MODEL:", best_model, "beats baseline!"))
  print(paste("Improvement:", comparison_df$Improvement_pct[1], "%"))
  print("Consider running full validation test with this configuration")
} else {
  print("\nNo VARIMA model with MA components beats the baseline VAR(2) sqrt")
  print("Consider other approaches or ensemble methods")
}