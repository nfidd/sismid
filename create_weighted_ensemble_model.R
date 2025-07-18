# Create Weighted Ensemble Model for Hub Submission
# Combines VAR(2) sqrt with individual ARIMA models using performance weights

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("hubEvals")
library("hubUtils")

set.seed(406)

# Load data
data(flu_data_hhs)

# Define transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")
hub_con <- connect_hub(hub_path)

# Define dates
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- seq(as.Date("2017-11-04"), as.Date("2019-05-18"), by = "week")

print("=== Creating Weighted Ensemble Model ===")

# Function to generate forecasts for a single origin date
generate_ensemble_forecasts <- function(origin_date, train_data, weights = NULL) {
  
  # Prepare data in different formats
  train_long <- train_data |>
    filter(origin_date <= !!origin_date) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week, key = location)
  
  train_wide <- train_data |>
    filter(origin_date <= !!origin_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    mutate(week = yearweek(origin_date)) |>
    as_tsibble(index = week)
  
  # Fit VAR(2) model
  var_model <- train_wide |>
    model(
      var2 = VAR(vars(
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
  
  # Fit individual ARIMA models
  arima_models <- train_long |>
    model(
      arima = ARIMA(my_sqrt(wili) ~ pdq(2,1,1) + PDQ(1,0,1))
    )
  
  # Generate forecasts
  var_fc <- var_model |> 
    forecast(h = 4) |>
    as_tibble() |>
    pivot_longer(cols = -c(.model, week), 
                 names_to = "location", 
                 values_to = ".distribution_var")
  
  arima_fc <- arima_models |> 
    forecast(h = 4) |>
    as_tibble() |>
    rename(.distribution_arima = wili)
  
  # Combine forecasts with weights
  combined_fc <- var_fc |>
    left_join(arima_fc, by = c("location", "week")) |>
    mutate(
      # Default weights if not provided
      var_weight = ifelse(is.null(weights), 0.6, weights$var),
      arima_weight = ifelse(is.null(weights), 0.4, weights$arima),
      
      # Weighted combination
      .distribution = purrr::map2(
        .distribution_var, .distribution_arima,
        function(d1, d2) {
          # Extract mean and variance
          m1 <- mean(d1)
          m2 <- mean(d2)
          v1 <- distributional::variance(d1)
          v2 <- distributional::variance(d2)
          
          # Weighted mean
          w1 <- var_weight[1]
          w2 <- arima_weight[1]
          mean_combined <- w1 * m1 + w2 * m2
          
          # Approximate combined variance
          var_combined <- w1^2 * v1 + w2^2 * v2
          
          # Return normal distribution
          distributional::dist_normal(mean_combined, sqrt(var_combined))
        }
      )
    )
  
  # Convert to quantile format
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  
  forecasts <- list()
  
  for (i in 1:nrow(combined_fc)) {
    row <- combined_fc[i, ]
    dist <- row$.distribution[[1]]
    
    # Extract quantiles
    q_values <- quantile(dist, probs = quantile_levels)
    
    for (j in seq_along(quantile_levels)) {
      forecasts[[length(forecasts) + 1]] <- data.frame(
        origin_date = origin_date,
        target = "wk inc flu hosp",
        target_end_date = as.Date(row$week),
        location = row$location,
        output_type = "quantile",
        output_type_id = quantile_levels[j],
        value = max(inv_sqrt(q_values[j]), 0)  # Transform back and ensure non-negative
      )
    }
  }
  
  return(forecasts)
}

# Step 1: Determine optimal weights using validation data
print("\nStep 1: Determining optimal weights using validation data...")

# Test different weight combinations
weight_options <- list(
  w1 = list(var = 0.7, arima = 0.3),
  w2 = list(var = 0.6, arima = 0.4),
  w3 = list(var = 0.5, arima = 0.5),
  w4 = list(var = 0.4, arima = 0.6),
  w5 = list(var = 0.3, arima = 0.7)
)

weight_results <- list()

# Test on subset of validation dates
test_val_dates <- validation_dates[seq(1, length(validation_dates), by = 4)][1:5]

for (weight_name in names(weight_options)) {
  weights <- weight_options[[weight_name]]
  print(paste("Testing weights:", weight_name, 
              "- VAR:", weights$var, "ARIMA:", weights$arima))
  
  all_forecasts <- list()
  
  for (origin_date in test_val_dates) {
    tryCatch({
      forecasts <- generate_ensemble_forecasts(origin_date, flu_data_hhs, weights)
      all_forecasts <- c(all_forecasts, forecasts)
    }, error = function(e) {
      print(paste("Error at", origin_date, ":", e$message))
    })
  }
  
  if (length(all_forecasts) > 0) {
    forecasts_df <- do.call(rbind, all_forecasts)
    
    # Get oracle data
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(test_val_dates) + 7,
             target_end_date <= max(test_val_dates) + 28) |>
      collect()
    
    # Score
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    weight_results[[weight_name]] <- mean_wis
    print(paste("  Mean WIS:", round(mean_wis, 4)))
  }
}

# Find best weights
best_weight_name <- names(which.min(weight_results))
best_weights <- weight_options[[best_weight_name]]
best_wis <- weight_results[[best_weight_name]]

print(paste("\nBest weights:", best_weight_name))
print(paste("VAR weight:", best_weights$var, "ARIMA weight:", best_weights$arima))
print(paste("Validation WIS:", round(best_wis, 4)))

# Step 2: Generate full validation forecasts with best weights
print("\nStep 2: Generating full validation forecasts with optimal weights...")

validation_forecasts <- list()

for (i in seq_along(validation_dates)) {
  origin_date <- validation_dates[i]
  
  tryCatch({
    forecasts <- generate_ensemble_forecasts(origin_date, flu_data_hhs, best_weights)
    validation_forecasts <- c(validation_forecasts, forecasts)
    
    if (i %% 10 == 0) {
      print(paste("Completed", i, "of", length(validation_dates), "validation forecasts"))
    }
  }, error = function(e) {
    print(paste("Error at", origin_date, ":", e$message))
  })
}

# Evaluate validation performance
if (length(validation_forecasts) > 0) {
  val_forecasts_df <- do.call(rbind, validation_forecasts)
  
  oracle_output <- connect_target_oracle_output(hub_path) |>
    filter(target_end_date >= min(validation_dates) + 7,
           target_end_date <= max(validation_dates) + 28) |>
    collect()
  
  val_scores <- score_model_out(
    val_forecasts_df,
    oracle_output,
    metrics = c("wis", "ae", "coverage")
  )
  
  val_mean_wis <- mean(val_scores$wis, na.rm = TRUE)
  
  print("\n=== Validation Results ===")
  print(paste("Mean WIS:", round(val_mean_wis, 4)))
  print(paste("Median WIS:", round(median(val_scores$wis, na.rm = TRUE), 4)))
  print(paste("Mean AE:", round(mean(val_scores$ae, na.rm = TRUE), 4)))
  
  # Compare to baseline
  improvement <- (0.208 - val_mean_wis) / 0.208 * 100
  print(paste("\nComparison to VAR(2) sqrt baseline (0.208):"))
  print(paste("Improvement:", round(improvement, 1), "%"))
  
  if (val_mean_wis < 0.208) {
    print("SUCCESS: Weighted ensemble beats baseline!")
  }
}

# Step 3: Generate test phase forecasts if validation successful
if (val_mean_wis < 0.208) {
  print("\nStep 3: Generating test phase forecasts for submission...")
  
  test_forecasts <- list()
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      forecasts <- generate_ensemble_forecasts(origin_date, flu_data_hhs, best_weights)
      test_forecasts <- c(test_forecasts, forecasts)
      
      if (i %% 10 == 0) {
        print(paste("Completed", i, "of", length(test_dates), "test forecasts"))
      }
    }, error = function(e) {
      print(paste("Error at", origin_date, ":", e$message))
    })
  }
  
  # Format for submission
  if (length(test_forecasts) > 0) {
    submission_df <- do.call(rbind, test_forecasts) |>
      # Add model metadata
      mutate(model_id = "team1-ensemble1") |>
      # Reorder columns for hub format
      select(model_id, origin_date, target, target_end_date, 
             location, output_type, output_type_id, value)
    
    # Save submission
    filename <- paste0("weighted_ensemble_submission_", Sys.Date(), ".csv")
    write.csv(submission_df, filename, row.names = FALSE)
    
    print(paste("\nSubmission saved to:", filename))
    print(paste("Total forecasts:", nrow(submission_df)))
    print(paste("Unique origin dates:", n_distinct(submission_df$origin_date)))
    print(paste("Unique locations:", n_distinct(submission_df$location)))
  }
}

# Save results summary
results_summary <- list(
  best_weights = best_weights,
  validation_wis = val_mean_wis,
  baseline_wis = 0.208,
  improvement = improvement,
  weight_comparison = weight_results
)

saveRDS(results_summary, "weighted_ensemble_results.rds")

print("\n=== Weighted Ensemble Model Complete ===")
print(paste("Best configuration: VAR weight =", best_weights$var, 
            ", ARIMA weight =", best_weights$arima))
print(paste("Validation WIS:", round(val_mean_wis, 4)))
print(paste("Improvement over baseline:", round(improvement, 1), "%"))