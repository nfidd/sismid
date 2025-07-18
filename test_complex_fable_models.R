# Test Complex Fable Models - Alternative approaches
# Focus on models that can actually be implemented in fable

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

log_transform <- function(x) log(x + 0.01)
inv_log <- function(x) exp(x) - 0.01
my_log <- new_transformation(log_transform, inv_log)

# Validation setup - test on subset first
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- validation_dates[1:10]

print("=== Testing Complex Fable Models ===")
print(paste("Testing", length(test_dates), "validation dates"))

# Store results
model_results <- list()

# Test different model approaches
model_configs <- list(
  "VAR_fourier_K3" = "VAR(2) with Fourier K=3",
  "VAR_fourier_trend" = "VAR(2) with Fourier K=2 + trend",
  "ARIMA_per_location" = "ARIMA per location with cross-correlation",
  "Dynamic_regression" = "Dynamic regression with exogenous variables",
  "Mixed_ensemble" = "Ensemble of VAR + individual ARIMAs"
)

for (model_name in names(model_configs)) {
  print(paste("\n--- Testing", model_configs[[model_name]], "---"))
  
  model_forecasts <- list()
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Prepare data in wide format for VAR
      train_data_wide <- flu_data_hhs |>
        filter(origin_date <= !!origin_date) |>
        select(origin_date, location, wili) |>
        pivot_wider(names_from = location, values_from = wili) |>
        mutate(week = yearweek(origin_date)) |>
        as_tsibble(index = week) |>
        fill_gaps()
      
      # Also prepare long format for individual models
      train_data_long <- flu_data_hhs |>
        filter(origin_date <= !!origin_date) |>
        mutate(week = yearweek(origin_date)) |>
        as_tsibble(index = week, key = location)
      
      if (model_name == "VAR_fourier_K3") {
        # VAR with higher order Fourier terms
        model_fit <- train_data_wide |>
          model(
            var_model = VAR(vars(
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
            ) ~ AR(2) + fourier(K = 3))
          )
        
      } else if (model_name == "VAR_fourier_trend") {
        # VAR with Fourier and trend
        model_fit <- train_data_wide |>
          model(
            var_model = VAR(vars(
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
            ) ~ AR(2) + trend() + fourier(K = 2))
          )
        
      } else if (model_name == "ARIMA_per_location") {
        # Individual ARIMA models per location with seasonal components
        model_fit <- train_data_long |>
          model(
            arima_model = ARIMA(my_sqrt(wili) ~ pdq(2,1,2) + PDQ(1,0,1))
          )
        
      } else if (model_name == "Dynamic_regression") {
        # Dynamic regression using national trend as exogenous variable
        # First get national data
        national_data <- train_data_wide |>
          select(week, `US National`) |>
          rename(national_wili = `US National`)
        
        # Join with long data
        train_data_dynamic <- train_data_long |>
          left_join(national_data, by = "week") |>
          filter(location != "US National")  # Don't model national against itself
        
        # Fit dynamic regression model
        model_fit <- train_data_dynamic |>
          model(
            dynamic_reg = ARIMA(my_sqrt(wili) ~ my_sqrt(national_wili) + pdq(2,0,1))
          )
        
        # Also need national model
        national_model <- train_data_long |>
          filter(location == "US National") |>
          model(
            national = ARIMA(my_sqrt(wili) ~ pdq(2,1,2) + PDQ(1,0,1))
          )
        
      } else if (model_name == "Mixed_ensemble") {
        # Ensemble approach: VAR for capturing cross-correlations + individual ARIMAs
        var_model <- train_data_wide |>
          model(
            var_part = VAR(vars(
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
        
        arima_model <- train_data_long |>
          model(
            arima_part = ARIMA(my_sqrt(wili))
          )
      }
      
      # Generate forecasts based on model type
      if (model_name == "Dynamic_regression") {
        # Special handling for dynamic regression
        # First forecast national
        nat_fc <- national_model |> forecast(h = 4)
        
        # Create future national values for dynamic regression
        future_national <- nat_fc |>
          as_tibble() |>
          select(week, wili) |>
          rename(national_wili = wili)
        
        # Forecast regions using future national values
        fc_regions <- model_fit |>
          forecast(new_data = future_national)
        
        # Combine forecasts
        fc <- bind_rows(fc_regions, nat_fc)
        
      } else if (model_name == "Mixed_ensemble") {
        # Generate forecasts from both models
        var_fc <- var_model |> forecast(h = 4)
        arima_fc <- arima_model |> forecast(h = 4)
        
        # Convert to consistent format and average
        var_tibble <- var_fc |>
          as_tibble() |>
          pivot_longer(cols = c(`HHS Region 1`:`US National`),
                      names_to = "location", 
                      values_to = ".distribution_var")
        
        arima_tibble <- arima_fc |>
          as_tibble() |>
          rename(.distribution_arima = wili)
        
        # Combine by averaging distributions
        fc <- var_tibble |>
          left_join(arima_tibble, by = c("week", "location")) |>
          mutate(
            # Simple average of the two models
            wili = map2(.distribution_var, .distribution_arima, 
                       ~(mean(.x) + mean(.y)) / 2),
            .distribution = .distribution_var  # Use VAR distribution for now
          ) |>
          select(week, location, wili, .distribution)
        
      } else {
        # Standard forecast
        fc <- model_fit |> forecast(h = 4)
      }
      
      # Convert to hub format
      quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
      
      fc_hub <- fc |>
        as_tibble() |>
        mutate(
          origin_date = origin_date,
          target = "wk inc flu hosp",
          target_end_date = as.Date(week)
        )
      
      # Extract quantiles
      for (j in 1:nrow(fc_hub)) {
        row <- fc_hub[j, ]
        
        if (".distribution" %in% names(row) && !is.null(row$.distribution[[1]])) {
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
              value = max(q_values[k], 0)
            )
          }
        }
      }
      
      if (i %% 3 == 0) {
        print(paste("Completed", i, "of", length(test_dates), "forecasts"))
      }
      
    }, error = function(e) {
      print(paste("Error at", origin_date, ":", e$message))
    })
  }
  
  # Calculate WIS
  if (length(model_forecasts) > 0) {
    forecasts_df <- do.call(rbind, model_forecasts)
    
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
    
    model_results[[model_name]] <- list(
      description = model_configs[[model_name]],
      mean_wis = mean_wis,
      n_forecasts = nrow(forecasts_df)
    )
    
    print(paste("Mean WIS:", round(mean_wis, 4)))
  }
}

# Summary
print("\n=== Complex Fable Model Comparison ===")
print("Baseline VAR(2) sqrt: WIS = 0.208")
print(paste(rep("=", 50), collapse = ""))

# Sort by WIS
if (length(model_results) > 0) {
  sorted_results <- model_results[order(sapply(model_results, function(x) x$mean_wis))]
  
  for (name in names(sorted_results)) {
    result <- sorted_results[[name]]
    improvement <- (0.208 - result$mean_wis) / 0.208 * 100
    
    print(sprintf("%-30s WIS: %.4f  Change: %+.1f%%", 
                  result$description, 
                  result$mean_wis, 
                  improvement))
  }
  
  # Best model
  best_name <- names(sorted_results)[1]
  if (sorted_results[[best_name]]$mean_wis < 0.208) {
    print(paste("\nBEST MODEL:", sorted_results[[best_name]]$description))
    print("This beats the baseline! Consider full validation.")
  }
}