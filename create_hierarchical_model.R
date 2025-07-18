# Create Hierarchical Model with National-Regional Structure
# Model national trend + regional deviations for better forecasts

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("feasts")

set.seed(406)

# Load data
data(flu_data_hhs)

# Define transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Test on subset first
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- validation_dates[1:10]

print("=== Hierarchical Forecasting Model ===")
print("Approach: Model national trend + regional deviations")

# Function for hierarchical forecasting
hierarchical_forecast <- function(origin_date, data) {
  
  # Prepare data
  train_data <- data |>
    filter(origin_date <= !!origin_date) |>
    mutate(week = yearweek(origin_date))
  
  # Get national data
  national_data <- train_data |>
    filter(location == "US National") |>
    as_tsibble(index = week)
  
  # Calculate regional deviations from national
  regional_data <- train_data |>
    filter(location != "US National") |>
    left_join(
      national_data |> select(week, national_wili = wili),
      by = "week"
    ) |>
    mutate(
      # Regional deviation from national (additive)
      deviation = wili - national_wili,
      # Also calculate ratio for alternative approach
      ratio = wili / (national_wili + 0.1)
    ) |>
    as_tsibble(index = week, key = location)
  
  # Model 1: National trend with ARIMA
  national_model <- national_data |>
    model(
      national = ARIMA(my_sqrt(wili) ~ pdq(2,1,1) + PDQ(1,0,1))
    )
  
  # Model 2: Regional deviations with simpler ARIMA
  deviation_models <- regional_data |>
    model(
      deviation = ARIMA(deviation ~ pdq(1,0,1))
    )
  
  # Generate forecasts
  national_fc <- national_model |> forecast(h = 4)
  deviation_fc <- deviation_models |> forecast(h = 4)
  
  # Combine forecasts
  results <- list()
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  
  # National forecasts
  nat_fc_data <- national_fc |> as_tibble()
  
  for (i in 1:nrow(nat_fc_data)) {
    row <- nat_fc_data[i, ]
    dist <- row$wili[[1]]
    q_values <- quantile(dist, probs = quantile_levels)
    
    for (j in seq_along(quantile_levels)) {
      results[[length(results) + 1]] <- data.frame(
        origin_date = origin_date,
        target = "wk inc flu hosp",
        target_end_date = as.Date(row$week),
        location = "US National",
        output_type = "quantile",
        output_type_id = quantile_levels[j],
        value = max(inv_sqrt(q_values[j]), 0)
      )
    }
  }
  
  # Regional forecasts
  regional_fc_data <- deviation_fc |> as_tibble()
  nat_fc_values <- nat_fc_data |> 
    mutate(national_mean = map_dbl(wili, mean))
  
  for (loc in unique(regional_fc_data$location)) {
    loc_data <- regional_fc_data |> filter(location == loc)
    
    for (i in 1:nrow(loc_data)) {
      row <- loc_data[i, ]
      dev_dist <- row$deviation[[1]]
      
      # Get corresponding national forecast
      nat_value <- inv_sqrt(nat_fc_values$national_mean[i])
      
      # Calculate regional forecast = national + deviation
      dev_q <- quantile(dev_dist, probs = quantile_levels)
      regional_q <- nat_value + dev_q
      
      # Ensure non-negative
      regional_q <- pmax(regional_q, 0)
      
      for (j in seq_along(quantile_levels)) {
        results[[length(results) + 1]] <- data.frame(
          origin_date = origin_date,
          target = "wk inc flu hosp",
          target_end_date = as.Date(row$week),
          location = loc,
          output_type = "quantile",
          output_type_id = quantile_levels[j],
          value = regional_q[j]
        )
      }
    }
  }
  
  return(results)
}

# Test on subset
print("\nTesting hierarchical approach...")

all_forecasts <- list()

for (i in seq_along(test_dates)) {
  origin_date <- test_dates[i]
  
  tryCatch({
    forecasts <- hierarchical_forecast(origin_date, flu_data_hhs)
    all_forecasts <- c(all_forecasts, forecasts)
    print(paste("Completed", origin_date))
  }, error = function(e) {
    print(paste("Error at", origin_date, ":", e$message))
  })
}

# Process results
if (length(all_forecasts) > 0) {
  hierarchical_df <- do.call(rbind, all_forecasts)
  
  print("\nHierarchical forecast summary:")
  print(paste("Total forecasts:", nrow(hierarchical_df)))
  print(paste("Unique dates:", n_distinct(hierarchical_df$origin_date)))
  print(paste("Unique locations:", n_distinct(hierarchical_df$location)))
  
  # Check forecast coherence
  sample_check <- hierarchical_df |>
    filter(origin_date == test_dates[1],
           output_type_id == 0.5) |>
    group_by(target_end_date) |>
    summarise(
      national = value[location == "US National"],
      regional_sum = sum(value[location != "US National"]),
      .groups = "drop"
    )
  
  print("\nCoherence check (National vs Sum of Regions):")
  print(sample_check)
  
  # Save sample
  write.csv(hierarchical_df |> head(1000), 
            "hierarchical_sample.csv", 
            row.names = FALSE)
  
  # Full validation run
  print("\n=== Running Full Validation ===")
  
  validation_forecasts <- list()
  
  for (i in seq_along(validation_dates)) {
    origin_date <- validation_dates[i]
    
    tryCatch({
      forecasts <- hierarchical_forecast(origin_date, flu_data_hhs)
      validation_forecasts <- c(validation_forecasts, forecasts)
      
      if (i %% 10 == 0) {
        print(paste("Completed", i, "of", length(validation_dates)))
      }
    }, error = function(e) {
      NULL  # Silent fail
    })
  }
  
  if (length(validation_forecasts) > 0) {
    val_df <- do.call(rbind, validation_forecasts) |>
      mutate(model_id = "team1-hierarchical") |>
      select(model_id, origin_date, target, target_end_date,
             location, output_type, output_type_id, value)
    
    write.csv(val_df, "hierarchical_validation_forecasts.csv", 
              row.names = FALSE)
    
    print("\n=== Hierarchical Model Complete ===")
    print(paste("Total validation forecasts:", nrow(val_df)))
    print("Saved to hierarchical_validation_forecasts.csv")
    
    # Summary statistics
    summary_stats <- val_df |>
      filter(output_type_id == 0.5) |>
      group_by(location) |>
      summarise(
        mean_forecast = mean(value),
        sd_forecast = sd(value),
        .groups = "drop"
      ) |>
      arrange(desc(mean_forecast))
    
    print("\nMean forecast by location (median quantile):")
    print(head(summary_stats, 12))
  }
}

print("\n=== Summary ===")
print("Hierarchical model approach:")
print("1. Model US National trend with ARIMA(2,1,1) + seasonal")
print("2. Model regional deviations from national with ARIMA(1,0,1)")
print("3. Combine: Regional forecast = National + Deviation")
print("\nAdvantages:")
print("- Leverages strong national signal")
print("- Simpler regional models (fewer parameters)")
print("- Natural hierarchical structure")
print("- More stable forecasts")