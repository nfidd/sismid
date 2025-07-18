# Create Practical Submission
# VAR(2) with sqrt transformation + small seasonal adjustment

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("lubridate")

set.seed(406)

# Load data
data(flu_data_hhs)

# Define transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Hub paths
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Get all origin dates for submission
origin_dates <- hub_path |> 
  hubUtils::read_config("tasks") |> 
  hubUtils::get_round_ids()

# We'll use test phase dates
test_dates <- origin_dates[origin_dates >= as.Date("2017-11-04")]

print("=== Creating Practical Submission ===")
print(paste("Number of forecast dates:", length(test_dates)))

# Store all forecasts
all_forecasts <- list()

# Generate forecasts for each origin date
for (i in seq_along(test_dates)) {
  origin_date <- test_dates[i]
  
  tryCatch({
    # Prepare data in wide format
    train_data <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      select(origin_date, location, wili) |>
      pivot_wider(names_from = location, values_from = wili) |>
      mutate(week = yearweek(origin_date)) |>
      select(-origin_date) |>
      as_tsibble(index = week) |>
      fill_gaps()
    
    # Add seasonal indicator
    train_data <- train_data |>
      mutate(
        week_num = week(week),
        # Simple seasonal adjustment - higher weight during typical flu season
        seasonal_weight = case_when(
          week_num >= 40 | week_num <= 20 ~ 1.1,  # Oct-Apr
          TRUE ~ 0.9  # May-Sep
        )
      )
    
    # Fit VAR(2) model
    model_fit <- train_data |>
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
    
    # Generate forecasts
    forecasts <- model_fit |> forecast(h = 4)
    
    # Get seasonal weights for forecast period
    forecast_weeks <- week(forecasts$week)
    seasonal_adj <- case_when(
      forecast_weeks >= 40 | forecast_weeks <= 20 ~ 1.05,
      TRUE ~ 0.95
    )
    
    # Convert to long format
    fc_long <- forecasts |>
      as_tibble() |>
      pivot_longer(
        cols = c(`HHS Region 1`:`US National`),
        names_to = "location",
        values_to = ".distribution"
      )
    
    # Extract quantiles
    quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    
    for (j in 1:nrow(fc_long)) {
      row <- fc_long[j, ]
      dist <- row$.distribution[[1]]
      
      # Get quantiles
      q_values <- quantile(dist, probs = quantile_levels)
      
      # Transform back and apply seasonal adjustment
      adj_idx <- which(fc_long$week == row$week)[1]
      q_values_adj <- inv_sqrt(q_values) * seasonal_adj[adj_idx]
      q_values_adj <- pmax(q_values_adj, 0)  # Ensure non-negative
      
      # Store forecasts
      for (k in seq_along(quantile_levels)) {
        all_forecasts[[length(all_forecasts) + 1]] <- data.frame(
          model_id = "nfidd-var2seasonal",
          origin_date = origin_date,
          target = "wk inc flu hosp",
          target_end_date = as.Date(row$week),
          location = row$location,
          output_type = "quantile",
          output_type_id = quantile_levels[k],
          value = q_values_adj[k]
        )
      }
    }
    
    if (i %% 10 == 0) {
      print(paste("Completed", i, "of", length(test_dates), "forecasts"))
    }
    
  }, error = function(e) {
    print(paste("Error at", origin_date, ":", e$message))
  })
}

# Combine all forecasts
if (length(all_forecasts) > 0) {
  submission_df <- do.call(rbind, all_forecasts)
  
  # Convert location names to hub format
  location_map <- data.frame(
    location = c("US National", paste("HHS Region", 1:10)),
    hub_location = c("US", paste0("0", 1:9), "10")
  )
  
  submission_df <- submission_df |>
    left_join(location_map, by = "location") |>
    mutate(location = hub_location) |>
    select(-hub_location)
  
  # Validate format
  print("\nSubmission summary:")
  print(paste("Total rows:", nrow(submission_df)))
  print(paste("Unique origin dates:", n_distinct(submission_df$origin_date)))
  print(paste("Unique locations:", n_distinct(submission_df$location)))
  print(paste("Unique quantiles:", n_distinct(submission_df$output_type_id)))
  
  # Check a sample
  sample_check <- submission_df |>
    filter(origin_date == test_dates[1],
           location == "US",
           output_type_id == 0.5) |>
    arrange(target_end_date)
  
  print("\nSample forecasts (US National, median):")
  print(sample_check)
  
  # Save submission
  filename <- paste0("VAR2_seasonal_submission_", Sys.Date(), ".csv")
  write.csv(submission_df, filename, row.names = FALSE)
  
  print(paste("\nSubmission saved to:", filename))
  print("Model: VAR(2) with sqrt transformation + seasonal adjustment")
  print("Ready for hub submission!")
  
} else {
  print("No forecasts generated!")
}

print("\n=== Submission Complete ===")
print("Model details:")
print("- Base: VAR(2) with sqrt transformation")
print("- Enhancement: Seasonal adjustment (5% boost Oct-Apr)")
print("- Simple, stable, and builds on proven approach")