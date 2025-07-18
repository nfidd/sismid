# Generate full test phase forecasts for VAR(2) log model

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")

set.seed(406)

# Load data
data(flu_data_hhs)

# Define log transformation
log_transform <- function(x) log(x + 0.1)
inv_log <- function(x) exp(x) - 0.1
my_log <- new_transformation(log_transform, inv_log)

# Get all test phase dates
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter for test phase dates (after 2017-05-06)
test_dates <- origin_dates[which(as.Date(origin_dates) > as.Date("2017-05-06"))]

print(paste("Generating forecasts for", length(test_dates), "test phase dates"))
print(paste("First date:", min(test_dates)))
print(paste("Last date:", max(test_dates)))

# Function to generate forecast for a single date
generate_single_forecast <- function(origin_date) {
  
  # Prepare data up to origin date
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= origin_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Fit VAR(2) log model
  model_fit <- flu_data_wide |>
    model(
      var2_log = VAR(vars(
        `HHS Region 1` = my_log(`HHS Region 1`),
        `HHS Region 2` = my_log(`HHS Region 2`),
        `HHS Region 3` = my_log(`HHS Region 3`),
        `HHS Region 4` = my_log(`HHS Region 4`),
        `HHS Region 5` = my_log(`HHS Region 5`),
        `HHS Region 6` = my_log(`HHS Region 6`),
        `HHS Region 7` = my_log(`HHS Region 7`),
        `HHS Region 8` = my_log(`HHS Region 8`),
        `HHS Region 9` = my_log(`HHS Region 9`),
        `HHS Region 10` = my_log(`HHS Region 10`),
        `US National` = my_log(`US National`)
      ) ~ AR(2))
    )
  
  # Generate forecasts
  forecast_samples <- model_fit |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  # Convert to long format
  location_cols <- c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                     "HHS Region 4", "HHS Region 5", "HHS Region 6",
                     "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                     "HHS Region 10", "US National")
  
  forecast_long <- forecast_samples |>
    pivot_longer(
      cols = all_of(location_cols),
      names_to = "location",
      values_to = "value"
    ) |>
    rename(target_end_date = origin_date) |>
    mutate(
      model_id = "sismid-var2-log",
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
  
  # Calculate quantiles
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- forecast_long |>
    group_by(model_id, location, origin_date, horizon, target_end_date) |>
    reframe(
      quantile_level = quantile_levels,
      value = quantile(value, quantile_levels, na.rm = TRUE)
    ) |>
    mutate(
      output_type = "quantile",
      output_type_id = quantile_level
    ) |>
    select(-quantile_level)
  
  return(forecast_quantiles)
}

# Generate forecasts for first 10 test dates as a sample
sample_dates <- test_dates[1:10]

print("=== Generating Sample Forecasts ===")

all_forecasts <- list()

for (i in seq_along(sample_dates)) {
  date <- sample_dates[i]
  print(paste("Processing date", i, "of", length(sample_dates), ":", date))
  
  tryCatch({
    fc <- generate_single_forecast(as.Date(date))
    all_forecasts[[as.character(date)]] <- fc
    print(paste("Generated", nrow(fc), "forecast rows"))
  }, error = function(e) {
    print(paste("Error for date", date, ":", e$message))
  })
}

# Combine all forecasts
if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  print(paste("Total forecasts generated:", nrow(combined_forecasts)))
  print(paste("Unique origin dates:", length(unique(combined_forecasts$origin_date))))
  print(paste("Unique locations:", length(unique(combined_forecasts$location))))
  
  # Save forecasts
  write.csv(combined_forecasts, "var2_log_sample_forecasts.csv", row.names = FALSE)
  
  # Quick evaluation using hubEvals
  library("hubData")
  library("hubEvals")
  
  # Get oracle output
  oracle_output <- connect_target_oracle_output(hub_path) |>
    collect()
  
  # Evaluate performance
  if (nrow(combined_forecasts) > 0) {
    model_scores <- score_model_out(combined_forecasts, oracle_output)
    
    print("\\n=== VAR(2) LOG MODEL PERFORMANCE ===")
    print(model_scores)
    
    # Save evaluation results
    write.csv(model_scores, "var2_log_wis_evaluation.csv", row.names = FALSE)
    
    # Compare to current model
    current_wis <- 0.466
    new_wis <- model_scores$wis[1]
    
    print(paste("\\nCurrent model WIS:", current_wis))
    print(paste("VAR(2) log model WIS:", round(new_wis, 3)))
    
    if (new_wis < current_wis) {
      improvement <- (current_wis - new_wis) / current_wis * 100
      print(paste("IMPROVEMENT:", round(improvement, 1), "% better!"))
      print("** This model should be used for full test phase submission **")
    } else {
      decline <- (new_wis - current_wis) / current_wis * 100
      print(paste("DECLINE:", round(decline, 1), "% worse"))
      print("Current model is still better")
    }
  }
  
} else {
  print("No forecasts generated successfully")
}