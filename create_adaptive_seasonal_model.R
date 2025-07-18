# Create Adaptive Seasonal VAR Model
# Based on our best VAR(2) sqrt but with adaptive seasonal adjustment
# Instead of fixed 2% boost, adjust based on recent flu activity

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("lubridate")
library("hubUtils")

set.seed(406)
data(flu_data_hhs)

# Define sqrt transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

print("=== Creating Adaptive Seasonal VAR Model ===")
print("Innovation: Dynamic seasonal adjustment based on recent flu activity")
print("Building on proven VAR(2) sqrt foundation")

# Get all dates from hub
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter to dates with sufficient data
min_data_date <- as.Date("2015-10-24") + 20 * 7
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("Generating forecasts for", length(valid_dates), "dates"))

# Location mapping
locs <- c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Function to calculate adaptive seasonal factor
calculate_adaptive_factor <- function(recent_data, target_week) {
  # Calculate recent flu activity level (last 4 weeks)
  recent_avg <- recent_data |>
    tail(4) |>
    summarise(avg_wili = mean(wili, na.rm = TRUE)) |>
    pull(avg_wili)
  
  # Calculate historical average for same period
  historical_avg <- recent_data |>
    mutate(week_num = week(origin_date)) |>
    filter(week_num %in% (target_week + (-2:2))) |>
    summarise(hist_avg = mean(wili, na.rm = TRUE)) |>
    pull(hist_avg)
  
  # Base seasonal factor
  base_factor <- if (target_week >= 40 | target_week <= 20) {
    1.02  # Flu season
  } else {
    0.98  # Summer
  }
  
  # Adaptive adjustment based on recent vs historical
  if (length(recent_avg) > 0 && length(historical_avg) > 0 && historical_avg > 0) {
    activity_ratio <- recent_avg / historical_avg
    # Moderate the adjustment (cap between 0.95 and 1.05)
    activity_adjustment <- pmin(pmax(activity_ratio, 0.95), 1.05)
    return(base_factor * activity_adjustment)
  } else {
    return(base_factor)
  }
}

# Generate forecasts function
generate_adaptive_var_forecast <- function(forecast_date) {
  # Prepare data
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Check for sufficient data
  if (nrow(flu_data_wide) < 20) {
    return(NULL)
  }
  
  # Fit VAR(2) sqrt model (our proven best)
  model_fit <- flu_data_wide |>
    model(
      var2_sqrt = VAR(vars(
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
  
  # Generate base forecasts
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
    left_join(loc_df, by = c("location" = "location")) |>
    mutate(
      origin_date = forecast_date,
      target = "ili perc"
    )
  
  # Apply adaptive seasonal adjustment
  adjusted_forecasts <- forecast_long |>
    group_by(location, target_end_date) |>
    mutate(
      target_week = week(target_end_date),
      # Calculate location-specific adaptive factor
      adaptive_factor = {
        recent_data <- flu_data_hhs |>
          filter(location == first(location), origin_date <= forecast_date)
        calculate_adaptive_factor(recent_data, first(target_week))
      },
      # Apply adjustment
      value = value * adaptive_factor,
      value = pmax(value, 0)  # Ensure non-negative
    ) |>
    ungroup() |>
    select(-target_week, -adaptive_factor)
  
  # Calculate quantiles
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- adjusted_forecasts |>
    group_by(location, origin_date, horizon, target_end_date) |>
    reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
    mutate(
      output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
      target = "ili perc",
      output_type = "quantile"
    ) |>
    select(-quantile)
  
  return(forecast_quantiles)
}

# Generate all forecasts
print("\n=== Generating Forecasts ===")
all_forecasts <- list()
successful <- 0

for (i in seq_along(valid_dates)) {
  date <- valid_dates[i]
  
  if (i %% 20 == 0) {
    print(paste("Progress:", i, "of", length(valid_dates)))
  }
  
  tryCatch({
    forecast <- generate_adaptive_var_forecast(as.Date(date))
    if (!is.null(forecast)) {
      all_forecasts[[length(all_forecasts) + 1]] <- forecast
      successful <- successful + 1
    }
  }, error = function(e) {
    # Silent skip
  })
}

print(paste("\nGenerated forecasts for", successful, "dates"))

if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  # Create output directory
  model_folder <- file.path(hub_path, "model-output", "sismid-adaptiveclaude")
  if (!dir.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Create metadata
  metadata_path <- file.path(hub_path, "model-metadata", "sismid-adaptiveclaude.yml")
  metadata_text <- c(
    'team_abbr: "sismid"',
    'model_abbr: "adaptiveclaude"',
    'designated_model: false',
    'model_details: "VAR(2) sqrt with adaptive seasonal adjustment developed with Claude. Dynamically adjusts seasonal factors based on recent flu activity vs historical patterns."'
  )
  writeLines(metadata_text, metadata_path)
  
  # Save individual forecast files
  print("\nSaving forecast files...")
  
  forecast_groups <- combined_forecasts |>
    group_by(origin_date) |>
    group_split()
  
  for (group in forecast_groups) {
    origin_date <- group$origin_date[1]
    filename <- paste0(origin_date, "-sismid-adaptiveclaude.csv")
    write_csv(group, file.path(model_folder, filename))
  }
  
  print(paste("Saved", length(forecast_groups), "forecast files"))
  print(paste("Model directory:", model_folder))
  
  print("\n=== Model Ready for Submission ===")
  print("Model: sismid-adaptiveclaude")
  print("Innovation: Adaptive seasonal adjustment based on recent flu activity")
  print("Foundation: Proven VAR(2) sqrt model")
}