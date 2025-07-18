# Generate ALL VAR(2) sqrt forecasts - complete data period
# Based on the working test phase script

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")

set.seed(406)
data(flu_data_hhs)

# Location mapping
locs <- c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Define sqrt transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

print("=== Generating ALL VAR(2) sqrt Forecasts ===")

# Get ALL dates from hub config
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

print(paste("Total dates available:", length(origin_dates)))
print(paste("Date range:", min(origin_dates), "to", max(origin_dates)))

# Filter to dates we have enough data for (need at least 20 weeks)
min_data_date <- as.Date("2015-10-24") + 20 * 7  # First date + 20 weeks
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("\nDates with sufficient data:", length(valid_dates)))
print(paste("First valid date:", min(valid_dates)))
print(paste("Last valid date:", max(valid_dates)))

# Function to generate forecast for a single date
generate_var2_sqrt_forecast <- function(forecast_date) {
  
  # Prepare data up to origin date
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Check if we have enough data
  if (nrow(flu_data_wide) < 20) {
    return(NULL)
  }
  
  # Fit VAR(2) sqrt model
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
  
  # Generate forecasts
  forecast_samples <- model_fit |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  # Convert to hub format
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
      origin_date = forecast_date,  # Use the forecast date as origin
      target = "ili perc",
      model_id = "sismid-var2sqrt"
    )
  
  # Calculate quantiles
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- forecast_long |>
    group_by(model_id, location, origin_date, horizon, target_end_date) |>
    reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
    mutate(
      output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
      target = "ili perc",
      output_type = "quantile"
    ) |>
    select(-quantile)
  
  return(forecast_quantiles)
}

# Generate forecasts for all valid dates
print("\n=== Processing All Dates ===")

all_forecasts <- list()
successful_dates <- 0

for (i in seq_along(valid_dates)) {
  date <- valid_dates[i]
  
  if (i %% 10 == 0) {
    print(paste("Processing", i, "of", length(valid_dates), "dates..."))
  }
  
  tryCatch({
    forecast <- generate_var2_sqrt_forecast(as.Date(date))
    if (!is.null(forecast)) {
      all_forecasts[[as.character(date)]] <- forecast
      successful_dates <- successful_dates + 1
    }
  }, error = function(e) {
    print(paste("Error for date", date, ":", e$message))
  })
}

print(paste("\nSuccessfully generated forecasts for", successful_dates, "dates"))

# Combine all forecasts and save
if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  print(paste("Total forecast rows:", nrow(combined_forecasts)))
  print(paste("Unique origin dates:", length(unique(combined_forecasts$origin_date))))
  
  # Create model output directory
  model_folder <- file.path(hub_path, "model-output", "sismid-var2sqrt")
  
  if (!dir.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Save individual forecast files
  print("\nSaving individual forecast files...")
  
  forecast_groups <- combined_forecasts |>
    group_by(origin_date) |>
    group_split()
  
  saved_files <- 0
  for (group in forecast_groups) {
    origin_date <- group$origin_date[1]
    filename <- paste0(origin_date, "-sismid-var2sqrt.csv")
    
    # Remove model_id column for hub format
    group_clean <- group |>
      select(-model_id)
    
    write.csv(group_clean, file.path(model_folder, filename), row.names = FALSE)
    saved_files <- saved_files + 1
  }
  
  print(paste("Saved", saved_files, "forecast files to:", model_folder))
  
  # List date ranges
  dates_created <- sort(unique(combined_forecasts$origin_date))
  print(paste("\nForecast date range:", min(dates_created), "to", max(dates_created)))
  
} else {
  print("No forecasts generated successfully")
}