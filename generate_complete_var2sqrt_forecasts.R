# Generate VAR2sqrt forecasts for complete data period
# This covers the validation phase and earlier periods (2015-10-24 to 2017-10-14)

library("dplyr")
library("readr")
library("lubridate")
library("tsibble")
library("fable")
library("tidyr")
library("here")

# Set hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Load flu data
flu_data <- read_csv(file.path(hub_path, "target-data", "time-series.csv"))

# Get all missing dates (before our current submissions)
current_var2_files <- list.files(
  file.path(hub_path, "model-output", "sismid-var2sqrt"), 
  pattern = "*.csv"
)
current_dates <- as.Date(sub("-sismid-var2sqrt.csv", "", current_var2_files))
min_current_date <- min(current_dates)

# Get all forecast dates from other models (e.g., delphi)
all_dates <- list.files(
  file.path(hub_path, "model-output", "delphi-epicast"), 
  pattern = "*.csv"
)
all_dates <- as.Date(sub("-delphi-epicast.csv", "", all_dates))

# Find missing dates
missing_dates <- all_dates[all_dates < min_current_date]

# Filter to only dates where we have enough historical data
# Need at least 20 weeks of data, data starts from 2015-10-24
min_forecast_date <- as.Date("2015-10-24") + weeks(20)
missing_dates <- missing_dates[missing_dates >= min_forecast_date]

print(paste("Number of missing dates:", length(missing_dates)))
print(paste("Date range:", min(missing_dates), "to", max(missing_dates)))

# Function to create VAR2sqrt forecasts for a given origin date
create_var2sqrt_forecast <- function(origin_date, flu_data) {
  # Filter data up to origin date
  train_data <- flu_data |>
    filter(target_end_date < origin_date) |>
    mutate(
      week_ending = as.Date(target_end_date),
      ili = sqrt(observation)
    ) |>
    as_tsibble(index = week_ending, key = location) |>
    filter(!is.na(ili))
  
  # Check if we have enough data
  if (nrow(train_data) == 0) {
    return(NULL)
  }
  
  min_weeks_by_location <- train_data |>
    group_by(location) |>
    summarise(n_weeks = n()) |>
    pull(n_weeks)
  
  if (length(min_weeks_by_location) == 0 || min(min_weeks_by_location) < 20) {
    # Need at least 20 weeks of data for VAR(2) model
    return(NULL)
  }
  
  # Also check total observations
  if (nrow(train_data) < 220) {  # 11 locations * 20 weeks minimum
    return(NULL)
  }
  
  # Fit VAR(2) model
  model_fit <- train_data |>
    model(var_model = VAR(ili ~ AR(2)))
  
  # Generate forecasts
  forecasts <- model_fit |>
    forecast(h = 4) |>
    mutate(
      # Transform back from sqrt
      .mean = .mean^2,
      .distribution = distributional::dist_transformed(
        .distribution, 
        transform = list(forward = function(x) x^2, inverse = sqrt),
        inverse_dist = .distribution
      )
    )
  
  # Extract quantiles
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  
  forecast_df <- forecasts |>
    hilo(quantiles * 100) |>
    unpack_hilo(names(quantiles)) |>
    select(location, week_ending, .mean, all_of(names(quantiles))) |>
    pivot_longer(
      cols = c(.mean, all_of(names(quantiles))),
      names_to = "output_type_id",
      values_to = "value"
    ) |>
    mutate(
      origin_date = origin_date,
      target = "wk inc flu hosp",
      horizon = as.numeric(difftime(week_ending, origin_date, units = "weeks")),
      target_end_date = week_ending,
      output_type = case_when(
        output_type_id == ".mean" ~ "mean",
        TRUE ~ "quantile"
      ),
      output_type_id = case_when(
        output_type_id == ".mean" ~ NA_character_,
        TRUE ~ as.character(readr::parse_number(output_type_id) / 100)
      )
    ) |>
    select(
      origin_date, target, horizon, target_end_date, 
      location, output_type, output_type_id, value
    )
  
  return(forecast_df)
}

# Generate forecasts for all missing dates
all_forecasts <- list()
for (i in seq_along(missing_dates)) {
  date <- missing_dates[i]
  
  if (i %% 10 == 1) {
    print(paste("Processing date", i, "of", length(missing_dates), ":", date))
  }
  
  tryCatch({
    forecast <- create_var2sqrt_forecast(date, flu_data)
    if (!is.null(forecast)) {
      all_forecasts[[length(all_forecasts) + 1]] <- forecast
    }
  }, error = function(e) {
    print(paste("Error for date", date, ":", e$message))
  })
}

# Combine and save forecasts
if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  # Save each origin date as separate file
  for (date in unique(combined_forecasts$origin_date)) {
    date_df <- combined_forecasts |>
      filter(origin_date == date) |>
      select(-origin_date)  # Remove origin_date as it's in filename
    
    filename <- paste0(date, "-sismid-var2sqrt.csv")
    filepath <- file.path(hub_path, "model-output", "sismid-var2sqrt", filename)
    
    write_csv(date_df, filepath)
  }
  
  print(paste("Created", length(unique(combined_forecasts$origin_date)), 
              "forecast files"))
} else {
  print("No forecasts were generated")
}