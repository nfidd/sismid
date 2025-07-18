# Create VAR(2) with log transformation and trimmed outliers
# Innovation: Pre-process data to handle outliers before modeling

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")

set.seed(406)
data(flu_data_hhs)

print("=== Creating VAR(2) Log with Trimmed Outliers Model ===")
print("Innovation: Outlier trimming before log transformation")
print("Hypothesis: Reducing extreme values improves forecast stability")

# Get all dates
hub_path <- here::here("sismid-ili-forecasting-sandbox")
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter to valid dates
min_data_date <- as.Date("2015-10-24") + 20 * 7
valid_dates <- origin_dates[as.Date(origin_dates) >= min_data_date]

print(paste("Total dates to process:", length(valid_dates)))

# Location mapping
location_formal_names <- c("US National", paste("HHS Region", 1:10))

# Function to trim outliers using IQR method
trim_outliers <- function(x, k = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - k * iqr
  upper <- q3 + k * iqr
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)
}

# Generate forecasts
generate_var2_logtrim_forecast <- function(forecast_date) {
  # Prepare data
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= forecast_date) |>
    select(origin_date, location, wili) |>
    # Trim outliers by location
    group_by(location) |>
    mutate(wili_trimmed = trim_outliers(wili)) |>
    ungroup() |>
    select(-wili) |>
    pivot_wider(names_from = location, values_from = wili_trimmed) |>
    as_tsibble(index = origin_date)
  
  if (nrow(flu_data_wide) < 20) {
    return(NULL)
  }
  
  # Define log transformation with small offset to handle zeros
  log_transform <- function(x) log(x + 0.01)
  inv_log <- function(x) exp(x) - 0.01
  my_log <- new_transformation(log_transform, inv_log)
  
  # Fit VAR(2) with log transformation
  tryCatch({
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
      mutate(horizon = row_number(), .by = c(.rep, .model)) |>
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
        origin_date = forecast_date,
        target = "ili perc",
        # Ensure non-negative
        value = pmax(value, 0)
      )
    
    # Calculate quantiles
    quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
    
    forecast_quantiles <- forecast_long |>
      group_by(location, origin_date, horizon, target_end_date) |>
      reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
      mutate(
        output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
        target = "ili perc",
        output_type = "quantile"
      ) |>
      select(-quantile)
    
    return(forecast_quantiles)
    
  }, error = function(e) {
    return(NULL)
  })
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
  
  forecast <- generate_var2_logtrim_forecast(as.Date(date))
  if (!is.null(forecast)) {
    all_forecasts[[length(all_forecasts) + 1]] <- forecast
    successful <- successful + 1
  }
}

print(paste("\nGenerated forecasts for", successful, "dates"))

if (length(all_forecasts) > 0) {
  combined_forecasts <- bind_rows(all_forecasts)
  
  # Create output directory
  model_folder <- file.path(hub_path, "model-output", "sismid-logtrimclaude")
  if (!dir.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Create metadata
  metadata_path <- file.path(hub_path, "model-metadata", "sismid-logtrimclaude.yml")
  metadata_text <- c(
    'team_abbr: "sismid"',
    'model_abbr: "logtrimclaude"',
    'designated_model: false',
    'model_details: "VAR(2) with log transformation and outlier trimming developed with Claude. Pre-processes data to handle extreme values before modeling."'
  )
  writeLines(metadata_text, metadata_path)
  
  # Save forecast files
  print("\nSaving forecast files...")
  
  forecast_groups <- combined_forecasts |>
    group_by(origin_date) |>
    group_split()
  
  for (group in forecast_groups) {
    origin_date <- group$origin_date[1]
    filename <- paste0(origin_date, "-sismid-logtrimclaude.csv")
    write_csv(group, file.path(model_folder, filename))
  }
  
  print(paste("Saved", length(forecast_groups), "forecast files"))
  print("\n=== Model Ready for Submission ===")
  print("Model: sismid-logtrimclaude")
}