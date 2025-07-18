# Generate Test Phase Forecasts for Seasonal VAR Model
# Updated model: VAR(2) + season(52) with sqrt transformation

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("stringr")

# Set seed
set.seed(406)

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Load data
data(flu_data_hhs)

# Location information
locs <- c("nat", 
          "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Define transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Quantile levels
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

# Get origin dates
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Filter TEST PHASE dates (seasons 3-5)
test_origin_dates <- origin_dates[which(as.Date(origin_dates) > as.Date("2017-05-06"))]

print(paste("Generating seasonal VAR test forecasts for", length(test_origin_dates), "origin dates"))
print(paste("Date range:", min(test_origin_dates), "to", max(test_origin_dates)))

# Prepare data in wide format for VAR
flu_data_wide <- flu_data_hhs |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

# Function to generate seasonal VAR forecast for one origin date
generate_seasonal_var_forecast <- function(data, forecast_date, quantile_levels, model_id) {
  
  # Filter training data
  train_data <- data |>
    filter(origin_date <= forecast_date)
  
  # Fit seasonal VAR model with sqrt transformation
  var_model <- train_data |>
    model(
      var2_seasonal = VAR(vars(
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
      ) ~ AR(2) + season(period = 52))
    )
  
  # Generate bootstrap samples
  var_samples <- var_model |>
    generate(h = 4, times = 100, bootstrap = TRUE)
  
  # Process samples (exclude innovation columns)
  special_cols <- c(".model", ".rep", "origin_date")
  innov_cols <- names(var_samples)[grep("^\\.innov", names(var_samples))]
  all_special_cols <- c(special_cols, innov_cols)
  location_cols <- names(var_samples)[!names(var_samples) %in% all_special_cols]
  
  # Convert to long format
  var_long <- var_samples |>
    as_tibble() |>
    select(all_of(c(".model", ".rep", "origin_date", location_cols))) |>
    pivot_longer(cols = all_of(location_cols),
                 names_to = "location",
                 values_to = "value") |>
    group_by(.rep, location) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    rename(target_end_date = origin_date) |>
    mutate(origin_date = forecast_date)
  
  # Calculate quantiles
  var_quantiles <- var_long |>
    group_by(location, origin_date, horizon, target_end_date) |>
    reframe(
      quantile = quantile_levels,
      value = quantile(value, quantile_levels)
    ) |>
    mutate(
      output_type_id = quantile,
      target = "ili perc",
      output_type = "quantile",
      model_id = model_id
    ) |>
    left_join(loc_df, by = "location")
  
  return(var_quantiles)
}

# Generate all test phase forecasts
all_test_forecasts <- list()
successful <- 0

for (i in seq_along(test_origin_dates)) {
  origin_date <- test_origin_dates[i]
  
  tryCatch({
    forecast_result <- generate_seasonal_var_forecast(
      flu_data_wide,
      as.Date(origin_date),
      quantile_levels,
      "sismid-var2-seasonal"
    )
    
    all_test_forecasts[[length(all_test_forecasts) + 1]] <- forecast_result
    successful <- successful + 1
    
    if (i %% 10 == 0) {
      print(paste("Completed", i, "of", length(test_origin_dates), "forecasts"))
    }
    
  }, error = function(e) {
    print(paste("Error with origin date", origin_date, ":", e$message))
  })
}

# Combine all forecasts
test_forecasts_combined <- bind_rows(all_test_forecasts)

print(paste("\nGenerated", nrow(test_forecasts_combined), "forecast rows from", 
            successful, "successful forecasts"))

# Create model metadata
this_model_id <- "sismid-var2-seasonal"

metadata_filepath <- file.path(
  hub_path,
  "model-metadata", 
  paste0(this_model_id, ".yml"))

metadata_text <- c(
  "team_abbr: \"sismid\"",
  "model_abbr: \"var2-seasonal\"",
  "model_name: \"Seasonal Vector Autoregression with Square Root Transform\"",
  "designated_model: true",
  "methods: \"VAR(2) model with square root transformation and annual seasonality (period=52). Captures both cross-regional dependencies and seasonal patterns in flu transmission. All 11 time series (10 HHS regions + national) modeled jointly. Square root transformation stabilizes variance. Two-week lag structure (AR2) captures autocorrelation. Annual seasonality captures yearly flu patterns. Bootstrap method used for prediction intervals. Substantially improved AICc over non-seasonal baseline.\"",
  "ensemble_of_models: false",
  "ensemble_of_hub_models: false"
)

writeLines(metadata_text, metadata_filepath)
print(paste("\nModel metadata saved to:", metadata_filepath))

# Write out forecast files
if (nrow(test_forecasts_combined) > 0) {
  # Group by origin date
  groups <- test_forecasts_combined |> 
    group_by(model_id, origin_date) |> 
    group_split()
  
  print(paste("\nWriting", length(groups), "forecast files..."))
  
  # Create model output directory
  model_folder <- file.path(hub_path, "model-output", this_model_id)
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  # Save each forecast file
  for (i in seq_along(groups)) {
    group_df <- groups[[i]]
    this_origin_date <- group_df$origin_date[1]
    
    # Select required columns in correct order
    output_df <- group_df |>
      select(location, target, horizon, target_end_date, 
             output_type, output_type_id, value) |>
      arrange(location, target, horizon, output_type_id)
    
    # Create filename
    filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
    filepath <- file.path(model_folder, filename)
    
    # Write CSV
    write.csv(output_df, filepath, row.names = FALSE)
    
    if (i %% 20 == 0) {
      print(paste("Written", i, "of", length(groups), "files"))
    }
  }
  
  print(paste("\n✓ All forecast files written to:", model_folder))
  
  # List files for verification
  forecast_files <- list.files(model_folder, pattern = "\\.csv$")
  print(paste("Total files created:", length(forecast_files)))
  print("Sample files:")
  print(head(forecast_files, 5))
  
} else {
  print("ERROR: No forecasts were generated!")
}

# Create submission summary
summary_text <- paste0(
  "# Seasonal VAR Model Test Phase Forecasts\n\n",
  "**Model ID**: ", this_model_id, "\n",
  "**Model**: VAR(2) + season(52) with sqrt transformation\n",
  "**Origin Dates**: ", length(test_origin_dates), "\n",
  "**Successful Forecasts**: ", successful, "\n",
  "**Total Forecast Rows**: ", nrow(test_forecasts_combined), "\n",
  "**Date Range**: ", min(test_origin_dates), " to ", max(test_origin_dates), "\n\n",
  "## Model Improvement:\n",
  "- **Seasonal VAR AICc**: -26,426\n",
  "- **Basic VAR AICc**: -15,254\n",
  "- **Improvement**: 11,172 AICc points\n\n",
  "## Key Features:\n",
  "- Annual seasonality (period=52) captures yearly flu patterns\n",
  "- Cross-regional dependencies via VAR framework\n",
  "- Square root transformation for variance stabilization\n",
  "- Bootstrap forecasting for proper uncertainty quantification\n\n",
  "## Next Steps:\n",
  "1. Copy forecast files to hub fork\n",
  "2. Commit and push changes\n",
  "3. Create pull request\n"
)

writeLines(summary_text, "seasonal_var_test_summary.txt")
print("\nSeasonal VAR submission summary saved to: seasonal_var_test_summary.txt")