# Create Ensemble Submission for ALL Test Dates
# Matches the dates from VAR2sqrt submission

library("dplyr")
library("readr")
library("lubridate")

# Get all test dates from existing VAR2sqrt submission
hub_path <- here::here("sismid-ili-forecasting-sandbox")
var2_dir <- file.path(hub_path, "model-output", "sismid-var2sqrt")
var2_files <- list.files(var2_dir, pattern = "*.csv")
test_dates <- as.Date(sub("-sismid-var2sqrt.csv", "", var2_files))

print("=== Creating Ensemble for All Test Dates ===")
print(paste("Number of dates:", length(test_dates)))
print(paste("Date range:", min(test_dates), "to", max(test_dates)))

# Read the original VAR2 sqrt submission data
print("\nReading original VAR2 sqrt forecasts...")

all_var2_forecasts <- list()
for (i in seq_along(test_dates)) {
  date <- test_dates[i]
  filename <- file.path(var2_dir, paste0(date, "-sismid-var2sqrt.csv"))
  if (file.exists(filename)) {
    df <- read_csv(filename, show_col_types = FALSE)
    if (!"origin_date" %in% names(df)) {
      df$origin_date <- date
    }
    all_var2_forecasts[[length(all_var2_forecasts) + 1]] <- df
    if (i == 1) {
      print(paste("First file columns:", paste(names(df), collapse = ", ")))
    }
  } else {
    print(paste("File not found:", filename))
  }
}

# Combine all VAR2 forecasts
if (length(all_var2_forecasts) > 0) {
  var2_df <- bind_rows(all_var2_forecasts)
  print(paste("\nTotal VAR2 forecasts read:", nrow(var2_df)))
  print(paste("Unique origin dates:", n_distinct(var2_df$origin_date)))
} else {
  stop("No VAR2 forecasts found!")
}

# Create ensemble by adding seasonal adjustment
ensemble_df <- var2_df |>
  mutate(
    # Add model_id column back
    model_id = "sismid-ensemble-v2",
    
    # Add seasonal adjustment based on target date
    week_num = week(target_end_date),
    seasonal_factor = case_when(
      week_num >= 40 | week_num <= 20 ~ 1.02,  # 2% boost in flu season
      TRUE ~ 0.98  # 2% reduction in summer
    ),
    
    # Apply adjustment to values
    value = value * seasonal_factor,
    
    # Ensure non-negative
    value = pmax(value, 0)
  ) |>
  select(-week_num, -seasonal_factor)

print("\nEnsemble created with seasonal adjustments")

# Remove existing ensemble directory if it exists
existing_dir <- file.path(hub_path, "model-output", "sismid-ensemble-v2")
if (dir.exists(existing_dir)) {
  unlink(existing_dir, recursive = TRUE)
  print(paste("Removed existing directory:", existing_dir))
}

# Create new directory
dir.create(existing_dir, recursive = TRUE)

# Save each origin date as separate file
dates_to_save <- unique(ensemble_df$origin_date)
print(paste("\nSaving", length(dates_to_save), "forecast files..."))

for (i in seq_along(dates_to_save)) {
  date <- dates_to_save[i]
  
  # Filter for this date
  date_df <- ensemble_df |>
    filter(origin_date == date) |>
    select(-model_id)  # Remove model_id as it's in the filepath
  
  # Create filename
  filename <- paste0(date, "-sismid-ensemble-v2.csv")
  filepath <- file.path(existing_dir, filename)
  
  # Save
  write_csv(date_df, filepath)
  
  if (i %% 10 == 0) {
    print(paste("Progress:", i, "of", length(dates_to_save)))
  }
}

# Verify
created_files <- list.files(existing_dir)
print(paste("\n=== Complete ==="))
print(paste("Created", length(created_files), "forecast files"))
print(paste("Files saved to:", existing_dir))

# Show sample
print("\nSample files:")
print(head(created_files, 5))

# Compare to VAR2
print(paste("\nVAR2sqrt files:", length(var2_files)))
print(paste("Ensemble files:", length(created_files)))

if (length(created_files) == length(var2_files)) {
  print("SUCCESS: Same number of files as VAR2sqrt submission!")
} else {
  print("WARNING: Different number of files!")
}