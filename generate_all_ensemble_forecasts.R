# Generate ALL Ensemble forecasts - complete data period
# Based on VAR2sqrt with seasonal adjustment

library("dplyr")
library("readr")
library("lubridate")

# Set hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Read all VAR2sqrt forecasts
var2_dir <- file.path(hub_path, "model-output", "sismid-var2sqrt")
var2_files <- list.files(var2_dir, pattern = "*.csv", full.names = TRUE)

print(paste("=== Creating Ensemble for ALL dates ==="))
print(paste("Found", length(var2_files), "VAR2sqrt files"))

# Read all VAR2sqrt forecasts
all_var2_forecasts <- list()
for (i in seq_along(var2_files)) {
  file <- var2_files[i]
  df <- read_csv(file, show_col_types = FALSE)
  all_var2_forecasts[[i]] <- df
}

# Combine all VAR2 forecasts
var2_df <- bind_rows(all_var2_forecasts)
print(paste("Total VAR2 forecasts read:", nrow(var2_df)))
print(paste("Unique origin dates:", n_distinct(var2_df$origin_date)))

# Create ensemble by adding seasonal adjustment
ensemble_df <- var2_df |>
  mutate(
    # Add seasonal adjustment based on target date
    week_num = week(target_end_date),
    seasonal_factor = case_when(
      week_num >= 40 | week_num <= 20 ~ 1.02,  # 2% boost in flu season (Oct-Apr)
      TRUE ~ 0.98  # 2% reduction in summer (May-Sep)
    ),
    
    # Apply adjustment to values
    value = value * seasonal_factor,
    
    # Ensure non-negative
    value = pmax(value, 0)
  ) |>
  select(-week_num, -seasonal_factor)

print("\nEnsemble created with seasonal adjustments")

# Create output directory
ensemble_dir <- file.path(hub_path, "model-output", "sismid-ensemblev2")
if (!dir.exists(ensemble_dir)) {
  dir.create(ensemble_dir, recursive = TRUE)
}

# Save each origin date as separate file
dates_to_save <- unique(ensemble_df$origin_date)
print(paste("\nSaving", length(dates_to_save), "forecast files..."))

for (i in seq_along(dates_to_save)) {
  date <- dates_to_save[i]
  
  # Filter for this date
  date_df <- ensemble_df |>
    filter(origin_date == date)
  
  # Create filename
  filename <- paste0(date, "-sismid-ensemblev2.csv")
  filepath <- file.path(ensemble_dir, filename)
  
  # Save
  write_csv(date_df, filepath)
  
  if (i %% 20 == 0) {
    print(paste("Progress:", i, "of", length(dates_to_save)))
  }
}

# Verify
created_files <- list.files(ensemble_dir)
print(paste("\n=== Complete ==="))
print(paste("Created", length(created_files), "ensemble forecast files"))
print(paste("Files saved to:", ensemble_dir))

# Compare to VAR2
print(paste("\nVAR2sqrt files:", length(var2_files)))
print(paste("Ensemble files:", length(created_files)))

if (length(created_files) == length(var2_files)) {
  print("SUCCESS: Same number of files as VAR2sqrt!")
} else {
  print("WARNING: Different number of files!")
}