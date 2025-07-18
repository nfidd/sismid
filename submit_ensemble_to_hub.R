# Submit Ensemble Model to Hub
# Following the same process as VAR model submission

library("dplyr")
library("readr")

# Read our ensemble submission
ensemble_df <- read_csv("ensemble_submission_final.csv")

print("=== Submitting Ensemble to Hub ===")
print(paste("Model ID:", unique(ensemble_df$model_id)))
print(paste("Total forecasts:", nrow(ensemble_df)))

# Hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Group by model_id and origin_date for saving
groups <- ensemble_df |>
  group_split(model_id, origin_date)

print(paste("\nNumber of forecast files to create:", length(groups)))

# Save each group as a separate CSV
for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  # Remove model_id from saved data, as it is implied from filepath
  group_df <- select(group_df, -model_id)
  
  # Path to the model output folder
  model_folder <- file.path(
    hub_path,
    "model-output", 
    this_model_id
  )
  
  # Filename format: YYYY-MM-DD-model_id.csv
  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  
  # Full path to the file
  results_path <- file.path(
    model_folder, 
    filename
  )
  
  # Create model directory if it doesn't exist
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
    print(paste("Created directory:", model_folder))
  }
  
  # Write the forecast file
  write.csv(group_df, file = results_path, row.names = FALSE)
  
  if (i == 1) {
    print(paste("\nFirst file saved:", results_path))
    print("Sample content:")
    print(head(group_df, 5))
  }
  
  if (i %% 10 == 0) {
    print(paste("Progress:", i, "of", length(groups), "files saved"))
  }
}

print(paste("\n=== Submission Complete ==="))
print(paste("All files saved to:", file.path(hub_path, "model-output", unique(ensemble_df$model_id))))
print("\nNext steps:")
print("1. Navigate to the hub directory")
print("2. git add model-output/sismid-ensemble-v2/")
print("3. git commit -m 'Add ensemble model submission'")
print("4. git push (if you have permissions)")

# List the created files
model_dir <- file.path(hub_path, "model-output", unique(ensemble_df$model_id))
if (dir.exists(model_dir)) {
  files <- list.files(model_dir)
  print(paste("\nCreated", length(files), "forecast files"))
  print("First few files:")
  print(head(files, 5))
}