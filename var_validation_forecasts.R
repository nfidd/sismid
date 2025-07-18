# VAR Model Validation Phase Forecasts

# Load required libraries
library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("ggplot2")
library("hubUtils")
library("hubEvals")
library("hubVis")
library("hubData")
library("hubEnsembles")

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

# Define transformations
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# Get origin dates
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

# Quantile levels
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

# Filter validation dates (seasons 1 and 2)
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

print(paste("Processing", length(validation_origin_dates), "validation forecasts"))

# Function to generate VAR forecasts for one origin date
generate_var_forecast_single <- function(data, origin_date, quantile_levels, model_id) {
  
  # Filter data up to origin date
  train_data <- data |>
    filter(origin_date <= !!origin_date) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Fit VAR model
  var_model <- train_data |>
    model(
      var2 = VAR(vars(
        `HHS Region 1` = my_fourth_root(`HHS Region 1`),
        `HHS Region 2` = my_fourth_root(`HHS Region 2`),
        `HHS Region 3` = my_fourth_root(`HHS Region 3`),
        `HHS Region 4` = my_fourth_root(`HHS Region 4`),
        `HHS Region 5` = my_fourth_root(`HHS Region 5`),
        `HHS Region 6` = my_fourth_root(`HHS Region 6`),
        `HHS Region 7` = my_fourth_root(`HHS Region 7`),
        `HHS Region 8` = my_fourth_root(`HHS Region 8`),
        `HHS Region 9` = my_fourth_root(`HHS Region 9`),
        `HHS Region 10` = my_fourth_root(`HHS Region 10`),
        `US National` = my_fourth_root(`US National`)
      ) ~ AR(2))
    )
  
  # Generate forecasts with bootstrap samples
  var_forecasts <- var_model |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    pivot_longer(cols = -c(origin_date, .model, .rep, .sim_id),
                 names_to = "location",
                 values_to = ".sim") |>
    group_by(.rep, location, origin_date) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble() |>
    rename(
      target_end_date = origin_date,
      value = .sim
    ) |>
    mutate(
      origin_date = !!origin_date,
      model_id = model_id
    ) |>
    left_join(loc_df, by = "location") |>
    group_by(model_id, location, origin_date, horizon, target_end_date) |>
    reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
    mutate(
      output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
      target = "ili perc",
      output_type = "quantile"
    ) |>
    select(-quantile)
  
  return(var_forecasts)
}

# Generate all validation forecasts
all_var_forecasts <- list()

# Process in batches to avoid memory issues
for (i in seq_along(validation_origin_dates)) {
  origin_date <- validation_origin_dates[i]
  
  tryCatch({
    forecast_result <- generate_var_forecast_single(
      flu_data_hhs,
      origin_date,
      quantile_levels,
      "sismid-var2"
    )
    
    all_var_forecasts[[i]] <- forecast_result
    print(paste("Completed forecast", i, "of", length(validation_origin_dates), 
                "- Origin date:", origin_date))
    
  }, error = function(e) {
    print(paste("Error with origin date", origin_date, ":", e$message))
  })
}

# Combine all forecasts
cv_forecasts_var <- bind_rows(all_var_forecasts)

print(paste("\nGenerated", nrow(cv_forecasts_var), "forecast rows"))

# Save model metadata for VAR model
this_model_id <- "sismid-var2"

metadata_filepath <- file.path(
  hub_path,
  "model-metadata", 
  paste0(this_model_id, ".yml"))

my_text <- c("team_abbr: \"sismid\"", 
             "model_abbr: \"var2\"", 
             "designated_model: true")

writeLines(my_text, metadata_filepath)

# Write out forecast files
groups <- cv_forecasts_var |> 
  group_by(model_id, target, location, origin_date, horizon, target_end_date) |> 
  group_split()

# Save each group as a separate CSV
for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  # remove model_id from saved data
  group_df <- select(group_df, -model_id)
  
  # path to the file
  model_folder <- file.path(
    hub_path,
    "model-output", 
    this_model_id)

  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  
  results_path <- file.path(
    model_folder, 
    filename)
  
  # create directory if needed
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  write.csv(group_df, file = results_path, row.names = FALSE)
}

print("\nVAR model validation forecasts complete!")

# Now evaluate the model
new_hub_con <- connect_hub(hub_path)

validation_forecasts <- new_hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

# Score models
model_scores <- hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
) |> 
  arrange(wis)

print("\n=== Model Performance Comparison ===")
print(knitr::kable(model_scores, digits = 2))