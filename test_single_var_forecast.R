# Test a single VAR forecast with detailed debugging

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("stringr")

# Load data
data(flu_data_hhs)

# Location information
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(
  region = c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
             "hhs6", "hhs7", "hhs8", "hhs9", "hhs10"),
  location = location_formal_names
)

# Define transformations
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# Test with a specific date
forecast_date <- as.Date("2016-01-02")

print(paste("Testing forecast for origin date:", forecast_date))

# Filter data
train_data <- flu_data_hhs |>
  filter(origin_date <= forecast_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

print(paste("Training data rows:", nrow(train_data)))
print(paste("Date range:", min(train_data$origin_date), "to", max(train_data$origin_date)))

# Fit VAR model
print("\nFitting VAR model...")
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

print("VAR model fitted!")

# Generate samples
print("\nGenerating bootstrap samples...")
var_samples <- var_model |>
  generate(h = 4, times = 10, bootstrap = TRUE)

print("Samples generated!")
print("Sample structure:")
print(str(var_samples))

# Check column names
print("\nColumn names in generated samples:")
print(names(var_samples))

# Find location columns - exclude special columns
special_cols <- c(".model", ".rep", "origin_date", ".sim_id")
# Check for .innov columns (they might have multiple columns)
innov_cols <- names(var_samples)[grep("^\\.innov", names(var_samples))]
all_special_cols <- c(special_cols, innov_cols)

location_cols <- names(var_samples)[!names(var_samples) %in% all_special_cols]
print(paste("\nLocation columns:", paste(location_cols, collapse = ", ")))

# Process samples
print("\nProcessing samples to long format...")
var_long <- var_samples |>
  as_tibble() |>  # Convert from tsibble to regular tibble
  select(all_of(c(".model", ".rep", "origin_date", location_cols))) |>
  pivot_longer(cols = all_of(location_cols),
               names_to = "location",
               values_to = "value") |>
  group_by(.rep, location) |>
  mutate(horizon = row_number()) |>
  ungroup() |>
  rename(target_end_date = origin_date) |>
  mutate(origin_date = forecast_date)

print("Long format data:")
print(var_long)

# Calculate quantiles
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

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
    model_id = "sismid-var2"
  ) |>
  left_join(loc_df, by = "location")

print("\nFinal quantile forecasts:")
print(var_quantiles)
print(paste("Total forecast rows:", nrow(var_quantiles)))