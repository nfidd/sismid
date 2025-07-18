# VAR Model with Proper Time-Series Cross-Validation
# Using fable's stretch_tsibble approach

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("stringr")
library("hubUtils")
library("hubEvals")
library("hubData")

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

# Quantile levels
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

# Get origin dates
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

print("=== VAR Model with Time-Series Cross-Validation ===")

# SEASON 1: Create time-series cross-validation dataset
# Need to pivot to wide format for VAR
flu_data_wide <- flu_data_hhs |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

# Create stretched tsibble for Season 1
flu_data_wide_tscv_season1 <- flu_data_wide |> 
  filter(origin_date <= "2016-05-07") |>  # last 2015/2016 date
  tsibble::stretch_tsibble(
    .init = 634,  # Same as ARIMA model
    .step = 1,
    .id = ".split"
  )

print(paste("Season 1 CV splits:", n_distinct(flu_data_wide_tscv_season1$.split)))

# Fit VAR models and generate forecasts for Season 1
print("\nFitting VAR models for Season 1...")
cv_var_season1 <- flu_data_wide_tscv_season1 |>
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
  ) |>
  generate(h = 4, times = 100, bootstrap = TRUE)

# Process Season 1 forecasts
# Extract location columns
special_cols <- c(".split", ".model", ".rep", "origin_date", ".sim_id")
innov_cols <- names(cv_var_season1)[grep("^\\.innov", names(cv_var_season1))]
all_special_cols <- c(special_cols, innov_cols)
location_cols <- names(cv_var_season1)[!names(cv_var_season1) %in% all_special_cols]

cv_forecasts_var_season1 <- cv_var_season1 |>
  as_tibble() |>
  select(all_of(c(".split", ".model", ".rep", "origin_date", location_cols))) |>
  pivot_longer(cols = all_of(location_cols),
               names_to = "location",
               values_to = ".sim") |>
  group_by(.split, .rep, location, .model) |>
  mutate(horizon = row_number()) |>
  ungroup() |>
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |>
  left_join(loc_df, by = "location") |>
  mutate(origin_date = target_end_date - horizon * 7L) |>
  group_by(model_id, location, origin_date, horizon, target_end_date) |>
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
  mutate(
    output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
    target = "ili perc",
    output_type = "quantile",
    model_id = "sismid-var2"
  ) |>
  select(-quantile)

print(paste("Season 1 forecasts generated:", nrow(cv_forecasts_var_season1)))

# SEASON 2: Create stretched tsibble for Season 2
first_idx_season2 <- nrow(flu_data_wide |> 
  filter(origin_date <= "2016-10-22"))

flu_data_wide_tscv_season2 <- flu_data_wide |> 
  filter(origin_date <= "2017-05-06") |>  # last 2016/2017 date
  tsibble::stretch_tsibble(
    .init = first_idx_season2,
    .step = 1,
    .id = ".split"
  )

print(paste("\nSeason 2 CV splits:", n_distinct(flu_data_wide_tscv_season2$.split)))

# Fit VAR models and generate forecasts for Season 2
print("\nFitting VAR models for Season 2...")
cv_var_season2 <- flu_data_wide_tscv_season2 |>
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
  ) |>
  generate(h = 4, times = 100, bootstrap = TRUE)

# Process Season 2 forecasts
innov_cols_s2 <- names(cv_var_season2)[grep("^\\.innov", names(cv_var_season2))]
all_special_cols_s2 <- c(special_cols, innov_cols_s2)
location_cols_s2 <- names(cv_var_season2)[!names(cv_var_season2) %in% all_special_cols_s2]

cv_forecasts_var_season2 <- cv_var_season2 |>
  as_tibble() |>
  select(all_of(c(".split", ".model", ".rep", "origin_date", location_cols_s2))) |>
  pivot_longer(cols = all_of(location_cols_s2),
               names_to = "location",
               values_to = ".sim") |>
  group_by(.split, .rep, location, .model) |>
  mutate(horizon = row_number()) |>
  ungroup() |>
  rename(
    target_end_date = origin_date,
    value = .sim,
    model_id = .model
  ) |>
  left_join(loc_df, by = "location") |>
  mutate(origin_date = target_end_date - horizon * 7L) |>
  group_by(model_id, location, origin_date, horizon, target_end_date) |>
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
  mutate(
    output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
    target = "ili perc",
    output_type = "quantile",
    model_id = "sismid-var2"
  ) |>
  select(-quantile)

print(paste("Season 2 forecasts generated:", nrow(cv_forecasts_var_season2)))

# Combine and save forecasts
all_var_forecasts <- bind_rows(cv_forecasts_var_season1, cv_forecasts_var_season2)

# Save model metadata
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
groups <- all_var_forecasts |> 
  group_by(model_id, target, location, origin_date, horizon, target_end_date) |> 
  group_split()

print(paste("\nWriting", length(groups), "forecast files..."))

for (i in seq_along(groups)) {
  group_df <- groups[[i]]
  this_model_id <- group_df$model_id[1]
  this_origin_date <- group_df$origin_date[1]
  
  group_df <- select(group_df, -model_id)
  
  model_folder <- file.path(hub_path, "model-output", this_model_id)
  filename <- paste0(this_origin_date, "-", this_model_id, ".csv")
  results_path <- file.path(model_folder, filename)
  
  if (!file.exists(model_folder)) {
    dir.create(model_folder, recursive = TRUE)
  }
  
  write.csv(group_df, file = results_path, row.names = FALSE)
}

print("\n=== Evaluating VAR model performance ===")

# Connect to hub and evaluate
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

new_hub_con <- connect_hub(hub_path)
validation_forecasts <- new_hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

model_scores <- hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
) |> 
  arrange(wis)

print("\n=== Model Performance Comparison ===")
print(knitr::kable(model_scores, digits = 2))

# Save workspace
save.image("var_cross_validation_results.RData")