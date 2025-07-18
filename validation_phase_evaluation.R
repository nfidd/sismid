# Validation Phase Evaluation - Following Playground Approach

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("hubUtils")
library("hubData")
library("hubEvals")

set.seed(406)
data(flu_data_hhs)

# Location mapping (from playground)
locs <- c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Define transformations to test
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

fourth_root_transform <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root_transform, inv_fourth_root)

print("=== Validation Phase Evaluation ===")
print("Following playground approach with validation periods:")
print("- Season 1: 2015/2016 (up to 2016-05-07)")
print("- Season 2: 2016/2017 (up to 2017-05-06)")
print("- Test phase: After 2017-05-06")

# Create validation forecasts for Season 1
print("\n=== Season 1 Validation Forecasts ===")

# Season 1 time-series cross-validation
flu_data_hhs_tscv_season1 <- flu_data_hhs |> 
  filter(origin_date <= "2016-05-07") |> 
  tsibble::stretch_tsibble(
    .init = 634, # From playground
    .step = 1,
    .id = ".split"
  )

print(paste("Season 1 splits:", max(flu_data_hhs_tscv_season1$.split)))

# Test VAR(2) with sqrt transformation
print("Testing VAR(2) sqrt transformation...")

# Convert to wide format for VAR
flu_wide_tscv_s1 <- flu_data_hhs_tscv_season1 |>
  select(origin_date, location, wili, .split) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date, key = .split)

# Generate VAR(2) sqrt forecasts
var2_sqrt_forecasts_s1 <- flu_wide_tscv_s1 |>
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
  ) |>
  generate(h = 4, times = 100, bootstrap = TRUE) |>
  group_by(.split, .rep, .model) |>
  mutate(horizon = row_number()) |>
  ungroup() |>
  as_tibble()

print(paste("VAR(2) sqrt samples generated:", nrow(var2_sqrt_forecasts_s1)))

# Convert to hub format
location_cols <- c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                   "HHS Region 4", "HHS Region 5", "HHS Region 6",
                   "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                   "HHS Region 10", "US National")

var2_sqrt_long_s1 <- var2_sqrt_forecasts_s1 |>
  pivot_longer(
    cols = all_of(location_cols),
    names_to = "location",
    values_to = "value"
  ) |>
  rename(
    target_end_date = origin_date,
    model_id = .model
  ) |>
  left_join(loc_df) |>
  mutate(
    origin_date = target_end_date - horizon * 7L,
    target = "ili perc",
    model_id = "sismid-var2-sqrt"
  )

# Calculate quantiles
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

var2_sqrt_quantiles_s1 <- var2_sqrt_long_s1 |>
  group_by(model_id, location, origin_date, horizon, target_end_date) |>
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
  mutate(
    output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
    target = "ili perc",
    output_type = "quantile"
  ) |>
  select(-quantile)

print(paste("VAR(2) sqrt quantiles:", nrow(var2_sqrt_quantiles_s1)))

# Now test VAR(2) with fourth-root transformation
print("Testing VAR(2) fourth-root transformation...")

var2_fourth_forecasts_s1 <- flu_wide_tscv_s1 |>
  model(
    var2_fourth = VAR(vars(
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
  generate(h = 4, times = 100, bootstrap = TRUE) |>
  group_by(.split, .rep, .model) |>
  mutate(horizon = row_number()) |>
  ungroup() |>
  as_tibble()

print(paste("VAR(2) fourth-root samples generated:", nrow(var2_fourth_forecasts_s1)))

var2_fourth_long_s1 <- var2_fourth_forecasts_s1 |>
  pivot_longer(
    cols = all_of(location_cols),
    names_to = "location",
    values_to = "value"
  ) |>
  rename(
    target_end_date = origin_date,
    model_id = .model
  ) |>
  left_join(loc_df) |>
  mutate(
    origin_date = target_end_date - horizon * 7L,
    target = "ili perc",
    model_id = "sismid-var2-fourth"
  )

var2_fourth_quantiles_s1 <- var2_fourth_long_s1 |>
  group_by(model_id, location, origin_date, horizon, target_end_date) |>
  reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
  mutate(
    output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
    target = "ili perc",
    output_type = "quantile"
  ) |>
  select(-quantile)

print(paste("VAR(2) fourth-root quantiles:", nrow(var2_fourth_quantiles_s1)))

# Combine both models for evaluation
all_validation_forecasts_s1 <- bind_rows(var2_sqrt_quantiles_s1, var2_fourth_quantiles_s1)

print(paste("Total validation forecasts Season 1:", nrow(all_validation_forecasts_s1)))

# Save validation forecasts
write.csv(all_validation_forecasts_s1, "validation_forecasts_season1.csv", row.names = FALSE)

# Evaluate against oracle
print("\n=== Evaluating Season 1 Performance ===")

hub_path <- here::here("sismid-ili-forecasting-sandbox")
oracle_output <- connect_target_oracle_output(hub_path) |>
  collect()

# Filter validation forecasts for proper evaluation
validation_origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

season1_dates <- validation_origin_dates[which(as.Date(validation_origin_dates) <= as.Date("2016-05-07"))]

validation_forecasts_filtered <- all_validation_forecasts_s1 |>
  filter(origin_date %in% season1_dates)

print(paste("Filtered validation forecasts:", nrow(validation_forecasts_filtered)))

# Score models
if (nrow(validation_forecasts_filtered) > 0) {
  model_scores <- score_model_out(validation_forecasts_filtered, oracle_output)
  
  print("=== VALIDATION PHASE RESULTS ===")
  print(model_scores)
  
  # Save results
  write.csv(model_scores, "validation_phase_results.csv", row.names = FALSE)
  
  # Compare models
  best_model <- model_scores |> slice_min(wis, n = 1)
  print(paste("Best model in validation:", best_model$model_id, "with WIS:", round(best_model$wis, 3)))
  
  # Compare to current terrible model
  current_wis <- 0.466
  best_wis <- best_model$wis
  
  if (best_wis < current_wis) {
    improvement <- (current_wis - best_wis) / current_wis * 100
    print(paste("IMPROVEMENT over current model:", round(improvement, 1), "%"))
  }
  
} else {
  print("No validation forecasts to evaluate")
}

print("\n=== Summary ===")
print("Using proper playground validation approach:")
print("- VAR(2) sqrt vs VAR(2) fourth-root tested")
print("- Cross-validation on Season 1 (2015/2016)")
print("- Best model selected for test phase")
print("- Following Nick Reich's guidance on model selection")