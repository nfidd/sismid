# Test Advanced Seasonality Models for VAR
# Focus on proper seasonal specification for flu data

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("feasts")
library("ggplot2")

# Set seed for reproducibility
set.seed(406)

# Load data
data(flu_data_hhs)

# Prepare data in wide format with proper seasonal structure
flu_data_wide <- flu_data_hhs |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date) |>
  fill_gaps()

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Check seasonal patterns in the data
print("=== Seasonal Pattern Analysis ===")
seasonal_analysis <- flu_data_wide |>
  select(origin_date, `US National`) |>
  mutate(
    week = lubridate::week(origin_date),
    month = lubridate::month(origin_date),
    year = lubridate::year(origin_date)
  ) |>
  group_by(week) |>
  summarise(
    avg_wili = mean(`US National`, na.rm = TRUE),
    sd_wili = sd(`US National`, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) |>
  arrange(week)

print(paste("Weekly seasonality - weeks with highest average ILI:"))
print(seasonal_analysis |> slice_max(avg_wili, n = 5))

# Check if we have enough data for annual seasonality (need ~52 weeks)
total_weeks <- nrow(flu_data_wide)
print(paste("Total weeks of data:", total_weeks))

# Split data for testing
train_cutoff <- as.Date("2018-05-05")  # Keep more data for training
train_data <- flu_data_wide |>
  filter(origin_date <= train_cutoff)

test_data <- flu_data_wide |>
  filter(origin_date > train_cutoff)

print(paste("Training data:", nrow(train_data), "weeks"))
print(paste("Test data:", nrow(test_data), "weeks"))

# Test models with proper seasonal specification
print("\n=== Testing Models with Seasonality ===")

# Since we have weekly data, try period = 52 (weeks per year)
# But start with simpler seasonal patterns first

# Model 1: VAR(2) with sqrt (baseline)
model_var2_sqrt <- train_data |>
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
  )

# Model 2: VAR with seasonal dummies (quarterly)
model_var_season_dummy <- train_data |>
  model(
    var_season_dummy = VAR(vars(
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
    ) ~ AR(2) + season(period = 13))  # Quarterly seasonality (13 weeks)
  )

# Model 3: VAR with annual seasonality (if we have enough data)
model_var_annual <- train_data |>
  model(
    var_annual = VAR(vars(
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
    ) ~ AR(2) + season(period = 52))  # Annual seasonality
  )

# Model 4: VAR with trend
model_var_trend <- train_data |>
  model(
    var_trend = VAR(vars(
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
    ) ~ AR(2) + trend())
  )

# Extract model diagnostics
print("\n=== Model Performance Comparison ===")

models_list <- list(
  model_var2_sqrt,
  model_var_season_dummy,
  model_var_annual,
  model_var_trend
)

names(models_list) <- c("VAR(2) sqrt", "VAR(2) + season(13)", 
                        "VAR(2) + season(52)", "VAR(2) + trend")

aic_comparison <- data.frame(
  Model = character(),
  AICc = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(models_list)) {
  model_name <- names(models_list)[i]
  model_obj <- models_list[[i]]
  
  tryCatch({
    aic_val <- glance(model_obj)$AICc
    aic_comparison <- rbind(aic_comparison, 
                           data.frame(Model = model_name, 
                                    AICc = aic_val, 
                                    Status = "Success"))
  }, error = function(e) {
    print(paste("Error getting AICc for", model_name, ":", e$message))
    aic_comparison <<- rbind(aic_comparison, 
                           data.frame(Model = model_name, 
                                    AICc = NA, 
                                    Status = paste("Error:", e$message)))
  })
}

# Remove failed models and rank
aic_comparison_clean <- aic_comparison |>
  filter(Status == "Success") |>
  arrange(AICc) |>
  mutate(
    Delta_AICc = AICc - min(AICc),
    Rank = row_number()
  )

print("Full comparison (including errors):")
print(aic_comparison)

print("\nSuccessful models ranked by AICc:")
print(aic_comparison_clean)

# Test cross-validation performance
print("\n=== Cross-Validation Performance ===")

if (nrow(aic_comparison_clean) > 0) {
  best_model_name <- aic_comparison_clean$Model[1]
  best_model_idx <- which(names(models_list) == best_model_name)
  best_model <- models_list[[best_model_idx]]
  
  print(paste("Best model:", best_model_name))
  
  # Generate forecasts with bootstrap
  print("Generating bootstrap forecasts...")
  forecasts <- best_model |>
    generate(h = 4, times = 100, bootstrap = TRUE)
  
  print(paste("Generated", nrow(forecasts), "forecast samples"))
  
  # Calculate simple point forecasts for accuracy
  point_forecasts <- forecasts |>
    group_by(origin_date, .model) |>
    summarise(
      across(starts_with("HHS"), mean, .names = "{.col}"),
      across(starts_with("US"), mean, .names = "{.col}"),
      .groups = "drop"
    ) |>
    as_tsibble(index = origin_date)
  
  print(paste("Point forecasts calculated for", nrow(point_forecasts), "time periods"))
  
  # Save results
  results_summary <- list(
    seasonal_analysis = seasonal_analysis,
    model_comparison = aic_comparison,
    model_comparison_clean = aic_comparison_clean,
    best_model = best_model_name,
    forecast_samples = nrow(forecasts),
    training_weeks = nrow(train_data),
    test_weeks = nrow(test_data)
  )
  
  saveRDS(results_summary, "seasonality_test_results.rds")
  write.csv(aic_comparison, "seasonality_model_comparison.csv", row.names = FALSE)
  
  print("\n=== Results Summary ===")
  print(paste("Best model:", best_model_name))
  if (nrow(aic_comparison_clean) > 1) {
    improvement <- aic_comparison_clean$AICc[aic_comparison_clean$Model == "VAR(2) sqrt"] - 
                   min(aic_comparison_clean$AICc)
    print(paste("AICc improvement over baseline:", round(improvement, 2)))
  }
  
  print("\nFiles saved:")
  print("- seasonality_test_results.rds")
  print("- seasonality_model_comparison.csv")
  
} else {
  print("No models completed successfully!")
}