# Test Advanced VARIMA Models with Seasonality
# Exploring more complex models now that we have multiple seasons of data

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

# Check data range to understand seasonal patterns
flu_data_hhs |>
  summarise(
    start_date = min(origin_date),
    end_date = max(origin_date),
    n_weeks = n_distinct(origin_date),
    n_locations = n_distinct(location)
  ) |>
  print()

# Convert to years for seasonal analysis
flu_data_with_season <- flu_data_hhs |>
  mutate(
    year = lubridate::year(origin_date),
    week = lubridate::week(origin_date),
    season = case_when(
      (lubridate::month(origin_date) >= 7) ~ paste0(year, "/", year + 1),
      TRUE ~ paste0(year - 1, "/", year)
    )
  )

# Count how many complete seasons we have
season_summary <- flu_data_with_season |>
  group_by(season) |>
  summarise(
    n_weeks = n_distinct(origin_date),
    start_date = min(origin_date),
    end_date = max(origin_date),
    .groups = "drop"
  ) |>
  arrange(season)

print("Season summary:")
print(season_summary)

# Prepare data in wide format
flu_data_wide <- flu_data_hhs |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date) |>
  fill_gaps()

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

log_transform <- function(x) log(x + 0.01)
inv_log <- function(x) exp(x) - 0.01
my_log <- new_transformation(log_transform, inv_log)

# Split data for cross-validation (using only first 5 seasons for training)
train_cutoff <- as.Date("2019-05-04")  # End of season 4
train_data <- flu_data_wide |>
  filter(origin_date <= train_cutoff)

test_data <- flu_data_wide |>
  filter(origin_date > train_cutoff)

print(paste("Training data:", nrow(train_data), "weeks"))
print(paste("Test data:", nrow(test_data), "weeks"))

# Test different VARIMA configurations
print("\n=== Testing VARIMA Models ===")

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

# Model 2: VAR with Fourier seasonality (K=2)
model_var_fourier2 <- train_data |>
  model(
    var_fourier2 = VAR(vars(
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
    ) ~ AR(2) + fourier(K = 2))
  )

# Model 3: VAR with Fourier seasonality (K=3)
model_var_fourier3 <- train_data |>
  model(
    var_fourier3 = VAR(vars(
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
    ) ~ AR(2) + fourier(K = 3))
  )

# Model 4: VAR with trend and seasonality
model_var_trend_season <- train_data |>
  model(
    var_trend_season = VAR(vars(
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
    ) ~ AR(2) + trend() + fourier(K = 2))
  )

# Extract model diagnostics
print("\n=== Model Performance Comparison ===")

# Compare AICc values
models_list <- list(
  model_var2_sqrt,
  model_var_fourier2,
  model_var_fourier3,
  model_var_trend_season
)

names(models_list) <- c("VAR(2) sqrt", "VAR(2) + fourier(K=2)", 
                        "VAR(2) + fourier(K=3)", "VAR(2) + trend + fourier(K=2)")

aic_comparison <- data.frame(
  Model = character(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(models_list)) {
  model_name <- names(models_list)[i]
  model_obj <- models_list[[i]]
  
  tryCatch({
    aic_val <- glance(model_obj)$AICc
    aic_comparison <- rbind(aic_comparison, 
                           data.frame(Model = model_name, AICc = aic_val))
  }, error = function(e) {
    print(paste("Error getting AICc for", model_name, ":", e$message))
  })
}

aic_comparison <- aic_comparison |>
  arrange(AICc) |>
  mutate(
    Delta_AICc = AICc - min(AICc),
    Rank = row_number()
  )

print(aic_comparison)

# Generate forecasts for best model
best_model <- models_list[[which.min(aic_comparison$AICc)]]
best_model_name <- aic_comparison$Model[which.min(aic_comparison$AICc)]

print(paste("\nBest model:", best_model_name))

# Generate forecasts
forecast_horizon <- 4
n_samples <- 100

print("\nGenerating forecasts...")
forecasts <- best_model |>
  generate(h = forecast_horizon, times = n_samples, bootstrap = TRUE)

print(paste("Generated", nrow(forecasts), "forecast samples"))

# Calculate forecast accuracy on test data
test_forecasts <- best_model |>
  forecast(h = min(nrow(test_data), forecast_horizon))

# Simple accuracy check
accuracy_results <- test_forecasts |>
  accuracy(test_data)

print("\n=== Forecast Accuracy ===")
print(accuracy_results)

# Create summary
results_summary <- list(
  model_comparison = aic_comparison,
  best_model = best_model_name,
  forecast_samples = nrow(forecasts),
  accuracy = accuracy_results
)

print("\n=== Summary ===")
print(paste("Best performing model:", best_model_name))
print(paste("AICc improvement over baseline:", 
            round(aic_comparison$AICc[aic_comparison$Model == "VAR(2) sqrt"] - 
                  min(aic_comparison$AICc), 2)))

# Save results
saveRDS(results_summary, "advanced_varima_results.rds")
write.csv(aic_comparison, "advanced_varima_comparison.csv", row.names = FALSE)

print("\nResults saved to:")
print("- advanced_varima_results.rds")
print("- advanced_varima_comparison.csv")