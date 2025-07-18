# Explore Additional VAR Model Variations
# Given the huge seasonality improvement, test other enhancements

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("feasts")

# Set seed
set.seed(406)

# Load data
data(flu_data_hhs)

# Prepare data
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

# Split data
train_cutoff <- as.Date("2018-05-05")
train_data <- flu_data_wide |>
  filter(origin_date <= train_cutoff)

print(paste("Training data:", nrow(train_data), "weeks"))

# Current best model (baseline for comparison)
model_seasonal_baseline <- train_data |>
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

print("=== Testing Additional Model Variations ===")

# Variation 1: Different lag orders with seasonality
model_var1_seasonal <- train_data |>
  model(
    var1_seasonal = VAR(vars(
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
    ) ~ AR(1) + season(period = 52))
  )

model_var3_seasonal <- train_data |>
  model(
    var3_seasonal = VAR(vars(
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
    ) ~ AR(3) + season(period = 52))
  )

# Variation 2: Log transformation with seasonality
model_log_seasonal <- train_data |>
  model(
    log_seasonal = VAR(vars(
      `HHS Region 1` = my_log(`HHS Region 1`),
      `HHS Region 2` = my_log(`HHS Region 2`),
      `HHS Region 3` = my_log(`HHS Region 3`),
      `HHS Region 4` = my_log(`HHS Region 4`),
      `HHS Region 5` = my_log(`HHS Region 5`),
      `HHS Region 6` = my_log(`HHS Region 6`),
      `HHS Region 7` = my_log(`HHS Region 7`),
      `HHS Region 8` = my_log(`HHS Region 8`),
      `HHS Region 9` = my_log(`HHS Region 9`),
      `HHS Region 10` = my_log(`HHS Region 10`),
      `US National` = my_log(`US National`)
    ) ~ AR(2) + season(period = 52))
  )

# Variation 3: Add trend component
model_trend_seasonal <- train_data |>
  model(
    trend_seasonal = VAR(vars(
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
    ) ~ AR(2) + trend() + season(period = 52))
  )

# Variation 4: Try automatic lag selection
model_auto_lag <- train_data |>
  model(
    auto_lag = VAR(vars(
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
    ) ~ AR(0:4, ic = "aic") + season(period = 52))
  )

# Variation 5: Different seasonal periods (maybe bi-annual patterns?)
model_biannual <- train_data |>
  model(
    biannual = VAR(vars(
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
    ) ~ AR(2) + season(period = 26))  # 26 weeks = 6 months
  )

# Compare all models
models_list <- list(
  model_seasonal_baseline,
  model_var1_seasonal,
  model_var3_seasonal,
  model_log_seasonal,
  model_trend_seasonal,
  model_auto_lag,
  model_biannual
)

names(models_list) <- c(
  "VAR(2) + season(52)", 
  "VAR(1) + season(52)", 
  "VAR(3) + season(52)",
  "VAR(2) + season(52) + log", 
  "VAR(2) + trend + season(52)",
  "VAR(auto) + season(52)",
  "VAR(2) + season(26)"
)

# Extract AICc values
print("\n=== Model Performance Comparison ===")

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
    aic_comparison <- rbind(aic_comparison, 
                           data.frame(Model = model_name, 
                                    AICc = NA, 
                                    Status = paste("Error:", e$message)))
  })
}

# Clean and rank results
aic_comparison_clean <- aic_comparison |>
  filter(Status == "Success") |>
  arrange(AICc) |>
  mutate(
    Delta_AICc = AICc - min(AICc),
    Rank = row_number()
  )

print("All models (including errors):")
print(aic_comparison)

print("\nRanked successful models:")
print(aic_comparison_clean)

# Test forecast generation for best model
if (nrow(aic_comparison_clean) > 0) {
  best_model_name <- aic_comparison_clean$Model[1]
  best_model_idx <- which(names(models_list) == best_model_name)
  best_model <- models_list[[best_model_idx]]
  
  print(paste("\nBest model:", best_model_name))
  print(paste("AICc:", round(aic_comparison_clean$AICc[1], 2)))
  
  # Quick forecast test
  print("Testing forecast generation...")
  forecasts <- best_model |>
    generate(h = 4, times = 50, bootstrap = TRUE)
  
  print(paste("✓ Generated", nrow(forecasts), "forecast samples"))
  
  # Calculate improvement over original VAR(2) sqrt baseline
  original_var2_aic <- -15253.50  # From previous analysis
  improvement <- original_var2_aic - aic_comparison_clean$AICc[1]
  
  print(paste("Improvement over original VAR(2) sqrt:", round(improvement, 2), "AICc points"))
}

# Save results
results_summary <- list(
  comparison = aic_comparison,
  ranked_models = aic_comparison_clean,
  best_model = if(nrow(aic_comparison_clean) > 0) best_model_name else "None",
  recommendation = if(nrow(aic_comparison_clean) > 0) {
    paste("Best model is", best_model_name, "with AICc =", 
          round(aic_comparison_clean$AICc[1], 2))
  } else "No successful models"
)

saveRDS(results_summary, "additional_var_models_results.rds")
write.csv(aic_comparison, "additional_var_models_comparison.csv", row.names = FALSE)

print("\n=== Final Recommendation ===")
print(results_summary$recommendation)
print("\nFiles saved:")
print("- additional_var_models_results.rds")
print("- additional_var_models_comparison.csv")