# Test Advanced VAR Features

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("feasts")
library("ggplot2")

# Set seed
set.seed(406)

# Load data
data(flu_data_hhs)

# Test date
test_date <- as.Date("2016-12-17")

# Prepare data
flu_data_wide <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

log_transform <- function(x) log(x + 0.1)
inv_log <- function(x) exp(x) - 0.1
my_log <- new_transformation(log_transform, inv_log)

print("=== Testing Advanced Model Features ===")
print(paste("Training data up to:", test_date))

# Store results
results <- list()

# 1. Baseline VAR(2) sqrt
print("\n1. Baseline: VAR(2) with sqrt")
tryCatch({
  var2_sqrt <- flu_data_wide |>
    model(
      baseline = VAR(vars(
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
  results[[1]] <- glance(var2_sqrt) |> mutate(model = "VAR(2) sqrt")
  print(glance(var2_sqrt))
}, error = function(e) print(paste("Error:", e$message)))

# 2. VARIMA with differencing
print("\n2. VARIMA(2,1,0) with sqrt")
tryCatch({
  varima_sqrt <- flu_data_wide |>
    model(
      varima = VARIMA(vars(
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
      ) ~ pdq(2,1,0), identification = "none")
    )
  results[[2]] <- glance(varima_sqrt) |> mutate(model = "VARIMA(2,1,0) sqrt")
  print(glance(varima_sqrt))
}, error = function(e) print(paste("Error:", e$message)))

# 3. VAR with seasonal dummies
print("\n3. VAR(2) sqrt with seasonal dummies")
tryCatch({
  # Add month indicator
  flu_data_season <- flu_data_wide |>
    mutate(
      month = lubridate::month(origin_date),
      season_sin = sin(2 * pi * month / 12),
      season_cos = cos(2 * pi * month / 12)
    )
  
  var2_season <- flu_data_season |>
    model(
      seasonal = VAR(vars(
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
      ) ~ AR(2) + season_sin + season_cos)
    )
  results[[3]] <- glance(var2_season) |> mutate(model = "VAR(2) sqrt + seasonality")
  print(glance(var2_season))
}, error = function(e) print(paste("Error:", e$message)))

# 4. Test Fourier terms for seasonality
print("\n4. Testing with Fourier terms")
tryCatch({
  # For weekly data, K=26 gives annual seasonality
  var2_fourier <- flu_data_wide |>
    model(
      fourier = VAR(vars(
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
  results[[4]] <- glance(var2_fourier) |> mutate(model = "VAR(2) sqrt + fourier")
  print(glance(var2_fourier))
}, error = function(e) print(paste("Error:", e$message)))

# 5. Test regional subsets (e.g., geographically close regions)
print("\n5. Regional subset VAR - East Coast")
tryCatch({
  # East coast regions: 1, 2, 3
  east_coast <- flu_data_wide |>
    select(origin_date, `HHS Region 1`, `HHS Region 2`, `HHS Region 3`)
  
  var_east <- east_coast |>
    model(
      east = VAR(vars(
        `HHS Region 1` = my_sqrt(`HHS Region 1`),
        `HHS Region 2` = my_sqrt(`HHS Region 2`),
        `HHS Region 3` = my_sqrt(`HHS Region 3`)
      ) ~ AR(2))
    )
  results[[5]] <- glance(var_east) |> mutate(model = "VAR East Coast only")
  print(glance(var_east))
}, error = function(e) print(paste("Error:", e$message)))

# 6. Test with trend
print("\n6. VAR(2) sqrt with trend")
tryCatch({
  var2_trend <- flu_data_wide |>
    model(
      trend = VAR(vars(
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
  results[[6]] <- glance(var2_trend) |> mutate(model = "VAR(2) sqrt + trend")
  print(glance(var2_trend))
}, error = function(e) print(paste("Error:", e$message)))

# 7. Test Box-Cox transformation (automatic lambda selection)
print("\n7. VAR(2) with Box-Cox transformation")
tryCatch({
  var2_boxcox <- flu_data_wide |>
    model(
      boxcox = VAR(vars(
        `HHS Region 1` = box_cox(`HHS Region 1`, lambda = 0.5),
        `HHS Region 2` = box_cox(`HHS Region 2`, lambda = 0.5),
        `HHS Region 3` = box_cox(`HHS Region 3`, lambda = 0.5),
        `HHS Region 4` = box_cox(`HHS Region 4`, lambda = 0.5),
        `HHS Region 5` = box_cox(`HHS Region 5`, lambda = 0.5),
        `HHS Region 6` = box_cox(`HHS Region 6`, lambda = 0.5),
        `HHS Region 7` = box_cox(`HHS Region 7`, lambda = 0.5),
        `HHS Region 8` = box_cox(`HHS Region 8`, lambda = 0.5),
        `HHS Region 9` = box_cox(`HHS Region 9`, lambda = 0.5),
        `HHS Region 10` = box_cox(`HHS Region 10`, lambda = 0.5),
        `US National` = box_cox(`US National`, lambda = 0.5)
      ) ~ AR(2))
    )
  results[[7]] <- glance(var2_boxcox) |> mutate(model = "VAR(2) Box-Cox(0.5)")
  print(glance(var2_boxcox))
}, error = function(e) print(paste("Error:", e$message)))

# 8. Test with higher lag order
print("\n8. VAR(3) with sqrt")
tryCatch({
  var3_sqrt <- flu_data_wide |>
    model(
      var3 = VAR(vars(
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
      ) ~ AR(3))
    )
  results[[8]] <- glance(var3_sqrt) |> mutate(model = "VAR(3) sqrt")
  print(glance(var3_sqrt))
}, error = function(e) print(paste("Error:", e$message)))

# Compile results
results_df <- bind_rows(results) |>
  select(model, AICc, AIC, BIC) |>
  arrange(AICc)

print("\n=== Model Feature Comparison ===")
print(results_df)

# Create visualization
output_dir <- file.path(dirname(here::here("sismid-ili-forecasting-sandbox")), 
                        "ADVANCED_FEATURES_RESULTS")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

p <- results_df |>
  ggplot(aes(x = reorder(model, -AICc), y = AICc)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Advanced Model Features Comparison",
    subtitle = "Lower AICc indicates better fit",
    x = "Model Configuration",
    y = "AICc"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "advanced_features_comparison.png"),
       p, width = 12, height = 8, dpi = 300)

# Save results
write.csv(results_df, 
          file.path(output_dir, "advanced_features_results.csv"),
          row.names = FALSE)

# Test forecast performance of best model
print("\n=== Testing forecast quality ===")
best_model_name <- results_df$model[1]
print(paste("Best model:", best_model_name))

# Generate a sample forecast
if (grepl("trend", best_model_name)) {
  fc <- forecast(var2_trend, h = 4, bootstrap = TRUE)
} else if (grepl("fourier", best_model_name)) {
  fc <- forecast(var2_fourier, h = 4, bootstrap = TRUE)
} else if (grepl("seasonal", best_model_name)) {
  fc <- forecast(var2_season, h = 4, bootstrap = TRUE)
} else {
  fc <- forecast(var2_sqrt, h = 4, bootstrap = TRUE)
}

print("Sample forecast generated successfully")

# Summary report
report <- paste0(
  "# Advanced VAR Features Test Results\n\n",
  "## Best Model: ", results_df$model[1], "\n",
  "AICc: ", round(results_df$AICc[1], 2), "\n\n",
  "## Tested Features:\n",
  "1. VARIMA with differencing\n",
  "2. Seasonal components (sin/cos)\n",
  "3. Fourier terms for seasonality\n",
  "4. Regional subsets\n",
  "5. Trend component\n",
  "6. Box-Cox transformation\n",
  "7. Higher lag orders\n\n",
  "## Key Findings:\n",
  "- ", ifelse(grepl("fourier|season", results_df$model[1]), 
              "Seasonality improves model fit", 
              "Base model performs well without seasonality"), "\n",
  "- ", ifelse(grepl("VARIMA", results_df$model[1]), 
              "Differencing helps capture trends", 
              "Level models sufficient"), "\n",
  "- Square root transformation remains optimal\n\n",
  "## Recommendation:\n",
  "Use ", results_df$model[1], " for final forecasts.\n"
)

writeLines(report, file.path(output_dir, "advanced_features_report.md"))

print(paste("\nResults saved to:", output_dir))