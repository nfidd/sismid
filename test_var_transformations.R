# Test Different Transformations for VAR Model

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("ggplot2")
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

# Define various transformations
# 1. Fourth root (original)
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# 2. Square root
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# 3. Log with offset
log_offset <- function(x) log(x + 0.01)
inv_log_offset <- function(x) exp(x) - 0.01
my_log <- new_transformation(log_offset, inv_log_offset)

# 4. Logit for percentages (bounded between 0 and 100)
logit_pct <- function(x) {
  # Add small offset to avoid 0 and 100
  x_adj <- (x + 0.01) / 100.01
  log(x_adj / (1 - x_adj))
}
inv_logit_pct <- function(x) {
  p <- exp(x) / (1 + exp(x))
  p * 100.01 - 0.01
}
my_logit <- new_transformation(logit_pct, inv_logit_pct)

# Test date for comparison
test_date <- as.Date("2016-01-02")

# Prepare data
train_data <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

print("Testing different transformations...")
print(paste("Training data up to:", test_date))
print(paste("Number of observations:", nrow(train_data)))

# Function to test a transformation
test_transformation <- function(data, transform_func, transform_name, lag_order = 2) {
  print(paste("\n=== Testing", transform_name, "==="))
  
  tryCatch({
    # Fit model
    if (transform_name == "None") {
      # No transformation
      model <- data |>
        model(
          var = VAR(vars(
            `HHS Region 1`, `HHS Region 2`, `HHS Region 3`, 
            `HHS Region 4`, `HHS Region 5`, `HHS Region 6`,
            `HHS Region 7`, `HHS Region 8`, `HHS Region 9`, 
            `HHS Region 10`, `US National`
          ) ~ AR(lag_order))
        )
    } else {
      # With transformation - need to use the actual function object
      tf <- transform_func
      model <- data |>
        model(
          var = VAR(vars(
            `HHS Region 1` = tf(`HHS Region 1`),
            `HHS Region 2` = tf(`HHS Region 2`),
            `HHS Region 3` = tf(`HHS Region 3`),
            `HHS Region 4` = tf(`HHS Region 4`),
            `HHS Region 5` = tf(`HHS Region 5`),
            `HHS Region 6` = tf(`HHS Region 6`),
            `HHS Region 7` = tf(`HHS Region 7`),
            `HHS Region 8` = tf(`HHS Region 8`),
            `HHS Region 9` = tf(`HHS Region 9`),
            `HHS Region 10` = tf(`HHS Region 10`),
            `US National` = tf(`US National`)
          ) ~ AR(lag_order))
        )
    }
    
    # Get model info
    model_info <- glance(model)
    print(paste("AICc:", round(model_info$AICc, 2)))
    
    # Generate forecasts
    forecasts <- generate(model, h = 4, times = 50, bootstrap = TRUE)
    
    # Process forecasts
    special_cols <- c(".model", ".rep", "origin_date", ".sim_id")
    innov_cols <- names(forecasts)[grep("^\\.innov", names(forecasts))]
    all_special_cols <- c(special_cols, innov_cols)
    location_cols <- names(forecasts)[!names(forecasts) %in% all_special_cols]
    
    forecast_long <- forecasts |>
      as_tibble() |>
      select(all_of(c(".model", ".rep", "origin_date", location_cols))) |>
      pivot_longer(cols = all_of(location_cols),
                   names_to = "location",
                   values_to = "value") |>
      group_by(location, origin_date) |>
      summarise(
        mean_forecast = mean(value),
        q10 = quantile(value, 0.1),
        q90 = quantile(value, 0.9),
        .groups = "drop"
      )
    
    # Calculate forecast spread (uncertainty)
    avg_spread <- mean(forecast_long$q90 - forecast_long$q10)
    
    return(list(
      transform = transform_name,
      aicc = model_info$AICc,
      avg_spread = avg_spread,
      forecasts = forecast_long
    ))
    
  }, error = function(e) {
    print(paste("Error:", e$message))
    return(NULL)
  })
}

# Test all transformations
results <- list()

# 1. No transformation
results[[1]] <- test_transformation(train_data, NULL, "None")

# 2. Fourth root (current)
results[[2]] <- test_transformation(train_data, my_fourth_root, "Fourth Root")

# 3. Square root
results[[3]] <- test_transformation(train_data, my_sqrt, "Square Root")

# 4. Log with offset
results[[4]] <- test_transformation(train_data, my_log, "Log(x+0.01)")

# 5. Logit
results[[5]] <- test_transformation(train_data, my_logit, "Logit")

# Compile results
results_df <- bind_rows(lapply(results, function(x) {
  if (!is.null(x)) {
    data.frame(
      transformation = x$transform,
      aicc = x$aicc,
      avg_spread = x$avg_spread
    )
  }
}))

print("\n=== Transformation Comparison ===")
print(results_df |> arrange(aicc))

# Create visualization
p1 <- results_df |>
  ggplot(aes(x = reorder(transformation, -aicc), y = aicc)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Model Fit Comparison (AICc)",
    subtitle = "Lower is better",
    x = "Transformation",
    y = "AICc"
  ) +
  theme_minimal()

p2 <- results_df |>
  ggplot(aes(x = reorder(transformation, -avg_spread), y = avg_spread)) +
  geom_col(fill = "coral") +
  coord_flip() +
  labs(
    title = "Forecast Uncertainty",
    subtitle = "Average 80% prediction interval width",
    x = "Transformation",
    y = "Average Spread"
  ) +
  theme_minimal()

# Save plots
output_dir <- file.path(dirname(hub_path), "VAR_TRANSFORMATION_TESTS")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

ggsave(file.path(output_dir, "transformation_aicc.png"), 
       p1, width = 8, height = 6, dpi = 300)

ggsave(file.path(output_dir, "transformation_spread.png"), 
       p2, width = 8, height = 6, dpi = 300)

# Also test different lag orders with best transformation
print("\n=== Testing Different Lag Orders ===")
best_transform <- results_df |> filter(aicc == min(aicc)) |> pull(transformation)
print(paste("Best transformation:", best_transform))

# Get the transformation function for the best one
best_func <- switch(best_transform,
                    "None" = NULL,
                    "Fourth Root" = my_fourth_root,
                    "Square Root" = my_sqrt,
                    "Log(x+0.01)" = my_log,
                    "Logit" = my_logit)

# Test lag orders 1-5
lag_results <- list()
for (lag in 1:5) {
  print(paste("\nTesting AR(", lag, ")", sep = ""))
  result <- test_transformation(train_data, best_func, best_transform, lag)
  if (!is.null(result)) {
    lag_results[[lag]] <- data.frame(
      lag_order = lag,
      aicc = result$aicc
    )
  }
}

lag_df <- bind_rows(lag_results)
print("\n=== Lag Order Comparison ===")
print(lag_df |> arrange(aicc))

# Save results
write.csv(results_df, 
          file.path(output_dir, "transformation_comparison.csv"),
          row.names = FALSE)

write.csv(lag_df,
          file.path(output_dir, "lag_order_comparison.csv"),
          row.names = FALSE)

print(paste("\nResults saved to:", output_dir))