# Create VAR(2) log model forecasts for test phase - Nick Reich approved approach!

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")

set.seed(406)
data(flu_data_hhs)

# Location mapping
locs <- c("nat", "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Log transformation
log_transform <- function(x) log(x + 0.1)
inv_log <- function(x) exp(x) - 0.1
my_log <- new_transformation(log_transform, inv_log)

# Test on a few dates first
test_dates <- c("2018-12-01", "2019-01-05", "2019-02-02")

print("=== Creating VAR(2) Log Model Forecasts ===")
print("Simple multivariate model without forced seasonality")

all_forecasts <- list()

for (date in test_dates) {
  print(paste("Processing date:", date))
  
  # Prepare data
  flu_data_wide <- flu_data_hhs |>
    filter(origin_date <= as.Date(date)) |>
    select(origin_date, location, wili) |>
    pivot_wider(names_from = location, values_from = wili) |>
    as_tsibble(index = origin_date)
  
  # Fit VAR(2) log model
  model_fit <- flu_data_wide |>
    model(
      var2_log = VAR(vars(
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
      ) ~ AR(2))
    )
  
  print(paste("AICc:", round(glance(model_fit)$AICc, 2)))
  
  # Generate forecasts with bootstrap
  forecast_samples <- model_fit |>
    generate(h = 4, times = 100, bootstrap = TRUE) |>
    group_by(.rep, .model) |>
    mutate(horizon = row_number()) |>
    ungroup() |>
    as_tibble()
  
  # Convert to hub format
  location_cols <- c("HHS Region 1", "HHS Region 2", "HHS Region 3", 
                     "HHS Region 4", "HHS Region 5", "HHS Region 6",
                     "HHS Region 7", "HHS Region 8", "HHS Region 9", 
                     "HHS Region 10", "US National")
  
  forecast_long <- forecast_samples |>
    pivot_longer(
      cols = all_of(location_cols),
      names_to = "location",
      values_to = "value"
    ) |>
    rename(target_end_date = origin_date) |>
    mutate(
      model_id = "sismid-var2-log",
      target = "ili perc",
      origin_date = target_end_date - horizon * 7
    )
  
  # Calculate quantiles
  quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)
  
  forecast_quantiles <- forecast_long |>
    group_by(model_id, location, origin_date, horizon, target_end_date) |>
    reframe(
      quantile_level = quantile_levels,
      value = quantile(value, quantile_levels, na.rm = TRUE)
    ) |>
    mutate(
      output_type = "quantile",
      output_type_id = quantile_level
    ) |>
    select(-quantile_level)
  
  all_forecasts[[date]] <- forecast_quantiles
  
  print(paste("Generated", nrow(forecast_quantiles), "forecast rows"))
}

# Combine all forecasts
combined_forecasts <- bind_rows(all_forecasts)

print(paste("Total forecasts generated:", nrow(combined_forecasts)))
print(paste("Unique origin dates:", length(unique(combined_forecasts$origin_date))))

# Save forecasts
write.csv(combined_forecasts, "var2_log_forecasts.csv", row.names = FALSE)

# Compare performance metrics
print("\n=== Performance Comparison ===")
print("Model configurations tested:")
print("1. VAR(1) + season(52) [our current]: 0.466 WIS")
print("2. VAR(2) log [this model]: Testing now...")
print("3. VAR(2) log had AICc of -6,956 (much better than seasonal)")
print("")
print("Key insight: Nick Reich is right - simpler models often perform better!")
print("The multivariate VAR captures cross-regional dependencies without overfitting")

# Show forecast characteristics
forecast_summary <- combined_forecasts |>
  group_by(model_id, location) |>
  summarise(
    mean_forecast = mean(value, na.rm = TRUE),
    median_forecast = median(value, na.rm = TRUE),
    min_forecast = min(value, na.rm = TRUE),
    max_forecast = max(value, na.rm = TRUE),
    .groups = "drop"
  )

print("\nForecast summary by location:")
print(forecast_summary)

print("\n=== Next Steps ===")
print("1. Evaluate WIS performance of this model")
print("2. If better than 0.466, generate full test phase forecasts")
print("3. Submit to hub as sismid-var2-log model")
print("4. Thank Nick Reich for steering us away from over-complicated seasonality!")