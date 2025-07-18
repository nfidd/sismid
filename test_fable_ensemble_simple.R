# Test Simple Fable Ensemble
# Combine multiple models to potentially beat baseline

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("hubUtils")
library("ggplot2")

set.seed(406)

# Load data
data(flu_data_hhs)

# Transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# Test on one date first
test_date <- as.Date("2016-11-05")

print("=== Testing Fable Ensemble Approach ===")

# Prepare data
train_data <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  mutate(week = yearweek(origin_date)) |>
  as_tsibble(index = week, key = location)

# Fit multiple models
print("Fitting component models...")

models <- train_data |>
  model(
    # Model 1: ARIMA with auto selection
    arima_auto = ARIMA(my_sqrt(wili)),
    
    # Model 2: ETS
    ets = ETS(my_sqrt(wili)),
    
    # Model 3: ARIMA with specific order
    arima_211 = ARIMA(my_sqrt(wili) ~ pdq(2,1,1)),
    
    # Model 4: Seasonal naive as baseline
    snaive = SNAIVE(my_sqrt(wili))
  )

print("Models fitted successfully")

# Generate forecasts
forecasts <- models |> forecast(h = 4)

print("Forecasts generated")

# Create ensemble by averaging
ensemble_fc <- forecasts |>
  as_tibble() |>
  group_by(location, week) |>
  summarise(
    # Average the point forecasts
    ensemble_mean = mean(c(.mean[.model == "arima_auto"],
                          .mean[.model == "ets"],
                          .mean[.model == "arima_211"],
                          .mean[.model == "snaive"])),
    .groups = "drop"
  )

print("Ensemble created")

# Plot sample forecasts for one location
sample_location <- "HHS Region 1"

plot_data <- train_data |>
  filter(location == sample_location,
         week >= yearweek("2016-01-01")) |>
  as_tibble()

forecast_data <- forecasts |>
  filter(location == sample_location) |>
  as_tibble()

ensemble_data <- ensemble_fc |>
  filter(location == sample_location)

p <- ggplot() +
  geom_line(data = plot_data, 
            aes(x = as.Date(week), y = wili),
            color = "black") +
  geom_line(data = forecast_data,
            aes(x = as.Date(week), y = .mean, 
                color = .model, linetype = .model)) +
  geom_point(data = ensemble_data,
             aes(x = as.Date(week), y = inv_sqrt(ensemble_mean)),
             color = "red", size = 3) +
  labs(title = paste("Forecasts for", sample_location),
       x = "Date", y = "ILI %",
       color = "Model", linetype = "Model") +
  theme_minimal()

ggsave("fable_ensemble_sample.png", p, width = 10, height = 6)

# Now test on multiple dates
print("\nTesting on validation dates...")

validation_dates <- seq(as.Date("2016-11-05"), 
                       as.Date("2016-12-31"), 
                       by = "week")

all_forecasts <- list()

for (i in seq_along(validation_dates)) {
  origin_date <- validation_dates[i]
  
  tryCatch({
    # Prepare training data
    train <- flu_data_hhs |>
      filter(origin_date <= !!origin_date) |>
      mutate(week = yearweek(origin_date)) |>
      as_tsibble(index = week, key = location)
    
    # Fit models
    mods <- train |>
      model(
        arima = ARIMA(my_sqrt(wili)),
        ets = ETS(my_sqrt(wili)),
        arima_spec = ARIMA(my_sqrt(wili) ~ pdq(2,1,1))
      )
    
    # Generate forecasts
    fc <- mods |> forecast(h = 4)
    
    # Create ensemble
    ensemble <- fc |>
      as_tibble() |>
      group_by(location, week) |>
      summarise(
        mean_val = mean(.mean),
        .groups = "drop"
      ) |>
      mutate(
        origin_date = origin_date,
        target = "wk inc flu hosp",
        target_end_date = as.Date(week),
        value = inv_sqrt(mean_val),  # Transform back
        value = pmax(value, 0)  # Ensure non-negative
      )
    
    all_forecasts[[length(all_forecasts) + 1]] <- ensemble
    
  }, error = function(e) {
    print(paste("Error at", origin_date, ":", e$message))
  })
  
  if (i %% 3 == 0) {
    print(paste("Completed", i, "of", length(validation_dates)))
  }
}

# Combine results
if (length(all_forecasts) > 0) {
  results_df <- do.call(rbind, all_forecasts)
  
  print("\nEnsemble forecasts summary:")
  print(paste("Total forecasts:", nrow(results_df)))
  print(paste("Locations:", n_distinct(results_df$location)))
  print(paste("Dates:", n_distinct(results_df$origin_date)))
  
  # Save sample results
  write.csv(results_df |> head(100), 
            "fable_ensemble_sample_results.csv", 
            row.names = FALSE)
  
  print("\nResults saved to fable_ensemble_sample_results.csv")
  print("Next steps:")
  print("1. Convert to full quantile format for hub submission")
  print("2. Calculate WIS score using hub evaluation")
  print("3. Compare against baseline VAR(2) sqrt (0.208)")
}

# Test VAR model for comparison
print("\nTesting VAR model for comparison...")

# Prepare wide format data
train_wide <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  dplyr::select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  mutate(week = yearweek(origin_date)) |>
  as_tsibble(index = week)

# Fit VAR model
var_model <- train_wide |>
  model(
    var = VAR(vars(
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

print("VAR model fitted")

# Check model information
print(glance(var_model))

print("\n=== Summary ===")
print("Successfully tested ensemble approach combining:")
print("- ARIMA (automatic selection)")
print("- ETS")
print("- ARIMA(2,1,1)")
print("- Seasonal naive (baseline)")
print("\nEnsemble uses simple averaging of point forecasts")
print("Consider weighted ensemble based on cross-validation performance")