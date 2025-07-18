# VAR Model Development for ILI Forecasting

# Load required libraries
library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("ggplot2")
library("hubUtils")
library("hubEvals")
library("hubVis")
library("hubData")
library("hubEnsembles")

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

print("=== Building VAR Model ===")

# Prepare data for VAR - need wide format
# First, let's test with a smaller subset to ensure it works
test_data <- flu_data_hhs |>
  filter(origin_date >= "2015-01-01",
         origin_date <= "2015-10-24") |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

print("Data structure for VAR:")
print(test_data)

# Define transformations (same as ARIMA model for consistency)
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# Fit VAR model with fourth-root transformation
# We'll use a VAR(2) to match the ARIMA(2,1,0) lag structure
print("\n=== Fitting VAR(2) model ===")
fit_var2 <- test_data |>
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

print("VAR model fitted successfully!")

# Generate forecasts
print("\n=== Generating VAR forecasts ===")
var_forecast <- forecast(fit_var2, h = 4, bootstrap = TRUE)

# Plot forecasts
var_plot <- var_forecast |>
  autoplot(test_data |> filter(origin_date >= "2014-09-01")) +
  facet_wrap(~.var, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(title = "VAR(2) Model Forecasts",
       x = "Date",
       y = "% visits due to ILI")

ggsave("var_forecasts.png", var_plot, width = 12, height = 8)

# Check model performance on one forecast
print("\n=== Model summary ===")
print(fit_var2)

# Now let's prepare the full validation phase forecasts
print("\n=== Preparing validation phase forecasts ===")

# We'll need to adapt the time-series cross-validation approach for VAR
# First, let's define a function to generate VAR forecasts in hub format
generate_var_forecasts <- function(data, origin_dates, quantile_levels, model_id) {
  
  all_forecasts <- list()
  
  for (i in seq_along(origin_dates)) {
    origin_date <- origin_dates[i]
    
    # Filter data up to origin date
    train_data <- data |>
      filter(origin_date <= !!origin_date) |>
      select(origin_date, location, wili) |>
      pivot_wider(names_from = location, values_from = wili) |>
      as_tsibble(index = origin_date)
    
    # Fit VAR model
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
    
    # Generate forecasts with bootstrap samples
    var_forecasts <- var_model |>
      generate(h = 4, times = 100, bootstrap = TRUE) |>
      pivot_longer(cols = -c(origin_date, .model, .rep, .sim_id),
                   names_to = "location",
                   values_to = ".sim") |>
      group_by(.rep, location, origin_date) |>
      mutate(horizon = row_number()) |>
      ungroup() |>
      as_tibble() |>
      rename(
        target_end_date = origin_date,
        value = .sim
      ) |>
      mutate(
        origin_date = !!origin_date,
        model_id = model_id
      ) |>
      left_join(loc_df, by = "location") |>
      group_by(model_id, location, origin_date, horizon, target_end_date) |>
      reframe(tibble::enframe(quantile(value, quantile_levels), "quantile", "value")) |>
      mutate(
        output_type_id = as.numeric(stringr::str_remove(quantile, "%"))/100,
        target = "ili perc",
        output_type = "quantile"
      ) |>
      select(-quantile)
    
    all_forecasts[[i]] <- var_forecasts
    
    print(paste("Completed forecast for origin date:", origin_date))
  }
  
  bind_rows(all_forecasts)
}

print("\nVAR model development complete - ready for validation phase forecasting!")

# Save the workspace for next steps
save.image("var_model_workspace.RData")