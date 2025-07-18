# Run hub-playground.qmd code to understand the workflow
# This script runs all code from the hub-playground session

# Set options
options(width=100)

# Load libraries
library("nfidd")
library("dplyr")
library("ggplot2")
library("epidatr")
library("fable")
library("hubUtils")
library("hubEvals")
library("hubVis")
library("hubData")
library("hubEnsembles")
theme_set(theme_bw())

# Set seed for reproducibility
set.seed(406)

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")
hub_con <- connect_hub(hub_path)

# Look at one forecast to understand structure
example_forecast <- hub_con |> 
  filter(model_id == "delphi-epicast",
         origin_date == "2015-11-14") |> 
  collect_hub()

print("Example forecast structure:")
print(head(example_forecast))

# Define location information
locs <- c("nat", 
          "hhs1", "hhs2", "hhs3", "hhs4", "hhs5", 
          "hhs6", "hhs7", "hhs8", "hhs9", "hhs10")
location_formal_names <- c("US National", paste("HHS Region", 1:10))
loc_df <- data.frame(region = locs, location = location_formal_names)

# Load the saved flu data
data(flu_data_hhs)

print("Flu data structure:")
print(flu_data_hhs)

# Build ARIMA(2,1,0) with fourth-root transformation
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)

# Fit ARIMA model to data up to 2015-10-24
fit_arima210 <- flu_data_hhs |>
  filter(origin_date <= "2015-10-24") |> 
  model(ARIMA(my_fourth_root(wili) ~ pdq(2,1,0)))

# Generate forecasts
arima_forecast <- forecast(fit_arima210, h=4)

# Plot forecasts
forecast_plot <- autoplot(arima_forecast, 
                         flu_data_hhs |> 
                           filter(origin_date <= "2015-10-24",
                                  origin_date >= "2014-09-01")) +
  facet_wrap(.~location, scales = "free_y") +
  labs(y = "% of visits due to ILI", x = NULL)

print("ARIMA(2,1,0) model fitted for all locations")

# Inspect models for specific regions
print("\nModel for HHS Region 1:")
fit_arima210 |> 
  filter(location == "HHS Region 1") |> 
  report(model)

print("\nModel for HHS Region 9:")
fit_arima210 |> 
  filter(location == "HHS Region 9") |> 
  report(model)

# Get origin dates from hub
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()
print("\nValid origin dates:")
print(origin_dates)

# Define quantile levels for hub submission
quantile_levels <- c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)

print("\nSetup complete - ready to design VAR model")