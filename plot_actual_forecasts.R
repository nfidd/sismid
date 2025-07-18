# Plot the actual VAR(2) sqrt forecasts with prediction intervals

library("nfidd")
library("dplyr")
library("ggplot2")
library("tidyr")
library("hubData")

# Set up theme
theme_set(theme_bw())

# Load our forecasts
forecasts <- read.csv("var2_sqrt_corrected_forecasts.csv")

print("=== Plotting Actual VAR(2) sqrt Forecasts ===")

# Get target data for comparison
hub_path <- here::here("sismid-ili-forecasting-sandbox")
oracle_output <- connect_target_oracle_output(hub_path) |>
  collect() |>
  rename(observation = oracle_value) |>
  mutate(target_end_date = as.Date(target_end_date))

# Load historical data for context
data(flu_data_hhs)
historical_data <- flu_data_hhs |>
  mutate(origin_date = as.Date(origin_date)) |>
  filter(origin_date >= as.Date("2017-01-01"))

# Select a few recent forecast dates to show
sample_dates <- c("2018-12-01", "2019-12-07", "2020-01-04")

print(paste("Plotting forecasts for dates:", paste(sample_dates, collapse = ", ")))

# Process forecasts for plotting
plot_forecasts <- forecasts |>
  filter(origin_date %in% sample_dates) |>
  mutate(
    origin_date = as.Date(origin_date),
    target_end_date = as.Date(target_end_date),
    forecast_date = paste("Forecast from", origin_date)
  )

# Get prediction intervals
forecast_intervals <- plot_forecasts |>
  filter(output_type_id %in% c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)) |>
  select(origin_date, forecast_date, location, target_end_date, output_type_id, value) |>
  pivot_wider(names_from = output_type_id, values_from = value, names_prefix = "q")

# Get observations for comparison
observations <- oracle_output |>
  filter(target_end_date >= as.Date("2017-01-01")) |>
  select(location, target_end_date, observation)

# Create forecast plot for each location
for (loc in unique(forecast_intervals$location)) {
  
  # Filter data for this location
  loc_forecasts <- forecast_intervals |> filter(location == loc)
  loc_observations <- observations |> filter(location == loc)
  loc_historical <- historical_data |> filter(location == loc)
  
  # Create the plot
  p <- ggplot() +
    # Historical data (before forecasts)
    geom_line(data = loc_historical, 
              aes(x = origin_date, y = wili), 
              color = "black", alpha = 0.7, size = 0.5) +
    
    # 95% prediction interval
    geom_ribbon(data = loc_forecasts, 
                aes(x = target_end_date, ymin = q0.025, ymax = q0.975, fill = forecast_date), 
                alpha = 0.2) +
    
    # 80% prediction interval
    geom_ribbon(data = loc_forecasts, 
                aes(x = target_end_date, ymin = q0.1, ymax = q0.9, fill = forecast_date), 
                alpha = 0.3) +
    
    # 50% prediction interval
    geom_ribbon(data = loc_forecasts, 
                aes(x = target_end_date, ymin = q0.25, ymax = q0.75, fill = forecast_date), 
                alpha = 0.5) +
    
    # Median forecast
    geom_line(data = loc_forecasts, 
              aes(x = target_end_date, y = q0.5, color = forecast_date), 
              size = 1.2) +
    
    # Actual observations
    geom_point(data = loc_observations, 
               aes(x = target_end_date, y = observation), 
               color = "red", size = 1.5, alpha = 0.8) +
    
    # Forecast origin markers
    geom_vline(data = loc_forecasts |> distinct(origin_date, forecast_date), 
               aes(xintercept = origin_date, color = forecast_date), 
               linetype = "dashed", alpha = 0.7) +
    
    scale_fill_viridis_d(name = "Forecast Date") +
    scale_color_viridis_d(name = "Forecast Date") +
    
    labs(
      title = paste("VAR(2) sqrt Forecasts -", loc),
      subtitle = "Shaded areas: 50%, 80%, 95% prediction intervals | Red dots: Actual observations",
      x = "Date",
      y = "% ILI Visits"
    ) +
    
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  print(p)
  
  # Save individual plot
  filename <- paste0("var2_sqrt_forecast_", gsub(" ", "_", loc), ".png")
  ggsave(filename, p, width = 12, height = 8, dpi = 300)
}

# Create a combined plot for key regions
key_regions <- c("US National", "HHS Region 1", "HHS Region 4", "HHS Region 9")

combined_data <- forecast_intervals |> 
  filter(location %in% key_regions) |>
  left_join(observations, by = c("location", "target_end_date"))

p_combined <- ggplot(combined_data) +
  # 80% prediction interval
  geom_ribbon(aes(x = target_end_date, ymin = q0.1, ymax = q0.9, fill = forecast_date), 
              alpha = 0.3) +
  
  # 50% prediction interval
  geom_ribbon(aes(x = target_end_date, ymin = q0.25, ymax = q0.75, fill = forecast_date), 
              alpha = 0.5) +
  
  # Median forecast
  geom_line(aes(x = target_end_date, y = q0.5, color = forecast_date), 
            size = 1) +
  
  # Actual observations
  geom_point(aes(x = target_end_date, y = observation), 
             color = "red", size = 1.2) +
  
  # Forecast origin markers
  geom_vline(aes(xintercept = origin_date, color = forecast_date), 
             linetype = "dashed", alpha = 0.7) +
  
  facet_wrap(~ location, scales = "free_y", ncol = 2) +
  
  scale_fill_viridis_d(name = "Forecast Date") +
  scale_color_viridis_d(name = "Forecast Date") +
  
  labs(
    title = "VAR(2) sqrt Forecasts - Key Regions",
    subtitle = "Shaded areas: 50% & 80% prediction intervals | Red dots: Actual observations",
    x = "Date",
    y = "% ILI Visits"
  ) +
  
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )

print(p_combined)

# Save combined plot
ggsave("var2_sqrt_forecasts_combined.png", p_combined, width = 14, height = 10, dpi = 300)

print("\n=== Forecast Plots Created! ===")
print("Files saved:")
print("- Individual plots for each location")
print("- var2_sqrt_forecasts_combined.png (key regions)")
print("")
print("What the plots show:")
print("- Black lines: Historical ILI data")
print("- Colored ribbons: 50%, 80%, 95% prediction intervals")
print("- Colored lines: Median forecasts")
print("- Red dots: Actual observations")
print("- Dashed lines: When forecasts were made")
print("")
print("This shows how well our VAR(2) sqrt model actually predicts!")