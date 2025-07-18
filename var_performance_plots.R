# Create Performance Plots for VAR Model Comparison

library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("hubUtils")
library("hubEvals")
library("hubData")
library("hubVis")
library("patchwork")

# Set theme
theme_set(theme_bw())

# Hub configuration
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Get validation forecasts
origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()
validation_origin_dates <- origin_dates[which(as.Date(origin_dates) <= as.Date("2017-05-06"))]

# Connect to hub
hub_con <- connect_hub(hub_path)
validation_forecasts <- hub_con |> 
  filter(origin_date %in% validation_origin_dates) |> 
  collect_hub()

# Get oracle output
oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

# Score models
model_scores <- hubEvals::score_model_out(
  validation_forecasts,
  oracle_output
)

print("Creating performance visualizations...")
print("Available columns in model_scores:")
print(names(model_scores))

# 1. Overall WIS Comparison
p1 <- model_scores |>
  group_by(model_id) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    se_wis = sd(wis, na.rm = TRUE) / sqrt(n())
  ) |>
  ggplot(aes(x = reorder(model_id, mean_wis), y = mean_wis)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(
    aes(ymin = mean_wis - 1.96 * se_wis, 
        ymax = mean_wis + 1.96 * se_wis),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "Model Performance Comparison",
    subtitle = "Mean Weighted Interval Score (WIS) - Lower is Better",
    x = "Model",
    y = "Mean WIS",
    caption = "Error bars show 95% confidence intervals"
  ) +
  theme_minimal()

# 2. WIS by Horizon
p2 <- model_scores |>
  group_by(model_id, horizon) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(x = horizon, y = mean_wis, color = model_id)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Forecast Accuracy by Horizon",
    subtitle = "How accuracy degrades with forecast horizon",
    x = "Forecast Horizon (weeks)",
    y = "Mean WIS",
    color = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 3. WIS by Location
p3 <- model_scores |>
  filter(model_id %in% c("sismid-var2", "sismid-arima210", "delphi-epicast")) |>
  group_by(model_id, location) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(x = location, y = mean_wis, fill = model_id)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  labs(
    title = "Model Performance by Location",
    subtitle = "Comparing VAR, ARIMA, and Delphi models",
    x = "Location",
    y = "Mean WIS",
    fill = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 4. Calibration Plot (PIT histograms)
# Get forecast data for PIT calculation
pit_data <- validation_forecasts |>
  filter(model_id %in% c("sismid-var2", "sismid-arima210", "delphi-epicast"),
         output_type == "quantile") |>
  left_join(oracle_output, by = c("location", "target_end_date", "target", "horizon")) |>
  filter(!is.na(observation)) |>
  group_by(model_id, location, origin_date, horizon, target_end_date, observation) |>
  summarise(
    pit = mean(value <= observation),
    .groups = "drop"
  )

p4 <- pit_data |>
  ggplot(aes(x = pit)) +
  geom_histogram(bins = 10, fill = "coral", alpha = 0.8, color = "white") +
  geom_hline(yintercept = nrow(pit_data) / 10, linetype = "dashed", color = "red") +
  facet_wrap(~model_id, scales = "free_y") +
  labs(
    title = "Calibration Assessment (PIT Histograms)",
    subtitle = "Well-calibrated models should have uniform distribution",
    x = "Probability Integral Transform",
    y = "Count",
    caption = "Red line shows expected count for uniform distribution"
  ) +
  theme_minimal()

# 5. Time series plot showing forecasts vs actuals for one location
# Get a specific forecast example
example_date <- as.Date("2016-12-17")
example_location <- "US National"

forecast_example <- validation_forecasts |>
  filter(
    model_id %in% c("sismid-var2", "sismid-arima210", "delphi-epicast"),
    origin_date == example_date,
    location == example_location,
    output_type == "quantile"
  )

actual_data <- oracle_output |>
  filter(
    location == example_location,
    target_end_date >= example_date - 90,
    target_end_date <= example_date + 35
  )

p5 <- ggplot() +
  # Historical data
  geom_line(
    data = actual_data |> filter(target_end_date <= example_date),
    aes(x = target_end_date, y = observation),
    size = 1.2
  ) +
  geom_point(
    data = actual_data |> filter(target_end_date <= example_date),
    aes(x = target_end_date, y = observation),
    size = 2
  ) +
  # Future observations
  geom_point(
    data = actual_data |> filter(target_end_date > example_date),
    aes(x = target_end_date, y = observation),
    size = 3,
    shape = 17,
    color = "red"
  ) +
  # Forecast ribbons
  geom_ribbon(
    data = forecast_example |> 
      filter(output_type_id %in% c(0.1, 0.9)) |>
      select(model_id, target_end_date, output_type_id, value) |>
      pivot_wider(names_from = output_type_id, values_from = value),
    aes(x = target_end_date, ymin = `0.1`, ymax = `0.9`, fill = model_id),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = forecast_example |> 
      filter(output_type_id %in% c(0.25, 0.75)) |>
      select(model_id, target_end_date, output_type_id, value) |>
      pivot_wider(names_from = output_type_id, values_from = value),
    aes(x = target_end_date, ymin = `0.25`, ymax = `0.75`, fill = model_id),
    alpha = 0.5
  ) +
  # Median forecasts
  geom_line(
    data = forecast_example |> filter(output_type_id == 0.5),
    aes(x = target_end_date, y = value, color = model_id),
    size = 1.2
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = paste("Forecast Comparison -", example_location),
    subtitle = paste("Origin date:", example_date),
    x = "Date",
    y = "% ILI visits",
    fill = "Model",
    color = "Model",
    caption = "Red triangles show actual future observations"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 6. Decomposition of WIS
p6 <- model_scores |>
  filter(model_id %in% c("sismid-var2", "sismid-arima210", "delphi-epicast")) |>
  select(model_id, overprediction, underprediction, dispersion) |>
  pivot_longer(cols = c(overprediction, underprediction, dispersion),
               names_to = "component",
               values_to = "value") |>
  group_by(model_id, component) |>
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = model_id, y = mean_value, fill = component)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("overprediction" = "#e41a1c",
                               "underprediction" = "#377eb8",
                               "dispersion" = "#4daf4a")) +
  labs(
    title = "WIS Decomposition",
    subtitle = "Breaking down forecast errors into components",
    x = "Model",
    y = "Mean Score Component",
    fill = "Component"
  ) +
  theme_minimal() +
  coord_flip()

# Create output directory
output_dir <- file.path(dirname(hub_path), "VAR_MODEL_PERFORMANCE_PLOTS")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save all plots
print(paste("\nSaving plots to:", output_dir))

ggsave(file.path(output_dir, "1_overall_wis_comparison.png"), 
       p1, width = 10, height = 6, dpi = 300)

ggsave(file.path(output_dir, "2_wis_by_horizon.png"), 
       p2, width = 10, height = 6, dpi = 300)

ggsave(file.path(output_dir, "3_wis_by_location.png"), 
       p3, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "4_calibration_pit_histograms.png"), 
       p4, width = 12, height = 6, dpi = 300)

ggsave(file.path(output_dir, "5_forecast_example.png"), 
       p5, width = 12, height = 8, dpi = 300)

ggsave(file.path(output_dir, "6_wis_decomposition.png"), 
       p6, width = 10, height = 6, dpi = 300)

# Create a combined dashboard
dashboard <- (p1 | p2) / (p3 | p6) / p5 +
  plot_annotation(
    title = "VAR Model Performance Dashboard",
    subtitle = "Comprehensive comparison of forecasting models on ILI data",
    theme = theme(plot.title = element_text(size = 20, face = "bold"))
  )

ggsave(file.path(output_dir, "0_PERFORMANCE_DASHBOARD.png"), 
       dashboard, width = 20, height = 20, dpi = 300)

# Create summary statistics table
summary_stats <- model_scores |>
  group_by(model_id) |>
  summarise(
    mean_wis = round(mean(wis, na.rm = TRUE), 3),
    median_wis = round(median(wis, na.rm = TRUE), 3),
    mean_bias = round(mean(bias, na.rm = TRUE), 3),
    coverage_50 = round(mean(interval_coverage_50, na.rm = TRUE), 3),
    coverage_90 = round(mean(interval_coverage_90, na.rm = TRUE), 3),
    n_forecasts = n()
  ) |>
  arrange(mean_wis)

write.csv(summary_stats, 
          file.path(output_dir, "model_performance_summary.csv"),
          row.names = FALSE)

print("\n=== Model Performance Summary ===")
print(knitr::kable(summary_stats))

print(paste("\n✓ All plots saved to:", output_dir))
print("\nFiles created:")
print("- 0_PERFORMANCE_DASHBOARD.png (main overview)")
print("- 1_overall_wis_comparison.png")
print("- 2_wis_by_horizon.png")
print("- 3_wis_by_location.png")
print("- 4_calibration_pit_histograms.png")
print("- 5_forecast_example.png")
print("- 6_wis_decomposition.png")
print("- model_performance_summary.csv")