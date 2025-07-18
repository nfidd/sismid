# Properly test seasonality in flu data

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("feasts")
library("ggplot2")
library("lubridate")

# Set seed
set.seed(406)

# Load data
data(flu_data_hhs)

# First, let's visualize the seasonality
print("=== Examining Seasonality in Flu Data ===")

# Plot seasonal patterns
p1 <- flu_data_hhs |>
  filter(location == "US National") |>
  gg_season(wili, period = "year") +
  labs(title = "Seasonal Pattern in US National ILI%",
       y = "% ILI visits")

ggsave("flu_seasonal_pattern.png", p1, width = 10, height = 6)

# STL decomposition to see seasonal component
p2 <- flu_data_hhs |>
  filter(location == "US National") |>
  model(STL(wili ~ season(window = "periodic"))) |>
  components() |>
  autoplot() +
  labs(title = "STL Decomposition of US National ILI%")

ggsave("flu_stl_decomposition.png", p2, width = 10, height = 8)

# Now test models with proper seasonality
test_date <- as.Date("2016-12-17")

# Prepare data
flu_data_wide <- flu_data_hhs |>
  filter(origin_date <= test_date) |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

# Add proper seasonal indicators
flu_data_seasonal <- flu_data_wide |>
  mutate(
    # Week of year (1-52)
    week = week(origin_date),
    # Month for monthly seasonality
    month = month(origin_date),
    # Season indicators
    winter = ifelse(month %in% c(12, 1, 2), 1, 0),
    spring = ifelse(month %in% c(3, 4, 5), 1, 0),
    summer = ifelse(month %in% c(6, 7, 8), 1, 0),
    fall = ifelse(month %in% c(9, 10, 11), 1, 0),
    # Harmonic terms for smooth seasonality
    sin1 = sin(2 * pi * week / 52),
    cos1 = cos(2 * pi * week / 52),
    sin2 = sin(4 * pi * week / 52),
    cos2 = cos(4 * pi * week / 52),
    sin3 = sin(6 * pi * week / 52),
    cos3 = cos(6 * pi * week / 52)
  )

# Define transformations
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

print("\n=== Testing VAR Models with Seasonality ===")

# 1. Baseline VAR(2) sqrt
print("\n1. Baseline: VAR(2) sqrt (no seasonality)")
var_base <- flu_data_wide |>
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
base_glance <- glance(var_base)
print(base_glance)

# 2. VAR with harmonic seasonality (K=1)
print("\n2. VAR(2) sqrt + harmonic seasonality (K=1)")
var_harm1 <- flu_data_seasonal |>
  model(
    harm1 = VAR(vars(
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
    ) ~ AR(2) + sin1 + cos1)
  )
harm1_glance <- glance(var_harm1)
print(harm1_glance)

# 3. VAR with more harmonics (K=2)
print("\n3. VAR(2) sqrt + harmonic seasonality (K=2)")
var_harm2 <- flu_data_seasonal |>
  model(
    harm2 = VAR(vars(
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
    ) ~ AR(2) + sin1 + cos1 + sin2 + cos2)
  )
harm2_glance <- glance(var_harm2)
print(harm2_glance)

# 4. VAR with seasonal dummies
print("\n4. VAR(2) sqrt + seasonal dummies")
var_seasonal <- flu_data_seasonal |>
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
    ) ~ AR(2) + winter + spring + summer)  # fall is baseline
  )
seasonal_glance <- glance(var_seasonal)
print(seasonal_glance)

# 5. Test with higher AR to see if it captures seasonality
print("\n5. VAR(4) sqrt (higher lags to capture seasonality)")
var_ar4 <- flu_data_wide |>
  model(
    ar4 = VAR(vars(
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
    ) ~ AR(4))
  )
ar4_glance <- glance(var_ar4)
print(ar4_glance)

# Compare results
results <- bind_rows(
  base_glance |> mutate(model = "VAR(2) baseline"),
  harm1_glance |> mutate(model = "VAR(2) + harmonics K=1"),
  harm2_glance |> mutate(model = "VAR(2) + harmonics K=2"),
  seasonal_glance |> mutate(model = "VAR(2) + season dummies"),
  ar4_glance |> mutate(model = "VAR(4) baseline")
) |>
  select(model, AICc, AIC, BIC) |>
  arrange(AICc)

print("\n=== Model Comparison ===")
print(results)

# Calculate improvement
improvement <- (results$AICc[1] - base_glance$AICc) / abs(base_glance$AICc) * 100
print(paste("\nBest model improves AICc by", round(improvement, 2), "% over baseline"))

# Test forecast accuracy with seasonality
print("\n=== Testing Forecast Quality ===")

# Generate forecasts from best model
if (grepl("harmonic", results$model[1])) {
  if (grepl("K=1", results$model[1])) {
    best_model <- var_harm1
  } else {
    best_model <- var_harm2
  }
} else if (grepl("season", results$model[1])) {
  best_model <- var_seasonal
} else if (grepl("VAR\\(4\\)", results$model[1])) {
  best_model <- var_ar4
} else {
  best_model <- var_base
}

# Generate a test forecast
fc <- forecast(best_model, h = 4, bootstrap = TRUE)
print("Forecast generated successfully")

# Visualize improvement
p3 <- results |>
  ggplot(aes(x = reorder(model, -AICc), y = AICc)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  labs(title = "VAR Model Comparison with Seasonality",
       subtitle = "Lower AICc is better",
       x = "Model", y = "AICc") +
  theme_minimal()

ggsave("var_seasonality_comparison.png", p3, width = 10, height = 6)

# Save results
write.csv(results, "var_seasonality_results.csv", row.names = FALSE)

print("\n=== CONCLUSION ===")
print(paste("Best model:", results$model[1]))
if (results$model[1] != "VAR(2) baseline") {
  print("SEASONALITY DOES HELP! You were right to question this.")
  print("The model benefits from explicit seasonal components.")
} else {
  print("Baseline VAR(2) is sufficient - it captures seasonality through lags.")
}