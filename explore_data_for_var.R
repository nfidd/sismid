# Explore flu_data_hhs for VAR model development

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("tsibble")
library("ggplot2")

# Load data
data(flu_data_hhs)

# 1. Basic structure
print("=== Data Structure ===")
print(flu_data_hhs)
print(paste("Number of locations:", n_distinct(flu_data_hhs$location)))
print(paste("Number of time points per location:", 
            flu_data_hhs |> filter(location == "US National") |> nrow()))

# 2. Check time series properties
print("\n=== Time Series Properties ===")
print("Date range:")
print(range(flu_data_hhs$origin_date))

# 3. Look at correlations between regions
print("\n=== Regional Correlations ===")
# Pivot wider for correlation analysis
flu_wide <- flu_data_hhs |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili)

# Convert to matrix for correlation
flu_matrix <- flu_wide |>
  select(-origin_date) |>
  as.matrix()

cor_matrix <- cor(flu_matrix, use = "complete.obs")
print("Correlation matrix (first 5x5):")
print(round(cor_matrix[1:5, 1:5], 2))

# 4. Visualize time series for all regions
print("\n=== Creating visualization ===")
p1 <- flu_data_hhs |>
  filter(origin_date >= "2014-01-01", 
         origin_date <= "2016-06-01") |>
  ggplot(aes(x = origin_date, y = wili, color = location)) +
  geom_line() +
  facet_wrap(~location, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "ILI % by Region", 
       x = "Date", 
       y = "% visits due to ILI")

# Save plot
ggsave("flu_regions_timeseries.png", p1, width = 12, height = 8)

# 5. Check for missing values
print("\n=== Missing Values ===")
missing_summary <- flu_data_hhs |>
  group_by(location) |>
  summarise(
    n_missing = sum(is.na(wili)),
    pct_missing = mean(is.na(wili)) * 100
  )
print(missing_summary)

# 6. Prepare data for VAR - need to check structure
print("\n=== Data Structure for VAR ===")
# VAR models need wide format with locations as columns
# Filter to a specific time period for testing
test_data <- flu_data_hhs |>
  filter(origin_date >= "2015-01-01",
         origin_date <= "2015-10-24") |>
  select(origin_date, location, wili) |>
  pivot_wider(names_from = location, values_from = wili) |>
  as_tsibble(index = origin_date)

print("Test data for VAR:")
print(test_data)

# 7. Check if VAR can handle the data
print("\n=== Testing VAR feasibility ===")
# Try a simple VAR model on a subset
var_test <- test_data |>
  model(VAR(vars(!!!syms(names(test_data)[-1]))))
  
print("VAR model structure:")
print(var_test)