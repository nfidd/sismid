# Evaluate test phase WIS for sismid-var1-optimal model

library("nfidd")
library("dplyr")
library("hubUtils")
library("hubData")
library("hubEvals")

# Set hub path
hub_path <- here::here("sismid-ili-forecasting-sandbox")

# Connect to hub and get our model data
hub_con <- connect_hub(hub_path)
var1_forecasts <- hub_con |> 
  filter(model_id == "sismid-var1-optimal") |> 
  collect_hub()

# Calculate origin_date from target_end_date and horizon
var1_forecasts <- var1_forecasts |> 
  mutate(
    origin_date = target_end_date - horizon * 7
  )

# Check test phase origin dates
test_origin_dates <- hub_path |> 
  read_config("tasks") |> 
  get_round_ids()

test_dates <- test_origin_dates[which(as.Date(test_origin_dates) > as.Date("2017-05-06"))]

# Filter for test phase
var1_test_forecasts <- var1_forecasts |> 
  filter(origin_date %in% test_dates)

print(paste("VAR1 test phase forecasts:", nrow(var1_test_forecasts)))
print(paste("Test phase origin dates:", length(unique(var1_test_forecasts$origin_date))))

# Get oracle output
oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

# Score our model using hubEvals
if (nrow(var1_test_forecasts) > 0) {
  
  # Calculate WIS and other metrics
  model_scores <- score_model_out(
    var1_test_forecasts,
    oracle_output
  )
  
  print("")
  print("=== TEST PHASE EVALUATION RESULTS ===")
  print("Model: sismid-var1-optimal (VAR(1) + season(52))")
  print("")
  print(model_scores)
  
  # Save results
  write.csv(model_scores, "test_phase_wis_results.csv", row.names = FALSE)
  
  # Extract key metrics
  wis_score <- model_scores$wis[1]
  coverage_50 <- model_scores$interval_coverage_50[1] 
  coverage_90 <- model_scores$interval_coverage_90[1]
  
  print("")
  print("=== KEY PERFORMANCE METRICS ===")
  print(paste("WIS Score:", round(wis_score, 3)))
  print(paste("50% Coverage:", round(coverage_50 * 100, 1), "%"))
  print(paste("90% Coverage:", round(coverage_90 * 100, 1), "%"))
  
  # Performance context
  print("")
  print("=== PERFORMANCE CONTEXT ===")
  print("Expected comparison (from training analysis):")
  print("- delphi-epicast: 0.31 WIS")
  print("- sismid-var2-sqrt: 0.34 WIS")
  print("- sismid-arima210: 0.40 WIS")
  print("- hist-avg: 0.45 WIS")
  print("")
  
  if (wis_score < 0.35) {
    print("EXCELLENT! Competitive with top models!")
  } else if (wis_score < 0.45) {
    print("GOOD! Better than historical average!")
  } else if (wis_score < 0.60) {
    print("ACCEPTABLE performance")
  } else {
    print("Performance needs investigation")
  }
  
} else {
  print("No test forecasts found for sismid-var1-optimal model")
}