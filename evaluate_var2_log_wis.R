# Evaluate VAR(2) log model WIS performance

library("nfidd")
library("dplyr")
library("hubUtils")
library("hubData")
library("hubEvals")

# Load our generated forecasts
forecasts <- read.csv("var2_log_forecasts.csv")

print("=== Evaluating VAR(2) Log Model ===")
print(paste("Loaded", nrow(forecasts), "forecast rows"))

# Prepare forecasts in correct format for hubEvals
forecasts_formatted <- forecasts |>
  mutate(
    origin_date = as.Date(origin_date),
    target_end_date = as.Date(target_end_date)
  ) |>
  # Add missing columns that hubEvals expects
  select(model_id, location, origin_date, horizon, target_end_date, 
         target, output_type, output_type_id, value)

print("Forecast format:")
print(names(forecasts_formatted))

# Get hub path and oracle output
hub_path <- here::here("sismid-ili-forecasting-sandbox")
oracle_output <- connect_target_oracle_output(hub_path) |>
  collect()

print("Oracle output columns:")
print(names(oracle_output))

# The issue might be that oracle_output has extra columns
# Let's try filtering to match what hubEvals expects
oracle_clean <- oracle_output |>
  select(location, target_end_date, oracle_value) |>
  # Remove the target column that was causing issues
  mutate(
    output_type = "quantile",
    output_type_id = as.numeric(NA)
  )

print("Cleaned oracle output:")
print(head(oracle_clean))

# Try evaluation
tryCatch({
  model_scores <- score_model_out(forecasts_formatted, oracle_clean)
  
  print("=== VAR(2) LOG MODEL PERFORMANCE ===")
  print(model_scores)
  
  # Compare to current model
  current_wis <- 0.466
  new_wis <- model_scores$wis[1]
  
  print("")
  print(paste("Current VAR(1) + season(52) WIS:", current_wis))
  print(paste("New VAR(2) log WIS:", round(new_wis, 3)))
  
  if (new_wis < current_wis) {
    improvement <- (current_wis - new_wis) / current_wis * 100
    print(paste("*** IMPROVEMENT:", round(improvement, 1), "% better! ***"))
    print("Nick Reich was right - simpler models work better!")
  } else {
    decline <- (new_wis - current_wis) / current_wis * 100
    print(paste("Decline:", round(decline, 1), "% worse"))
  }
  
  # Save results
  write.csv(model_scores, "var2_log_final_evaluation.csv", row.names = FALSE)
  
}, error = function(e) {
  print(paste("Evaluation error:", e$message))
  print("Let's try a different approach...")
  
  # Alternative: Use the same approach as the working evaluation
  # Load existing hub forecasts to compare format
  hub_con <- connect_hub(hub_path)
  existing_forecasts <- hub_con |> 
    filter(model_id == "hist-avg") |> 
    head(100) |>
    collect_hub()
  
  print("Existing hub forecast format:")
  print(names(existing_forecasts))
  print("Sample:")
  print(head(existing_forecasts))
  
  # Try to match this format
  forecasts_hub_style <- forecasts_formatted |>
    # Add any missing columns that hub forecasts have
    mutate(
      # Fill in any missing standard columns
      horizon = ifelse(is.na(horizon), 1, horizon)
    )
  
  print("Trying hub-style evaluation...")
  
  # Try evaluation again
  model_scores_v2 <- score_model_out(forecasts_hub_style, oracle_output)
  
  print("=== VAR(2) LOG MODEL PERFORMANCE (V2) ===")
  print(model_scores_v2)
  
  current_wis <- 0.466
  new_wis <- model_scores_v2$wis[1]
  
  print("")
  print(paste("Current VAR(1) + season(52) WIS:", current_wis))
  print(paste("New VAR(2) log WIS:", round(new_wis, 3)))
  
  if (new_wis < current_wis) {
    improvement <- (current_wis - new_wis) / current_wis * 100
    print(paste("*** IMPROVEMENT:", round(improvement, 1), "% better! ***"))
    print("Simpler VAR(2) log model outperforms complex seasonal model!")
  } else {
    decline <- (new_wis - current_wis) / current_wis * 100
    print(paste("Decline:", round(decline, 1), "% worse"))
  }
  
  write.csv(model_scores_v2, "var2_log_final_evaluation.csv", row.names = FALSE)
})

print("\n=== Summary ===")
print("VAR(2) with log transformation:")
print("- Simple multivariate model")
print("- Captures cross-regional dependencies") 
print("- No forced seasonality")
print("- Much better AICc than seasonal model")
print("- Following Nick Reich's guidance on model simplicity")