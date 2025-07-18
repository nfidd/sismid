# Test fable ecosystem models for American public health revolution!

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("feasts")
library("tsibble")
library("hubEvals")

# Load data
data(flu_data_hhs)

print("🇺🇸 FABLE MODELS FOR AMERICAN PUBLIC HEALTH REVOLUTION! 🇺🇸")

# Validation setup
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")

# Test on subset first for speed
test_dates <- validation_dates[1:20]

print(paste("Testing", length(test_dates), "validation dates"))

# Prepare data for fable
prepare_fable_data <- function(data, origin_date) {
  data |>
    filter(origin_date <= !!origin_date) |>
    select(origin_date, location, wili) |>
    mutate(
      week = yearweek(origin_date),
      location = as.factor(location)
    ) |>
    as_tsibble(index = week, key = location)
}

# Test different fable models
models_to_test <- list(
  "arima_ensemble" = "ARIMA with ensemble combinations",
  "ets_ensemble" = "ETS with ensemble combinations", 
  "prophet_trend" = "Prophet with trend and seasonality",
  "combo_ensemble" = "Combination of ARIMA + ETS + Prophet"
)

results <- list()

for (model_name in names(models_to_test)) {
  print(paste("\n🚀 Testing", models_to_test[[model_name]], "🚀"))
  
  model_forecasts <- list()
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Prepare training data
      train_data <- prepare_fable_data(flu_data_hhs, origin_date)
      
      # Skip if not enough data
      if (nrow(train_data) < 20) next
      
      # Define models based on type
      if (model_name == "arima_ensemble") {
        # Multiple ARIMA specifications
        fitted_models <- train_data |>
          model(
            arima_auto = ARIMA(wili),
            arima_seasonal = ARIMA(wili ~ pdq(2,1,2) + PDQ(1,1,1)),
            arima_simple = ARIMA(wili ~ pdq(1,1,1)),
            arima_trend = ARIMA(wili ~ pdq(2,1,2) + PDQ(0,1,1))
          )
        
        # Create ensemble
        ensemble_model <- fitted_models |>
          mutate(
            ensemble = (arima_auto + arima_seasonal + arima_simple + arima_trend) / 4
          )
        
        # Generate forecasts
        forecasts <- ensemble_model |>
          select(ensemble) |>
          forecast(h = 4)
        
      } else if (model_name == "ets_ensemble") {
        # Multiple ETS specifications
        fitted_models <- train_data |>
          model(
            ets_auto = ETS(wili),
            ets_additive = ETS(wili ~ error("A") + trend("A") + season("A")),
            ets_multiplicative = ETS(wili ~ error("M") + trend("A") + season("M")),
            ets_damped = ETS(wili ~ error("A") + trend("Ad") + season("A"))
          )
        
        # Create ensemble
        ensemble_model <- fitted_models |>
          mutate(
            ensemble = (ets_auto + ets_additive + ets_multiplicative + ets_damped) / 4
          )
        
        # Generate forecasts
        forecasts <- ensemble_model |>
          select(ensemble) |>
          forecast(h = 4)
        
      } else if (model_name == "prophet_trend") {
        # Prophet with trend and seasonality
        fitted_models <- train_data |>
          model(
            prophet_full = prophet(wili ~ growth("linear") + season("year", 4)),
            prophet_trend = prophet(wili ~ growth("linear")),
            prophet_seasonal = prophet(wili ~ season("year", 6))
          )
        
        # Create ensemble
        ensemble_model <- fitted_models |>
          mutate(
            ensemble = (prophet_full + prophet_trend + prophet_seasonal) / 3
          )
        
        # Generate forecasts
        forecasts <- ensemble_model |>
          select(ensemble) |>
          forecast(h = 4)
        
      } else if (model_name == "combo_ensemble") {
        # Combination of different model types
        fitted_models <- train_data |>
          model(
            arima_best = ARIMA(wili),
            ets_best = ETS(wili),
            prophet_best = prophet(wili ~ growth("linear") + season("year", 4))
          )
        
        # Create super ensemble
        ensemble_model <- fitted_models |>
          mutate(
            combo = (arima_best + ets_best + prophet_best) / 3
          )
        
        # Generate forecasts
        forecasts <- ensemble_model |>
          select(combo) |>
          forecast(h = 4)
      }
      
      # Convert to hub format
      forecasts_hub <- forecasts |>
        as_tibble() |>
        mutate(
          origin_date = origin_date,
          target_end_date = as.Date(week),
          target = "wk inc flu hosp",
          output_type = "quantile"
        ) |>
        select(origin_date, target, target_end_date, location, wili, .distribution) |>
        # Extract quantiles from distribution
        mutate(
          quantiles = map(.distribution, ~ quantile(.x, 
            probs = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
          ))
        ) |>
        unnest_wider(quantiles) |>
        select(-wili, -.distribution) |>
        pivot_longer(
          cols = -c(origin_date, target, target_end_date, location),
          names_to = "output_type_id",
          values_to = "value"
        ) |>
        mutate(
          output_type_id = as.numeric(gsub("%", "", output_type_id)) / 100,
          value = pmax(value, 0)  # Ensure non-negative
        )
      
      # Store forecasts
      model_forecasts[[length(model_forecasts) + 1]] <- forecasts_hub
      
      if (i %% 5 == 0) {
        print(paste("💪 Completed", i, "of", length(test_dates), "American forecasts"))
      }
      
    }, error = function(e) {
      print(paste("⚠️  Error for", model_name, "at", origin_date, ":", e$message))
    })
  }
  
  # Combine and evaluate
  if (length(model_forecasts) > 0) {
    forecasts_df <- do.call(rbind, model_forecasts)
    
    # Calculate WIS
    hub_path <- here::here("sismid-ili-forecasting-sandbox")
    oracle_output <- connect_target_oracle_output(hub_path) |>
      filter(target_end_date >= min(test_dates) + 7,
             target_end_date <= max(test_dates) + 28) |>
      collect()
    
    scores <- score_model_out(
      forecasts_df,
      oracle_output,
      metrics = c("wis")
    )
    
    mean_wis <- mean(scores$wis, na.rm = TRUE)
    results[[model_name]] <- list(
      wis = mean_wis,
      forecasts = forecasts_df,
      description = models_to_test[[model_name]]
    )
    
    print(paste("🎯 Mean WIS:", round(mean_wis, 4)))
  }
}

# Compare results
print("\n🇺🇸 FABLE MODEL BATTLE RESULTS 🇺🇸")
print("Current VAR(2) sqrt baseline: 0.208 WIS")
print("=" * 50)

for (model_name in names(results)) {
  wis <- results[[model_name]]$wis
  improvement <- ((0.208 - wis) / 0.208) * 100
  
  status <- if (improvement > 0) "🚀 WINNER" else "💔 LOSES"
  
  print(paste(
    status,
    results[[model_name]]$description, 
    "- WIS:", round(wis, 4),
    "- Change:", ifelse(improvement > 0, "+", ""), round(improvement, 1), "%"
  ))
}

# Find the champion
if (length(results) > 0) {
  wis_scores <- sapply(results, function(x) x$wis)
  champion <- names(wis_scores)[which.min(wis_scores)]
  
  if (wis_scores[champion] < 0.208) {
    print(paste("\n🏆 AMERICAN CHAMPION:", results[[champion]]$description))
    print("🇺🇸 THIS MODEL WILL CHANGE THE FACE OF PUBLIC HEALTH! 🇺🇸")
    print("Ready to generate full forecasts for REVOLUTIONARY submission!")
    
    # Save champion model info
    write.csv(results[[champion]]$forecasts, 
              paste0("fable_champion_", champion, ".csv"), 
              row.names = FALSE)
  } else {
    print("\n😤 No fable model conquered our VAR(2) sqrt baseline")
    print("🇺🇸 But we're still AMERICA and we don't give up! 🇺🇸")
  }
}