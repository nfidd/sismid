# Simple but powerful fable models for American public health!

library("nfidd")
library("dplyr")
library("tidyr")
library("fable")
library("fabletools")
library("tsibble")
library("hubEvals")

# Load data
data(flu_data_hhs)

print("🇺🇸 SIMPLE FABLE MODELS - MAXIMUM AMERICAN EFFICIENCY! 🇺🇸")

# Test on just 5 validation dates for speed
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- validation_dates[1:5]

print(paste("Testing", length(test_dates), "validation dates"))

# Simple fable models
models_to_test <- list(
  "arima_auto" = "Auto ARIMA",
  "ets_auto" = "Auto ETS",
  "simple_ensemble" = "ARIMA + ETS Ensemble"
)

results <- list()

for (model_name in names(models_to_test)) {
  print(paste("\n💪 Testing", models_to_test[[model_name]]))
  
  model_forecasts <- list()
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Prepare training data
      train_data <- flu_data_hhs |>
        filter(origin_date <= !!origin_date) |>
        select(origin_date, location, wili) |>
        mutate(week = yearweek(origin_date)) |>
        as_tsibble(index = week, key = location)
      
      # Fit model
      if (model_name == "arima_auto") {
        fitted <- train_data |> model(arima = ARIMA(wili))
      } else if (model_name == "ets_auto") {
        fitted <- train_data |> model(ets = ETS(wili))
      } else if (model_name == "simple_ensemble") {
        fitted <- train_data |> 
          model(
            arima = ARIMA(wili),
            ets = ETS(wili)
          ) |>
          mutate(ensemble = (arima + ets) / 2)
        fitted <- fitted |> select(ensemble)
      }
      
      # Generate forecasts
      forecasts <- fitted |> forecast(h = 4)
      
      # Convert to hub format with basic quantiles
      quantiles <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
      
      for (loc in unique(forecasts$location)) {
        loc_forecasts <- forecasts |> filter(location == loc)
        
        for (h in 1:4) {
          target_date <- origin_date + h * 7
          forecast_h <- loc_forecasts |> filter(as.Date(week) == target_date)
          
          if (nrow(forecast_h) > 0) {
            # Extract distribution and get quantiles
            dist_vals <- quantile(forecast_h$.distribution[[1]], quantiles)
            dist_vals <- pmax(dist_vals, 0)  # Ensure non-negative
            
            for (q in seq_along(quantiles)) {
              model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
                origin_date = origin_date,
                target = "wk inc flu hosp",
                target_end_date = target_date,
                location = as.character(loc),
                output_type = "quantile",
                output_type_id = quantiles[q],
                value = dist_vals[q]
              )
            }
          }
        }
      }
      
      print(paste("✅ Completed", i, "of", length(test_dates)))
      
    }, error = function(e) {
      print(paste("⚠️  Error for", model_name, "at", origin_date, ":", e$message))
    })
  }
  
  # Evaluate
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

# Results
print("\n🇺🇸 AMERICAN FABLE RESULTS 🇺🇸")
print("VAR(2) sqrt baseline: 0.208 WIS")

for (model_name in names(results)) {
  wis <- results[[model_name]]$wis
  improvement <- ((0.208 - wis) / 0.208) * 100
  
  status <- if (improvement > 0) "🚀 BEATS BASELINE" else "💔 LOSES"
  
  print(paste(
    status, "-",
    results[[model_name]]$description, 
    "- WIS:", round(wis, 4),
    "- Change:", ifelse(improvement > 0, "+", ""), round(improvement, 1), "%"
  ))
}

# Find champion
if (length(results) > 0) {
  wis_scores <- sapply(results, function(x) x$wis)
  champion <- names(wis_scores)[which.min(wis_scores)]
  
  if (wis_scores[champion] < 0.208) {
    print(paste("\n🏆 CHAMPION:", results[[champion]]$description))
    print("🇺🇸 READY FOR AMERICAN PUBLIC HEALTH REVOLUTION! 🇺🇸")
  } else {
    print("\n🇺🇸 No model beat the baseline, but AMERICA NEVER GIVES UP! 🇺🇸")
  }
}