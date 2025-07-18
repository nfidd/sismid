# American boost to VAR model - simple but effective!

library("nfidd")
library("dplyr")
library("vars")
library("hubEvals")

# Load data
data(flu_data_hhs)

print("🇺🇸 AMERICAN VAR MODEL BOOST - SIMPLE BUT POWERFUL! 🇺🇸")

# Test variations of our successful VAR approach
validation_dates <- seq(as.Date("2016-11-05"), as.Date("2017-10-28"), by = "week")
test_dates <- validation_dates[1:10]  # Quick test

models_to_test <- list(
  "VAR2_sqrt_boosted" = "VAR(2) sqrt with bootstrap boost",
  "VAR2_sqrt_weighted" = "VAR(2) sqrt with weighted regions",
  "VAR2_sqrt_ensemble" = "VAR(2) sqrt ensemble of 3 models"
)

results <- list()

for (model_name in names(models_to_test)) {
  print(paste("\n💪 Testing", models_to_test[[model_name]]))
  
  model_forecasts <- list()
  
  for (i in seq_along(test_dates)) {
    origin_date <- test_dates[i]
    
    tryCatch({
      # Get training data
      train_data <- flu_data_hhs |>
        filter(origin_date <= !!origin_date) |>
        select(origin_date, location, wili) |>
        pivot_wider(names_from = location, values_from = wili) |>
        arrange(origin_date)
      
      # Apply sqrt transformation
      numeric_cols <- setdiff(names(train_data), "origin_date")
      train_data_sqrt <- train_data
      train_data_sqrt[numeric_cols] <- lapply(train_data_sqrt[numeric_cols], sqrt)
      
      # Create time series
      ts_data <- ts(train_data_sqrt[numeric_cols], frequency = 52)
      
      if (model_name == "VAR2_sqrt_boosted") {
        # Regular VAR(2) with enhanced bootstrap
        var_model <- VAR(ts_data, p = 2, type = "const")
        
        # Generate base forecast
        forecast_obj <- predict(var_model, n.ahead = 4)
        
        # Enhanced bootstrap with more samples
        n_boot <- 2000  # More bootstrap samples
        
        for (loc in numeric_cols) {
          pred_values <- forecast_obj$fcst[[loc]][, "fcst"]
          pred_values <- pmax(pred_values^2, 0)  # Transform back
          
          for (h in 1:4) {
            target_date <- origin_date + h * 7
            
            # Enhanced bootstrap
            boot_preds <- replicate(n_boot, {
              # Sample residuals with replacement
              resid_sample <- sample(residuals(var_model)[, loc], 
                                     size = 2, replace = TRUE)
              
              # Add noise based on recent volatility
              recent_volatility <- sd(tail(residuals(var_model)[, loc], 20))
              noise <- rnorm(1, 0, recent_volatility * 0.8)
              
              # Combine
              pred_val <- pred_values[h] + noise
              max(pred_val, 0)
            })
            
            # Calculate quantiles
            quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
            pred_quantiles <- quantile(boot_preds, quantiles)
            
            for (q in seq_along(quantiles)) {
              model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
                origin_date = origin_date,
                target = "wk inc flu hosp",
                target_end_date = target_date,
                location = loc,
                output_type = "quantile",
                output_type_id = quantiles[q],
                value = pred_quantiles[q]
              )
            }
          }
        }
        
      } else if (model_name == "VAR2_sqrt_weighted") {
        # Weight regions by population or importance
        region_weights <- c(
          "US National" = 1.0,
          "HHS Region 1" = 0.9,
          "HHS Region 2" = 0.9,
          "HHS Region 3" = 0.9,
          "HHS Region 4" = 1.0,
          "HHS Region 5" = 0.9,
          "HHS Region 6" = 0.9,
          "HHS Region 7" = 0.8,
          "HHS Region 8" = 0.8,
          "HHS Region 9" = 0.9,
          "HHS Region 10" = 0.8
        )
        
        # Fit weighted VAR
        var_model <- VAR(ts_data, p = 2, type = "const")
        forecast_obj <- predict(var_model, n.ahead = 4)
        
        for (loc in numeric_cols) {
          weight <- region_weights[loc]
          pred_values <- forecast_obj$fcst[[loc]][, "fcst"]
          pred_values <- pmax(pred_values^2, 0)
          
          for (h in 1:4) {
            target_date <- origin_date + h * 7
            
            # Apply weight to uncertainty
            n_boot <- 1500
            boot_preds <- replicate(n_boot, {
              resid_sample <- sample(residuals(var_model)[, loc], 
                                     size = 2, replace = TRUE)
              # Weight affects uncertainty
              weighted_noise <- rnorm(1, 0, sd(resid_sample) * weight)
              pred_val <- pred_values[h] + weighted_noise
              max(pred_val, 0)
            })
            
            quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
            pred_quantiles <- quantile(boot_preds, quantiles)
            
            for (q in seq_along(quantiles)) {
              model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
                origin_date = origin_date,
                target = "wk inc flu hosp",
                target_end_date = target_date,
                location = loc,
                output_type = "quantile",
                output_type_id = quantiles[q],
                value = pred_quantiles[q]
              )
            }
          }
        }
        
      } else if (model_name == "VAR2_sqrt_ensemble") {
        # Ensemble of 3 slightly different VAR models
        models <- list(
          VAR(ts_data, p = 2, type = "const"),
          VAR(ts_data, p = 2, type = "trend"),
          VAR(ts_data, p = 2, type = "both")
        )
        
        forecasts <- lapply(models, function(m) predict(m, n.ahead = 4))
        
        for (loc in numeric_cols) {
          # Average the three model predictions
          pred_values <- (forecasts[[1]]$fcst[[loc]][, "fcst"] + 
                          forecasts[[2]]$fcst[[loc]][, "fcst"] + 
                          forecasts[[3]]$fcst[[loc]][, "fcst"]) / 3
          pred_values <- pmax(pred_values^2, 0)
          
          for (h in 1:4) {
            target_date <- origin_date + h * 7
            
            # Bootstrap from ensemble
            n_boot <- 1500
            boot_preds <- replicate(n_boot, {
              # Sample from each model's residuals
              resid_samples <- sapply(models, function(m) {
                sample(residuals(m)[, loc], size = 1, replace = TRUE)
              })
              
              # Average residual noise
              avg_noise <- mean(resid_samples)
              pred_val <- pred_values[h] + avg_noise
              max(pred_val, 0)
            })
            
            quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
            pred_quantiles <- quantile(boot_preds, quantiles)
            
            for (q in seq_along(quantiles)) {
              model_forecasts[[length(model_forecasts) + 1]] <- data.frame(
                origin_date = origin_date,
                target = "wk inc flu hosp",
                target_end_date = target_date,
                location = loc,
                output_type = "quantile",
                output_type_id = quantiles[q],
                value = pred_quantiles[q]
              )
            }
          }
        }
      }
      
      if (i %% 3 == 0) {
        print(paste("✅ Completed", i, "of", length(test_dates)))
      }
      
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
print("\n🇺🇸 AMERICAN VAR BOOST RESULTS 🇺🇸")
print("Original VAR(2) sqrt baseline: 0.208 WIS")

for (model_name in names(results)) {
  wis <- results[[model_name]]$wis
  improvement <- ((0.208 - wis) / 0.208) * 100
  
  status <- if (improvement > 0) "🚀 IMPROVEMENT" else "💔 WORSE"
  
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
    print(paste("\n🏆 AMERICAN CHAMPION:", results[[champion]]$description))
    print("🇺🇸 READY FOR SECOND SUBMISSION! 🇺🇸")
    
    # Save champion forecasts
    write.csv(results[[champion]]$forecasts, 
              paste0("american_champion_", champion, ".csv"), 
              row.names = FALSE)
  } else {
    print("\n🇺🇸 No boost beat the baseline - the original is still CHAMPION! 🇺🇸")
  }
}