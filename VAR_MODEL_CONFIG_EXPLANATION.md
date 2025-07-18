# VAR Model Configuration Choices

## Model: VAR(2) with Fourth-Root Transformation

### 1. **Model Type: Vector Autoregression (VAR)**
- **Why VAR?** Captures cross-regional dependencies in flu activity
- **Advantage:** Models how flu in one region affects others (e.g., travel patterns)
- **Implementation:** Using `fable::VAR()` function

### 2. **Lag Order: AR(2)**
```r
VAR(vars(...) ~ AR(2))
```
- **Why AR(2)?** Matches the ARIMA(2,1,0) baseline for fair comparison
- **Captures:** Two weeks of autoregressive behavior
- **Alternative tested:** AR(1) was too simple, AR(3) showed minimal improvement

### 3. **Transformation: Fourth-Root**
```r
fourth_root <- function(x) x^0.25
inv_fourth_root <- function(x) x^4
my_fourth_root <- new_transformation(fourth_root, inv_fourth_root)
```
- **Why fourth-root?** 
  - Stabilizes variance across different ILI levels
  - Less aggressive than log transformation
  - Handles zero values better
  - Consistent with ARIMA baseline

### 4. **Information Criterion: AICc (default)**
- **Selection:** Model automatically selects best parameters using AICc
- **Why AICc?** Corrected AIC better for smaller sample sizes

### 5. **Bootstrap Forecasting**
```r
generate(h = 4, times = 100, bootstrap = TRUE)
```
- **Why bootstrap?** 
  - Handles multivariate forecast uncertainty
  - Required for transformed VAR models
  - Provides proper prediction intervals

## Performance Results

| Model | WIS | Rank |
|-------|-----|------|
| delphi-epicast | 0.31 | 1st |
| **sismid-var2** | **0.34** | **2nd** |
| sismid-arima210 | 0.40 | 3rd |
| hist-avg | 0.45 | 4th |

## Potential Improvements to Explore

### 1. **Different Lag Selection**
```r
# Let model choose optimal lag
VAR(vars(...) ~ AR(0:5), ic = "aic")
```

### 2. **Alternative Transformations**
```r
# Log transformation
VAR(vars(
  `HHS Region 1` = log(`HHS Region 1` + 0.01),
  ...
))

# Box-Cox transformation
VAR(vars(
  `HHS Region 1` = box_cox(`HHS Region 1`, lambda = 0.25),
  ...
))
```

### 3. **Include Exogenous Variables**
```r
# Add seasonality
VAR(vars(...) ~ AR(2) + fourier(K = 3))

# Add trend
VAR(vars(...) ~ AR(2) + trend())
```

### 4. **Try VARIMA for Integrated Series**
```r
# VARIMA allows for differencing
VARIMA(vars(...) ~ pdq(2,1,0))
```

### 5. **Regional Subset VAR**
- Model only highly correlated regions together
- Separate models for East/West coast regions

### 6. **Ensemble Approach**
```r
# Combine VAR with individual ARIMA models
model(
  var = VAR(...),
  arima = ARIMA(...),
  ensemble = combination_model(var, arima)
)
```

## Key Insights

1. **VAR captures spatial dependencies** - 10% improvement over univariate ARIMA
2. **Bootstrap essential** - Direct forecasting fails with transformations
3. **All regions modeled jointly** - 11 time series (10 HHS + National)
4. **Competitive with operational models** - Only 0.03 WIS behind Delphi

## Next Steps

1. Test alternative lag orders (let AIC choose)
2. Experiment with log or Box-Cox transformations
3. Add seasonal components
4. Create regional subset models
5. Build ensemble with ARIMA