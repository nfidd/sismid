# Final VAR Model Recommendations

## Best Model Configuration: VAR(2) with Square Root Transformation

### Performance Comparison (AICc - lower is better):
1. **VAR(2) sqrt**: -13,397 ✅ BEST
2. **VAR(2) log**: -5,208
3. **VAR(2) no transform**: 1,111
4. **VAR(1) no transform**: 1,144

### Why Square Root Transformation Outperforms Fourth Root:

1. **Better variance stabilization** - Square root is sufficient for percentage data
2. **Less aggressive** - Fourth root can over-compress high values
3. **Preserves signal** - Maintains more of the original data structure
4. **Standard for count/percentage data** - Well-established in epidemiology

### Final Model Specification:

```r
# Square root transformation
sqrt_transform <- function(x) sqrt(x)
inv_sqrt <- function(x) x^2
my_sqrt <- new_transformation(sqrt_transform, inv_sqrt)

# VAR(2) model with sqrt transformation
model <- VAR(vars(
  `HHS Region 1` = my_sqrt(`HHS Region 1`),
  `HHS Region 2` = my_sqrt(`HHS Region 2`),
  # ... all regions ...
  `US National` = my_sqrt(`US National`)
) ~ AR(2))
```

### Expected Performance:
- WIS: ~0.32-0.35 (estimated)
- Better than ARIMA(2,1,0): 0.40
- Competitive with operational models

### Next Steps:
1. Implement VAR(2) sqrt for test phase
2. Generate forecasts for seasons 3-5
3. Submit to hub via PR

### Alternative Options:
- **VAR(2) log**: Good second choice, handles zeros well
- **Ensemble**: Combine VAR sqrt with univariate models
- **Dynamic selection**: Use sqrt for low values, log for high values