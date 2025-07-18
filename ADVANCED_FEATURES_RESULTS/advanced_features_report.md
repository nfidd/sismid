# Advanced VAR Features Test Results

## Best Model: VAR(2) sqrt
AICc: -13397.48

## Tested Features:
1. VARIMA with differencing
2. Seasonal components (sin/cos)
3. Fourier terms for seasonality
4. Regional subsets
5. Trend component
6. Box-Cox transformation
7. Higher lag orders

## Key Findings:
- Base model performs well without seasonality
- Level models sufficient
- Square root transformation remains optimal

## Recommendation:
Use VAR(2) sqrt for final forecasts.

