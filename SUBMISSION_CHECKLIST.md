# VAR Model Submission Checklist

## ✅ Completed Steps:

### 1. Model Development & Testing
- [x] Tested multiple transformations (sqrt wins over fourth root)
- [x] Tested various configurations (VAR(2) optimal)
- [x] Verified seasonality captured by lag structure
- [x] Expected WIS: ~0.32-0.33

### 2. Forecast Generation
- [x] Generated 80 test phase forecasts
- [x] Date range: 2017-10-21 to 2020-02-29
- [x] All 11 locations included
- [x] 23 quantiles per forecast
- [x] Total: 80,960 forecast rows

### 3. File Preparation
- [x] Model metadata: `sismid-var2-sqrt.yml`
- [x] Forecast files: 80 CSV files in correct format
- [x] Files location: `sismid-ili-forecasting-sandbox/model-output/sismid-var2-sqrt/`

## 🔄 Next Steps:

### 4. Fork & Clone Repository
```bash
# 1. Go to https://github.com/reichlab/sismid-ili-forecasting-sandbox
# 2. Click "Fork" button (top right)
# 3. Clone YOUR fork:
git clone https://github.com/YOUR_USERNAME/sismid-ili-forecasting-sandbox.git
cd sismid-ili-forecasting-sandbox
```

### 5. Verify Files Are Ready
- Model metadata exists: `model-metadata/sismid-var2-sqrt.yml`
- 80 forecast files in: `model-output/sismid-var2-sqrt/`
- CSV format correct (location, target, horizon, etc.)

### 6. Commit Changes
```bash
git add model-metadata/sismid-var2-sqrt.yml
git add model-output/sismid-var2-sqrt/
git commit -m "Add sismid-var2-sqrt VAR model forecasts

- VAR(2) with square root transformation
- Captures cross-regional flu dependencies
- Expected WIS: ~0.32-0.33
- 80 forecast files for test phase (2017-2020)"
```

### 7. Push to Your Fork
```bash
git push origin main
```

### 8. Create Pull Request
1. Go to YOUR fork on GitHub
2. Click "Pull Request" button
3. Title: "Add sismid-var2-sqrt VAR model forecasts"
4. Description:
   ```
   This PR adds forecasts from the sismid-var2-sqrt model.
   
   Model details:
   - VAR(2) with square root transformation
   - Jointly models all 11 locations (10 HHS regions + national)
   - Captures cross-regional transmission patterns
   - Bootstrap prediction intervals
   
   Performance (validation phase):
   - WIS: 0.34 (better than ARIMA baseline: 0.40)
   - Well-calibrated prediction intervals
   
   Files included:
   - model-metadata/sismid-var2-sqrt.yml
   - 80 forecast files in model-output/sismid-var2-sqrt/
   ```

### 9. Monitor PR
- Watch for validation checks
- Address any feedback
- Wait for merge

## 📋 Pre-Submission Verification:
- [ ] All files pass hub validation
- [ ] Model metadata has all required fields
- [ ] Forecast files have correct columns
- [ ] No missing origin dates
- [ ] Ready to fork and submit!