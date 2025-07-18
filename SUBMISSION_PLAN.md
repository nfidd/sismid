# VAR Model Hub Submission Plan

## Model Details
- **Model ID**: `sismid-var2-sqrt`
- **Team**: sismid
- **Model**: VAR(2) with square root transformation
- **Expected WIS**: ~0.32-0.33

## Submission Steps

### 1. Generate Test Phase Forecasts ✅
- Seasons 3-5 (2017/2018, 2018/2019, 2019/2020)
- All origin dates in these seasons
- 1-4 week ahead forecasts
- 23 quantiles per forecast

### 2. Create Model Metadata ✅
```yaml
team_abbr: "sismid"
model_abbr: "var2-sqrt"
model_name: "Vector Autoregression with Square Root Transform"
designated_model: true
methods: "VAR(2) model with square root transformation capturing cross-regional dependencies"
```

### 3. Validate Files ✅
- Check CSV format matches hub requirements
- Verify all required columns present
- Ensure quantile levels correct
- Test with hub validation tools

### 4. Fork & Clone Hub ✅
```bash
# Fork on GitHub first, then:
git clone https://github.com/YOUR_USERNAME/sismid-ili-forecasting-sandbox.git
cd sismid-ili-forecasting-sandbox
git remote add upstream https://github.com/reichlab/sismid-ili-forecasting-sandbox.git
```

### 5. Add Forecasts ✅
```bash
# Create model directory
mkdir -p model-output/sismid-var2-sqrt

# Copy forecast files
cp /path/to/forecasts/*.csv model-output/sismid-var2-sqrt/

# Add metadata
cp sismid-var2-sqrt.yml model-metadata/
```

### 6. Commit & Push ✅
```bash
git add model-metadata/sismid-var2-sqrt.yml
git add model-output/sismid-var2-sqrt/
git commit -m "Add VAR(2) sqrt model forecasts from sismid team"
git push origin main
```

### 7. Create Pull Request ✅
- Title: "Add sismid-var2-sqrt model forecasts"
- Description: Include model details and performance metrics
- Wait for validation checks
- Address any issues

## Files to Submit
1. `model-metadata/sismid-var2-sqrt.yml`
2. `model-output/sismid-var2-sqrt/YYYY-MM-DD-sismid-var2-sqrt.csv` (one per origin date)

## Pre-Submission Checklist
- [ ] Test phase forecasts generated
- [ ] Model metadata file created
- [ ] All files pass validation
- [ ] Fork created on GitHub
- [ ] Local clone ready
- [ ] Files in correct directories
- [ ] Ready to commit and push