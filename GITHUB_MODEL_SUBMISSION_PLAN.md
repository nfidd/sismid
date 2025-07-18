# GitHub Model Submission Plan for SISMID ILI Forecasting Hub

## Overview
This plan outlines the complete workflow for submitting a custom model to the SISMID ILI Forecasting Hub using GitHub's fork-and-pull-request workflow.

## Phase 1: Repository Setup and Forking

### 1.1 Fork the Main Hub Repository
```bash
# Navigate to the main hub repository (to be confirmed)
# Fork via GitHub web interface or CLI
gh repo fork reichlab/sismid-ili-forecasting-hub --clone=true
cd sismid-ili-forecasting-hub
```

### 1.2 Set Up Local Development Environment
```bash
# Add upstream remote
git remote add upstream https://github.com/reichlab/sismid-ili-forecasting-hub.git

# Create a new branch for your model
git checkout -b add-custom-model-[MODEL_NAME]

# Install R dependencies
Rscript -e "pak::pak(c('hubUtils', 'hubData', 'hubEvals', 'arrow'))"
```

## Phase 2: Model Development and Integration

### 2.1 Create Model Metadata File
Create `/model-metadata/[team_abbr]-[model_abbr].yml`:

```yaml
team_abbr: "sismid"
model_abbr: "[your_model_name]"
designated_model: true
model_name: "[Full Model Name]"
model_version: "1.0"
model_contributors: [
  {
    "name": "Your Name",
    "affiliation": "Your Institution",
    "email": "your.email@institution.edu"
  }
]
website_url: "https://github.com/your-username/your-model-repo"
license: "CC-BY-4.0"
methods: "Brief description of your forecasting method"
data_inputs: "Data sources used (e.g., CDC ILI surveillance data)"
methods_long: |
  Detailed description of your model including:
  - Statistical/mathematical approach
  - Data preprocessing steps
  - Model fitting procedures
  - Uncertainty quantification methods
```

### 2.2 Implement Model Code
Create model implementation in `/model-code/[team_abbr]-[model_abbr]/`:

```r
# main.R - Entry point for model execution
library(hubUtils)
library(dplyr)
library(tidyr)

# Load model-specific functions
source("model_functions.R")

# Read hub configuration
hub_config <- hubUtils::read_config("../../hub-config/")

# Load target data
target_data <- read.csv("../../target-data/target-data.csv")

# Generate forecasts for all required combinations
forecasts <- generate_forecasts(
  target_data = target_data,
  hub_config = hub_config,
  model_name = "[team_abbr]-[model_abbr]"
)

# Save forecasts in hubverse format
write.csv(forecasts, "../../model-output/[team_abbr]-[model_abbr]/forecasts.csv", 
          row.names = FALSE)
```

### 2.3 Create Model Output Structure
```bash
# Create model output directory
mkdir -p model-output/[team_abbr]-[model_abbr]

# Generate forecasts for all required task combinations:
# - origin_date: All dates from 2015-2020 (from tasks.json)
# - horizon: 1, 2, 3, 4 weeks ahead
# - location: US National + 10 HHS regions
# - target: "ili perc" (Weighted ILI percentage)
# - output_type: 23 quantiles (0.01, 0.025, ..., 0.99)
```

## Phase 3: Forecast Generation and Validation

### 3.1 Required Forecast Format
Each forecast file must contain:
- `origin_date`: Forecast creation date (YYYY-MM-DD)
- `target`: "ili perc"
- `horizon`: 1, 2, 3, or 4
- `location`: One of 11 required locations
- `target_end_date`: Calculated from origin_date + horizon
- `output_type`: "quantile"
- `output_type_id`: Quantile level (0.01-0.99)
- `value`: Forecast value (non-negative double)

### 3.2 Validation Checks
```r
# Validate forecast format
library(hubValidations)

# Check file structure
hubValidations::validate_submission(
  hub_path = ".",
  file_path = "model-output/[team_abbr]-[model_abbr]/forecasts.csv"
)

# Check against hub configuration
hubValidations::validate_model_data(
  hub_path = ".",
  model_id = "[team_abbr]-[model_abbr]"
)
```

## Phase 4: GitHub Submission Process

### 4.1 Commit Changes
```bash
# Stage all model files
git add model-metadata/[team_abbr]-[model_abbr].yml
git add model-code/[team_abbr]-[model_abbr]/
git add model-output/[team_abbr]-[model_abbr]/

# Commit with descriptive message
git commit -m "Add [MODEL_NAME] forecasting model

- Implement [brief method description]
- Generate forecasts for 2015-2020 flu seasons
- Include metadata and validation checks
- Performance: WIS = [if known]

🤖 Generated with Claude Code

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### 4.2 Push and Create Pull Request
```bash
# Push to your fork
git push origin add-custom-model-[MODEL_NAME]

# Create pull request
gh pr create \
  --title "Add [MODEL_NAME] forecasting model" \
  --body "$(cat <<'EOF'
## Summary
- New forecasting model: [MODEL_NAME]
- Method: [Brief description]
- Coverage: All required task combinations (2015-2020)
- Validation: Passes hubValidations checks

## Model Details
- **Team**: [team_abbr]
- **Model**: [model_abbr]  
- **Approach**: [Statistical/ML method]
- **Data**: CDC ILI surveillance data
- **Horizons**: 1-4 weeks ahead
- **Locations**: US National + 10 HHS regions

## Validation Checklist
- [ ] Model metadata file created
- [ ] Forecast format matches hub requirements
- [ ] All required quantiles included (23 levels)
- [ ] Non-negative forecast values
- [ ] Passes hubValidations checks
- [ ] Code is documented and reproducible

## Performance
[If available: Expected WIS score, comparison to baselines]

🤖 Generated with Claude Code
EOF
)"
```

## Phase 5: Automated Validation and Review

### 5.1 GitHub Actions Workflow
The hub will automatically:
- Run validation checks on your submission
- Generate performance metrics
- Create visualisations of forecasts
- Update the leaderboard dashboard

### 5.2 Address Review Feedback
```bash
# Make requested changes
git add .
git commit -m "Address review feedback: [specific changes]"
git push origin add-custom-model-[MODEL_NAME]
```

## Required Files Summary

### Required Files and Structure

```
sismid-ili-forecasting-hub/
├── model-metadata/
│   └── [team_abbr]-[model_abbr].yml        # Model description and metadata
├── model-code/
│   └── [team_abbr]-[model_abbr]/
│       ├── main.R                          # Model execution script
│       ├── model_functions.R               # Model-specific functions
│       ├── README.md                       # Model documentation
│       └── requirements.txt                # Dependencies (if needed)
├── model-output/
│   └── [team_abbr]-[model_abbr]/
│       └── forecasts.csv                   # All forecasts in hubverse format
└── hub-config/
    └── tasks.json                          # Hub configuration (read-only)
```

## Hub Configuration Details

Based on the local sandbox hub configuration:

### Task Requirements
- **Origin Dates**: 2015-2020 flu seasons (specific dates defined in tasks.json)
- **Horizons**: 1, 2, 3, 4 weeks ahead
- **Locations**: 11 total
  - US National
  - HHS Region 1-10
- **Target**: "ili perc" (Weighted ILI percentage)
- **Output Type**: 23 quantiles (0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

### Existing Model Performance Benchmarks
- **delphi-epicast**: WIS = 0.31 (operational model)
- **sismid-var2**: WIS = 0.34 (VAR(2) with cross-regional dependencies)
- **sismid-arima210**: WIS = 0.40 (univariate ARIMA)
- **hist-avg**: WIS = 0.45 (baseline)

## Key Success Factors

1. **Format Compliance**: Ensure exact adherence to hubverse format specifications
2. **Complete Coverage**: Generate forecasts for all required task combinations
3. **Validation**: Pass all hubValidations checks before submission
4. **Documentation**: Provide clear model description and reproducible code
5. **Performance**: Aim for competitive WIS scores against existing models

## Timeline Estimate
- **Setup**: 30 minutes
- **Model Development**: 2-4 hours (depending on complexity)
- **Validation**: 30 minutes
- **GitHub Submission**: 15 minutes
- **Review Process**: 1-3 days

## Notes for Implementation

- Use the local sandbox hub (`sismid-ili-forecasting-sandbox/`) for development and testing
- Leverage existing Stan models in `/inst/stan/` if applicable
- Follow the hubverse format exactly as specified in `tasks.json`
- Consider using the VAR model implementation as a reference for time-series forecasting approaches
- Ensure all forecasts are non-negative (percentage values)