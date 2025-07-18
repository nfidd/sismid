# Hub Submission Validation Troubleshooting Guide

This guide documents the common validation issues encountered when submitting forecasts to modeling hubs and their solutions.

## Overview

When submitting forecasts to a modeling hub (like the sismid-ili-forecasting-sandbox), validation failures can occur due to metadata format issues, filename parsing problems, or invalid configuration values.

## Common Validation Issues and Solutions

### 1. Filename Parsing Errors

**Problem**: Error message like:
```
File name '2017-10-21-sismid-var2-sqrt' does not match expected pattern of [round_id]-[team_abbr]-[model_abbr]
```

**Root Cause**: The hub filename parser expects exactly 3 parts separated by hyphens: `[round_id]-[team_abbr]-[model_abbr]`. If your `model_abbr` contains hyphens, the parser sees more than 3 parts.

**Solution**:
1. Remove hyphens from your `model_abbr` in the metadata YAML file
2. Rename the metadata file to match the new abbreviation
3. Rename the model output directory
4. Rename all CSV files to use the new abbreviation

**Example Fix**:
```yaml
# Before (causes validation failure)
team_abbr: "sismid"
model_abbr: "var2-sqrt"  # ❌ Contains hyphen

# After (passes validation)
team_abbr: "sismid"
model_abbr: "var2sqrt"   # ✅ No hyphens
```

### 2. Invalid Round IDs

**Problem**: Error message like:
```
`round_id` must be one of "2015-10-24", "2015-10-31", ... NOT "2015-10-17"
```

**Root Cause**: Your submission contains files with round_ids (dates) that aren't in the hub's valid list.

**Solution**:
1. Check the hub's `hub-config/tasks.json` file for valid `origin_date` values
2. Remove any files with invalid round_ids
3. Ensure all your forecast files use only valid round_ids

**How to Check**:
```bash
# Find files with invalid round_ids
ls model-output/your-model/2015-10-17* 2>/dev/null

# Remove invalid files
rm model-output/your-model/2015-10-17-your-model.csv
```

### 3. Missing Required Locations

**Problem**: Validation fails because forecasts don't cover all required locations.

**Solution**: Ensure your forecasts include all required locations as specified in the hub configuration:
- US National
- HHS Region 1 through HHS Region 10

**Verification**:
```bash
# Check locations in a forecast file
cut -d',' -f1 model-output/your-model/2017-10-21-your-model.csv | sort | uniq
```

## Debugging Workflow

### Step 1: Check PR Validation Status
```bash
# View PR validation status
gh pr checks <PR_NUMBER> --repo <HUB_REPO>

# Get detailed error logs
gh run view <RUN_ID> --repo <HUB_REPO> --log | grep -A10 -B5 "Error\|Failed"
```

### Step 2: Common File Operations
```bash
# Rename metadata file
mv model-metadata/team-model-old.yml model-metadata/team-model-new.yml

# Rename model directory
mv model-output/team-model-old model-output/team-model-new

# Batch rename CSV files
cd model-output/team-model-new
for file in *-team-model-old.csv; do 
    mv "$file" "${file//-team-model-old/-team-model-new}"
done
```

### Step 3: Commit and Push Fixes
```bash
# Stage all changes
git add -A

# Commit with descriptive message
git commit -m "Fix validation: remove hyphens from model_abbr"

# Push to your fork
git push fork main
```

## Hub Configuration Reference

Key files to check in the hub repository:
- `hub-config/tasks.json` - Valid round_ids, locations, and other task parameters
- `hub-config/model-metadata-schema.json` - Required metadata fields
- `hub-config/admin.json` - General hub configuration

## Prevention Tips

1. **Model Naming**: Use simple names without hyphens, underscores, or special characters
2. **Date Validation**: Always check valid round_ids in the hub configuration before generating forecasts
3. **Location Coverage**: Ensure all required locations are included in your forecasts
4. **Local Testing**: Test your submission locally using hub validation tools before creating a PR

## Example Working Submission Structure

```
model-metadata/
├── sismid-var2sqrt.yml          # ✅ No hyphens in filename or model_abbr

model-output/
├── sismid-var2sqrt/             # ✅ Directory matches metadata
    ├── 2017-10-21-sismid-var2sqrt.csv  # ✅ Valid round_id, proper format
    ├── 2017-10-28-sismid-var2sqrt.csv  # ✅ All files follow same pattern
    └── ...
```

## When to Wait vs. Fix

- **Wait**: If validation is still running (status: "in_progress")
- **Fix**: If validation fails with specific error messages
- **Monitor**: Check validation status every 1-2 minutes during the run

This troubleshooting guide should help future agents quickly identify and resolve common hub submission validation issues.