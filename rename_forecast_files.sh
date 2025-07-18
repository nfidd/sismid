#!/bin/bash

# Rename VAR2sqrt files
cd sismid-ili-forecasting-sandbox/model-output/sismid-var2sqrtclaude
for file in *-sismid-var2sqrt.csv; do
    newname="${file//-sismid-var2sqrt.csv/-sismid-var2sqrtclaude.csv}"
    mv "$file" "$newname"
done

# Rename ensemble files  
cd ../sismid-ensemblev2claude
for file in *-sismid-ensemblev2.csv; do
    newname="${file//-sismid-ensemblev2.csv/-sismid-ensemblev2claude.csv}"
    mv "$file" "$newname"
done

echo "Files renamed successfully"