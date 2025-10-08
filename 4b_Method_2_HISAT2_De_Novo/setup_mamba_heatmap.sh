#!/bin/bash

# ===============================================
# Simplified Heatmap Environment Setup
# ===============================================

set -euo pipefail

ENV_NAME="Heatmap_ENV"
YAML_FILE="$(dirname "${BASH_SOURCE[0]}")/Heatmap_ENV_R.yml"

echo "Setting up Heatmap_ENV environment..."

# Check if conda/mamba is available
if command -v mamba >/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda >/dev/null; then
    CONDA_CMD="conda"
else
    echo "Error: Neither conda nor mamba found. Please install conda first." >&2
    exit 1
fi

echo "Using: $CONDA_CMD"

# Remove existing environment if it exists
if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
    echo "Removing existing $ENV_NAME environment..."
    $CONDA_CMD env remove -n "$ENV_NAME" -y
fi

# Create new environment from YAML
echo "Creating $ENV_NAME environment from YAML..."
$CONDA_CMD env create -f "$YAML_FILE"

# Activate environment and install additional R packages
echo "Installing additional R packages..."
eval "$($CONDA_CMD shell.bash hook)"
$CONDA_CMD activate "$ENV_NAME"

# Install required R packages that might not be available via conda
Rscript -e "
# Install packages that need special handling
packages_to_install <- c()

# Check and install ComplexHeatmap if not available
if (!require('ComplexHeatmap', quietly = TRUE)) {
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
        install.packages('BiocManager', repos = 'https://cloud.r-project.org/')
    }
    BiocManager::install('ComplexHeatmap', ask = FALSE, update = FALSE)
    packages_to_install <- c(packages_to_install, 'ComplexHeatmap')
}

# Check and install circlize if not available
if (!require('circlize', quietly = TRUE)) {
    install.packages('circlize', repos = 'https://cloud.r-project.org/')
    packages_to_install <- c(packages_to_install, 'circlize')
}

# Check and install grid (usually comes with base R)
if (!require('grid', quietly = TRUE)) {
    install.packages('grid', repos = 'https://cloud.r-project.org/')
    packages_to_install <- c(packages_to_install, 'grid')
}

if (length(packages_to_install) > 0) {
    cat('Successfully installed:', paste(packages_to_install, collapse = ', '), '\n')
} else {
    cat('All required packages are already available\n')
}

# Test that everything works
cat('Testing package loading...\n')
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(grid)
cat('All packages loaded successfully!\n')
"

echo ""
echo "==============================================="
echo "Setup complete!"
echo "To use the environment, run:"
echo "  conda activate $ENV_NAME"
echo "==============================================="
