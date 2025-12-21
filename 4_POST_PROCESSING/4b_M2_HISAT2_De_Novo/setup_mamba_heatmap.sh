#!/bin/bash

# ===============================================
# Simplified Heatmap Environment Setup
# ===============================================

set -euo pipefail

# Define paths relative to the script's location
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
ENV_NAME="Heatmap_ENV"
YAML_FILE="$SCRIPT_DIR/Heatmap_ENV_R.yml"
R_SCRIPT_FILE="$SCRIPT_DIR/install_heatmap_libraries.R" # <-- Added this

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

#$CONDA_CMD install r-base r-essentials -y

: << 'OFF'
# Remove existing environment if it exists
#if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
#    echo "Removing existing $ENV_NAME environment..."
#    $CONDA_CMD env remove -n "$ENV_NAME" -y
#fi
OFF

conda install -c conda-forge xorg-libxt xorg-libx11 xorg-libxrender xorg-libxext xorg-libsm xorg-libice -y

conda install -c conda-forge xorg-libxau xorg-libxdmcp xorg-libxpm -y

# Create new environment from YAML
#echo "Creating $ENV_NAME environment from YAML..."
#$CONDA_CMD env create -f "$YAML_FILE"

# Instead of creating a new environment, update the active one (or 'base' if none)
echo "Updating existing environment '$ENV_NAME' from YAML (no new env will be created)..."
#$CONDA_CMD env update -n "$ENV_NAME" -f "$YAML_FILE"

# Install additional R packages using 'conda run'
echo "Installing additional R packages..."
#$CONDA_CMD run -n "$ENV_NAME" Rscript "$R_SCRIPT_FILE"

#Rscript --no-save "$R_SCRIPT_FILE"


echo ""
echo "==============================================="
echo "Setup complete!"
echo "To use the environment, run:"
echo "  $CONDA_CMD activate $ENV_NAME"
echo "==============================================="
