#!/bin/bash

# Exit on error
set -e

# --- Conda Environment Setup ---
# This script sets up the conda environment with all necessary dependencies.
# It uses mamba for faster environment creation.

# Environment name
ENV_NAME="gea"

# Check if mamba is installed
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Please install mamba first."
    echo "e.g., conda install -c conda-forge mamba"
    exit 1
fi

# Create the conda environment
echo "Creating conda environment: $ENV_NAME"
mamba create -n $ENV_NAME -c conda-forge -c bioconda --yes \
    r-base \
    r-essentials \
    r-deseq2 \
    r-complexheatmap \
    r-tidyverse \
    r-rcolorbrewer \
    r-circlize \
    r-getopt

echo "Conda environment '$ENV_NAME' created successfully."
echo "To activate the environment, run: conda activate $ENV_NAME"
