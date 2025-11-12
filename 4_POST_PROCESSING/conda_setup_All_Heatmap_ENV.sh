#!/bin/bash

# ==============================================================================
# COMBINED CONDA/MAMBA ENVIRONMENT SETUP SCRIPT
# ==============================================================================
# Description: This script sets up a single conda environment named
#              'ALL_Heatmap_ENV' with all necessary dependencies from
#              multiple sub-project setup scripts.
#
# Environment Name: ALL_Heatmap_ENV
#
# Instructions:
# 1. Make sure you have Conda or Miniconda installed.
# 2. Run this script from your terminal: ./conda_setup_All_Heatmap_ENV.sh
# 3. Activate the environment: conda activate ALL_Heatmap_ENV
# ==============================================================================

# Stop on first error, and handle unbound variables and pipe failures
set -euo pipefail

# --- Configuration ---
#ENV_NAME="All_Heatmap_ENV"
ENV_NAME="geaheat"
MAMBA_FLAGS="--quiet --yes"

# List of packages to install
PACKAGES=(
    # Alignment and quantification tools
    "hisat2"
    "stringtie"
    "samtools"
    "salmon"
    "bowtie2"
    "rsem"
    "trinity"
    
    # R base and essentials
    "r-base=4.3"
    "r-essentials"
    
    # Bioconductor packages
    "bioconductor-deseq2"
    "bioconductor-complexheatmap"
    "bioconductor-ballgown"
    "bioconductor-tximport"
    "bioconductor-tximeta"
    "bioconductor-annotationdbi"
    "bioconductor-org.mm.eg.db"
    
    # CRAN packages for data manipulation
    "r-tidyverse"
    "r-dplyr"
    "r-tibble"
    "r-readr"
    "r-biocmanager"
    
    # Visualization packages
    "r-ggplot2"
    "r-rcolorbrewer"
    "r-circlize"
    "r-getopt"
)

# --- Main Execution ---

echo "================================================="
echo "Setting up Conda Environment: ${ENV_NAME}"
echo "================================================="

CONDA_CMD="conda"
: << 'OFF'
# Check if Mamba is available, otherwise use Conda
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found, using Conda. For faster installation, consider installing Mamba ('conda install -c conda-forge mamba')."
    CONDA_CMD="conda"
else
    CONDA_CMD="mamba"
fi
OFF

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Check if the environment already exists
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo "Environment '${ENV_NAME}' already exists. Updating..."
    "${CONDA_CMD}" install -n ${ENV_NAME} -c conda-forge -c bioconda ${MAMBA_FLAGS} "${PACKAGES[@]}"
else
    echo "Creating new Conda environment '${ENV_NAME}'..."
    "${CONDA_CMD}" create -n ${ENV_NAME} python=3.11 ${MAMBA_FLAGS}
    echo "Installing packages from conda-forge and bioconda..."
    # We specify channels to ensure correct package versions are found
    # Channel priority is important: conda-forge > bioconda > defaults
    "${CONDA_CMD}" install -n ${ENV_NAME} -c conda-forge -c bioconda ${MAMBA_FLAGS} "${PACKAGES[@]}"
fi

# Activate the environment
echo "Activating environment '${ENV_NAME}'..."
conda activate ${ENV_NAME}

#conda install -c conda-forge r-base

mamba install -c conda-forge r-essentials

# Install additional R packages via Rscript
echo "Installing additional R packages..."
if [ -f "install_R_packages.R" ]; then
    Rscript --no-save install_R_packages.R
else
    echo "Warning: install_R_packages.R not found. Skipping R package installation."
fi

echo "-------------------------------------------------"
echo "Environment setup complete and activated!"
echo "Current environment: ${ENV_NAME}"
echo "================================================="