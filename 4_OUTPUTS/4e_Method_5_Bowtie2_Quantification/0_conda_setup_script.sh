#!/bin/bash

# ===============================================
# CONDA ENVIRONMENT SETUP FOR METHOD 5 - BOWTIE2/RSEM HEATMAP PIPELINE
# ===============================================
# Sets up conda environment with all necessary dependencies
# Uses mamba for faster environment creation

set -e

ENV_NAME="ALL_Heatmap_ENV"

echo "================================================"
echo "METHOD 5 HEATMAP ENVIRONMENT SETUP"
echo "================================================"
echo ""
echo "Environment name: $ENV_NAME"
echo ""

# Check if mamba is available, fallback to conda
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
    echo "Using mamba for faster installation"
else
    INSTALLER="conda"
    echo "Using conda (mamba not found)"
    echo "Tip: Install mamba for faster environment creation:"
    echo "  conda install -c conda-forge mamba"
fi

echo ""
echo "Creating conda environment..."
echo ""

# Create the conda environment with all required packages
$INSTALLER create -n $ENV_NAME -c conda-forge -c bioconda --yes \
    r-base \
    r-essentials \
    bioconductor-deseq2 \
    bioconductor-tximport \
    bioconductor-complexheatmap \
    r-tidyverse \
    r-ggplot2 \
    r-dplyr \
    r-tibble \
    r-rcolorbrewer \
    r-circlize \
    r-grid \
    r-getopt

echo ""
echo "================================================"
echo "INSTALLATION COMPLETE"
echo "================================================"
echo ""
echo "Environment '$ENV_NAME' created successfully!"
echo ""
echo "To activate the environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To run the heatmap pipeline:"
echo "  cd 4_OUTPUTS/4e_Method_5_Bowtie2_Quantification"
echo "  bash 0_Heatmap_Wrapper.sh"
echo ""
echo "Required R packages installed:"
echo "  - DESeq2 (differential expression analysis)"
echo "  - tximport (RSEM data import)"
echo "  - ComplexHeatmap (heatmap generation)"
echo "  - ggplot2 (bar graphs)"
echo "  - circlize (color schemes)"
echo "  - tidyverse tools (data manipulation)"
echo ""
