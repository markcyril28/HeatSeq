#!/bin/bash

#===============================================================================
# CONDA ENVIRONMENT SETUP FOR METHOD 4 - SALMON SAF HEATMAP PIPELINE
#===============================================================================
# Purpose: Sets up conda environment with all necessary R and Bioconductor packages
# Features:
#   - Uses mamba for faster installation (falls back to conda if unavailable)
#   - Installs tximport for statistically sound Salmon data import
#   - Includes DESeq2 for normalization and ComplexHeatmap for visualization
# Usage: bash 0_conda_setup_script.sh
#===============================================================================

set -e  # Exit on any error

ENV_NAME="ALL_Heatmap_ENV"

echo "==============================================================================="
echo "METHOD 4: SALMON SAF QUANTIFICATION - ENVIRONMENT SETUP"
echo "==============================================================================="
echo ""
echo "Environment name: $ENV_NAME"
echo ""

# Check if mamba is available for faster package installation, fallback to conda
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
    echo "✓ Using mamba for faster installation"
else
    INSTALLER="conda"
    echo "⚠ Using conda (mamba not found)"
    echo "  Tip: Install mamba for significantly faster environment creation:"
    echo "    conda install -c conda-forge mamba"
fi

echo ""
echo "Creating conda environment with R/Bioconductor packages..."
echo "This may take several minutes..."
echo ""

# Create conda environment with all required packages for:
#   - Statistical analysis: DESeq2, tximport
#   - Visualization: ComplexHeatmap, ggplot2
#   - Data manipulation: tidyverse tools
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
echo "==============================================================================="
echo "✓ INSTALLATION COMPLETE"
echo "==============================================================================="
echo ""
echo "Environment '$ENV_NAME' created successfully!"
echo ""
echo "Next steps:"
echo "  1. Activate the environment:"
echo "       conda activate $ENV_NAME"
echo ""
echo "  2. Run the heatmap pipeline:"
echo "       cd 4_POST_PROCESSING/4d_Method_4_Salmon_Saf_Quantification"
echo "       bash 1_run_Heatmap_Wrapper.sh"
echo ""
echo "Installed packages:"
echo "  • DESeq2          - Normalization and differential expression"
echo "  • tximport        - Statistical import of Salmon quantification"
echo "  • ComplexHeatmap  - Advanced heatmap visualization"
echo "  • ggplot2         - Bar graph generation"
echo "  • circlize        - Color gradient schemes"
echo "  • tidyverse       - Data manipulation and processing"
echo ""
echo "==============================================================================="
echo ""
