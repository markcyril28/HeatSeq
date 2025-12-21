#!/bin/bash

# ===============================================
# CONDA ENVIRONMENT SETUP FOR METHOD 5 - BOWTIE2/RSEM HEATMAP PIPELINE
# ===============================================

set -e

ENV_NAME="gea"

echo "================================================"
echo "METHOD 5 HEATMAP ENVIRONMENT SETUP"
echo "================================================"

# Use mamba if available, fallback to conda
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
else
    INSTALLER="conda"
fi

# Create conda environment with required packages
$INSTALLER create -n $ENV_NAME -c conda-forge -c bioconda --yes \
    r-base r-essentials \
    bioconductor-deseq2 bioconductor-tximport bioconductor-complexheatmap \
    r-tidyverse r-ggplot2 r-dplyr r-tibble r-ggrepel \
    r-rcolorbrewer r-circlize r-grid r-getopt \
    r-wgcna r-dynamictreecut r-fastcluster

echo "================================================"
echo "Environment '$ENV_NAME' created successfully!"
echo "Packages installed:"
echo "  - Bioconductor: DESeq2, tximport, ComplexHeatmap"
echo "  - Visualization: tidyverse, ggplot2, RColorBrewer"
echo "  - Coexpression: WGCNA, dynamicTreeCut, fastcluster"
echo "Activate with: conda activate $ENV_NAME"
echo "================================================"
