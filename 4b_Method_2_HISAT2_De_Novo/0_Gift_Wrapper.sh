#!/bin/bash 

: << 'DESCRIPTIONS'
Changes conda environment midway through the HISAT2 De Novo pipeline to run StringTie and then create matrices and heatmaps.
CONDA ENV NAME: Heatmap_ENV

Run this script in sequence:
1. bash a_stringtie_matrix_builder_m2_v4.sh
2. Rscript b_make_heatmap_of_matrices_v4.R
DESCRIPTIONS

# Initialize conda/mamba if not already done
if command -v mamba &> /dev/null; then
    echo "Initializing mamba shell..."
    eval "$(mamba shell hook --shell bash)"
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    echo "Initializing conda shell..."
    eval "$(conda shell.bash hook)"
    CONDA_CMD="conda"
else
    echo "Error: Neither mamba nor conda found in PATH"
    exit 1
fi

# Activate the conda environment for the remaining pipeline steps
echo "Activating environment: Heatmap_ENV"
$CONDA_CMD activate Heatmap_ENV

# Verify environment activation
if [ "$CONDA_DEFAULT_ENV" = "Heatmap_ENV" ]; then
    echo "Successfully activated environment: $CONDA_DEFAULT_ENV"
else
    echo "Error: Failed to activate Heatmap_ENV environment"
    echo "Current environment: $CONDA_DEFAULT_ENV"
    exit 1
fi

#Rscript --verbose install_heatmap_libraries.R
#bash setup_mamba_Heatmap_ENV.sh

#bash a_stringtie_matrix_builder_m2_v4.sh

#Rscript --verbose b_make_heatmap_of_matrices_v4.R
