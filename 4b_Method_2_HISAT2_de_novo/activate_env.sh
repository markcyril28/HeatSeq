#!/bin/bash

# ===============================================
# Quick Environment Activation Script
# ===============================================
# This script activates the GEA_ENV environment
# and provides quick access to pipeline commands

ENV_NAME="GEA_ENV"

# Check if environment exists
if ! conda env list | grep -q "^$ENV_NAME "; then
    echo "Error: Environment '$ENV_NAME' not found!"
    echo "Please run 'bash setup_conda.sh' first to create the environment."
    exit 1
fi

# Activate environment
echo "Activating conda environment: $ENV_NAME"
conda activate "$ENV_NAME"

# Show environment info
echo "Environment activated successfully!"
echo ""
echo "Available commands:"
echo "  bash a_stringtie_method_2.2.sh    # Run StringTie pipeline"
echo "  Rscript b_heatmap_DeSeq2_v2.R      # Generate heatmaps"
echo "  conda deactivate                   # Exit environment"
echo ""

# Optional: Start an interactive shell in the environment
if [[ "${1:-}" == "--shell" ]]; then
    echo "Starting interactive shell in $ENV_NAME environment..."
    bash --rcfile <(echo "PS1='($ENV_NAME) \u@\h:\w\$ '")
fi