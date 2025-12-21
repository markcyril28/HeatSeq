#!/bin/bash

#===============================================================================
# METHOD 4: SALMON SAF QUANTIFICATION - CONFIGURATION FILE
#===============================================================================
# Purpose: Centralized configuration for all Method 4 scripts
# Defines:
#   - Directory paths for input/output
#   - Master reference name
#   - Output structure organization
# Note: Keep synchronized with modules/methods.sh for consistency
#===============================================================================

# Base directory (parent of current module directory)
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input directories
MATRICES_DIR="$BASE_DIR/6_matrices_from_Salmon"      # tximport output matrices
GENE_GROUPS_DIR="$BASE_DIR/A_GeneGroups_InputList" # Gene lists for subsetting

# Output directories (organized by visualization type)
HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/I_Basic_Heatmap"
CV_HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/II_Heatmap_with_CV"
BAR_GRAPH_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/III_Bar_Graphs"
CONSOLIDATED_HEATMAPS_DIR="$BASE_DIR/7_Heatmap_Outputs"

# Create output directories if they don't exist
mkdir -p "$HEATMAP_OUT_DIR" "$CV_HEATMAP_OUT_DIR" "$BAR_GRAPH_OUT_DIR" 

# Master reference identifier (must match Salmon output directory name)
MASTER_REFERENCE="All_Smel_Genes"
