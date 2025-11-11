#!/bin/bash

# ===============================================
# METHOD 4 - SALMON SAF QUANTIFICATION
# Configuration File
# ===============================================

# Base directories - sync with modules/methods.sh
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRICES_DIR="$BASE_DIR/6_matrices_from_Salmon"
GENE_GROUPS_DIR="$BASE_DIR/a_gene_groups_input_list"

# Output directories - standardized naming
HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/I_Basic_Heatmap"
CV_HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/II_Heatmap_with_CV"
BAR_GRAPH_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/III_Bar_Graphs"

# Consolidated output directory
CONSOLIDATED_HEATMAPS_DIR="$BASE_DIR/ALL_HEATMAPS_CONSOLIDATED"

# Create output directories
#mkdir -p "$HEATMAP_OUT_DIR" "$CV_HEATMAP_OUT_DIR" "$BAR_GRAPH_OUT_DIR" "$CONSOLIDATED_HEATMAPS_DIR"

# Master reference name
MASTER_REFERENCE="All_Smel_Genes"
