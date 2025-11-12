#!/bin/bash

# ===============================================
# METHOD 5 - BOWTIE2 + RSEM QUANTIFICATION
# Configuration File
# ===============================================

# Base directories - sync with modules/methods.sh
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRICES_DIR="$BASE_DIR/4_matrices"
GENE_GROUPS_DIR="$BASE_DIR/a_gene_groups_input_list"

# Output directories - standardized naming
HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/I_Basic_Heatmap"
CV_HEATMAP_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/II_Heatmap_with_CV"
BAR_GRAPH_OUT_DIR="$BASE_DIR/7_Heatmap_Outputs/III_Bar_Graphs"

# Consolidated output directory
CONSOLIDATED_HEATMAPS_DIR="$BASE_DIR/ALL_HEATMAPS_CONSOLIDATED"

# Master reference name
MASTER_REFERENCE="All_Smel_Genes"