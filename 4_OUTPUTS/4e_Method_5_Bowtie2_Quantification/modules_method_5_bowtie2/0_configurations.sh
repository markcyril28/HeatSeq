#!/bin/bash

# ===============================================
# METHOD 5 - BOWTIE2 + RSEM QUANTIFICATION
# Configuration File
# ===============================================

# Base directories
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRICES_DIR="$BASE_DIR/4_matrices"
GENE_GROUPS_DIR="$BASE_DIR/2_gene_groups"


# Output directories
HEATMAP_OUT_DIR="$BASE_DIR/6_Basic_Heatmap_Visualizations"
CV_HEATMAP_OUT_DIR="$BASE_DIR/7_CV_Heatmap_Visualizations"
BAR_GRAPH_OUT_DIR="$BASE_DIR/8_Bar_Graph_Visualizations"

# Consolidated output directory for all heatmaps
CONSOLIDATED_HEATMAPS_DIR="$BASE_DIR/ALL_HEATMAPS_CONSOLIDATED"

# Create output directories
mkdir -p "$HEATMAP_OUT_DIR" "$CV_HEATMAP_OUT_DIR" "$BAR_GRAPH_OUT_DIR" "$CONSOLIDATED_HEATMAPS_DIR"

# Master reference name
MASTER_REFERENCE="All_Smel_Genes"

