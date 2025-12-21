#!/bin/bash

# ===============================================
# METHOD 5 - BOWTIE2 + RSEM QUANTIFICATION
# Configuration File
# ===============================================

# Base directories - sync with modules/methods.sh
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRICES_DIR="$BASE_DIR/6_matrices"

# Master reference name
#MASTER_REFERENCE="All_Smel_Genes"
MASTER_REFERENCE="Eggplant_V4.1_transcripts.function"

# Gene groups directory - points to subdirectory based on MASTER_REFERENCE
GENE_GROUPS_DIR="$BASE_DIR/A_GeneGroups_InputList/for_${MASTER_REFERENCE}"

# Output directories - standardized naming
HEATMAP_OUT_DIR="$BASE_DIR/7_Figures_Outputs/I_Basic_Heatmap"
CV_HEATMAP_OUT_DIR="$BASE_DIR/7_Figures_Outputs/II_Heatmap_with_CV"
BAR_GRAPH_OUT_DIR="$BASE_DIR/7_Figures_Outputs/III_Bar_Graphs"
WGCNA_OUT_DIR="$BASE_DIR/7_Figures_Outputs/IV_Coexpression_WGCNA"

# Consolidated output directory
CONSOLIDATED_HEATMAPS_DIR="$BASE_DIR/7_Figures_Outputs"