#!/bin/bash

#===============================================================================
# METHOD 5 (BOWTIE2/RSEM) - HEATMAP WRAPPER SCRIPT
#===============================================================================
# Processes RSEM quantification data through tximport and generates heatmaps
# 
# Features:
#   - Gene group selection (comment/uncomment in GENE_GROUPS array)
#   - SRR sample filtering (comment/uncomment in SRR_SAMPLES array)
#   - Toggle individual analysis types (matrix creation, heatmaps, bar graphs)
#   - Overwrite control for existing files
#
# Usage: bash 1_run_Heatmap_Wrapper.sh
#===============================================================================

set -euo pipefail

# Source configurations and logging
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

#==================================================================================

# Gene groups to process (toggle by commenting/uncommenting)
GENE_GROUPS=(
    # Control Gene Groups
    #"Best_Cell_Cycle_Associated_Control_Genes"
    #"Best_Control_Genes"

    # Individual Gene Groups
    #"SmelDMPs"
    #"SmelDMPs_with_18s_rRNA"
    "SmelDMPs_with_18s_rRNA_PSBMB"
    #"SmelGIFs"
    #"SmelGRFs"
    #"Selected_GRF_GIF_Genes"
    #"Selected_GRF_GIF_Genes_vAll_GIFs"
    #"Selected_GRF_GIF_Genes_vTwo_GIFs"

    # Combined Gene Groups
    #"SmelGIF_with_Cell_Cycle_Control_genes"
    #"SmelGRF_with_Cell_Cycle_Control_genes"
    #"SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# SRR samples to include (leave empty to use all available samples)
# Comment out samples to exclude them from analysis
SRR_SAMPLES=(
    #"SRR3884675"    # Roots
    #"SRR3884690"    # Stems
    #"SRR3884689"    # Leaves
    #"SRR3884684"    # Senescent_leaves
    "SRR3884686"    # Buds
    "SRR3884687"    # Opened_Buds
    "SRR3884597"    # Flowers
    "SRR3884631"    # Fruits
    #"SRR3884608"    # Fruits_1cm
    #"SRR3884620"    # Fruits_Stage_1
    #"SRR3884642"    # Fruits_Skin_Stage_2
    #"SRR3884653"    # Fruits_Flesh_Stage_2
    #"SRR3884664"    # Fruits_Calyx_Stage_2
    #"SRR3884680"    # Fruits_Skin_Stage_3
    #"SRR3884681"    # Fruits_Flesh_Stage_3
    #"SRR3884678"    # Fruits_peduncles
    #"SRR3884685"    # Radicles
    #"SRR3884677"    # Cotyledons
    #"SRR3884679"    # Pistils
)

# Analysis toggles
RUN_MATRIX_CREATION=false
RUN_BASIC_HEATMAP=true
RUN_CV_HEATMAP=true
RUN_BAR_GRAPHS=true

# Output control
OVERWRITE_EXISTING=true  # Set to false to skip files that already exist

# Log management
CLEAR_LOGS_ON_RUN=true  # Set to true to clear all previous logs before each run


#==================================================================================


source "b_modules_for_Method_5/0_configurations.sh"
source "b_modules_for_Method_5/1_logging.sh"

setup_logging "$CLEAR_LOGS_ON_RUN"

log_step "METHOD 5 HEATMAP GENERATION PIPELINE"
log_info "Starting post-processing of Bowtie2/RSEM quantification data"

# Export SRR samples list if specified
if [ ${#SRR_SAMPLES[@]} -gt 0 ]; then
    log_info "SRR sample selection enabled: ${#SRR_SAMPLES[@]} samples specified"
    log_info "Selected samples: ${SRR_SAMPLES[*]}"
    export SRR_FILTER_LIST="${SRR_SAMPLES[*]}"
else
    log_info "Processing all available SRR samples (no filter applied)"
    export SRR_FILTER_LIST=""
fi

# Helper function to write config files for R scripts
write_r_config() {
    printf "%s\n" "${GENE_GROUPS[@]}" > "b_modules_for_Method_5/.gene_groups_temp.txt"
    echo "$OVERWRITE_EXISTING" > "b_modules_for_Method_5/.overwrite_temp.txt"
}

# ===============================================================================
# STEP 1: Process RSEM data with tximport
# ===============================================================================


if [ "$RUN_MATRIX_CREATION" = true ]; then
    log_step "Processing RSEM data with tximport"
    log_info "Using tximport for statistically sound RSEM quantification handling"
    log_info "Running a_tximport_rsem_to_matrices_m5.R"

    Rscript --no-save b_modules_for_Method_5/a_tximport_rsem_to_matrices_m5.R
    
    if [ $? -ne 0 ]; then
        log_error "tximport processing failed. Cannot proceed with heatmap generation."
        exit 1
    fi
    
    log_info "tximport processing complete - matrices ready for visualization"
else
    log_info "Matrix creation skipped (toggle RUN_MATRIX_CREATION in config)"
fi

# ===============================================================================
# STEP 2: Generate Basic Heatmaps
# ===============================================================================

if [ "$RUN_BASIC_HEATMAP" = true ]; then
    log_step "Generating Basic Heatmaps"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save b_modules_for_Method_5/3_make_heatmap_of_matrices.R
    log_info "Basic heatmap generation complete - Output: $HEATMAP_OUT_DIR"
else
    log_info "Basic heatmap generation skipped (toggle RUN_BASIC_HEATMAP)"
fi

# ===============================================================================
# STEP 3: Generate Heatmaps with Coefficient of Variation
# ===============================================================================

if [ "$RUN_CV_HEATMAP" = true ]; then
    log_step "Generating Heatmaps with CV"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save b_modules_for_Method_5/4_make_heatmap_with_CV_of_matrices.R
    log_info "CV heatmap generation complete - Output: $CV_HEATMAP_OUT_DIR"
else
    log_info "CV heatmap generation skipped (toggle RUN_CV_HEATMAP)"
fi

# ===============================================================================
# STEP 4: Generate Bar Graphs
# ===============================================================================

if [ "$RUN_BAR_GRAPHS" = true ]; then
    log_step "Generating Bar Graphs"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save b_modules_for_Method_5/5_make_BarGraph_of_matrices.R
    log_info "Bar graph generation complete - Output: $BAR_GRAPH_OUT_DIR"
else
    log_info "Bar graph generation skipped (toggle RUN_BAR_GRAPHS)"
fi

# ===============================================================================
# FINAL SUMMARY
# ===============================================================================

log_step "PIPELINE COMPLETE"
log_info "All requested visualizations have been generated"
log_info "Gene groups processed: ${GENE_GROUPS[*]}"
log_info "Check output directories for results"
log_info "Log file: $LOG_FILE"

# Count total files generated
TOTAL_FILES=0
if [ -d "$CONSOLIDATED_HEATMAPS_DIR" ]; then
    TOTAL_FILES=$(find "$CONSOLIDATED_HEATMAPS_DIR" -type f \( -name "*.png" -o -name "*.pdf" \) 2>/dev/null | wc -l)
fi

echo ""
echo "==========================================="
echo "METHOD 5 HEATMAP PIPELINE - COMPLETE"
echo "==========================================="
echo "Outputs:"
echo "  • ALL VISUALIZATIONS: $CONSOLIDATED_HEATMAPS_DIR"
if [ "$RUN_BASIC_HEATMAP" = true ]; then
    echo "    - Basic Heatmaps"
fi
if [ "$RUN_CV_HEATMAP" = true ]; then
    echo "    - CV Heatmaps"
fi
if [ "$RUN_BAR_GRAPHS" = true ]; then
    echo "    - Bar Graphs"
fi
echo ""
echo "  • Total files generated: $TOTAL_FILES"
echo "  • Log: $LOG_FILE"
echo "==========================================="
