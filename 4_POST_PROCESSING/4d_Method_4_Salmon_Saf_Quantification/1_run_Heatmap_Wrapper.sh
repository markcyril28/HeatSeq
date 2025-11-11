#!/bin/bash

#===============================================================================
# METHOD 4 (SALMON SAF) - HEATMAP WRAPPER SCRIPT
#===============================================================================
# Processes Salmon quantification data through tximport/DESeq2 and generates heatmaps
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
    #"SmelDMPs_with_18s_rRNA_PSBMB"
    #"SmelGIFs"
    #"SmelGRFs"
    #"Selected_GRF_GIF_Genes"
    "Selected_GRF_GIF_Genes_vAll_GIFs"
    #"Selected_GRF_GIF_Genes_vTwo_GIFs"

    # Combined Gene Groups
    #"SmelGIF_with_Cell_Cycle_Control_genes"
    #"SmelGRF_with_Cell_Cycle_Control_genes"
    #"SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Analysis toggles
RUN_MATRIX_CREATION=true
RUN_BASIC_HEATMAP=true
RUN_CV_HEATMAP=true
RUN_BAR_GRAPHS=true

# Output control
OVERWRITE_EXISTING=true  # Set to false to skip files that already exist

# Log management
CLEAR_LOGS_ON_RUN=true  # Set to true to clear all previous logs before each run


#==================================================================================


source "b_modules_for_Method_4/0_configurations.sh"
source "b_modules_for_Method_4/1_logging.sh"

# Setup logging
setup_logging "$CLEAR_LOGS_ON_RUN"

log_step "METHOD 4 HEATMAP GENERATION PIPELINE"
log_info "Starting post-processing of Salmon quantification data"

# ===============================================================================
# STEP 1: Process Salmon data with tximport
# ===============================================================================


if [ "$RUN_MATRIX_CREATION" = true ]; then
    log_step "Processing Salmon data with tximport"
    log_info "Using tximport for statistically sound Salmon quantification handling"
    log_info "Running a_tximport_salmon_to_matrices_m4.R"

    Rscript --no-save b_modules_for_Method_4/a_tximport_salmon_to_matrices_m4.R
    
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
    
    # Write gene groups to temporary file for R scripts to read
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    # Write overwrite setting
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"
    log_info "Running b_modules_for_Method_4/3_make_heatmap_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/3_make_heatmap_of_matrices.R

    log_info "Basic heatmap generation complete"
    log_info "Output: $HEATMAP_OUT_DIR"
else
    log_info "Basic heatmap generation skipped (toggle RUN_BASIC_HEATMAP in config)"
fi

# ===============================================================================
# STEP 3: Generate Heatmaps with Coefficient of Variation
# ===============================================================================

if [ "$RUN_CV_HEATMAP" = true ]; then
    log_step "Generating Heatmaps with CV"
    
    # Write gene groups to temporary file for R scripts to read
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    # Write overwrite setting
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"
    log_info "Running b_modules_for_Method_4/4_make_heatmap_with_CV_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/4_make_heatmap_with_CV_of_matrices.R

    log_info "CV heatmap generation complete"
    log_info "Output: $CV_HEATMAP_OUT_DIR"
else
    log_info "CV heatmap generation skipped (toggle RUN_CV_HEATMAP in config)"
fi

# ===============================================================================
# STEP 4: Generate Bar Graphs
# ===============================================================================

if [ "$RUN_BAR_GRAPHS" = true ]; then
    log_step "Generating Bar Graphs"
    
    # Write gene groups to temporary file for R scripts to read
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    # Write overwrite setting
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"
    log_info "Running b_modules_for_Method_4/5_make_BarGraph_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/5_make_BarGraph_of_matrices.R

    log_info "Bar graph generation complete"
    log_info "Output: $BAR_GRAPH_OUT_DIR"
else
    log_info "Bar graph generation skipped (toggle RUN_BAR_GRAPHS in config)"
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
echo "METHOD 4 HEATMAP PIPELINE - COMPLETE"
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
