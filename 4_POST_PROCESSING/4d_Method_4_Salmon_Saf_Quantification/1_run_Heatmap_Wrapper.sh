#!/bin/bash

#===============================================================================
# METHOD 4: SALMON SAF QUANTIFICATION - HEATMAP GENERATION WRAPPER
#===============================================================================
# Purpose: Orchestrates complete visualization pipeline for Salmon quantification
# Pipeline Steps:
#   1. Import Salmon data with tximport (statistically sound)
#   2. Generate basic heatmaps (multiple normalization schemes)
#   3. Generate CV-annotated heatmaps (biological variability)
#   4. Generate individual gene bar graphs
# Usage: bash 1_run_Heatmap_Wrapper.sh
# Requirements: conda environment ALL_Heatmap_ENV must be activated
#===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Navigate to script directory for relative path resolution
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

#===============================================================================
# USER CONFIGURATION
#===============================================================================

# Gene groups to process (uncomment desired groups)
GENE_GROUPS=(
    # Control gene sets for normalization validation
    #"Best_Cell_Cycle_Associated_Control_Genes"
    #"Best_Control_Genes"

    # Individual gene families
    #"SmelDMPs"                              # DORMANCY-ASSOCIATED MADS-box genes
    #"SmelDMPs_with_18s_rRNA"               # DMPs with ribosomal control
    "SmelDMPs_with_18s_rRNA_PSBMB"         # DMPs with photosystem control
    #"SmelGIFs"                              # GROWTH INHIBITING FACTORs
    #"SmelGRFs"                              # GROWTH REGULATING FACTORs
    #"Selected_GRF_GIF_Genes"               # Curated GRF-GIF subset
    "Selected_GRF_GIF_Genes_vAll_GIFs"     # All GIFs with selected GRFs
    #"Selected_GRF_GIF_Genes_vTwo_GIFs"     # Two GIFs with selected GRFs

    # Combined gene sets with controls
    #"SmelGIF_with_Cell_Cycle_Control_genes"
    #"SmelGRF_with_Cell_Cycle_Control_genes"
    "SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Pipeline step toggles (enable/disable specific analyses)
RUN_MATRIX_CREATION=true   # Import Salmon data via tximport (required for first run)
RUN_BASIC_HEATMAP=true     # Generate standard heatmaps
RUN_CV_HEATMAP=true       # Generate heatmaps with coefficient of variation
RUN_BAR_GRAPHS=false       # Generate individual gene expression bar graphs

# File handling options
OVERWRITE_EXISTING=true    # true: Regenerate all files | false: Skip existing files

# Log file management
CLEAR_LOGS_ON_RUN=true     # true: Fresh logs each run | false: Append to existing logs

# Sample IDs for matrix creation (must match Salmon output directories)
# These samples will be used by tximport to create expression matrices
SRR_SAMPLES=(
    # Roots
    "SRR3884675"   # Roots_1 (PRJNA328564)
    #"SRR20722229"  # Roots_2 (SAMN28540077)
    #"SRR31755282"  # Roots_3 (SAMN28540068)
    # Stems
    "SRR3884690"   # Stems_1 (PRJNA328564)
    #"SRR20722227"  # Stems_2 (SAMN28540077)
    #"SRR20722384"  # Stems_3 (SAMN28540068)
    # Leaves
    "SRR3884689"   # Leaves_1 (PRJNA328564)
    #"SRR20722230"  # Leaves_2 (SAMN28540077)
    #"SRR20722386"  # Leaves_3 (SAMN28540068)
    #"SRR3884684"   # Senescent_leaves (PRJNA328564)
    # Buds
    "SRR3884686"   # Buds_1 (PRJNA328564)
    #"SRR21010466"  # Buds_2 (SAMN28540077)
    #"SRR20722297"  # Buds_3 (SAMN28540068)
    # Opened Buds
    "SRR3884687"   # Opened_Buds_1 (PRJNA328564)
    # Flowers
    "SRR3884597"   # Flowers_1 (PRJNA328564)
    #"SRR20722234"  # Flowers_2 (SAMN28540077)
    #"SRR23909863"  # Flowers_3 (SAMN28540068)
    # Fruits
    "SRR3884608"   # Fruits_1cm (PRJNA328564)
    "SRR3884631"   # Fruits_1 (PRJNA328564)
    #"SRR2072232"   # Fruits_2 (SAMN28540077)
    "SRR20722387"  # Fruits_3 (SAMN28540068)
    #"SRR3884620"   # Fruits_Stage_1 (PRJNA328564)
    #"SRR3884642"   # Fruits_Skin_Stage_2 (PRJNA328564)
    #"SRR3884653"   # Fruits_Flesh_Stage_2 (PRJNA328564)
    #"SRR3884664"   # Fruits_Calyx_Stage_2 (PRJNA328564)
    #"SRR3884680"   # Fruits_Skin_Stage_3 (PRJNA328564)
    #"SRR3884681"   # Fruits_Flesh_Stage_3 (PRJNA328564)
    "SRR3884678"   # Fruits_peduncle (PRJNA328564)
    # Other organs
    "SRR3884685"   # Radicles (PRJNA328564)
    "SRR3884677"   # Cotyledons (PRJNA328564)
    "SRR3884679"   # Pistils (PRJNA328564)
)

# Export for use by R scripts
export SRR_SAMPLES

#===============================================================================
# INTERNAL CONFIGURATION (Do not modify)
#===============================================================================

# Generate unique run identifier for this execution
export RUN_ID="$(date +%Y%m%d_%H%M%S)"

# Load configuration and logging modules
source "b_modules_for_Method_4/0_configurations.sh"
source "b_modules_for_Method_4/1_logging.sh"

# Initialize logging system
setup_logging "$CLEAR_LOGS_ON_RUN"

log_step "METHOD 4: SALMON SAF QUANTIFICATION - VISUALIZATION PIPELINE"
log_info "Run ID: $RUN_ID"
log_info "Processing Salmon quantification data for heatmap generation"
log_info "Gene groups: ${#GENE_GROUPS[@]} selected"

#===============================================================================
# STEP 1: IMPORT SALMON DATA WITH TXIMPORT
#===============================================================================
# Imports Salmon quantification files using tximport package
# - Handles transcript-to-gene mapping
# - Applies appropriate length scaling
# - Generates DESeq2-compatible count matrices
# - Creates both gene-level and isoform-level matrices

if [ "$RUN_MATRIX_CREATION" = true ]; then
    log_step "STEP 1: Importing Salmon quantification with tximport"
    log_info "Approach: Statistically sound import preserving transcript-level uncertainty"
    log_info "Output: Gene-level and isoform-level count matrices"
    log_info "Executing: a_tximport_salmon_to_matrices_m4.R"

    Rscript --no-save b_modules_for_Method_4/a_tximport_salmon_to_matrices_m4.R
    
    if [ $? -ne 0 ]; then
        log_error "tximport processing failed - check Salmon output files"
        log_error "Verify quant.sf files exist in 5_Salmon_Quant_WD/All_Smel_Genes/"
        exit 1
    fi
    
    log_info "✓ tximport processing complete"
    log_info "✓ Count matrices generated and ready for visualization"
else
    log_info "⊘ Matrix creation skipped (RUN_MATRIX_CREATION=false)"
    log_info "  Using existing matrices from previous run"
fi

#===============================================================================
# STEP 2: GENERATE BASIC HEATMAPS
#===============================================================================
# Creates standard heatmaps with multiple options:
# - Multiple normalization schemes (raw, log2, z-score, etc.)
# - Original and transposed orientations
# - Sorted by organ or by expression level

if [ "$RUN_BASIC_HEATMAP" = true ]; then
    log_step "STEP 2: Generating Basic Heatmaps"
    
    # Pass configuration to R scripts via temporary files
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Gene groups to process: ${GENE_GROUPS[*]}"
    log_info "Overwrite mode: $OVERWRITE_EXISTING"
    log_info "Executing: 3_make_heatmap_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/3_make_heatmap_of_matrices.R

    log_info "✓ Basic heatmap generation complete"
    log_info "  Output directory: $HEATMAP_OUT_DIR"
else
    log_info "⊘ Basic heatmap generation skipped (RUN_BASIC_HEATMAP=false)"
fi

#===============================================================================
# STEP 3: GENERATE CV-ANNOTATED HEATMAPS
#===============================================================================
# Creates heatmaps with Coefficient of Variation annotations
# - CV calculated on raw data for biological interpretation
# - Heatmap displays normalized values
# - Useful for assessing gene expression variability across samples

if [ "$RUN_CV_HEATMAP" = true ]; then
    log_step "STEP 3: Generating Heatmaps with Coefficient of Variation"
    
    # Pass configuration to R scripts
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Gene groups to process: ${GENE_GROUPS[*]}"
    log_info "CV calculation: On raw, untransformed counts"
    log_info "Executing: 4_make_heatmap_with_CV_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/4_make_heatmap_with_CV_of_matrices.R

    log_info "✓ CV heatmap generation complete"
    log_info "  Output directory: $CV_HEATMAP_OUT_DIR"
else
    log_info "⊘ CV heatmap generation skipped (RUN_CV_HEATMAP=false)"
fi

#===============================================================================
# STEP 4: GENERATE INDIVIDUAL GENE BAR GRAPHS
#===============================================================================
# Creates separate bar graph for each gene showing expression across organs
# - One PNG file per gene
# - Useful for detailed single-gene analysis
# - Sorted by organ or expression level

if [ "$RUN_BAR_GRAPHS" = true ]; then
    log_step "STEP 4: Generating Individual Gene Bar Graphs"
    
    # Pass configuration to R scripts
    GENE_GROUPS_FILE="b_modules_for_Method_4/.gene_groups_temp.txt"
    printf "%s\n" "${GENE_GROUPS[@]}" > "$GENE_GROUPS_FILE"
    
    OVERWRITE_FILE="b_modules_for_Method_4/.overwrite_temp.txt"
    echo "$OVERWRITE_EXISTING" > "$OVERWRITE_FILE"
    
    log_info "Gene groups to process: ${GENE_GROUPS[*]}"
    log_info "Output: One bar graph per gene"
    log_info "Executing: 5_make_BarGraph_of_matrices.R"

    Rscript --no-save b_modules_for_Method_4/5_make_BarGraph_of_matrices.R

    log_info "✓ Bar graph generation complete"
    log_info "  Output directory: $BAR_GRAPH_OUT_DIR"
else
    log_info "⊘ Bar graph generation skipped (RUN_BAR_GRAPHS=false)"
fi

#===============================================================================
# PIPELINE COMPLETION SUMMARY
#===============================================================================

log_step "PIPELINE EXECUTION COMPLETE"
log_info "All requested analyses have been completed successfully"
log_info "Gene groups processed: ${GENE_GROUPS[*]}"
log_info "Results available in output directories"
log_info "Detailed log: $LOG_FILE"

# Calculate total output files generated
TOTAL_FILES=0
if [ -d "$CONSOLIDATED_HEATMAPS_DIR" ]; then
    TOTAL_FILES=$(find "$CONSOLIDATED_HEATMAPS_DIR" -type f \( -name "*.png" -o -name "*.pdf" \) 2>/dev/null | wc -l)
fi

echo ""
echo "==============================================================================="
echo "✓ METHOD 4: SALMON SAF QUANTIFICATION - PIPELINE COMPLETE"
echo "==============================================================================="
echo ""
echo "Results Summary:"
echo "  • Main output directory: $CONSOLIDATED_HEATMAPS_DIR"
echo ""
if [ "$RUN_BASIC_HEATMAP" = true ]; then
    echo "  ✓ Basic Heatmaps         → I_Basic_Heatmap/"
fi
if [ "$RUN_CV_HEATMAP" = true ]; then
    echo "  ✓ CV-Annotated Heatmaps  → II_Heatmap_with_CV/"
fi
if [ "$RUN_BAR_GRAPHS" = true ]; then
    echo "  ✓ Gene Bar Graphs        → III_Bar_Graphs/"
fi
echo ""
echo "Statistics:"
echo "  • Gene groups processed: ${#GENE_GROUPS[@]}"
echo "  • Total files generated: $TOTAL_FILES"
echo "  • Detailed log file: $LOG_FILE"
echo ""
echo "Next steps:"
echo "  • Review visualizations in output directories"
echo "  • Check log file for detailed processing information"
echo "  • Adjust gene groups or toggles in this script for different analyses"
echo ""
echo "==============================================================================="
