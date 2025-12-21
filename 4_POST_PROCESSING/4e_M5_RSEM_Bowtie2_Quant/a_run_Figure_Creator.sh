#!/bin/bash

#===============================================================================
# METHOD 5 (BOWTIE2/RSEM) - HEATMAP WRAPPER SCRIPT
#===============================================================================
# Processes RSEM quantification data through tximport and generates heatmaps
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

#==================================================================================
# CONFIGURATION
#==================================================================================

# Master Reference to use 
MASTER_REFERENCE="Eggplant_V4.1_transcripts.function"
#MASTER_REFERENCE="All_Smel_Genes"

# Genes of interest to extract for Figure Creation and Analyses
GENE_GROUPS=(
    #"SmelGRF-GIFs"
    #"SmelDMPs"
    "SmelDMP10_200"
)

# NOTE: Sample selection is controlled by SAMPLE_LABELS in B_M5_modules/0_shared_config.R
# Edit SAMPLE_LABELS there to add/remove samples and set organ names

# Analysis toggles
ANALYSIS_STEPS=(
    "MATRIX_CREATION"
    #"BASIC_HEATMAP"
    #"CV_HEATMAP"
    #"BAR_GRAPHS"
    "COEXPRESSION_ANALYSIS"
    #"DIFFERENTIAL_EXPRESSION"
    #"GENE_SET_ENRICHMENT"
    #"DIMENSIONALITY_REDUCTION"
    #"SAMPLE_CORRELATION"
    #"TISSUE_SPECIFICITY"
)

OVERWRITE_EXISTING=true
CLEAR_LOGS_ON_RUN=true

#==================================================================================

source "B_M5_modules/0_configurations.sh"
source "B_M5_modules/1_logging.sh"

setup_logging "$CLEAR_LOGS_ON_RUN"

log_step "METHOD 5 HEATMAP GENERATION PIPELINE"
log_info "Starting post-processing of Bowtie2/RSEM quantification data"
log_info "Sample selection controlled by SAMPLE_LABELS in B_M5_modules/0_shared_config.R"

# Helper function to write config files for R scripts
write_r_config() {
    printf "%s\n" "${GENE_GROUPS[@]}" > "B_M5_modules/.gene_groups_temp.txt"
    echo "$MASTER_REFERENCE" > "B_M5_modules/.master_reference_temp.txt"
    echo "$OVERWRITE_EXISTING" > "B_M5_modules/.overwrite_temp.txt"
}

# Helper function to check if step is enabled
is_step_enabled() {
    local step="$1"
    for s in "${ANALYSIS_STEPS[@]}"; do
        [[ "$s" == "$step" ]] && return 0
    done
    return 1
}

# ===============================================================================
# STEP 1: Process RSEM data with tximport
# ===============================================================================

if is_step_enabled "MATRIX_CREATION"; then
    log_step "Processing RSEM data with tximport"
    log_info "Using tximport for statistically sound RSEM quantification handling"
    log_info "Running a_tximport_rsem_to_matrices_m5.R"

    Rscript --no-save B_M5_modules/a_tximport_rsem_to_matrices_m5.R
    
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

if is_step_enabled "BASIC_HEATMAP"; then
    log_step "Generating Basic Heatmaps"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save B_M5_modules/3_make_heatmap_of_matrices.R
    log_info "Basic heatmap generation complete - Output: $HEATMAP_OUT_DIR"
else
    log_info "Basic heatmap generation skipped (toggle RUN_BASIC_HEATMAP)"
fi

# ===============================================================================
# STEP 3: Generate Heatmaps with Coefficient of Variation
# ===============================================================================

if is_step_enabled "CV_HEATMAP"; then
    log_step "Generating Heatmaps with CV"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save B_M5_modules/4_make_heatmap_with_CV_of_matrices.R
    log_info "CV heatmap generation complete - Output: $CV_HEATMAP_OUT_DIR"
else
    log_info "CV heatmap generation skipped (toggle RUN_CV_HEATMAP)"
fi

# ===============================================================================
# STEP 4: Generate Bar Graphs
# ===============================================================================

if is_step_enabled "BAR_GRAPHS"; then
    log_step "Generating Bar Graphs"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save B_M5_modules/5_make_BarGraph_of_matrices.R
    log_info "Bar graph generation complete - Output: $BAR_GRAPH_OUT_DIR"
else
    log_info "Bar graph generation skipped (toggle RUN_BAR_GRAPHS)"
fi

# ===============================================================================
# STEP 5: WGCNA Coexpression Analysis
# ===============================================================================

if is_step_enabled "COEXPRESSION_ANALYSIS"; then
    log_step "Running WGCNA Coexpression Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    log_info "Overwrite existing files: $OVERWRITE_EXISTING"

    Rscript --no-save B_M5_modules/6_coexpression_wgcna.R
    log_info "WGCNA analysis complete - Output: $WGCNA_OUT_DIR"
else
    log_info "WGCNA coexpression analysis skipped (toggle RUN_COEXPRESSION_ANALYSIS)"
fi

# ===============================================================================
# STEP 6: Differential Expression Analysis (DESeq2)
# ===============================================================================

if is_step_enabled "DIFFERENTIAL_EXPRESSION"; then
    log_step "Running Differential Expression Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    Rscript --no-save B_M5_modules/7_differential_expression.R
    log_info "DEA complete - Output: 7_Figures_Outputs/V_Differential_Expression"
else
    log_info "Differential expression analysis skipped"
fi

# ===============================================================================
# STEP 7: Gene Set Enrichment Analysis
# ===============================================================================

if is_step_enabled "GENE_SET_ENRICHMENT"; then
    log_step "Running Gene Set Enrichment Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    Rscript --no-save B_M5_modules/8_gene_set_enrichment.R
    log_info "GSEA complete - Output: 7_Figures_Outputs/VI_Gene_Set_Enrichment"
else
    log_info "Gene set enrichment analysis skipped"
fi

# ===============================================================================
# STEP 8: Dimensionality Reduction (PCA, t-SNE, UMAP)
# ===============================================================================

if is_step_enabled "DIMENSIONALITY_REDUCTION"; then
    log_step "Running Dimensionality Reduction Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    Rscript --no-save B_M5_modules/9_pca_dimensionality_reduction.R
    log_info "Dim reduction complete - Output: 7_Figures_Outputs/VII_Dimensionality_Reduction"
else
    log_info "Dimensionality reduction skipped"
fi

# ===============================================================================
# STEP 9: Sample Correlation & Clustering
# ===============================================================================

if is_step_enabled "SAMPLE_CORRELATION"; then
    log_step "Running Sample Correlation Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    Rscript --no-save B_M5_modules/10_sample_correlation_clustering.R
    log_info "Correlation analysis complete - Output: 7_Figures_Outputs/VIII_Sample_Correlation"
else
    log_info "Sample correlation analysis skipped"
fi

# ===============================================================================
# STEP 10: Tissue Specificity Analysis
# ===============================================================================

if is_step_enabled "TISSUE_SPECIFICITY"; then
    log_step "Running Tissue Specificity Analysis"
    write_r_config
    
    log_info "Processing gene groups: ${GENE_GROUPS[*]}"
    Rscript --no-save B_M5_modules/11_tissue_specificity.R
    log_info "Tissue specificity complete - Output: 7_Figures_Outputs/IX_Tissue_Specificity"
else
    log_info "Tissue specificity analysis skipped"
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
if [ -d "$HEATMAP_OUT_DIR" ] || [ -d "$CV_HEATMAP_OUT_DIR" ] || [ -d "$BAR_GRAPH_OUT_DIR" ] || [ -d "$WGCNA_OUT_DIR" ]; then
    TOTAL_FILES=$(find "$HEATMAP_OUT_DIR" "$CV_HEATMAP_OUT_DIR" "$BAR_GRAPH_OUT_DIR" "$WGCNA_OUT_DIR" -type f \( -name "*.png" -o -name "*.tsv" \) 2>/dev/null | wc -l)
fi

echo ""
echo "==========================================="
echo "METHOD 5 HEATMAP PIPELINE - COMPLETE"
echo "==========================================="
echo "Outputs: 7_Figures_Outputs/"
is_step_enabled "BASIC_HEATMAP" && echo "  - Basic Heatmaps"
is_step_enabled "CV_HEATMAP" && echo "  - CV Heatmaps"
is_step_enabled "BAR_GRAPHS" && echo "  - Bar Graphs"
is_step_enabled "COEXPRESSION_ANALYSIS" && echo "  - WGCNA Coexpression Analysis"
is_step_enabled "DIFFERENTIAL_EXPRESSION" && echo "  - Differential Expression (DESeq2)"
is_step_enabled "GENE_SET_ENRICHMENT" && echo "  - Gene Set Enrichment"
is_step_enabled "DIMENSIONALITY_REDUCTION" && echo "  - PCA/t-SNE/UMAP"
is_step_enabled "SAMPLE_CORRELATION" && echo "  - Sample Correlation"
is_step_enabled "TISSUE_SPECIFICITY" && echo "  - Tissue Specificity"
echo "Total files generated: $TOTAL_FILES"
echo "Log: $LOG_FILE"
echo "==========================================="
