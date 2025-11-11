#!/bin/bash

#===============================================================================
# METHOD 3 (TRINITY DE NOVO) - HEATMAP WRAPPER SCRIPT
#===============================================================================
# Generates heatmaps from StringTie quantification data (Trinity de novo assembly)
# Usage: bash 1_run_Heatmap_Wrapper.sh
#===============================================================================

set -euo pipefail

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
OVERWRITE_EXISTING=true

#==================================================================================

FASTA_TAG="All_Smel_Genes"

echo "======================================="
echo "METHOD 3 HEATMAP GENERATION PIPELINE"
echo "======================================="
echo "Trinity De Novo + StringTie Quantification"
echo ""

# Check for gene count matrix
MATRIX_PATH="5_stringtie_WD/${FASTA_TAG}/deseq2_input/gene_count_matrix.csv"
if [ ! -f "$MATRIX_PATH" ]; then
    echo "ERROR: Matrix not found at $MATRIX_PATH"
    echo "Run the Trinity de novo alignment pipeline first (0_run_Method_3_Pipeline.sh)"
    exit 1
fi

echo "Found matrix: $MATRIX_PATH"
echo ""

#==================================================================================
# STEP 1: Create normalized matrices with DESeq2
#==================================================================================

if [ "$RUN_MATRIX_CREATION" = true ]; then
    echo "STEP 1: Creating normalized expression matrices"
    echo "Running StringTie to DESeq2 matrix conversion..."
    
    Rscript --no-save b_modules_for_Method_3/d_stringtie_to_deseq2_m3.R
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Matrix creation failed"
        exit 1
    fi
    
    echo "Matrix creation complete"
    echo ""
else
    echo "STEP 1: Matrix creation SKIPPED"
    echo ""
fi

#==================================================================================
# STEP 2: Generate heatmaps for each gene group
#==================================================================================

process_gene_group() {
    local gene_group=$1
    local gene_group_file="a_gene_groups_input_list/${gene_group}.txt"
    
    if [ ! -f "$gene_group_file" ]; then
        echo "WARNING: Gene group file not found: $gene_group_file"
        return
    fi
    
    echo "Processing: $gene_group"
    
    # Basic heatmap
    if [ "$RUN_BASIC_HEATMAP" = true ]; then
        echo "  - Generating basic heatmap..."
        Rscript --vanilla b_modules_for_Method_3/b_make_heatmap_of_matrices.R \
            -d "$SCRIPT_DIR" \
            -f "$FASTA_TAG" \
            -c "counts" \
            -g "$gene_group_file"
    fi
    
    # CV heatmap
    if [ "$RUN_CV_HEATMAP" = true ]; then
        echo "  - Generating CV heatmap..."
        Rscript --vanilla b_modules_for_Method_3/c_make_heatmap_with_CV_of_matrices.R \
            -d "$SCRIPT_DIR" \
            -f "$FASTA_TAG" \
            -c "counts" \
            -g "$gene_group_file"
    fi
    
    # Bar graphs
    if [ "$RUN_BAR_GRAPHS" = true ]; then
        echo "  - Generating bar graphs..."
        Rscript --vanilla b_modules_for_Method_3/d_make_BarGraph_of_matrices.R \
            -d "$SCRIPT_DIR" \
            -f "$FASTA_TAG" \
            -c "counts" \
            -g "$gene_group_file"
    fi
    
    echo ""
}

echo "STEP 2: Generating visualizations for gene groups"
echo ""

for gene_group in "${GENE_GROUPS[@]}"; do
    process_gene_group "$gene_group"
done

#==================================================================================
# COMPLETION
#==================================================================================

echo "======================================="
echo "HEATMAP GENERATION COMPLETE"
echo "======================================="
echo "Outputs:"
echo "  - Basic Heatmaps:   7_Heatmap_Outputs/I_Basic_Heatmap/"
echo "  - CV Heatmaps:      7_Heatmap_Outputs/II_Heatmap_with_CV/"
echo "  - Bar Graphs:       7_Heatmap_Outputs/III_Bar_Graphs/"
echo "======================================="
