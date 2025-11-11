#!/bin/bash
# ==============================================================================
# HEATMAP WRAPPER FOR TRINITY + SALMON TXIMPORT RESULTS
# ==============================================================================
# Description: Generate visualizations from Trinity + Salmon + tximport results
# Author: Mark Cyril R. Mercado
# Version: v2 (Updated for Salmon/tximport workflow)
# Date: November 2025
# ==============================================================================

set -euo pipefail

# Navigate to Method 3 directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Configuration
BASE_DIR="$(pwd)"
FASTA_TAG="All_Smel_Genes"  # Update based on your FASTA file
GENE_GROUPS_DIR="$BASE_DIR/a_gene_groups_input_list"
MATRIX_DIR="$BASE_DIR/6_matrices_from_stringtie"
VISUALIZATION_SCRIPT="$BASE_DIR/b_modules_for_Method_3/e_visualize_tximport_results.R"

echo "========================================================================"
echo "TRINITY + SALMON VISUALIZATION WRAPPER"
echo "========================================================================"
echo "Base directory: $BASE_DIR"
echo "FASTA tag: $FASTA_TAG"
echo "Gene groups: $GENE_GROUPS_DIR"
echo ""

# Check if tximport was run
if [[ ! -f "$MATRIX_DIR/deseq2_dataset_trinity.rds" ]]; then
    echo "ERROR: DESeq2 dataset not found!"
    echo "Please run tximport first:"
    echo "  Rscript $MATRIX_DIR/run_tximport_trinity_salmon.R"
    exit 1
fi

# Check if visualization script exists
if [[ ! -f "$VISUALIZATION_SCRIPT" ]]; then
    echo "ERROR: Visualization script not found: $VISUALIZATION_SCRIPT"
    echo "Please ensure e_visualize_tximport_results.R exists"
    exit 1
fi

# Check for gene group files
if [[ ! -d "$GENE_GROUPS_DIR" ]]; then
    echo "ERROR: Gene groups directory not found: $GENE_GROUPS_DIR"
    exit 1
fi

GENE_GROUP_FILES=("$GENE_GROUPS_DIR"/*.txt)
if [[ ${#GENE_GROUP_FILES[@]} -eq 0 ]] || [[ ! -f "${GENE_GROUP_FILES[0]}" ]]; then
    echo "ERROR: No gene group files (*.txt) found in: $GENE_GROUPS_DIR"
    exit 1
fi

echo "Found ${#GENE_GROUP_FILES[@]} gene group file(s)"
echo ""

# Process each gene group
for gene_group_file in "${GENE_GROUP_FILES[@]}"; do
    [[ ! -f "$gene_group_file" ]] && continue
    
    group_name=$(basename "$gene_group_file" .txt)
    echo "------------------------------------------------------------------------"
    echo "Processing gene group: $group_name"
    echo "------------------------------------------------------------------------"
    
    # Run visualization script
    if Rscript --vanilla "$VISUALIZATION_SCRIPT" \
        -d "$BASE_DIR" \
        -f "$FASTA_TAG" \
        -g "$gene_group_file" \
        -t "all"; then
        echo "✓ Successfully generated visualizations for $group_name"
    else
        echo "✗ Failed to generate visualizations for $group_name"
    fi
    echo ""
done

echo "========================================================================"
echo "VISUALIZATION COMPLETED"
echo "========================================================================"
echo "Check outputs in:"
echo "  - Basic heatmaps: $BASE_DIR/7_Heatmap_Outputs/I_Basic_Heatmap"
echo "  - CV heatmaps: $BASE_DIR/7_Heatmap_Outputs/II_Heatmap_with_CV"
echo "  - Bar graphs: $BASE_DIR/7_Heatmap_Outputs/III_Bar_Graphs"
echo "========================================================================"

exit 0
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
