#!/bin/bash
set -euo pipefail

# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
OUT_DIR="5_stringtie/method_2.2_count_matrices"
mkdir -p "$OUT_DIR"

# Which column to extract from the *_gene_abundances_de_novo_v2.tsv
# Check your files: typical columns include "Gene.ID", "Gene.Name", "Coverage", "FPKM", "TPM"
# Adjust COLNUM to the *raw counts* column (example: 2=Gene.ID, 7=Coverage, etc.)
GENE_COL=1     # Gene.ID column
COUNT_COL=7    # Coverage column (adjust if needed)

# ===============================================
# FUNCTIONS
# ===============================================

merge_group_counts() {
    local group_path="$1"
    local group=$(basename "$group_path")
    local group_dir="$group_path"

    local tmpdir
    tmpdir=$(mktemp -d)

    # Collect abundance files
    mapfile -t files < <(find "$group_dir" -type f -name "*gene_abundances_de_novo_v2.tsv" | sort)
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No files found in $group_dir"
        rm -r "$tmpdir"
        return
    fi

    # Extract Gene.ID column from the first file
    cut -f"$GENE_COL" "${files[0]}" > "$tmpdir/gene_ids.txt"

    # For each sample file, extract counts
    for f in "${files[@]}"; do
        sample=$(basename "$f" | cut -d_ -f1)  # e.g., SRR3884597
        echo "  Adding sample: $sample"
        cut -f"$COUNT_COL" "$f" > "$tmpdir/${sample}.txt"
    done

    # Paste all together
    {
        printf "GeneID"
        for f in "${files[@]}"; do
            sample=$(basename "$f" | cut -d_ -f1)
            printf "\t%s" "$sample"
        done
        printf "\n"

        paste "$tmpdir/gene_ids.txt" "$tmpdir"/*.txt
    } > "$OUT_DIR/${group}_counts.tsv"

    echo "Matrix saved: $OUT_DIR/${group}_counts.tsv"
    rm -r "$tmpdir"
}

# ===============================================
# MAIN
# ===============================================

# Find directories containing target files (excluding 6_Visualization)
# Use a simpler approach to avoid path issues
find "$BASE_DIR" -mindepth 1 -maxdepth 3 -type d -not -path "*/6_Visualization*" | while read -r group_path; do
    # Check if this directory contains target files
    if find "$group_path" -name "*gene_abundances_de_novo_v2.tsv" -type f | head -1 | grep -q .; then
        group=$(basename "$group_path")
        echo "Found group directory: $group_path"
        echo "Processing group: $group"
        merge_group_counts "$group_path"
    fi
done

echo "All groups processed. Count matrices are in: $OUT_DIR"

# To ready and run DESeq2 in R:
#conda install -c conda-forge -c bioconda r-base r-tidyverse bioconductor-deseq2 r-pheatmap

# Check if R script exists before running
R_SCRIPT="b_heatmap_DeSeq2_v2.R"
if [[ -f "$R_SCRIPT" ]]; then
    echo "Running R script: $R_SCRIPT"
    Rscript "$R_SCRIPT"
else
    echo "Warning: R script '$R_SCRIPT' not found in current directory"
fi