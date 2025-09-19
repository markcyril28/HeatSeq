#!/bin/bash
set -euo pipefail

# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
OUT_DIR="03_stringtie/method_2.2_count_matrices"
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
    local group="$1"
    local group_dir="$BASE_DIR/$group"

    echo "Processing group: $group"
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

# Loop through each group directory under BASE_DIR
for group in $(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n"); do
    merge_group_counts "$group"
done

echo "All groups processed. Count matrices are in: $OUT_DIR"

conda install -c conda-forge -c bioconda r-base r-tidyverse bioconductor-deseq2 r-pheatmap
