#!/bin/bash
# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
OUT_DIR="5_stringtie/a_METHOD2_RESULTS_matrices"
mkdir -p "$OUT_DIR"

INPUTS_DIR="5_stringtie/a_Method2_RESULTS"

# Fasta Group 
Fasta_Groups=(
	TEST
	SmelDMP_CDS_Control_Best
	SmelGIF_with_Best_Control_Cyclo
	SmelGRF_with_Best_Control_Cyclo
	SmelGRF-GIF_with_Best_Control_Cyclo
)

SRR_LIST_PRJNA328564=(
	# List of SRR sample IDs to process
	SRR3884631 # Fruits 6 cm
	#SRR3884664 # Fruits Calyx Stage 2
	#SRR3884653 # Fruits Flesh Stage 2
	SRR3884677 # Cotyledons
	SRR3884679 # Pistils
	SRR3884597 # Flowers
	SRR3884686 # Buds 0.7 cm
	SRR3884687 # Buds, Opened Buds
	SRR3884689 # Leaves
	SRR3884690 # Stems
	SRR3884685 # Radicles
	SRR3884675 # Roots
)

: << 'SCRATCH'
REF_TSVs_v2=(
    "$INPUTS_DIR/TEST/SRR3884686/SRR3884686_TEST_gene_abundances_de_novo_v2.tsv"
    "$INPUTS_DIR/SmelDMP_CDS_Control_Best/SRR3884686/SRR3884686_SmelDMP_CDS_Control_Best_gene_abundances_de_novo_v2.tsv"
    "$INPUTS_DIR/SmelGIF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGIF_with_Best_Control_Cyclo_gene_abundances_de_novo_v2.tsv"
    "$INPUTS_DIR/SmelGRF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGRF_with_Best_Control_Cyclo_gene_abundances_de_novo_v2.tsv"
    "$INPUTS_DIR/SmelGRF-GIF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGRF-GIF_with_Best_Control_Cyclo_gene_abundances_de_novo_v2.tsv"
)

REF_TSVs_v1=(
    "$INPUTS_DIR/TEST/SRR3884686/SRR3884686_TEST_gene_abundances_de_novo_v1.tsv"
    "$INPUTS_DIR/SmelDMP_CDS_Control_Best/SRR3884686/SRR3884686_SmelDMP_CDS_Control_Best_gene_abundances_de_novo_v1.tsv"
    "$INPUTS_DIR/SmelGIF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGIF_with_Best_Control_Cyclo_gene_abundances_de_novo_v1.tsv"
    "$INPUTS_DIR/SmelGRF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGRF_with_Best_Control_Cyclo_gene_abundances_de_novo_v1.tsv"
    "$INPUTS_DIR/SmelGRF-GIF_with_Best_Control_Cyclo/SRR3884686/SRR3884686_SmelGRF-GIF_with_Best_Control_Cyclo_gene_abundances_de_novo_v1.tsv"
)

SCRATCH

# Check your files: typical columns include "Gene.ID", "Gene.Name", "Coverage", "FPKM", "TPM"
# Adjust COLNUM to the *raw counts* column (example: 2=Gene.ID, 7=Coverage, etc.)
GENE_COL=1          # Gene.ID column
COVERAGE_COL=7      # Coverage column (adjust if needed)
FPKM_COL=8          # FPKM column (adjust if needed)
TPM_COL=9           # TPM column (adjust if needed)

# ===============================================
# FUNCTIONS
# ===============================================

merge_group_counts() {
    local gene_group_path="$1"
    local REF_TSV="$2"
    local version="$3"

    local group_name=$(basename "$gene_group_path")

    # Create output directory for this group and version
    mkdir -p "$OUT_DIR/$group_name/$version"

    local tmpdir
    tmpdir=$(mktemp -d)

    # Collect abundance files
    mapfile -t files < <(find "$gene_group_path" -type f -name "*gene_abundances_de_novo_${version}.tsv" | sort)
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No files found in $gene_group_path"
        rm -r "$tmpdir"
        return
    fi

    # Extract Gene.ID column from the first file
    cut -f"$GENE_COL" "${REF_TSV}" > "$tmpdir/gene_ids.txt"

    for count in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

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
        } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts.tsv"

        echo "Matrix saved: $OUT_DIR/$group_name/$version/${group_name}_${count}_counts.tsv"
        
        # Clean up temp files for this count type
        #rm -f "$tmpdir"/*.txt
    done
    
    #rm -r "$tmpdir"
}

# ===============================================
# MAIN
# ===============================================

for version in v1 v2; do
    for Gene_group in "${Fasta_Groups[@]}"; do
        REF_TSV="$INPUTS_DIR/$Gene_group/SRR3884686/SRR3884686_${Gene_group}_gene_abundances_de_novo_${version}.tsv"
        echo "Checking for reference TSV: $REF_TSV"
        Gene_group_path="$INPUTS_DIR/$Gene_group"
        
        echo "Merging counts for GENE_GROUP: $Gene_group, VERSION: $version"
        merge_group_counts "$Gene_group_path" "$REF_TSV" "$version"
    done
done

echo "All groups processed. Count matrices are in: $OUT_DIR"
