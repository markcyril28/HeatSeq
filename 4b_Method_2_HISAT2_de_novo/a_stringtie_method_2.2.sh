#!/bin/bash
# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
OUT_DIR="5_stringtie/a_METHOD2_RESULTS_matrices"
mkdir -p "$OUT_DIR"
rm -rf "$OUT_DIR"/*
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
	SRR3884687 # Buds, Opened Buds
	SRR3884686 # Buds 0.7 cm
	SRR3884689 # Leaves
	SRR3884690 # Stems
	SRR3884685 # Radicles
	SRR3884675 # Roots
)

# SRR to Organ mapping
declare -A SRR_TO_ORGAN=(
	["SRR3884631"]="Fruits_6cm"
	["SRR3884677"]="Cotyledons"
	["SRR3884679"]="Pistils"
	["SRR3884597"]="Flowers"
	["SRR3884687"]="Buds_Opened"
	["SRR3884686"]="Buds_0.7cm"
	["SRR3884689"]="Leaves"
	["SRR3884690"]="Stems"
	["SRR3884685"]="Radicles"
	["SRR3884675"]="Roots"
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
GENENAME_COL=3      # Gene.Name column - verify this is correct!
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

    # Extract Gene.ID and Gene.Name columns from the reference file (skip header)
    tail -n +2 "${REF_TSV}" | cut -f"$GENE_COL" > "$tmpdir/gene_ids.txt"
    tail -n +2 "${REF_TSV}" | cut -f"$GENENAME_COL" > "$tmpdir/gene_names.txt"

    # Debug: Check if gene_names.txt has content
    echo "Gene names extracted: $(wc -l < "$tmpdir/gene_names.txt") lines"
    echo "First few gene names: $(head -3 "$tmpdir/gene_names.txt")"

    for count in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        # Create array to store sample file paths for paste command
        sample_files=()

        # For each sample file, extract counts (skip header)
        for f in "${files[@]}"; do
            sample=$(basename "$f" | cut -d_ -f1)  # e.g., SRR3884597
            echo "  Adding sample: $sample"
            tail -n +2 "$f" | cut -f"$COUNT_COL" > "$tmpdir/${sample}.txt"
            sample_files+=("$tmpdir/${sample}.txt")
        done

        # Create Gene ID matrix with SRR headers
        {
            printf "GeneID"
            for f in "${files[@]}"; do
                sample=$(basename "$f" | cut -d_ -f1)
                printf "\t%s" "$sample"
            done
            printf "\n"

            paste "$tmpdir/gene_ids.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
        } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneID_SRR.tsv"

        # Create Gene ID matrix with Organ headers
        {
            printf "GeneID"
            for f in "${files[@]}"; do
                sample=$(basename "$f" | cut -d_ -f1)
                organ="${SRR_TO_ORGAN[$sample]}"
                printf "\t%s" "$organ"
            done
            printf "\n"

            paste "$tmpdir/gene_ids.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
        } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneID_Organ.tsv"

        # Create Gene Name matrix with SRR headers
        {
            printf "GeneName"
            for f in "${files[@]}"; do
                sample=$(basename "$f" | cut -d_ -f1)
                printf "\t%s" "$sample"
            done
            printf "\n"

            paste "$tmpdir/gene_names.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
        } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_SRR.tsv"

        # Create Gene Name matrix with Organ headers
        {
            printf "GeneName"
            for f in "${files[@]}"; do
                sample=$(basename "$f" | cut -d_ -f1)
                organ="${SRR_TO_ORGAN[$sample]}"
                printf "\t%s" "$organ"
            done
            printf "\n"

            paste "$tmpdir/gene_names.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
        } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_Organ.tsv"

        echo "Matrix saved: $OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneID_SRR.tsv"
        echo "Matrix saved: $OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneID_Organ.tsv"
        echo "Matrix saved: $OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_SRR.tsv"
        echo "Matrix saved: $OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_Organ.tsv"
        
        # Clean up temp files for this count type
        rm -f "${sample_files[@]}"
    done
    
    rm -r "$tmpdir"
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
