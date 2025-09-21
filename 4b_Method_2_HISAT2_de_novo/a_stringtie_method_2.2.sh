#!/bin/bash
# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
OUT_DIR="5_stringtie/a_METHOD2_RESULTS_matrices"
mkdir -p "$OUT_DIR"
rm -rf "$OUT_DIR"/*
INPUTS_DIR="5_stringtie/a_Method2_RESULTS"

Fasta_Groups=(
	#TEST
	#SmelDMP_CDS_Control_Best
	SmelGIF_with_Best_Control_Cyclo
	SmelGRF_with_Best_Control_Cyclo
	SmelGRF-GIF_with_Best_Control_Cyclo
)

SRR_LIST_PRJNA328564=(
	SRR3884631 # Fruits 6 cm
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

GENENAME_COL=3      # Gene.Name column: tsv column header is reference but this is actually the geneName
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

    # Collect abundance files in the specified order
    files=()
    for srr in "${SRR_LIST_PRJNA328564[@]}"; do
        file_path=$(find "$gene_group_path" -type f -name "${srr}_*gene_abundances_de_novo_${version}.tsv")
        if [[ -n "$file_path" ]]; then
            files+=("$file_path")
        fi
    done
    
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No files found in $gene_group_path"
        rm -r "$tmpdir"
        return
    fi

    # Extract geneName from reference column (column 3) of the reference file (skip header)
    tail -n +2 "${REF_TSV}" | cut -f1 > "$tmpdir/gene_names.txt"

    # Debug: Check if gene_names.txt has content
    echo "Gene names extracted: $(wc -l < "$tmpdir/gene_names.txt") lines"
    echo "First few gene names: $(head -5 "$tmpdir/gene_names.txt")"

    for count in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        # Create array to store sample file paths for paste command
        sample_files=()

        # For each sample in the specified order, extract counts (skip header)
        for srr in "${SRR_LIST_PRJNA328564[@]}"; do
            # Find the corresponding file for this SRR
            sample_file=""
            for f in "${files[@]}"; do
                if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                    sample_file="$f"
                    break
                fi
            done
            
            if [[ -n "$sample_file" ]]; then
                tail -n +2 "$sample_file" | cut -f"$COUNT_COL" > "$tmpdir/${srr}.txt"
                sample_files+=("$tmpdir/${srr}.txt")
            fi
        done

    # Create matrix: gene names + SRR sample columns
    {
            # Print header: GeneName and SRR IDs
            printf "GeneName"
            for srr in "${SRR_LIST_PRJNA328564[@]}"; do
                # Check if this SRR has a corresponding file
                # Add SRR ID to header if sample file exists
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        printf "\t%s" "$srr"
                        break
                    fi
                done
            done
            printf "\n"

            # Paste gene names and sample counts, fill missing values with 0
            paste "$tmpdir/gene_names.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
    } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_SRR.tsv"

    # Create matrix: gene names + Organ columns
    {
            # Print header: GeneName and Organ names
            printf "GeneName"
            for srr in "${SRR_LIST_PRJNA328564[@]}"; do
                # Check if this SRR has a corresponding file
                # Add Organ name to header if sample file exists
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        organ="${SRR_TO_ORGAN[$srr]}"
                        printf "\t%s" "$organ"
                        break
                    fi
                done
            done
            printf "\n"

            # Paste gene names and sample counts, fill missing values with 0
            paste "$tmpdir/gene_names.txt" "${sample_files[@]}" | sed 's/\t\t/\t0\t/g; s/\t$/\t0/; s/^\t/0\t/'
    } > "$OUT_DIR/$group_name/$version/${group_name}_${count}_counts_geneName_Organ.tsv"

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
        REF_TSV="5_stringtie/${Gene_group}_REF_BOILERPLATE_TSV.tsv"
        echo "Checking for reference TSV: $REF_TSV"
        Gene_group_path="$INPUTS_DIR/$Gene_group"
        
        echo "Merging counts for GENE_GROUP: $Gene_group, VERSION: $version"
        merge_group_counts "$Gene_group_path" "$REF_TSV" "$version"
    done
done

echo "All groups processed. Count matrices are in: $OUT_DIR"
