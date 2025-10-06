#!/bin/bash

: << 'INSTRUCTIONS'
# ===============================================
# a_stringtie_method_2.2.sh
#
# Description:
# This script generates gene expression count matrices for multiple gene groups and sample sets.
# Step-by-step:
# 1. For each gene group, locate abundance files for all samples.
# 2. Extract gene names from a reference file.
# 3. For each count type (coverage, FPKM, TPM):
#    a. Extract gene name and count columns from each sample file.
#    b. Build matrices with gene names as rows and samples/organs as columns using a Python script.
# 4. Output matrices to a results directory.
INSTRUCTIONS

# ===============================================
# CONFIGURATION
# ===============================================
BASE_DIR="$PWD"
INPUTS_DIR="5_stringtie_WD/a_Method_2_RAW_RESULTs"
OUT_DIR="5_stringtie_WD/b_Method_2_COUNT_MATRICES"
mkdir -p "$OUT_DIR"

QUERY_AGAINST_ALL_SMELGENES=FALSE
FASTA_Groups=(
	#TEST
	SmelDMP_CDS_Control_Best
	#SmelGIF_with_Best_Control_Cyclo
	#SmelGRF_with_Best_Control_Cyclo
	#SmelGRF-GIF_with_Best_Control_Cyclo
    #SmelGIF_with_Cell_Cycle_Control_genes
    #SmelGRF_with_Cell_Cycle_Control_genes
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

GENENAME_COL=3      # Gene name column (TSV column header is 'reference', but this is actually the gene name)
COVERAGE_COL=7      # Coverage column (adjust if needed)
FPKM_COL=8          # FPKM column (adjust if needed)
TPM_COL=9           # TPM column (adjust if needed)

# ===============================================
# FUNCTIONS
# ===============================================

merge_group_counts() {
    local gene_group="$1"
    local gene_group_path="$2"
    local ref_tsv="$3"
    local group_name=$(basename "$gene_group_path")

    # Create output directory for this group
    mkdir -p "$OUT_DIR/$group_name"

    local tmpdir
    tmpdir=$(mktemp -d)

    # Collect abundance files in the specified order
    files=()
    for srr in "${SRR_LIST_PRJNA328564[@]}"; do
        if QUERY_AGAINST_ALL_SMELGENES=TRUE; then
            gene_group="All_SmelGenes"
            file_path="All_SmelGenes/$srr/${srr}_All_SmelGenes_gene_abundances_de_novo.tsv"
        else
            file_path=${gene_group_path}/$srr/${srr}_${gene_group}_gene_abundances_de_novo.tsv
        fi
        if [[ -n "$file_path" ]]; then
            files+=("$file_path")
        fi
    done
    
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No files found in $gene_group_path."
        rm -r "$tmpdir"
        return
    fi

    # Extract gene names from the reference column (column 3) of the reference file (skip header)
    tail -n +2 "${ref_tsv}" | cut -f1 > "$tmpdir/gene_names.txt"

    # Debug: Check if gene_names.txt has content
    echo "Gene names extracted: $(wc -l < "$tmpdir/gene_names.txt") lines."
    #echo "Last of gene names:"
    #head -5 "$tmpdir/gene_names.txt"
    #tail -n 5 "$tmpdir/gene_names.txt"

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
                # Extract gene name (column 3) and count column (tab-separated)
                tail -n +2 "$sample_file" | cut -f"$GENENAME_COL","$COUNT_COL" > "$tmpdir/${srr}.txt"
                sample_files+=("$tmpdir/${srr}.txt")
                #tail $tmpdir/${srr}.txt
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

            # Match gene names and append sample counts to the correct row
            python3 "$(dirname "$0")/matrix_builder.py" "$tmpdir/gene_names.txt" ${sample_files[*]}
        } > "$OUT_DIR/$group_name/${group_name}_${count}_counts_geneName_SRR.tsv"

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

            # Match gene names and append sample counts to the correct row
            python3 "$(dirname "$0")/matrix_builder.py" "$tmpdir/gene_names.txt" ${sample_files[*]}
        } > "$OUT_DIR/$group_name/${group_name}_${count}_counts_geneName_Organ.tsv"

        # Clean up temporary files for this count type
        rm -f "${sample_files[@]}"
    done
    
    rm -r "$tmpdir"
}

# ===============================================
# MAIN
# ===============================================
 
for Gene_Group in "${FASTA_Groups[@]}"; do
    REF_TSV="5_stringtie_WD/0_Ref_Boilerplate_TSVs/${Gene_Group}_REF_BOILERPLATE_TSV.tsv"
    echo -e "\nChecking for reference TSV: $REF_TSV"
    if [[ ! -f "$REF_TSV" ]]; then 
        echo "Reference TSV file not found: $REF_TSV. Skipping group $Gene_Group."
        continue
    else
        echo "Found reference TSV: $REF_TSV" 
    fi
    Gene_Group_PATH="$INPUTS_DIR/$Gene_Group"
    
    echo "Merging counts for gene group: $Gene_Group"
    merge_group_counts "$Gene_Group" "$Gene_Group_PATH" "$REF_TSV"
done

echo "All groups processed. Count matrices are located in: $OUT_DIR"
