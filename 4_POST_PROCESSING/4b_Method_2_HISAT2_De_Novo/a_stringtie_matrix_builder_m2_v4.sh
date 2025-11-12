#!/bin/bash

# ===============================================
# StringTie Matrix Builder v4 - Method 2 HISAT2 De Novo
# ===============================================
# Description:
# Generates gene expression count matrices for multiple gene groups and sample sets.
# Process:
# 1. For each gene group, locate abundance files for all samples
# 2. Extract gene names from reference TSV files
# 3. Build matrices (coverage, FPKM, TPM) with genes as rows, samples/organs as columns
# 4. Output matrices to results directory
# ===============================================

# ===============================================
# CONFIGURATION
# ===============================================
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Gene groups to process
Gene_Groups_Boilerplates=(
    # Available gene groups based on TSV files
    "Best_Cell_Cycle_Associated_Control_Genes"
    "Best_Control_Genes"
    "SmelDMPs"
    "SmelDMPs_with_18s_rRNA"
    "SmelDMPs_with_18s_rRNA_PSBMB"
    "SmelGIFs"
    "SmelGRFs"
    "Selected_GRF_GIF_Genes_v1"
    "Selected_GRF_GIF_Genes_vAll_GIFs"
    "Selected_GRF_GIF_Genes_vTwo_GIFs"
    "SmelGIF_with_Cell_Cycle_Control_genes"
    "SmelGRF_with_Cell_Cycle_Control_genes"
    "SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Directories
BASE_DIR="$PWD"
INPUTS_DIR="5_stringtie_WD/a_Method_2_RAW_RESULTs"
OUT_DIR="5_stringtie_WD/b_Method_2_COUNT_MATRICES"
LOG_FILE="$OUT_DIR/logs/matrix_builder_$(date +%Y%m%d_%H%M%S).log"

# Create output directory and initialize logging
mkdir -p "$OUT_DIR/logs"
#exec 1> >(tee -a "$LOG_FILE")
#exec 2> >(tee -a "$LOG_FILE" >&2)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting StringTie Matrix Builder v4"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Working directory: $BASE_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input directory: $INPUTS_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file: $LOG_FILE"

# Query configuration
QUERY_AGAINST_MASTER_REFERENCE="TRUE"
#MASTER_REFERENCE="SmelGRF-GIF_with_Best_Control_Cyclo"
MASTER_REFERENCE="All_Smel_Genes"
MASTER_SUFFIX="_from_${MASTER_REFERENCE}"  # Suffix added when QUERY_AGAINST_MASTER_REFERENCE=TRUE

# Sample lists from different projects
SRR_LIST_PRJNA328564=(
    SRR3884685  # Radicles
    SRR3884677  # Cotyledons  
    SRR3884675  # Roots  
    SRR3884690  # Stems  
    SRR3884689  # Leaves_vegetative 
    SRR3884684  # Senescent_leaves 
    SRR3884686  # Buds_0.7cm  
    SRR3884687  # Opened_Buds  
    SRR3884597  # Flowers 
    SRR3884679  # Pistils 
    SRR3884608  # Fruits_1cm 
    #SRR3884620  # Fruits_Stage_1 
    SRR3884631  # Fruits_6cm 
    #SRR3884642  # Fruits_Skin_Stage_2 
    #SRR3884653  # Fruits_Flesh_Stage_2 
    #SRR3884664  # Fruits_Calyx_Stage_2 
    #SRR3884680  # Fruits_Skin_Stage_3 
    #SRR3884681  # Fruits_Flesh_Stage_3 
    SRR3884678  # Fruits_peduncle 
)

SRR_LIST_SAMN28540077=(
    SRR20722229  # Roots
    SRR20722227  # Stems
    SRR20722230  # Leaves
    SRR21010466  # Buds
    SRR20722234  # Flowers
    SRR2072232   # Fruits
)

# Source: SAMN28540068 project samples
SRR_LIST_SAMN28540068=(
    SRR31755282  # Roots
    SRR20722384  # Stems
    SRR20722386  # Leaves
    SRR20722297  # Buds
    SRR23909863  # Flowers
    SRR20722387  # Fruits
)


# Combined list of all samples across projects
SRR_LIST_COMBINED=(
    "${SRR_LIST_PRJNA328564[@]}" 
    #"${SRR_LIST_SAMN28540077[@]}" 
    #"${SRR_LIST_SAMN28540068[@]}"
)

SRR_LIST_COMBINED_OFF=(
    # Roots
    "SRR3884675"   # Roots_1 (PRJNA328564)
    "SRR20722229"  # Roots_2 (SAMN28540077)
    "SRR31755282"  # Roots_3 (SAMN28540068)

    # Stems
    "SRR3884690"   # Stems_1 (PRJNA328564)
    "SRR20722227"  # Stems_2 (SAMN28540077)
    "SRR20722384"  # Stems_3 (SAMN28540068)

    # Leaves
    "SRR3884689"   # Leaves_1 (PRJNA328564)
    "SRR20722230"  # Leaves_2 (SAMN28540077)
    "SRR20722386"  # Leaves_3 (SAMN28540068)
    "SRR3884684"   # Senescent_leaves (PRJNA328564)

    # Buds
    "SRR3884686"   # Buds_1 (PRJNA328564)
    "SRR21010466"  # Buds_2 (SAMN28540077)
    "SRR20722297"  # Buds_3 (SAMN28540068)

    # Opened Buds
    "SRR3884687"   # Opened_Buds_1 (PRJNA328564)

    # Flowers
    "SRR3884597"   # Flowers_1 (PRJNA328564)
    "SRR20722234"  # Flowers_2 (SAMN28540077)
    "SRR23909863"  # Flowers_3 (SAMN28540068)

    # Fruits
    "SRR3884631"   # Fruits_1 (PRJNA328564)
    "SRR2072232"   # Fruits_2 (SAMN28540077)
    "SRR20722387"  # Fruits_3 (SAMN28540068)
    "SRR3884608"   # Fruits_1cm (PRJNA328564)
    "SRR3884620"   # Fruits_Stage_1 (PRJNA328564)
    "SRR3884642"   # Fruits_Skin_Stage_2 (PRJNA328564)
    "SRR3884653"   # Fruits_Flesh_Stage_2 (PRJNA328564)
    "SRR3884664"   # Fruits_Calyx_Stage_2 (PRJNA328564)
    "SRR3884680"   # Fruits_Skin_Stage_3 (PRJNA328564)
    "SRR3884681"   # Fruits_Flesh_Stage_3 (PRJNA328564)
    "SRR3884678"   # Fruits_peduncle (PRJNA328564)

    # Other organs
    "SRR3884685"   # Radicles (PRJNA328564)
    "SRR3884677"   # Cotyledons (PRJNA328564)
    "SRR3884679"   # Pistils (PRJNA328564)
)

# Mapping from SRR IDs to organ names for matrix headers
declare -A SRR_TO_ORGAN=(
    # Roots
    ["SRR3884675"]="Roots"      # PRJNA328564
    ["SRR20722229"]="Roots_2"     # SAMN28540077
    ["SRR31755282"]="Roots_3"     # SAMN28540068

    # Stems
    ["SRR3884690"]="Stems"      # PRJNA328564
    ["SRR20722227"]="Stems_2"     # SAMN28540077
    ["SRR20722384"]="Stems_3"     # SAMN28540068

    # Leaves
    ["SRR3884689"]="Leaves"     # PRJNA328564
    ["SRR20722230"]="Leaves_2"    # SAMN28540077
    ["SRR20722386"]="Leaves_3"    # SAMN28540068
    ["SRR3884684"]="Senescent_leaves" # PRJNA328564

    # Buds
    ["SRR3884686"]="Buds"       # PRJNA328564
    ["SRR21010466"]="Buds_2"      # SAMN28540077
    ["SRR20722297"]="Buds_3"      # SAMN28540068

    # Opened Buds
    ["SRR3884687"]="Opened_Buds" # PRJNA328564

    # Flowers
    ["SRR3884597"]="Flowers"    # PRJNA328564
    ["SRR20722234"]="Flowers_2"   # SAMN28540077
    ["SRR23909863"]="Flowers_3"   # SAMN28540068

    # Fruits
    ["SRR3884631"]="Fruits"	  # PRJNA328564
    ["SRR2072232"]="Fruits_2"     # SAMN28540077
    ["SRR20722387"]="Fruits_3"    # SAMN28540068
    ["SRR3884608"]="Fruits_1cm"   # PRJNA328564
    ["SRR3884620"]="Fruits_Stage_1" # PRJNA328564
    ["SRR3884642"]="Fruits_Skin_Stage_2" # PRJNA328564
    ["SRR3884653"]="Fruits_Flesh_Stage_2" # PRJNA328564
    ["SRR3884664"]="Fruits_Calyx_Stage_2" # PRJNA328564
    ["SRR3884680"]="Fruits_Skin_Stage_3" # PRJNA328564
    ["SRR3884681"]="Fruits_Flesh_Stage_3" # PRJNA328564
    ["SRR3884678"]="Fruits_peduncle" # PRJNA328564

    # Other organs
    ["SRR3884685"]="Radicles"     # PRJNA328564
    ["SRR3884677"]="Cotyledons"   # PRJNA328564
    ["SRR3884679"]="Pistils"      # PRJNA328564
)

# StringTie abundance file column indices (1-based)
GENENAME_COL=3      # Gene name column ('reference' header contains gene names)
COVERAGE_COL=7      # Coverage values
FPKM_COL=8          # FPKM values  
TPM_COL=9           # TPM values

# ===============================================
# FUNCTIONS
# ===============================================

# Function: merge_group_counts
# Purpose: Process abundance files for a gene group and create count matrices
# Args: $1=gene_group_name $2=gene_group_path $3=reference_tsv_file
merge_group_counts() {
    local gene_group="$1"
    local gene_group_path="$2"
    local ref_tsv="$3"
    local group_name=$(basename "$gene_group_path")

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $group_name"
    
    # Create output directory for this group
    mkdir -p "$OUT_DIR/$group_name"

    local tmpdir
    tmpdir=$(mktemp -d)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Created temporary directory: $tmpdir"

    # Collect abundance files in the specified sample order
    files=()
    local files_found=0
    for srr in "${SRR_LIST_COMBINED[@]}"; do
        if [[ "$QUERY_AGAINST_MASTER_REFERENCE" == "TRUE" ]]; then
            #local gene_group="All_SmelGenes"
            local gene_group=$MASTER_REFERENCE
            local gene_group_path="$INPUTS_DIR/$gene_group"
        fi 
        file_path="$gene_group_path/$srr/${srr}_${gene_group}_gene_abundances_de_novo.tsv"
        if [[ -f "$file_path" ]]; then
            files+=("$file_path")
            ((files_found++))
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: File not found: $file_path"
        fi
    done
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $files_found abundance files for $group_name"
    
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: No abundance files found for $group_name"
        rm -r "$tmpdir"
        return 1
    fi

    # Extract gene names from the reference column (column 3) of the reference file (skip header)
    tail -n +2 "${ref_tsv}" | cut -f1 > "$tmpdir/gene_names.txt"

    # Debug: Check if gene_names.txt has content
    echo "Gene names extracted: $(wc -l < "$tmpdir/gene_names.txt") lines."
    #echo "Last of gene names:"
    #head -5 "$tmpdir/gene_names.txt"
    #tail -n 5 "$tmpdir/gene_names.txt"

    for count_type in coverage fpkm tpm; do
        local COUNT_COL_VAR="${count_type^^}_COL"
        local COUNT_COL="${!COUNT_COL_VAR}"

        # Extract count data from each sample file
        sample_files=()
        local samples_processed=0
        
        for srr in "${SRR_LIST_COMBINED[@]}"; do
            # Find the corresponding abundance file for this SRR
            sample_file=""
            for f in "${files[@]}"; do
                if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                    sample_file="$f"
                    break
                fi
            done
            
            if [[ -n "$sample_file" ]]; then
                # Extract gene name (column 3) and count column (tab-separated) from the sample file
                tail -n +2 "$sample_file" | cut -f"$GENENAME_COL","$COUNT_COL" > "$tmpdir/${srr}.txt"
                sample_files+=("$tmpdir/${srr}.txt")
                ((samples_processed++))
            fi
        done
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processed $samples_processed samples for $count_type"

        # GENERATING THE OUTPUT FILENAMES
        if [[ "$QUERY_AGAINST_MASTER_REFERENCE" == "TRUE" ]]; then
            local output_geneName_SRR_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR${MASTER_SUFFIX}.tsv"
        else
            local output_geneName_SRR_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_SRR.tsv"
        fi
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating SRR matrix: $(basename "$output_geneName_SRR_tsv")"
        
        # Create matrix: gene names + SRR sample columns
        {
            # Print header: GeneName and SRR IDs
            printf "GeneName"
            for srr in "${SRR_LIST_COMBINED[@]}"; do
                # Add SRR ID to header if corresponding file exists
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        printf "\t%s" "$srr"
                        break
                    fi
                done
            done
            printf "\n"

            # Build matrix using Python script
            python3 "$(dirname "$0")/matrix_builder.py" "$tmpdir/gene_names.txt" "${sample_files[@]}"
        } > "$output_geneName_SRR_tsv"
        
        # Generate organ matrix
        if [[ "$QUERY_AGAINST_MASTER_REFERENCE" == "TRUE" ]]; then
            local output_geneName_Organ_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ${MASTER_SUFFIX}.tsv"
        else
            local output_geneName_Organ_tsv="$OUT_DIR/$group_name/${group_name}_${count_type}_counts_geneName_Organ.tsv"
        fi
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating Organ matrix: $(basename "$output_geneName_Organ_tsv")"
        
        # Create matrix: gene names + Organ columns
        {
            # Print header: GeneName and Organ names
            printf "GeneName"
            for srr in "${SRR_LIST_COMBINED[@]}"; do
                # Add organ name to header if corresponding file exists
                for f in "${files[@]}"; do
                    if [[ "$(basename "$f")" == "${srr}_"* ]]; then
                        local organ="${SRR_TO_ORGAN[$srr]:-Unknown}"
                        printf "\t%s" "$organ"
                        break
                    fi
                done
            done
            printf "\n"

            # Build matrix using Python script
            python3 "$(dirname "$0")/matrix_builder.py" "$tmpdir/gene_names.txt" "${sample_files[@]}"
        } > "$output_geneName_Organ_tsv"

        # Clean up temporary files for this count type
        rm -f "${sample_files[@]}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed $count_type matrix generation"
    done
    
    # Cleanup temporary directory
    rm -r "$tmpdir"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed processing for $group_name"
}

# ===============================================
# MAIN EXECUTION
# ===============================================
 
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting count matrix generation for ${#Gene_Groups_Boilerplates[@]} gene groups"

# Process each gene group
for Gene_Group_Boilerplate in "${Gene_Groups_Boilerplates[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing gene group: $Gene_Group_Boilerplate"
    
    # Check for reference TSV file
    REF_TSV="5_stringtie_WD/0_Ref_Boilerplate_TSVs/${Gene_Group_Boilerplate}.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Looking for reference TSV: $REF_TSV"
    
    if [[ ! -f "$REF_TSV" ]]; then 
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: Reference TSV not found, skipping $Gene_Group_Boilerplate"
        continue
    fi
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found reference TSV with $(tail -n +2 "$REF_TSV" | wc -l) genes"
    
    # Set input path and process
    Gene_Group_Boilerplate_PATH="$INPUTS_DIR/$Gene_Group_Boilerplate"
    
    if merge_group_counts "$Gene_Group_Boilerplate" "$Gene_Group_Boilerplate_PATH" "$REF_TSV"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Successfully processed $Gene_Group_Boilerplate"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed to process $Gene_Group_Boilerplate"
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ========================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Matrix generation completed"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file: $LOG_FILE"
