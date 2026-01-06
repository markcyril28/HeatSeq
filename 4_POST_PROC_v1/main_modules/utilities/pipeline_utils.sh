#!/bin/bash
#===============================================================================
# PIPELINE UTILITIES - RNA-SEQ ANALYSIS
# Shared functions for post-processing scripts
# Note: Logging functions are provided by logging_utils.sh
#===============================================================================

#===============================================================================
# CSV PARSING
#===============================================================================

# Parse CSV and extract SRR entries as SRR_ID:Organ format
# Usage: mapfile -t SRR_LIST < <(parse_srr_csv "path/to/file.csv")
parse_srr_csv() {
    local csv_file="$1"
    [[ ! -f "$csv_file" ]] && { echo "Warning: CSV not found: $csv_file" >&2; return; }
    
    while IFS=',' read -r srr_id organ notes || [[ -n "$srr_id" ]]; do
        [[ "$srr_id" =~ ^#.*$ || "$srr_id" == "SRR_ID" || -z "$srr_id" ]] && continue
        echo "${srr_id}:${organ}"
    done < "$csv_file"
}

#===============================================================================
# ANALYSIS SCRIPT MAPPING
#===============================================================================

# Map analysis name to R script filename
get_analysis_script() {
    local -A scripts=(
        ["Matrix_Creation"]="3_Matrix_Creation.R"
        ["Basic_Heatmap"]="4_Basic_Heatmap.R"
        ["Heatmap_with_CV"]="5_Heatmap_with_CV.R"
        ["BarGraph"]="6_BarGraph.R"
        ["Coexpression_using_WGCNA"]="7_Coexpression_WGCNA.R"
        ["Differential_Expression"]="8_Differential_Expression.R"
        ["Gene_Set_Enrichment"]="9_Gene_Set_Enrichment.R"
        ["PCA_Dimensionality_Reduction"]="10_PCA_Dimensionality_Reduction.R"
        ["Sample_Correlation_Clustering"]="11_Sample_Correlation_Clustering.R"
        ["Tissue_Specificity"]="12_Tissue_Specificity.R"
        # Data import scripts
        ["Tximport_Salmon"]="a_tximport_salmon_to_matrices.R"
        ["Tximport_RSEM"]="a_tximport_rsem_to_matrices.R"
    )
    echo "${scripts[$1]:-}"
}

# Get method-specific preprocessing script
# Usage: get_preprocessing_script "method_name"
# Returns the script that converts raw quant output to count matrices
get_preprocessing_script() {
    local method=$1
    case "$method" in
        "M1_HISAT2_RefGuided")
            # Reference-guided also uses StringTie for quantification
            echo "stringtie_matrix_builder.sh"
            ;;
        "M2_HISAT2_DeNovo")
            echo "stringtie_matrix_builder.sh"
            ;;
        "M3_STAR_Align")
            # STAR alignment - quantification via featureCounts or RSEM
            echo ""  # TODO: implement star_matrix_builder.sh
            ;;
        "M4_Salmon_Saf")
            echo "a_tximport_salmon_to_matrices.R"
            ;;
        "M5_RSEM_Bowtie2")
            echo "a_tximport_rsem_to_matrices.R"
            ;;
        *)
            echo ""
            ;;
    esac
}

#===============================================================================
# METHOD ANALYSIS RUNNER
#===============================================================================

# Run analysis for a single method
# Usage: run_method_analysis "method_name" "master_reference"
run_method_analysis() {
    local method=$1 master_ref=$2
    local method_dir="$SCRIPT_DIR/$method"
    
    [[ ! -d "$method_dir" ]] && { log_error "Method directory not found: $method_dir"; return 1; }
    
    log_step "Processing Method: $method"
    pushd "$method_dir" > /dev/null
    
    export CURRENT_METHOD="$method" MASTER_REFERENCE="$master_ref"
    
    # Export GENE_GROUPS_DIR as absolute path for R scripts
    export GENE_GROUPS_DIR="$SCRIPT_DIR/main_modules/gene_groups_csv"
    
    # Setup temp config files for R scripts
    local modules_dir="B_${method#4}_modules"
    [[ -d "$modules_dir" ]] || modules_dir="."
    
    printf '%s\n' "${GENE_GROUPS[@]}" > "$modules_dir/.gene_groups_temp.txt"
    echo "$master_ref" > "$modules_dir/.master_reference_temp.txt"
    echo "TRUE" > "$modules_dir/.overwrite_temp.txt"
    
    # Rebuild GENE_GROUPS array from exported string in parallel mode
    [[ -n "${GENE_GROUPS_STR:-}" ]] && IFS=' ' read -ra GENE_GROUPS <<< "$GENE_GROUPS_STR"
    
    # Rebuild ANALYSES array from exported string in parallel mode
    [[ -n "${ANALYSES_STR:-}" ]] && IFS=' ' read -ra ANALYSES <<< "$ANALYSES_STR"
    
    # Run method-specific preprocessing if needed
    local preprocess_script=$(get_preprocessing_script "$method")
    if [[ -n "$preprocess_script" ]]; then
        local preprocess_path=""
        # Check in utilities directory first
        if [[ -f "$UTILITIES_DIR/$preprocess_script" ]]; then
            preprocess_path="$UTILITIES_DIR/$preprocess_script"
        elif [[ -f "$ANALYSIS_MODULES_DIR/$preprocess_script" ]]; then
            preprocess_path="$ANALYSIS_MODULES_DIR/$preprocess_script"
        elif [[ -f "$modules_dir/$preprocess_script" ]]; then
            preprocess_path="$modules_dir/$preprocess_script"
        fi
        
        if [[ -n "$preprocess_path" ]]; then
            log_info "Running preprocessing: $preprocess_script"
            if [[ "$preprocess_script" == *.R ]]; then
                run_with_error_capture Rscript "$preprocess_path" || log_error "Failed: preprocessing"
            else
                run_with_error_capture bash "$preprocess_path" || log_error "Failed: preprocessing"
            fi
        fi
    fi
    
    # Run each enabled analysis
    for analysis in "${ANALYSES[@]}"; do
        [[ -z "$analysis" ]] && continue
        
        # Skip tximport analyses if not applicable to method
        if [[ "$analysis" == "Tximport_Salmon" && ! "$method" =~ Salmon|M4 ]]; then
            continue
        fi
        if [[ "$analysis" == "Tximport_RSEM" && ! "$method" =~ RSEM|M5 ]]; then
            continue
        fi
        
        local script=$(get_analysis_script "$analysis")
        local script_path="$ANALYSIS_MODULES_DIR/$script"
        
        if [[ -f "$script_path" ]]; then
            log_info "Running: $analysis"
            run_with_error_capture Rscript "$script_path" || log_error "Failed: $analysis"
        elif [[ -f "$modules_dir/${script%.R}.R" ]]; then
            log_info "Running method script: $analysis"
            run_with_error_capture Rscript "$modules_dir/${script%.R}.R" || log_error "Failed: $analysis"
        else
            log_warn "Script not found for: $analysis"
        fi
    done
    
    popd > /dev/null
    log_info "Method $method complete"
}

#===============================================================================
# EXPORT FUNCTIONS FOR PARALLEL
#===============================================================================

export_utils_for_parallel() {
    # Export logging functions (from logging_utils.sh)
    export -f log log_info log_warn log_error log_step timestamp
    # Export error capture (from logging_utils.sh)
    export -f run_with_error_capture capture_stderr_errors 2>/dev/null || true
    # Export pipeline functions
    export -f run_method_analysis get_analysis_script get_preprocessing_script parse_srr_csv
}
