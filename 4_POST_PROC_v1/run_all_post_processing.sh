#!/bin/bash
#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RNA-SEQ ANALYSIS PIPELINE
# Runs selected analyses across all configured methods and gene groups
#===============================================================================

#set -euo pipefail

#===============================================================================
# DIRECTORY PATHS and SOURCE UTILITIES
#===============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ANALYSIS_MODULES_DIR="$SCRIPT_DIR/main_modules/analysis_modules"
GENE_GROUPS_DIR="$SCRIPT_DIR/main_modules/gene_groups_csv"
SRR_CSV_DIR="$SCRIPT_DIR/main_modules/SRR_csv"
UTILITIES_DIR="$SCRIPT_DIR/main_modules/utilities"
LOGGING_UTILS="$BASE_DIR/modules/logging/logging_utils.sh"

# Source logging utilities first (provides log_info, log_step, log_error, etc.)
source "$LOGGING_UTILS"
# Source pipeline utilities (provides parse_srr_csv, run_method_analysis, etc.)
source "$UTILITIES_DIR/pipeline_utils.sh"

#===============================================================================
# GLOBAL SETTINGS
#===============================================================================

# System resources (24GB RAM available - prioritize accuracy over memory)
THREADS=12
ENABLE_GPU="TRUE"
ENABLE_GNU_PARALLEL="FALSE"
DESIRED_CPU_PER_JOB=4
AVAILABLE_RAM_GB=18

# Logging and output
CLEAR_LOGS="TRUE"
CLEAR_OUTPUT_FOLDER="TRUE"

# Export RAM info for R scripts
export AVAILABLE_RAM_GB

#===============================================================================
# DATASETS (comment/uncomment to enable/disable)
#===============================================================================

# Master reference (uncomment ONE)
MASTER_REFERENCES=(
    #"All_Smel_Genes"
    "Eggplant_V4.1_transcripts.function"
)

# Extract first enabled reference
for ref in "${MASTER_REFERENCES[@]}"; do
    MASTER_REFERENCE="$ref"
    break
done

# Gene groups to analyze
GENE_GROUPS=(
    "SmelDMPs"
    "SmelGRF-GIFs"
    #"Selected_GRF_GIF_Genes_vAll_GIFs"
)

# Each name corresponds to a CSV file in main_modules/SRR_csv/
# CSV format: SRR_ID,Organ,Notes
SRR_DATASETS=(
    "PRJNA328564"      # Main Dataset - Eggplant tissue atlas (PRJNA328564)
    "SAMN28540077"     # Chinese Dataset 1 - replicability
    "SAMN28540068"     # Chinese Dataset 2 - replicability
    #"PRJNA865018"     # SET_1: Good Dataset for SmelDMP GEA
    #"PRJNA941250"     # SET_2: Good Dataset for SmelDMP GEA
    #"OTHER_SRR_LIST"  # Other miscellaneous samples
)

#===============================================================================
# CONFIGURATION (comment/uncomment to enable/disable)
#===============================================================================

# ALignment Methods
METHODS=(
    #"M1_HISAT2_RefGuided"
    #"M2_HISAT2_DeNovo"
    #"M3_STAR_Align"
    "M4_Salmon_Saf"
    "M5_RSEM_Bowtie2"
)

# Analysis modules (comment/uncomment to enable/disable)
# Note: Tximport scripts are auto-selected based on method type
ANALYSES=(
    "Tximport_Salmon"          # For Salmon-based methods (M4) - auto-run based on method
    "Tximport_RSEM"            # For RSEM-based methods (M5) - auto-run based on method
    "Matrix_Creation"
    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"
    "Coexpression_using_WGCNA"
    #"Differential_Expression"
    #"Gene_Set_Enrichment"
    #"PCA_Dimensionality_Reduction"
    #"Sample_Correlation_Clustering"
    #"Tissue_Specificity"
)


#===============================================================================
# INITIALIZATION
#===============================================================================

# Build combined SRR list from selected datasets
SRR_COMBINED_LIST=()
for dataset in "${SRR_DATASETS[@]}"; do
    csv_file="$SRR_CSV_DIR/${dataset}.csv"
    if [[ -f "$csv_file" ]]; then
        mapfile -t -O "${#SRR_COMBINED_LIST[@]}" SRR_COMBINED_LIST < <(parse_srr_csv "$csv_file")
    else
        echo "Warning: SRR CSV not found: $csv_file"
    fi
done

# Validate arrays (after SRR list is built)
[[ -z "$MASTER_REFERENCE" ]] && { echo "ERROR: No master reference configured"; exit 1; }
[[ ${#METHODS[@]} -eq 0 ]] && { echo "ERROR: No methods configured"; exit 1; }
[[ ${#GENE_GROUPS[@]} -eq 0 ]] && { echo "ERROR: No gene groups configured"; exit 1; }
[[ ${#SRR_COMBINED_LIST[@]} -eq 0 ]] && { echo "ERROR: No SRR samples loaded from datasets"; exit 1; }

# Calculate parallel jobs
JOBS=1
[[ "$ENABLE_GNU_PARALLEL" == "TRUE" ]] && JOBS=$((THREADS / DESIRED_CPU_PER_JOB)) && [[ $JOBS -lt 1 ]] && JOBS=1

# Initialize conda environment
eval "$(conda shell.bash hook)"
conda activate gea 2>/dev/null || echo "Warning: conda env 'gea' not found, using current env"

# Source optional shared utilities
[[ -f "$BASE_DIR/modules/config.sh" ]] && source "$BASE_DIR/modules/config.sh"
[[ -f "$BASE_DIR/modules/shared_utils.sh" ]] && source "$BASE_DIR/modules/shared_utils.sh"

# Setup logging directories with ABSOLUTE paths (critical for subprocesses that change directories)
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
ERROR_WARN_DIR="$SCRIPT_DIR/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$SCRIPT_DIR/logs/software_catalogs"
GPU_LOG_DIR="$SCRIPT_DIR/logs/gpu_log"
export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR ERROR_WARN_DIR SOFTWARE_CATALOG_DIR GPU_LOG_DIR

# Initialize logging system
setup_logging "$CLEAR_LOGS"

# Export actual log FILE paths after setup_logging has set them (for subprocesses)
export LOG_FILE TIME_FILE SPACE_FILE SPACE_TIME_FILE ERROR_WARN_FILE SOFTWARE_FILE GPU_LOG_FILE

# Clear output folders if requested (before processing)
if [[ "$CLEAR_OUTPUT_FOLDER" == "TRUE" ]]; then
    log_info "Clearing output folders..."
    for method in "${METHODS[@]}"; do
        output_dir="$SCRIPT_DIR/$method/7_Figures_Outputs"
        [[ -d "$output_dir" ]] && rm -rf "$output_dir"/* 2>/dev/null || true
    done
fi

# Export environment variables for R scripts and subprocesses
export THREADS ENABLE_GPU ANALYSIS_MODULES_DIR GENE_GROUPS_DIR UTILITIES_DIR SRR_CSV_DIR

#===============================================================================
# MAIN EXECUTION
#===============================================================================

log_step "Starting Post-Processing Pipeline"
log_info "Threads: $THREADS | GPU: $ENABLE_GPU | Parallel: $ENABLE_GNU_PARALLEL (Jobs: $JOBS)"
log_info "Master Reference: $MASTER_REFERENCE"
log_info "Methods: ${METHODS[*]}"
log_info "Gene Groups: ${GENE_GROUPS[*]}"
log_info "SRR Datasets: ${SRR_DATASETS[*]} (${#SRR_COMBINED_LIST[@]} samples)"
log_info "Analyses: ${ANALYSES[*]}"
log_info "Note: Method-specific preprocessing auto-detected (HISAT2->StringTie, Salmon->tximport, RSEM->tximport)"

# Run methods (parallel or sequential)
if [[ "$ENABLE_GNU_PARALLEL" == "TRUE" && $JOBS -gt 1 ]] && command -v parallel &>/dev/null; then
    log_info "Using GNU Parallel with $JOBS jobs"
    
    # Export for parallel subshells
    export_utils_for_parallel
    export SCRIPT_DIR LOG_FILE ERROR_WARN_FILE RUN_ID ANALYSIS_MODULES_DIR GENE_GROUPS_DIR
    export MASTER_REFERENCE THREADS ENABLE_GPU
    export GENE_GROUPS_STR="${GENE_GROUPS[*]}" ANALYSES_STR="${ANALYSES[*]}"
    
    printf '%s\n' "${METHODS[@]}" | parallel -j "$JOBS" run_method_analysis {} "$MASTER_REFERENCE"
else
    log_info "Running methods sequentially"
    for method in "${METHODS[@]}"; do
        run_method_analysis "$method" "$MASTER_REFERENCE"
    done
fi

#===============================================================================
# SUMMARY
#===============================================================================

log_step "Post-Processing Complete"
log_info "Methods processed: ${#METHODS[@]}"
log_info "Log file: $LOG_FILE"
log_info "Time metrics: $TIME_FILE"

if [[ -f "$ERROR_WARN_FILE" && -s "$ERROR_WARN_FILE" ]]; then
    log_info "Errors/Warnings: $(wc -l < "$ERROR_WARN_FILE") (see $ERROR_WARN_FILE)"
else
    log_info "No errors encountered"
fi

echo -e "\n========================================"
echo "Pipeline completed at $(date)"
echo "========================================"