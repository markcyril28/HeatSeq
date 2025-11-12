#!/bin/bash

#===============================================================================
# METHOD 4: SALMON SAF QUANTIFICATION PIPELINE - MAIN SCRIPT
#===============================================================================
# Purpose: Executes complete Salmon quantification workflow
# Steps:
#   1. Index creation with decoy-aware Salmon
#   2. Quantification of RNA-seq samples
#   3. Output ready for tximport/DESeq2 processing
# Usage: bash 0_run_Method_4_Pipeline.sh
# Note: Ensure Salmon is installed and accessible in PATH
#===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Generate unique run identifier if not already set
export RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"

# Source centralized logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Configure log directories for Method 4 (local to this method)
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
LOG_FILE="$LOG_DIR/Method_4_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/Method_4_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/Method_4_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/Method_4_${RUN_ID}_combined_metrics.csv"

# Initialize logging system
setup_logging

# Source main pipeline functions from centralized methods module
source "$BASE_DIR/modules/methods.sh"

#===============================================================================
# CONFIGURATION
#===============================================================================

# Input files
FASTA_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes.fasta"
GENOME_FILE="$BASE_DIR/0_INPUT_FASTAs/TEST.fa"  # TODO: Update path if needed

# RNA-seq samples to quantify
SRR_LIST=(
    "SRR3884597"  # Flowers
    "SRR3884684"  # Senescent leaves
    "SRR3884686"  # Buds
    "SRR3884687"  # Opened buds
)

#===============================================================================
# PIPELINE EXECUTION
#===============================================================================

log_step "Initiating Method 4: Salmon SAF Quantification Pipeline"
log_info "Processing ${#SRR_LIST[@]} RNA-seq samples"
log_info "Reference: $(basename "$FASTA_FILE")"
log_info "Genome: $(basename "$GENOME_FILE")"

# Validate input files exist
if [ ! -f "$FASTA_FILE" ]; then
    log_error "FASTA file not found: $FASTA_FILE"
    exit 1
fi

if [ ! -f "$GENOME_FILE" ]; then
    log_error "Genome file not found: $GENOME_FILE"
    exit 1
fi

log_info "Input files validated successfully"

# Execute Salmon SAF quantification pipeline
# This function handles indexing and quantification
log_info "Starting Salmon quantification..."
salmon_saf_pipeline \
    --FASTA "$FASTA_FILE" \
    --GENOME "$GENOME_FILE" \
    --RNASEQ_LIST "${SRR_LIST[@]}"

log_step "Method 4 pipeline completed successfully"
log_info "Quantification files ready for tximport processing"
log_info "Output location: $SCRIPT_DIR/5_Salmon_Quant_WD/"
log_info "Next step: Run 1_run_Heatmap_Wrapper.sh for visualization"
