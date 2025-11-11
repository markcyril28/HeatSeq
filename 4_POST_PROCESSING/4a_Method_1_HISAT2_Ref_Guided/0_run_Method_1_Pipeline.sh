#!/bin/bash
# Method 1: HISAT2 Reference-Guided Pipeline Main Script
# Runs complete pipeline from alignment to visualization

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Source logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories to use local logs folder
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_files"
LOG_FILE="$LOG_DIR/Method_1_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/Method_1_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"

# Setup logging
setup_logging

# Source main pipeline functions
source "$BASE_DIR/modules/methods.sh"

# Configuration
FASTA_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes_Full_Name.fasta"
GTF_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes_transcript_based.gtf"

SRR_LIST=(
    "SRR3884597"
    "SRR3884684"
    "SRR3884686"
    "SRR3884687"
)

log_step "Starting Method 1: HISAT2 Reference-Guided Pipeline"

# Run HISAT2 reference-guided pipeline
hisat2_ref_guided_pipeline \
    --FASTA "$FASTA_FILE" \
    --GTF "$GTF_FILE" \
    --RNASEQ_LIST "${SRR_LIST[@]}"

log_info "Method 1 pipeline complete"
