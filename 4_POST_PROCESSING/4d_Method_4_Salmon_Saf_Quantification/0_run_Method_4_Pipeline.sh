#!/bin/bash
# Method 4: Salmon SAF Quantification Pipeline Main Script
# Runs complete pipeline from indexing to quantification

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Source logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories to use local logs folder
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
LOG_FILE="$LOG_DIR/Method_4_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/Method_4_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/Method_4_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/Method_4_${RUN_ID}_combined_metrics.csv"

# Setup logging
setup_logging

# Source main pipeline functions
source "$BASE_DIR/modules/methods.sh"

# Configuration
FASTA_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes.fasta"
GENOME_FILE="$BASE_DIR/0_INPUT_FASTAs/Eggplant_genome.fa"  # Update with actual genome file

SRR_LIST=(
    "SRR3884597"
    "SRR3884684"
    "SRR3884686"
    "SRR3884687"
)

log_step "Starting Method 4: Salmon SAF Quantification Pipeline"

# Run Salmon SAF pipeline
salmon_saf_pipeline \
    --FASTA "$FASTA_FILE" \
    --GENOME "$GENOME_FILE" \
    --RNASEQ_LIST "${SRR_LIST[@]}"

log_info "Method 4 pipeline complete"
