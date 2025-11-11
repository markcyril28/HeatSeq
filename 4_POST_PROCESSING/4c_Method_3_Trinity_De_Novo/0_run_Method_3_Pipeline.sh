#!/bin/bash
# Method 3: Trinity De Novo Pipeline Main Script
# Runs complete pipeline from assembly to visualization

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Source logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories to use local logs folder
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_files"
LOG_FILE="$LOG_DIR/Method_3_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/Method_3_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"

# Setup logging
setup_logging

# Source main pipeline functions
source "$BASE_DIR/modules/methods.sh"

# Configuration
FASTA_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes.fasta"

SRR_LIST=(
    "SRR3884597"
    "SRR3884684"
    "SRR3884686"
    "SRR3884687"
)

log_step "Starting Method 3: Trinity De Novo Pipeline"

# Run Trinity de novo pipeline
trinity_de_novo_alignment_pipeline \
    --FASTA "$FASTA_FILE" \
    --RNASEQ_LIST "${SRR_LIST[@]}"

log_info "Method 3 pipeline complete"
