#!/bin/bash
# Method 5: Bowtie2/RSEM Quantification Pipeline Main Script
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
LOG_FILE="$LOG_DIR/Method_5_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/Method_5_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/Method_5_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/Method_5_${RUN_ID}_combined_metrics.csv"

# Setup logging
setup_logging

# Source main pipeline functions
source "$BASE_DIR/modules/methods.sh"

# Configuration
FASTA_FILE="$BASE_DIR/0_INPUT_FASTAs/All_Smel_Genes.fasta"

SRR_LIST=(
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
    #"SRR3884620"   # Fruits_Stage_1 (PRJNA328564)
    #"SRR3884642"   # Fruits_Skin_Stage_2 (PRJNA328564)
    #"SRR3884653"   # Fruits_Flesh_Stage_2 (PRJNA328564)
    #"SRR3884664"   # Fruits_Calyx_Stage_2 (PRJNA328564)
    #"SRR3884680"   # Fruits_Skin_Stage_3 (PRJNA328564)
    #"SRR3884681"   # Fruits_Flesh_Stage_3 (PRJNA328564)
    #"SRR3884678"   # Fruits_peduncle (PRJNA328564)
    # Other organs
    "SRR3884685"   # Radicles (PRJNA328564)
    "SRR3884677"   # Cotyledons (PRJNA328564)
)

log_step "Starting Method 5: Bowtie2/RSEM Quantification Pipeline"

# Run Bowtie2/RSEM pipeline
bowtie2_rsem_pipeline \
    --FASTA "$FASTA_FILE" \
    --RNASEQ_LIST "${SRR_LIST[@]}"

log_info "Method 5 pipeline complete"
