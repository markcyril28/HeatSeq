#!/bin/bash

# ==============================================================================
# MASTER SCRIPT TO RUN ALL HEATSEQ ANALYSIS METHODS
# ==============================================================================
# Description: This script serves as a master controller to run the analysis
#              pipelines for different methods. Each method can be toggled
#              on or off using the configuration variables below.
#
# Author: Mark Cyril R. Mercado
# Version: v2
# Date: November 2025
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Initialize conda for bash
eval "$(conda shell.bash hook)"

conda activate GEA_ENV

# Source logging utilities from parent directory
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories to use local logs folder
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
LOG_FILE="$LOG_DIR/run_all_methods_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/run_all_methods_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/run_all_methods_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/run_all_methods_${RUN_ID}_combined_metrics.csv"

# Setup logging
setup_logging

# Toggle which methods to run
RUN_METHOD_1_HISAT2_REF_GUIDED=true
RUN_METHOD_3_TRINITY_DE_NOVO=true
RUN_METHOD_4_SALMON_SAF=true
RUN_METHOD_5_BOWTIE2_RSEM=true

log_step "Starting All Methods"

if [ "$RUN_METHOD_1_HISAT2_REF_GUIDED" = true ]; then
    log_step "RUNNING METHOD 1: HISAT2 REFERENCE-GUIDED"
    pushd 4a_Method_1_HISAT2_Ref_Guided
    run_with_space_time_log bash 0_run_Method_1_Pipeline.sh
    popd
    log_info "METHOD 1 COMPLETE"
fi

if [ "$RUN_METHOD_3_TRINITY_DE_NOVO" = true ]; then
    log_step "RUNNING METHOD 3: TRINITY DE NOVO"
    pushd 4c_Method_3_Trinity_De_Novo
    run_with_space_time_log bash 0_run_Method_3_Pipeline.sh
    popd
    log_info "METHOD 3 COMPLETE"
fi

if [ "$RUN_METHOD_4_SALMON_SAF" = true ]; then
    log_step "RUNNING METHOD 4: SALMON SAF QUANTIFICATION"
    pushd 4d_Method_4_Salmon_Saf_Quantification
    run_with_space_time_log bash 0_run_Method_4_Pipeline.sh
    popd
    log_info "METHOD 4 COMPLETE"
fi

if [ "$RUN_METHOD_5_BOWTIE2_RSEM" = true ]; then
    log_step "RUNNING METHOD 5: BOWTIE2/RSEM QUANTIFICATION"
    pushd 4e_Method_5_Bowtie2_Quantification
    run_with_space_time_log bash 0_run_Method_5_Pipeline.sh
    popd
    log_info "METHOD 5 COMPLETE"
fi

log_step "All Methods Complete"