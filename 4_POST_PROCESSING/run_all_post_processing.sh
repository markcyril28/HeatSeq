#!/bin/bash

#===============================================================================
# MASTER POST-PROCESSING SCRIPT - RUN ALL HEATMAP GENERATION
#===============================================================================
# Runs heatmap generation for all analysis methods
# Assumes alignment/quantification steps are already complete
#
# NOTE: This script coordinates running multiple method post-processing scripts.
#       To configure gene groups and analysis options for each method individually,
#       edit the respective 1_run_Heatmap_Wrapper.sh in each method folder:
#         - 4a_Method_1_HISAT2_Ref_Guided/1_run_Heatmap_Wrapper.sh
#         - 4c_Method_3_Trinity_De_Novo/1_run_Heatmap_Wrapper.sh
#         - 4d_Method_4_Salmon_Saf_Quantification/1_run_Heatmap_Wrapper.sh
#         - 4e_Method_5_Bowtie2_Quantification/1_run_Heatmap_Wrapper.sh
#
#       Each wrapper script can be run independently within its own folder.
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

#==================================================================================
# CONFIGURATION - TOGGLE METHODS TO RUN
#==================================================================================

# Methods to process (comment out methods you want to skip)
METHODS_TO_RUN=(
    "METHOD_1"  # HISAT2 Reference-Guided
    #"METHOD_2"  # HISAT2 De Novo (archived - do not use)
    #"METHOD_3"  # Trinity De Novo
    "METHOD_4"  # Salmon SAF Quantification
    "METHOD_5"  # Bowtie2/RSEM Quantification
)

# NOTE: Gene groups, analysis toggles, and output controls are configured
#       in each method's individual 1_run_Heatmap_Wrapper.sh script

#==================================================================================

# Initialize conda for bash
eval "$(conda shell.bash hook)"

conda activate GEA_ENV

source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
LOG_FILE="$LOG_DIR/run_all_post_processing_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/run_all_post_processing_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/run_all_post_processing_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/run_all_post_processing_${RUN_ID}_combined_metrics.csv"

setup_logging

log_step "Starting Post-Processing for All Methods"
log_info "Methods to run: ${METHODS_TO_RUN[*]}"

# Helper function to check if method should run
should_run_method() {
    local method=$1
    for m in "${METHODS_TO_RUN[@]}"; do
        if [ "$m" = "$method" ]; then
            return 0
        fi
    done
    return 1
}

if should_run_method "METHOD_1"; then
    log_step "RUNNING METHOD 1 POST-PROCESSING"
    pushd 4a_Method_1_HISAT2_Ref_Guided > /dev/null
    run_with_space_time_log bash 1_run_Heatmap_Wrapper.sh
    popd > /dev/null
    log_info "METHOD 1 POST-PROCESSING COMPLETE"
else
    log_info "METHOD 1 POST-PROCESSING SKIPPED"
fi

if should_run_method "METHOD_3"; then
    log_step "RUNNING METHOD 3 POST-PROCESSING"
    pushd 4c_Method_3_Trinity_De_Novo > /dev/null
    run_with_space_time_log bash 1_run_Heatmap_Wrapper.sh
    popd > /dev/null
    log_info "METHOD 3 POST-PROCESSING COMPLETE"
else
    log_info "METHOD 3 POST-PROCESSING SKIPPED"
fi

if should_run_method "METHOD_4"; then
    log_step "RUNNING METHOD 4 POST-PROCESSING"
    pushd 4d_Method_4_Salmon_Saf_Quantification > /dev/null
    run_with_space_time_log bash 1_run_Heatmap_Wrapper.sh
    popd > /dev/null
    log_info "METHOD 4 POST-PROCESSING COMPLETE"
else
    log_info "METHOD 4 POST-PROCESSING SKIPPED"
fi

if should_run_method "METHOD_5"; then
    log_step "RUNNING METHOD 5 POST-PROCESSING"
    pushd 4e_Method_5_Bowtie2_Quantification > /dev/null
    run_with_space_time_log bash 1_run_Heatmap_Wrapper.sh
    popd > /dev/null
    log_info "METHOD 5 POST-PROCESSING COMPLETE"
else
    log_info "METHOD 5 POST-PROCESSING SKIPPED"
fi

log_step "All Post-Processing Complete"
log_info "Methods processed: ${METHODS_TO_RUN[*]}" 

