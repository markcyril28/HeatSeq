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
    "METHOD_2"  # HISAT2 De Novo
    #"METHOD_3"  # Trinity De Novo
    "METHOD_4"  # Salmon SAF Quantification
    "METHOD_5"  # Bowtie2/RSEM Quantification
)

# NOTE: Gene groups, analysis toggles, and output controls are configured
#       in each method's individual 1_run_Heatmap_Wrapper.sh script

#==================================================================================
# PARALLEL PROCESSING CONFIGURATION
#==================================================================================

# Number of threads for multi-threaded operations
THREADS=32

# Number of parallel jobs (set to 1 to disable parallel processing)
JOBS=4

# Use GNU Parallel for method execution (true/false)
USE_PARALLEL=true

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
log_info "Threads: $THREADS | Jobs: $JOBS | Parallel: $USE_PARALLEL"

# Export variables for child processes
export THREADS
export JOBS

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

# Function to run a single method
run_method() {
    local method_num=$1
    local method_name=$2
    local method_dir=$3
    
    if should_run_method "METHOD_$method_num"; then
        log_step "RUNNING METHOD $method_num POST-PROCESSING ($method_name)"
        pushd "$method_dir" > /dev/null
        run_with_space_time_log bash 1_run_Heatmap_Wrapper.sh
        popd > /dev/null
        log_info "METHOD $method_num POST-PROCESSING COMPLETE"
    else
        log_info "METHOD $method_num POST-PROCESSING SKIPPED"
    fi
}

# Run methods based on parallel configuration
if [ "$USE_PARALLEL" = true ] && [ "$JOBS" -gt 1 ] && command -v parallel &> /dev/null; then
    log_info "Using GNU Parallel with $JOBS jobs"
    
    # Build list of methods to run
    METHODS_LIST=()
    should_run_method "METHOD_1" && METHODS_LIST+=("1:HISAT2_Ref_Guided:4a_Method_1_HISAT2_Ref_Guided")
    should_run_method "METHOD_3" && METHODS_LIST+=("3:Trinity_De_Novo:4c_Method_3_Trinity_De_Novo")
    should_run_method "METHOD_4" && METHODS_LIST+=("4:Salmon_Saf:4d_Method_4_Salmon_Saf_Quantification")
    should_run_method "METHOD_5" && METHODS_LIST+=("5:Bowtie2:4e_Method_5_Bowtie2_Quantification")
    
    # Export variables and functions for parallel
    export BASE_DIR
    export LOG_DIR TIME_DIR SPACE_DIR SPACE_TIME_DIR
    export LOG_FILE TIME_FILE TIME_TEMP SPACE_FILE SPACE_TIME_FILE
    export RUN_ID
    export -f run_method
    export -f should_run_method
    export -f log_step
    export -f log_info
    export -f log_error
    export -f run_with_space_time_log
    
    # Run in parallel
    printf '%s\n' "${METHODS_LIST[@]}" | parallel -j "$JOBS" --colsep ':' run_method {1} {2} {3}
else
    log_info "Running methods sequentially"
    
    run_method 1 "HISAT2_Ref_Guided" "4a_Method_1_HISAT2_Ref_Guided"
    run_method 3 "Trinity_De_Novo" "4c_Method_3_Trinity_De_Novo"
    run_method 4 "Salmon_Saf" "4d_Method_4_Salmon_Saf_Quantification"
    run_method 5 "Bowtie2" "4e_Method_5_Bowtie2_Quantification"
fi

log_step "All Post-Processing Complete"
log_info "Methods processed: ${METHODS_TO_RUN[*]}"