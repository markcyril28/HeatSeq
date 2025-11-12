#!/bin/bash 

set -euo pipefail

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate GEA_ENV

# Source logging utilities
source "modules/logging_utils.sh"

# Setup logging system
setup_logging

log_info "Preparing pipeline execution environment"

# Chmod and Converts the file to Unix line endings
chmod +x ./*.sh
dos2unix ./*.sh
find . -type f \( -name "*.fa" -o -name "*.fasta"\) -exec dos2unix {} + 2>/dev/null || true
find . -type f \( -name "*.txt" -o -name "*.sh" -o -name "*.py" -o -name "*.R" \) -exec dos2unix {} + 2>/dev/null || true
find . -type f \( -name "*.sh" -o -name "*.py" -o -name "*.pl" \) -exec chmod +x {} + 2>/dev/null || true


# Run the main pipeline script and log resource usage
log_step "Executing main GEA pipeline script"
run_with_space_time_log ./b_GEA_script_v11.sh

log_step "Creating compressed archive for folders: 4_POST_PROCESSING, 4_OUTPUTS, and logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
tar -czf "CMSC244_${TIMESTAMP}.tar.gz" 4_POST_PROCESSING 4_OUTPUTS logs
log_info "Archive created: CMSC244_${TIMESTAMP}.tar.gz"

log_info "Pipeline execution completed successfully"
