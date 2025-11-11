#!/bin/bash 

set -euo pipefail

# Source logging utilities
source "modules/logging_utils.sh"

# Setup logging system
setup_logging

log_info "Preparing pipeline execution environment"

# Chmod and Converts the file to Unix line endings
chmod +x ./*.sh
dos2unix ./*.sh

# Run the main pipeline script and log resource usage
log_step "Executing main GEA pipeline script"
run_with_space_time_log ./b_GEA_script_v12.sh

log_step "Creating compressed archive of outputs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
tar -czf "4_OUTPUTS_${TIMESTAMP}.tar.gz" 4_OUTPUTS
log_info "Archive created: 4_OUTPUTS_${TIMESTAMP}.tar.gz"