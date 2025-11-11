#!/bin/bash

# ===============================================
# LOGGING SYSTEM FOR METHOD 3
# Sources centralized logging_utils.sh
# ===============================================

# Get the base directory (two levels up from this module)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Source centralized logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories for Method 3 local logs
METHOD_3_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="${LOG_DIR:-$METHOD_3_DIR/logs/log_files}"
TIME_DIR="${TIME_DIR:-$METHOD_3_DIR/logs/time_files}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/Method_3_${RUN_ID}_full_log.log}"
TIME_FILE="${TIME_FILE:-$TIME_DIR/Method_3_${RUN_ID}_time_metrics.csv}"
TIME_TEMP="${TIME_TEMP:-$TIME_DIR/.time_temp_${RUN_ID}.txt}"

# Legacy function aliases for backward compatibility
log_header() { log_step "$@"; }
log() { log_info "$@"; }
log_error_and_exit() { log_error "$@"; exit 1; }