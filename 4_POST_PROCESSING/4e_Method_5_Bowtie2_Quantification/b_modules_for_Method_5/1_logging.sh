#!/bin/bash

# ===============================================
# LOGGING SYSTEM FOR METHOD 5
# Sources centralized logging_utils.sh
# ===============================================

# Get the base directory (two levels up from this module)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Source centralized logging utilities
source "$BASE_DIR/modules/logging_utils.sh"

# Override log directories for Method 5 local logs
METHOD_5_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="${LOG_DIR:-$METHOD_5_DIR/logs/log_files}"
TIME_DIR="${TIME_DIR:-$METHOD_5_DIR/logs/time_logs}"
SPACE_DIR="${SPACE_DIR:-$METHOD_5_DIR/logs/space_logs}"
SPACE_TIME_DIR="${SPACE_TIME_DIR:-$METHOD_5_DIR/logs/space_time_logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/Method_5_${RUN_ID}_full_log.log}"
TIME_FILE="${TIME_FILE:-$TIME_DIR/Method_5_${RUN_ID}_time_metrics.csv}"
TIME_TEMP="${TIME_TEMP:-$TIME_DIR/.time_temp_${RUN_ID}.txt}"
SPACE_FILE="${SPACE_FILE:-$SPACE_DIR/Method_5_${RUN_ID}_space_metrics.csv}"
SPACE_TIME_FILE="${SPACE_TIME_FILE:-$SPACE_TIME_DIR/Method_5_${RUN_ID}_combined_metrics.csv}"
