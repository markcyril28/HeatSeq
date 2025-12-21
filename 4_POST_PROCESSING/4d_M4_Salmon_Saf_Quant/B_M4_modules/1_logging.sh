#!/bin/bash

#===============================================================================
# METHOD 4: SALMON SAF QUANTIFICATION - LOGGING CONFIGURATION
#===============================================================================
# Purpose: Configures logging system for Method 4
# Features:
#   - Sources centralized logging utilities
#   - Overrides log directories to use local Method 4 logs
#   - Maintains separate log categories (full, time, space, combined)
# Dependencies: Requires ../../modules/logging_utils.sh
#===============================================================================

# Determine base directory (three levels up from this module)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Source centralized logging utilities from main pipeline
source "$BASE_DIR/modules/logging_utils.sh"

# Generate unique run identifier if not already set
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"

# Configure Method 4-specific log directories (override global paths)
METHOD_4_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export LOG_DIR="$METHOD_4_DIR/logs/log_files"             # Main execution logs
export TIME_DIR="$METHOD_4_DIR/logs/time_logs"            # Performance timing
export SPACE_DIR="$METHOD_4_DIR/logs/space_logs"          # Disk usage tracking
export SPACE_TIME_DIR="$METHOD_4_DIR/logs/space_time_logs" # Combined metrics

# Define log file paths with unique run identifiers
export LOG_FILE="$LOG_DIR/Method_4_${RUN_ID}_full_log.log"
export TIME_FILE="$TIME_DIR/Method_4_${RUN_ID}_time_metrics.csv"
export TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
export SPACE_FILE="$SPACE_DIR/Method_4_${RUN_ID}_space_metrics.csv"
export SPACE_TIME_FILE="$SPACE_TIME_DIR/Method_4_${RUN_ID}_combined_metrics.csv"
