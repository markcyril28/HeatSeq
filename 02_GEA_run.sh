#!/bin/bash 

# Set variables
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_script_log.log}"

mkdir -p "$LOG_DIR"
# Forces removal and ignores errors if files do not exist.
rm -rf "$LOG_DIR"/*.log

chmod +x ./02_GEA_script_v9.sh

# Run the main pipeline script and log resource usage
/usr/bin/time -v ./02_GEA_script_v9.sh >> "$LOG_FILE" 2>&1
