#!/bin/bash 

# Set variables
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_script_log.log}"

mkdir -p "$LOG_DIR"
# Forces removal and ignores errors if files do not exist.
rm -rf "$LOG_DIR"/*

# Chmod and Converts the file to Unix line endings
chmod +x ./*.sh
dos2unix ./*.sh

# Run the main pipeline script and log resource usage
/usr/bin/time -v ./b_GEA_script_v9_HPC.sh >> "$LOG_FILE" 2>&1

/usr/bin/time -v ./4b_Method_2_HISAT2_De_Novo/a_stringtie_method_2.2_HPC.sh >> "$LOG_FILE" 2>&1

#/usr/bin/time -v ./4b_Method_2_HISAT2_De_Novo/b_heatmap_METHOD2_RESULTS_matrices_HPC.R >> "$LOG_FILE" 2>&1
