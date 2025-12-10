#!/bin/bash

# ==============================================================================
# LOGGING UTILITIES
# ==============================================================================
# Four-version logging system for comprehensive pipeline tracking:
# 1. Full logs (logs/log_files/*.log) - Complete execution output
# 2. Time logs (logs/time_logs/*.csv) - Time/CPU/memory metrics only
# 3. Space logs (logs/space_logs/*.csv) - File/directory size metrics only
# 4. Combined logs (logs/space_time_logs/*.csv) - Time + space metrics together
# ==============================================================================

set -euo pipefail

# ==============================================================================
# LOGGING CONFIGURATION
# ==============================================================================

# Force correct paths regardless of environment variables
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="logs/log_files"
TIME_DIR="logs/time_logs"
SPACE_DIR="logs/space_logs"
SPACE_TIME_DIR="logs/space_time_logs"
LOG_FILE="$LOG_DIR/pipeline_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/pipeline_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/pipeline_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/pipeline_${RUN_ID}_combined_metrics.csv"

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; }
log_error() { log ERROR "$@"; }
log_step() { log INFO "=============== $* ==============="; }

setup_logging() {
	# Set up logging and output redirection with dual-format support
	# Usage: setup_logging [clear_logs_flag]
	# clear_logs_flag: "true" to clear existing logs, anything else to keep them
	local clear_logs="${1:-false}"
	
	# Skip if already initialized
	if [[ "${LOGGING_INITIALIZED:-}" == "true" ]]; then
		log_info "Logging already initialized, skipping setup"
		return 0
	fi
	
	keep_bam_global="${keep_bam_global:-n}"

	# Create all log directories - ensure they exist before any file operations
	mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" || {
		echo "ERROR: Failed to create logging directories" >&2
		return 1
	}
	
	# Clear previous logs if requested
	if [[ "$clear_logs" == "true" ]]; then
		rm -f "$LOG_DIR"/*.log 2>/dev/null || true
		rm -f "$TIME_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_TIME_DIR"/*.csv 2>/dev/null || true
		echo "Previous logs cleared"
	fi
	
	# Initialize CSV header for time metrics if file doesn't exist
	if [[ ! -f "$TIME_FILE" ]]; then
		echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Exit_Status" > "$TIME_FILE" || {
			echo "ERROR: Failed to create TIME_FILE at $TIME_FILE" >&2
			return 1
		}
	fi
	
	# Initialize CSV header for space metrics if file doesn't exist
	if [[ ! -f "$SPACE_FILE" ]]; then
		echo "Timestamp,Type,Path,Size_KB,Size_MB,Size_GB,File_Count,Description" > "$SPACE_FILE" || {
			echo "ERROR: Failed to create SPACE_FILE at $SPACE_FILE" >&2
			return 1
		}
	fi
	
	# Initialize CSV header for combined space+time metrics if file doesn't exist
	if [[ ! -f "$SPACE_TIME_FILE" ]]; then
		echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Input_Size_MB,Output_Size_MB,Exit_Status" > "$SPACE_TIME_FILE" || {
			echo "ERROR: Failed to create SPACE_TIME_FILE at $SPACE_TIME_FILE" >&2
			return 1
		}
	fi
	
	log_choice="${log_choice:-1}"
	if [[ "$log_choice" == "2" ]]; then
		exec >"$LOG_FILE" 2>&1
	else
		exec > >(tee -a "$LOG_FILE") 2>&1
	fi
	
	export LOGGING_INITIALIZED="true"
	log_info "Logging to: $LOG_FILE"
	log_info "Time metrics to: $TIME_FILE"
	log_info "Space metrics to: $SPACE_FILE"
	log_info "Combined metrics to: $SPACE_TIME_FILE"
}

# Error handling and cleanup traps
trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
trap 'log_info "Script finished. See log: $LOG_FILE"; log_info "Time metrics: $TIME_FILE"' EXIT

run_with_space_time_log() {
	# Run a command and log resource usage (tracks time and memory)
	# Logs full verbose output to log file and extracts metrics to CSV
	# Optionally logs space metrics if input/output paths provided
	# Usage: run_with_space_time_log [--input PATH] [--output PATH] COMMAND...
	
	local input_path=""
	local output_path=""
	
	# Parse optional space tracking arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--input)
				input_path="$2"
				shift 2
				;;
			--output)
				output_path="$2"
				shift 2
				;;
			*)
				break
				;;
		esac
	done
	
	local cmd_string="$*"
	local cmd_short="$cmd_string"  # Full command without truncation
	local start_ts="$(timestamp)"
	
	# Measure input size before running command
	local input_size_mb="0"
	if [[ -n "$input_path" && -e "$input_path" ]]; then
		local input_kb=$(du -sk "$input_path" 2>/dev/null | awk '{print $1}')
		input_size_mb=$(echo "scale=2; $input_kb / 1024" | bc)
	fi
	
	# Ensure TIME_DIR exists before running command
	mkdir -p "$TIME_DIR" || {
		log_error "Failed to create TIME_DIR: $TIME_DIR"
		return 1
	}
	
	# Run command with time and capture output to temp file
	local exit_code=0
	/usr/bin/time -v "$@" >> "$LOG_FILE" 2> "$TIME_TEMP" || exit_code=$?
	
	# Also append full time output to log file
	cat "$TIME_TEMP" >> "$LOG_FILE" 2>&1
	
	# Extract key metrics from time output
	local elapsed_time=$(grep "Elapsed (wall clock)" "$TIME_TEMP" | awk '{print $NF}' | awk -F: '{if (NF==3) print ($1*3600)+($2*60)+$3; else if (NF==2) print ($1*60)+$2; else print $1}')
	local cpu_percent=$(grep "Percent of CPU" "$TIME_TEMP" | awk '{print $NF}' | tr -d '%')
	local max_rss=$(grep "Maximum resident set size" "$TIME_TEMP" | awk '{print $NF}')
	local user_time=$(grep "User time" "$TIME_TEMP" | awk '{print $NF}')
	local system_time=$(grep "System time" "$TIME_TEMP" | awk '{print $NF}')
	
	# Measure output size after running command
	local output_size_mb="0"
	if [[ -n "$output_path" && -e "$output_path" ]]; then
		local output_kb=$(du -sk "$output_path" 2>/dev/null | awk '{print $1}')
		output_size_mb=$(echo "scale=2; $output_kb / 1024" | bc)
	fi
	
	# Append to time-only CSV
	echo "${start_ts},\"${cmd_short}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${exit_code}" >> "$TIME_FILE"
	
	# Append to combined space+time CSV
	echo "${start_ts},\"${cmd_short}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${input_size_mb},${output_size_mb},${exit_code}" >> "$SPACE_TIME_FILE"
	
	# Clean up temp file
	rm -f "$TIME_TEMP"
	
	return $exit_code
}

# ==============================================================================
# SPACE LOGGING FUNCTIONS
# ==============================================================================

log_file_size() {
	# Log size of a single file or directory
	local file_path="$1"
	local description="${2:-}"
	local type="FILE"
	
	if [[ ! -e "$file_path" ]]; then
		log_warn "Path does not exist: $file_path"
		return 1
	fi
	
	if [[ -d "$file_path" ]]; then
		type="DIR"
	fi
	
	# Get size in KB using du
	local size_kb=$(du -sk "$file_path" 2>/dev/null | awk '{print $1}')
	local size_mb=$(echo "scale=2; $size_kb / 1024" | bc)
	local size_gb=$(echo "scale=2; $size_kb / 1048576" | bc)
	
	# Count files if directory
	local file_count="-"
	if [[ -d "$file_path" ]]; then
		file_count=$(find "$file_path" -type f 2>/dev/null | wc -l)
	fi
	
	local ts="$(timestamp)"
	echo "${ts},${type},\"${file_path}\",${size_kb},${size_mb},${size_gb},${file_count},\"${description}\"" >> "$SPACE_FILE"
	log_info "Space logged: $file_path = ${size_mb}MB"
}

log_input_output_size() {
	# Log sizes of input and output files/directories
	local input_path="$1"
	local output_path="$2"
	local step_description="${3:-}"
	
	log_info "Logging space for: $step_description"
	
	if [[ -e "$input_path" ]]; then
		log_file_size "$input_path" "${step_description} - INPUT"
	else
		log_warn "Input path not found: $input_path"
	fi
	
	if [[ -e "$output_path" ]]; then
		log_file_size "$output_path" "${step_description} - OUTPUT"
	else
		log_warn "Output path not found: $output_path"
	fi
}

log_disk_usage() {
	# Log overall disk usage for the workspace
	local workspace_path="${1:-.}"
	local description="${2:-Workspace disk usage}"
	
	log_file_size "$workspace_path" "$description"
}

# ==============================================================================
# SOFTWARE CATALOG FUNCTION
# ==============================================================================

log_software_catalog() {
	# Log all software versions used during pipeline execution
	# Temporarily disable strict error handling for version checks
	set +e
	
	log_step "SOFTWARE CATALOG"
	echo ""
	echo "================================================================================"
	echo "SOFTWARE VERSIONS AND ENVIRONMENT INFORMATION"
	echo "Generated: $(timestamp)"
	echo "================================================================================"
	echo ""

	# System info
	echo "--- SYSTEM ---"
	echo "OS: $(uname -s -r -m 2>/dev/null || echo 'N/A')"
	echo "Hostname: $(hostname 2>/dev/null || echo 'N/A')"
	echo "User: $(whoami 2>/dev/null || echo 'N/A')"
	echo ""

	# Conda environment
	echo "--- CONDA ENVIRONMENT ---"
	if command -v conda &>/dev/null; then
		echo "Conda: $(conda --version 2>/dev/null || echo 'N/A')"
		echo "Active env: ${CONDA_DEFAULT_ENV:-N/A}"
	else
		echo "Conda: Not available"
	fi
	echo ""

	# Core bioinformatics tools
	echo "--- ALIGNMENT & ASSEMBLY ---"
	if command -v hisat2 &>/dev/null; then echo "HISAT2: $(hisat2 --version 2>&1 | head -1)"; else echo "HISAT2: Not installed"; fi
	if command -v bowtie2 &>/dev/null; then echo "Bowtie2: $(bowtie2 --version 2>&1 | head -1)"; else echo "Bowtie2: Not installed"; fi
	if command -v salmon &>/dev/null; then echo "Salmon: $(salmon --version 2>&1 | head -1)"; else echo "Salmon: Not installed"; fi
	if command -v Trinity &>/dev/null; then echo "Trinity: $(Trinity --version 2>&1 | head -1)"; else echo "Trinity: Not installed"; fi
	echo ""

	# Quantification tools
	echo "--- QUANTIFICATION ---"
	if command -v stringtie &>/dev/null; then echo "StringTie: $(stringtie --version 2>&1)"; else echo "StringTie: Not installed"; fi
	if command -v rsem-calculate-expression &>/dev/null; then echo "RSEM: $(rsem-calculate-expression --version 2>&1 | head -1)"; else echo "RSEM: Not installed"; fi
	echo ""

	# Trimming tools
	echo "--- READ PROCESSING ---"
	if command -v trim_galore &>/dev/null; then echo "Trim Galore: $(trim_galore --version 2>&1 | grep -i version | head -1)"; else echo "Trim Galore: Not installed"; fi
	if command -v trimmomatic &>/dev/null; then echo "Trimmomatic: $(trimmomatic -version 2>&1 | head -1)"; else echo "Trimmomatic: Not installed"; fi
	if command -v cutadapt &>/dev/null; then echo "Cutadapt: $(cutadapt --version 2>&1)"; else echo "Cutadapt: Not installed"; fi
	echo ""

	# QC tools
	echo "--- QUALITY CONTROL ---"
	if command -v fastqc &>/dev/null; then echo "FastQC: $(fastqc --version 2>&1)"; else echo "FastQC: Not installed"; fi
	if command -v multiqc &>/dev/null; then echo "MultiQC: $(multiqc --version 2>&1)"; else echo "MultiQC: Not installed"; fi
	echo ""

	# SAM/BAM tools
	echo "--- SAM/BAM UTILITIES ---"
	if command -v samtools &>/dev/null; then echo "SAMtools: $(samtools --version 2>&1 | head -1)"; else echo "SAMtools: Not installed"; fi
	echo ""

	# SRA tools
	echo "--- SRA TOOLS ---"
	if command -v prefetch &>/dev/null; then echo "SRA-tools (prefetch): $(prefetch --version 2>&1 | head -2 | tail -1)"; else echo "SRA-tools: Not installed"; fi
	if command -v fasterq-dump &>/dev/null; then echo "fasterq-dump: available"; else echo "fasterq-dump: Not installed"; fi
	echo ""

	# Utilities
	echo "--- UTILITIES ---"
	if command -v parallel &>/dev/null; then echo "GNU Parallel: $(parallel --version 2>&1 | head -1)"; else echo "GNU Parallel: Not installed"; fi
	if command -v pigz &>/dev/null; then echo "pigz: $(pigz --version 2>&1)"; else echo "pigz: Not installed"; fi
	if command -v gzip &>/dev/null; then echo "gzip: $(gzip --version 2>&1 | head -1)"; else echo "gzip: Not installed"; fi
	echo ""

	# R if available
	echo "--- R ENVIRONMENT ---"
	if command -v R &>/dev/null; then
		echo "R: $(R --version 2>&1 | head -1)"
	else
		echo "R: Not installed"
	fi
	echo ""

	# Python if available
	echo "--- PYTHON ENVIRONMENT ---"
	if command -v python &>/dev/null; then
		echo "Python: $(python --version 2>&1)"
	elif command -v python3 &>/dev/null; then
		echo "Python: $(python3 --version 2>&1)"
	else
		echo "Python: Not installed"
	fi
	echo ""

	echo "================================================================================"
	echo "END OF SOFTWARE CATALOG"
	echo "================================================================================"
	
	# Re-enable strict error handling
	set -e
}
