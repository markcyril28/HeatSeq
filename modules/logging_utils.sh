#!/bin/bash

# ==============================================================================
# LOGGING UTILITIES
# ==============================================================================
# Five-version logging system for comprehensive pipeline tracking:
# 1. Full logs (logs/log_files/*.log) - Complete execution output
# 2. Time logs (logs/time_logs/*.csv) - Time/CPU/memory metrics only
# 3. Space logs (logs/space_logs/*.csv) - File/directory size metrics only
# 4. Combined logs (logs/space_time_logs/*.csv) - Time + space metrics together
# 5. Error/Warning logs (logs/error_warn_logs/*.log) - Errors and warnings only
# 6. Software catalog (logs/software_catalogs/*.csv) - Software versions used
# ==============================================================================
# ERROR CAPTURE:
# - Monitors for: error, exception, fatal, failed, command not found, 
#   no such file, cannot find, not installed, permission denied, traceback
# - Captures from: stderr, stdout, time output, and exit codes
# - Use run_with_error_capture() for simple commands
# - Use run_with_space_time_log() for resource-intensive commands
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
ERROR_WARN_DIR="logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="logs/software_catalogs"
LOG_FILE="$LOG_DIR/pipeline_${RUN_ID}_full_log.log"
TIME_FILE="$TIME_DIR/pipeline_${RUN_ID}_time_metrics.csv"
TIME_TEMP="$TIME_DIR/.time_temp_${RUN_ID}.txt"
SPACE_FILE="$SPACE_DIR/pipeline_${RUN_ID}_space_metrics.csv"
SPACE_TIME_FILE="$SPACE_TIME_DIR/pipeline_${RUN_ID}_combined_metrics.csv"
ERROR_WARN_FILE="$ERROR_WARN_DIR/pipeline_${RUN_ID}_errors_warnings.log"
SOFTWARE_FILE="$SOFTWARE_CATALOG_DIR/software_catalog_${RUN_ID}.csv"

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; [[ -n "${ERROR_WARN_FILE:-}" ]] && printf '[%s] [WARN] %s\n' "$(timestamp)" "$*" >> "$ERROR_WARN_FILE"; }
log_error() { log ERROR "$@"; [[ -n "${ERROR_WARN_FILE:-}" ]] && printf '[%s] [ERROR] %s\n' "$(timestamp)" "$*" >> "$ERROR_WARN_FILE"; }
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
	mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" "$ERROR_WARN_DIR" "$SOFTWARE_CATALOG_DIR" || {
		echo "ERROR: Failed to create logging directories" >&2
		return 1
	}
	
	# Clear previous logs if requested
	if [[ "$clear_logs" == "true" ]]; then
		rm -f "$LOG_DIR"/*.log 2>/dev/null || true
		rm -f "$TIME_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_TIME_DIR"/*.csv 2>/dev/null || true
		rm -f "$ERROR_WARN_DIR"/*.log 2>/dev/null || true
		rm -f "$SOFTWARE_CATALOG_DIR"/*.csv 2>/dev/null || true
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
	
	# Initialize combined error/warning log file
	if [[ ! -f "$ERROR_WARN_FILE" ]]; then
		touch "$ERROR_WARN_FILE" || {
			echo "ERROR: Failed to create ERROR_WARN_FILE at $ERROR_WARN_FILE" >&2
			return 1
		}
	fi
	
	# Initialize CSV header for software catalog if file doesn't exist
	if [[ ! -f "$SOFTWARE_FILE" ]]; then
		echo "Software/Tool,Version" > "$SOFTWARE_FILE" || {
			echo "ERROR: Failed to create SOFTWARE_FILE at $SOFTWARE_FILE" >&2
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
	log_info "Errors & Warnings to: $ERROR_WARN_FILE"
	log_info "Software catalog to: $SOFTWARE_FILE"
}

# Error handling and cleanup traps
trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
trap 'log_info "Script finished. See log: $LOG_FILE"; log_info "Time metrics: $TIME_FILE"; log_info "Errors & Warnings: $ERROR_WARN_FILE"' EXIT

capture_stderr_errors() {
	# Monitor stderr/stdout stream and capture errors to error log
	# Usage: command 2>&1 | capture_stderr_errors
	while IFS= read -r line; do
		echo "$line"
		if echo "$line" | grep -qiE 'error|exception|fatal|failed|command not found|no such file|cannot find|not installed|permission denied|traceback'; then
			printf '[%s] [ERROR] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
		if echo "$line" | grep -qiE 'warning|warn'; then
			printf '[%s] [WARN] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
	done
}

run_with_error_capture() {
	# Simple wrapper to run commands with error capture (without time logging)
	# Usage: run_with_error_capture COMMAND...
	local cmd_string="$*"
	local exit_code=0
	
	# Run command and capture output, monitoring for errors
	"$@" 2>&1 | capture_stderr_errors || exit_code=$?
	
	# Log explicit failure
	if [[ $exit_code -ne 0 ]]; then
		log_error "Command failed (exit=$exit_code): $cmd_string"
	fi
	
	return $exit_code
}

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
	
	# Capture errors/exceptions to error log
	if [[ $exit_code -ne 0 ]] || grep -qiE 'exception|error|fatal|failed|traceback|command not found|no such file|cannot find|not installed' "$TIME_TEMP" 2>/dev/null; then
		{
			printf '[%s] [ERROR] Command failed (exit=%d): %s\n' "$(timestamp)" "$exit_code" "$cmd_short"
			grep -iE 'exception|error|fatal|failed|traceback|filenotfound|no such file|command not found|cannot find|not installed|permission denied|access denied' "$TIME_TEMP" 2>/dev/null || true
		} >> "$ERROR_WARN_FILE"
	fi
	
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

log_software_to_csv() {
	# Helper function to write software version to CSV
	local tool_name="$1"
	local version="$2"
	echo "\"${tool_name}\",\"${version}\"" >> "$SOFTWARE_FILE"
}

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
	local os_info="$(uname -s -r -m 2>/dev/null || echo 'N/A')"
	echo "OS: $os_info"
	log_software_to_csv "OS" "$os_info"
	echo "Hostname: $(hostname 2>/dev/null || echo 'N/A')"
	echo "User: $(whoami 2>/dev/null || echo 'N/A')"
	echo ""

	# Conda environment
	echo "--- CONDA ENVIRONMENT ---"
	if command -v conda &>/dev/null; then
		local conda_ver="$(conda --version 2>/dev/null || echo 'N/A')"
		echo "Conda: $conda_ver"
		log_software_to_csv "Conda" "$conda_ver"
		echo "Active env: ${CONDA_DEFAULT_ENV:-N/A}"
	else
		echo "Conda: Not available"
		log_software_to_csv "Conda" "Not available"
	fi
	echo ""

	# Core bioinformatics tools
	echo "--- ALIGNMENT & ASSEMBLY ---"
	local ver
	if command -v hisat2 &>/dev/null; then ver="$(hisat2 --version 2>&1 | head -1)"; echo "HISAT2: $ver"; log_software_to_csv "HISAT2" "$ver"; else echo "HISAT2: Not installed"; log_software_to_csv "HISAT2" "Not installed"; fi
	if command -v bowtie2 &>/dev/null; then ver="$(bowtie2 --version 2>&1 | head -1)"; echo "Bowtie2: $ver"; log_software_to_csv "Bowtie2" "$ver"; else echo "Bowtie2: Not installed"; log_software_to_csv "Bowtie2" "Not installed"; fi
	if command -v salmon &>/dev/null; then ver="$(salmon --version 2>&1 | head -1)"; echo "Salmon: $ver"; log_software_to_csv "Salmon" "$ver"; else echo "Salmon: Not installed"; log_software_to_csv "Salmon" "Not installed"; fi
	if command -v STAR &>/dev/null; then ver="$(STAR --version 2>&1 | head -1)"; echo "STAR: $ver"; log_software_to_csv "STAR" "$ver"; else echo "STAR: Not installed"; log_software_to_csv "STAR" "Not installed"; fi
	echo ""

	# Quantification tools
	echo "--- QUANTIFICATION ---"
	if command -v stringtie &>/dev/null; then ver="$(stringtie --version 2>&1)"; echo "StringTie: $ver"; log_software_to_csv "StringTie" "$ver"; else echo "StringTie: Not installed"; log_software_to_csv "StringTie" "Not installed"; fi
	if command -v rsem-calculate-expression &>/dev/null; then ver="$(rsem-calculate-expression --version 2>&1 | head -1)"; echo "RSEM: $ver"; log_software_to_csv "RSEM" "$ver"; else echo "RSEM: Not installed"; log_software_to_csv "RSEM" "Not installed"; fi
	echo ""

	# Trimming tools
	echo "--- READ PROCESSING ---"
	if command -v trim_galore &>/dev/null; then ver="$(trim_galore --version 2>&1 | grep -i version | head -1)"; echo "Trim Galore: $ver"; log_software_to_csv "Trim Galore" "$ver"; else echo "Trim Galore: Not installed"; log_software_to_csv "Trim Galore" "Not installed"; fi
	if command -v trimmomatic &>/dev/null; then ver="$(trimmomatic -version 2>&1 | head -1)"; echo "Trimmomatic: $ver"; log_software_to_csv "Trimmomatic" "$ver"; else echo "Trimmomatic: Not installed"; log_software_to_csv "Trimmomatic" "Not installed"; fi
	if command -v cutadapt &>/dev/null; then ver="$(cutadapt --version 2>&1)"; echo "Cutadapt: $ver"; log_software_to_csv "Cutadapt" "$ver"; else echo "Cutadapt: Not installed"; log_software_to_csv "Cutadapt" "Not installed"; fi
	echo ""

	# QC tools
	echo "--- QUALITY CONTROL ---"
	if command -v fastqc &>/dev/null; then ver="$(fastqc --version 2>&1)"; echo "FastQC: $ver"; log_software_to_csv "FastQC" "$ver"; else echo "FastQC: Not installed"; log_software_to_csv "FastQC" "Not installed"; fi
	if command -v multiqc &>/dev/null; then ver="$(multiqc --version 2>&1)"; echo "MultiQC: $ver"; log_software_to_csv "MultiQC" "$ver"; else echo "MultiQC: Not installed"; log_software_to_csv "MultiQC" "Not installed"; fi
	echo ""

	# SAM/BAM tools
	echo "--- SAM/BAM UTILITIES ---"
	if command -v samtools &>/dev/null; then ver="$(samtools --version 2>&1 | head -1)"; echo "SAMtools: $ver"; log_software_to_csv "SAMtools" "$ver"; else echo "SAMtools: Not installed"; log_software_to_csv "SAMtools" "Not installed"; fi
	echo ""

	# SRA tools
	echo "--- SRA TOOLS ---"
	if command -v prefetch &>/dev/null; then ver="$(prefetch --version 2>&1 | head -2 | tail -1)"; echo "SRA-tools (prefetch): $ver"; log_software_to_csv "SRA-tools" "$ver"; else echo "SRA-tools: Not installed"; log_software_to_csv "SRA-tools" "Not installed"; fi
	if command -v fasterq-dump &>/dev/null; then echo "fasterq-dump: available"; log_software_to_csv "fasterq-dump" "available"; else echo "fasterq-dump: Not installed"; log_software_to_csv "fasterq-dump" "Not installed"; fi
	echo ""

	# Utilities
	echo "--- UTILITIES ---"
	if command -v parallel &>/dev/null; then ver="$(parallel --version 2>&1 | head -1)"; echo "GNU Parallel: $ver"; log_software_to_csv "GNU Parallel" "$ver"; else echo "GNU Parallel: Not installed"; log_software_to_csv "GNU Parallel" "Not installed"; fi
	if command -v pigz &>/dev/null; then ver="$(pigz --version 2>&1)"; echo "pigz: $ver"; log_software_to_csv "pigz" "$ver"; else echo "pigz: Not installed"; log_software_to_csv "pigz" "Not installed"; fi
	if command -v gzip &>/dev/null; then ver="$(gzip --version 2>&1 | head -1)"; echo "gzip: $ver"; log_software_to_csv "gzip" "$ver"; else echo "gzip: Not installed"; log_software_to_csv "gzip" "Not installed"; fi
	echo ""

	# R if available
	echo "--- R ENVIRONMENT ---"
	if command -v R &>/dev/null; then
		ver="$(R --version 2>&1 | head -1)"
		echo "R: $ver"
		log_software_to_csv "R" "$ver"
	else
		echo "R: Not installed"
		log_software_to_csv "R" "Not installed"
	fi
	echo ""

	# Python if available
	echo "--- PYTHON ENVIRONMENT ---"
	if command -v python &>/dev/null; then
		ver="$(python --version 2>&1)"
		echo "Python: $ver"
		log_software_to_csv "Python" "$ver"
	elif command -v python3 &>/dev/null; then
		ver="$(python3 --version 2>&1)"
		echo "Python: $ver"
		log_software_to_csv "Python" "$ver"
	else
		echo "Python: Not installed"
		log_software_to_csv "Python" "Not installed"
	fi
	echo ""

	echo "================================================================================"
	echo "END OF SOFTWARE CATALOG"
	echo "================================================================================"
	log_info "Software catalog CSV saved to: $SOFTWARE_FILE"
	
	# Re-enable strict error handling
	set -e
}
