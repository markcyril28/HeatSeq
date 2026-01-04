#!/bin/bash
# ==============================================================================
# QUALITY CONTROL FUNCTIONS
# ==============================================================================
# FastQC and MultiQC quality control utilities
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${QC_SOURCED:-}" == "true" ]] && return 0
export QC_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_preproc.sh"

# ==============================================================================
# QUALITY CONTROL FUNCTIONS
# ==============================================================================

run_quality_control() {
	local SRR="$1"
	local RAW_DIR="$RAW_DIR_ROOT/$SRR"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	
	log_step "QC: $SRR"
	
	if command -v fastqc >/dev/null 2>&1; then
		# QC for raw files
		if [[ -d "$RAW_DIR" ]]; then
			local raw_outdir="$FASTQC_ROOT/${SRR}_raw"
			mkdir -p "$raw_outdir"
			if ! compgen -G "$raw_outdir/*_fastqc.html" >/dev/null; then
				log_info "Running FastQC on raw files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$raw_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || log_warn "FastQC failed for raw $SRR"
			fi
		fi
		
		# QC for trimmed files
		if [[ -d "$TrimGalore_DIR" ]]; then
			local trim_outdir="$FASTQC_ROOT/${SRR}_trimmed"
			mkdir -p "$trim_outdir"
			if ! compgen -G "$trim_outdir/*_fastqc.html" >/dev/null; then
				log_info "Running FastQC on trimmed files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-2}" -o "$trim_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || log_warn "FastQC failed for trimmed $SRR"
			fi
		fi
	else
		log_warn "FastQC not found. Skipping QC for $SRR."
	fi
	
	# Run MultiQC to aggregate results
	run_multiqc
}

run_multiqc() {
	if command -v multiqc >/dev/null 2>&1; then
		run_with_space_time_log multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || true
	else
		log_warn "MultiQC not found. Skipping aggregation."
	fi
}

run_quality_control_all() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for QC"; return 1; }
	
	log_step "Running Quality Control for ${#SRR_LIST[@]} samples"
	
	if should_use_parallel; then
		log_info "Running FastQC with GNU Parallel (JOBS=${JOBS:-2})"
		run_quality_control_parallel "${SRR_LIST[@]}"
	else
		log_info "Running FastQC sequentially (USE_GNU_PARALLEL=${USE_GNU_PARALLEL:-FALSE})"
		for SRR in "${SRR_LIST[@]}"; do
			run_quality_control "$SRR"
		done
	fi
	
	# Run MultiQC once at the end to aggregate all results
	run_multiqc
	
	log_info "Quality control completed for all samples."
}

# ==============================================================================
# PARALLEL QC FUNCTIONS
# ==============================================================================

run_quality_control_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && { log_error "No SRR IDs provided for parallel QC"; return 1; }
	
	# Export required variables and functions for parallel execution
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	export RAW_DIR_ROOT TRIM_DIR_ROOT FASTQC_ROOT THREADS_PER_JOB
	export -f log_info log_warn log_error log_step run_with_space_time_log 2>/dev/null || true
	
	_qc_worker() {
		local SRR="$1"
		
		# Activate conda environment in subshell
		if [[ -n "$CONDA_PREFIX" ]]; then
			source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
			conda activate "$CONDA_DEFAULT_ENV" 2>/dev/null || true
		fi
		
		local RAW_DIR="$RAW_DIR_ROOT/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		
		# QC for raw files
		if [[ -d "$RAW_DIR" ]]; then
			local raw_outdir="$FASTQC_ROOT/${SRR}_raw"
			mkdir -p "$raw_outdir"
			if ! compgen -G "$raw_outdir/*_fastqc.html" >/dev/null; then
				echo "[INFO] Running FastQC on raw files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$raw_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || echo "[WARN] FastQC failed for raw $SRR"
			fi
		fi
		
		# QC for trimmed files
		if [[ -d "$TrimGalore_DIR" ]]; then
			local trim_outdir="$FASTQC_ROOT/${SRR}_trimmed"
			mkdir -p "$trim_outdir"
			if ! compgen -G "$trim_outdir/*_fastqc.html" >/dev/null; then
				echo "[INFO] Running FastQC on trimmed files for $SRR"
				fastqc -t "${THREADS_PER_JOB:-2}" -o "$trim_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || echo "[WARN] FastQC failed for trimmed $SRR"
			fi
		fi
	}
	export -f _qc_worker
	
	# Run FastQC in parallel for all samples
	printf "%s\n" "${SRR_LIST[@]}" | parallel --env PATH --env CONDA_PREFIX -j "${JOBS:-2}" _qc_worker {}
}

# ==============================================================================
# QC SUMMARY FUNCTIONS
# ==============================================================================

generate_qc_summary() {
	local output_file="${1:-$FASTQC_ROOT/qc_summary.txt}"
	
	log_step "Generating QC Summary"
	
	{
		echo "=========================================="
		echo "FastQC Summary Report"
		echo "Generated: $(date)"
		echo "=========================================="
		echo ""
		
		# Find all fastqc_data.txt files and extract key metrics
		for summary_file in "$FASTQC_ROOT"/*/*_fastqc/summary.txt; do
			if [[ -f "$summary_file" ]]; then
				local sample_dir=$(dirname "$summary_file")
				local sample_name=$(basename "$sample_dir" | sed 's/_fastqc$//')
				echo "Sample: $sample_name"
				cat "$summary_file"
				echo ""
			fi
		done
	} > "$output_file"
	
	log_info "QC summary saved to: $output_file"
}
