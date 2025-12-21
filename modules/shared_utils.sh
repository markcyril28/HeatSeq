#!/bin/bash
# ==============================================================================
# SHARED UTILITY FUNCTIONS
# ==============================================================================
# Common helper functions used across all pipeline scripts
# Sourced by: pipeline_utils.sh, methods/methods.sh
# ==============================================================================

set -euo pipefail

# Guard against double-sourcing
[[ "${SHARED_UTILS_SOURCED:-}" == "true" ]] && return 0
export SHARED_UTILS_SOURCED="true"

# ==============================================================================
# FASTQ FILE DETECTION
# ==============================================================================
# Find trimmed FASTQ files for a given SRR ID
# Sets: trimmed1, trimmed2 (empty if single-end)
find_trimmed_fastq() {
	local SRR="$1"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	trimmed1="" trimmed2=""
	
	# Paired-end patterns (ordered by priority)
	if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
		trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
	elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
		trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
	elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1; then
		local files1=("$TrimGalore_DIR"/${SRR}*val_1.*) files2=("$TrimGalore_DIR"/${SRR}*val_2.*)
		trimmed1="${files1[0]}"
		[[ -f "${files2[0]:-}" ]] && trimmed2="${files2[0]}"
	# Single-end patterns
	elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq"
	elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq.gz" ]]; then
		trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq.gz"
	elif compgen -G "$TrimGalore_DIR/${SRR}*trimmed.fq*" >/dev/null 2>&1; then
		local files=("$TrimGalore_DIR"/${SRR}*trimmed.fq*)
		trimmed1="${files[0]}"
	fi
}

# Find raw FASTQ files for a given SRR ID
# Sets: raw1, raw2 (empty if single-end or not found)
find_raw_fastq() {
	local SRR="$1"
	local raw_dir="$RAW_DIR_ROOT/$SRR"
	raw1="" raw2=""
	
	if [[ -f "$raw_dir/${SRR}_1.fastq" && -f "$raw_dir/${SRR}_2.fastq" ]]; then
		raw1="$raw_dir/${SRR}_1.fastq"
		raw2="$raw_dir/${SRR}_2.fastq"
	elif [[ -f "$raw_dir/${SRR}_1.fastq.gz" && -f "$raw_dir/${SRR}_2.fastq.gz" ]]; then
		raw1="$raw_dir/${SRR}_1.fastq.gz"
		raw2="$raw_dir/${SRR}_2.fastq.gz"
	fi
}

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================
# Validate count matrix quality
validate_count_matrix() {
	local matrix="$1"
	local matrix_type="${2:-gene}"
	local min_samples="${3:-2}"
	
	[[ ! -f "$matrix" ]] && { log_error "Matrix not found: $matrix"; return 1; }
	
	log_step "Validating $matrix_type matrix: $(basename "$matrix")"
	
	local delim="\t"
	[[ "$matrix" == *.csv ]] && delim=","
	
	local header=$(head -n1 "$matrix")
	local num_samples=$(echo "$header" | awk -F"$delim" '{print NF-1}')
	
	[[ $num_samples -lt $min_samples ]] && { log_error "Insufficient samples: $num_samples (need ≥$min_samples)"; return 1; }
	
	local total_rows=$(tail -n +2 "$matrix" | wc -l)
	local unique_ids=$(tail -n +2 "$matrix" | cut -d"${delim:0:1}" -f1 | sort -u | wc -l)
	
	[[ $total_rows -ne $unique_ids ]] && { log_error "Duplicate ${matrix_type} IDs detected!"; return 1; }
	
	log_info "[VALIDATION] Passed - Samples: $num_samples, ${matrix_type^}s: $unique_ids"
	return 0
}

# Detect read length from FASTQ file
detect_read_length() {
	local fastq="$1"
	local default_length="${2:-150}"
	
	[[ ! -f "$fastq" ]] && { echo "$default_length"; return 1; }
	
	local decompress_cmd="cat"
	case "$fastq" in
		*.gz) decompress_cmd="zcat" ;;
		*.bz2) decompress_cmd="bzcat" ;;
	esac
	
	local avg_length=$($decompress_cmd "$fastq" 2>/dev/null | \
		awk 'NR%4==2 {sum+=length($0); count++} count==1000 {print int(sum/count); exit}')
	
	if [[ -z "$avg_length" || $avg_length -lt 50 || $avg_length -gt 300 ]]; then
		echo "$default_length"
		return 1
	fi
	echo "$avg_length"
}

# ==============================================================================
# METADATA FUNCTIONS
# ==============================================================================
# Load sample metadata from file
load_sample_metadata() {
	local metadata_file="${1:-sample_conditions.txt}"
	local -n metadata_array=$2
	
	[[ ! -f "$metadata_file" ]] && { log_warn "Metadata file not found: $metadata_file"; return 1; }
	
	local line_count=$(wc -l < "$metadata_file")
	[[ $line_count -lt 2 ]] && { log_error "Metadata file must contain at least 2 samples"; return 1; }
	
	while IFS=$'\t' read -r srr condition batch; do
		[[ -z "$srr" || "$srr" =~ ^# ]] && continue
		metadata_array["${srr}_condition"]="$condition"
		metadata_array["${srr}_batch"]="$batch"
	done < "$metadata_file"
	
	log_info "Loaded metadata: $(( ${#metadata_array[@]} / 2 )) samples"
	return 0
}

# Create sample metadata CSV for DESeq2
create_sample_metadata() {
	local metadata_file="$1"
	local -a sample_list=("${!2}")
	local metadata_source="${3:-sample_conditions.txt}"
	
	local delim=","
	[[ "$metadata_file" == *.tsv ]] && delim=$'\t'
	
	declare -A sample_metadata
	local has_external=false
	load_sample_metadata "$metadata_source" sample_metadata 2>/dev/null && has_external=true
	
	echo -e "sample${delim}condition${delim}batch" > "$metadata_file"
	for SRR in "${sample_list[@]}"; do
		if [[ "$has_external" == "true" ]]; then
			local condition="${sample_metadata[${SRR}_condition]:-unknown}"
			local batch="${sample_metadata[${SRR}_batch]:-1}"
		else
			local condition="condition_${SRR}"
			local batch="1"
		fi
		echo -e "$SRR${delim}$condition${delim}$batch" >> "$metadata_file"
	done
	
	[[ "$has_external" == "false" ]] && \
		log_warn "Sample conditions need manual specification in: $metadata_file"
}

# Generate tximport R script
generate_tximport_script() {
	local method="$1"
	local quant_dir="$2"
	local output_script="$3"
	local metadata_file="$4"
	local template="modules/tximport.R"
	
	[[ ! -f "$template" ]] && { log_error "Template not found: $template"; return 1; }
	
	cp "$template" "$output_script"
	sed -i "s|METHOD_PLACEHOLDER|$method|g" "$output_script"
	sed -i "s|QUANT_DIR_PLACEHOLDER|$quant_dir|g" "$output_script"
	sed -i "s|METADATA_FILE_PLACEHOLDER|$metadata_file|g" "$output_script"
	chmod +x "$output_script"
	log_info "[TXIMPORT] Generated: $output_script"
}

# ==============================================================================
# TRIMMING VERIFICATION
# ==============================================================================
# Verify trimming success and optionally cleanup raw files
verify_trimming_and_cleanup() {
	local SRR="$1"
	local trimmed1="$2"
	local trimmed2="$3"
	local raw1="${4:-}"
	local raw2="${5:-}"
	
	local success=false
	if { [[ -f "$trimmed1" && -s "$trimmed1" ]] || [[ -f "${trimmed1}.gz" && -s "${trimmed1}.gz" ]]; } && \
	   { [[ -z "$trimmed2" ]] || [[ -f "$trimmed2" && -s "$trimmed2" ]] || [[ -f "${trimmed2}.gz" && -s "${trimmed2}.gz" ]]; }; then
		success=true
		log_info "Trimming completed for $SRR"
		
		if [[ "$DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING" == "TRUE" && -n "$raw1" ]]; then
			log_info "Cleaning up raw files for $SRR..."
			rm -f "$raw1" "$raw2"
			rm -rf "$RAW_DIR_ROOT/$SRR/$SRR"
			rmdir "$RAW_DIR_ROOT/$SRR" 2>/dev/null || true
		fi
	else
		log_warn "Trimming may have failed for $SRR - keeping raw files"
	fi
	
	[[ "$success" == "true" ]]
}
