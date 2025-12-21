#!/bin/bash
# ==============================================================================
# PIPELINE UTILITIES
# ==============================================================================
# Core pipeline functions for data download, trimming, and QC
# ==============================================================================

set -euo pipefail

# ==============================================================================
# TRIMMING CONFIGURATION
# ==============================================================================
HEADCROP_BASES=10
TAILCROP_BASES=0
DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING="${DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING:-FALSE}"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/logging_utils.sh"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/shared_utils.sh"

# Initialize directories
init_directories

# ==============================================================================
# TOOL INSTALLATION
# ==============================================================================
mamba_install() {
	log_info "Installing prerequisites via mamba..."
	command -v mamba >/dev/null 2>&1 || { log_error "Mamba not found."; return 1; }
	mamba install -c conda-forge -c bioconda \
		aria2 parallel-fastq-dump sra-tools \
		hisat2 stringtie samtools bowtie2 rsem salmon star trim-galore trimmomatic \
		fastqc multiqc parallel -y
	log_info "Prerequisites installation completed."
}

# ==============================================================================
# QUALITY CONTROL
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
				run_with_space_time_log fastqc -t "${THREADS:-1}" -o "$raw_outdir" \
					"$RAW_DIR"/${SRR}*.fastq* 2>/dev/null || log_warn "FastQC failed for raw $SRR"
			fi
		fi
		
		# QC for trimmed files
		if [[ -d "$TrimGalore_DIR" ]]; then
			local trim_outdir="$FASTQC_ROOT/${SRR}_trimmed"
			mkdir -p "$trim_outdir"
			if ! compgen -G "$trim_outdir/*_fastqc.html" >/dev/null; then
				log_info "Running FastQC on trimmed files for $SRR"
				run_with_space_time_log fastqc -t "${THREADS:-1}" -o "$trim_outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || log_warn "FastQC failed for trimmed $SRR"
			fi
		fi
	fi
	
	command -v multiqc >/dev/null 2>&1 && \
		run_with_space_time_log multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || true
}

gzip_trimmed_fastq_files() {
	log_info "Compressing trimmed FASTQ files in $TRIM_DIR_ROOT..."
	find "$TRIM_DIR_ROOT" -type f -name "*.fq" -print0 | \
		xargs -0 -P "${JOBS:-2}" -I {} gzip {} 2>/dev/null || true
	log_info "Compression completed."
}

# ==============================================================================
# DOWNLOAD FUNCTIONS
# ==============================================================================
download_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		[[ -n "$raw1" ]] && { log_info "Raw files for $SRR exist. Skipping."; continue; }
		
		log_info "Downloading $SRR..."
		run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir"
		run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
			"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir"
		
		[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
	done
	log_info "All downloads completed."
}

# ==============================================================================
# TRIMMING FUNCTIONS
# ==============================================================================
# Core trimming logic used by all trim functions
_trim_single_srr() {
	local SRR="$1"
	local raw_dir="$RAW_DIR_ROOT/$SRR"
	local trim_dir="$TRIM_DIR_ROOT/$SRR"
	local trimmed1="$trim_dir/${SRR}_1_val_1.fq"
	local trimmed2="$trim_dir/${SRR}_2_val_2.fq"
	
	mkdir -p "$trim_dir"
	
	# Skip if already trimmed
	find_trimmed_fastq "$SRR"
	[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; return 0; }
	
	# Find raw files
	find_raw_fastq "$SRR"
	[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; return 1; }
	
	log_info "Trimming $SRR with TrimGalore..."
	run_with_space_time_log trim_galore --cores "$THREADS" \
		--paired "$raw1" "$raw2" --output_dir "$trim_dir"
	
	# Trimmomatic HEADCROP
	log_info "Applying HEADCROP:${HEADCROP_BASES} for $SRR..."
	local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
	local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
	[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz"
	[[ -f "${tg_r2}.gz" ]] && gunzip "${tg_r2}.gz"
	
	local tmp_r1="$trim_dir/${SRR}_1_headcrop.fq"
	local tmp_r2="$trim_dir/${SRR}_2_headcrop.fq"
	run_with_space_time_log trimmomatic PE -threads "$THREADS" \
		"$tg_r1" "$tg_r2" "$tmp_r1" /dev/null "$tmp_r2" /dev/null \
		HEADCROP:${HEADCROP_BASES}
	mv "$tmp_r1" "$tg_r1"
	mv "$tmp_r2" "$tg_r2"
	
	verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
}

trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	local total=${#SRR_LIST[@]} current=0
	for SRR in "${SRR_LIST[@]}"; do
		((current++))
		log_info "[$current/$total] Processing $SRR..."
		_trim_single_srr "$SRR" || log_warn "Trimming failed for $SRR"
	done
	gzip_trimmed_fastq_files
	log_info "All trimming completed."
}

trim_srrs_trimmomatic() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		local out1="$trim_dir/${SRR}_1_val_1.fq"
		local out2="$trim_dir/${SRR}_2_val_2.fq"
		
		mkdir -p "$trim_dir"
		find_trimmed_fastq "$SRR"
		[[ -n "${trimmed1:-}" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; continue; }
		
		log_info "Trimming $SRR with Trimmomatic..."
		run_with_space_time_log trimmomatic PE -threads "$THREADS" \
			"$raw1" "$raw2" "$out1" /dev/null "$out2" /dev/null \
			ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
			HEADCROP:${HEADCROP_BASES} SLIDINGWINDOW:4:20 MINLEN:36
		
		# Trim last N bases using cutadapt (relative 3' trimming)
		log_info "Applying TAILCROP:${TAILCROP_BASES} for $SRR..."
		run_with_error_capture cutadapt -u -${TAILCROP_BASES} -U -${TAILCROP_BASES} -o "${out1}.tmp" -p "${out2}.tmp" "$out1" "$out2"
		mv "${out1}.tmp" "$out1"
		mv "${out2}.tmp" "$out2"
		
		verify_trimming_and_cleanup "$SRR" "$out1" "$out2" "$raw1" "$raw2"
	done
	gzip_trimmed_fastq_files
}

download_and_trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	for SRR in "${SRR_LIST[@]}"; do
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed files for $SRR exist. Skipping."; continue; }
		
		find_raw_fastq "$SRR"
		if [[ -z "$raw1" ]]; then
			log_info "Downloading $SRR..."
			run_with_space_time_log prefetch "$SRR" --output-directory "$raw_dir"
			run_with_space_time_log fasterq-dump --split-files --threads "$THREADS" \
				"$raw_dir/$SRR/$SRR.sra" -O "$raw_dir"
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
			find_raw_fastq "$SRR"
		fi
		
		[[ -z "$raw1" ]] && { log_warn "Raw FASTQ not found for $SRR"; continue; }
		_trim_single_srr "$SRR"
	done
	gzip_trimmed_fastq_files
}

# ==============================================================================
# PARALLEL PROCESSING
# ==============================================================================
download_and_trim_srrs_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS HEADCROP_BASES DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING LOG_FILE
	export -f timestamp log log_info log_warn log_error run_with_space_time_log
	export -f find_trimmed_fastq find_raw_fastq verify_trimming_and_cleanup
	
	_parallel_worker() {
		local SRR="$1"
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		# Check if done
		find_trimmed_fastq "$SRR"
		[[ -n "$trimmed1" ]] && { log_info "Trimmed $SRR exists. Skipping."; return 0; }
		
		# Download if needed
		find_raw_fastq "$SRR"
		if [[ -z "$raw1" ]]; then
			prefetch "$SRR" --output-directory "$raw_dir" || return 1
			fasterq-dump --split-files --threads 2 "$raw_dir/$SRR/$SRR.sra" -O "$raw_dir" || return 1
			[[ -f "$raw_dir/${SRR}_1.fastq" ]] && gzip "$raw_dir/${SRR}_1.fastq" "$raw_dir/${SRR}_2.fastq"
			find_raw_fastq "$SRR"
		fi
		
		[[ -z "$raw1" ]] && { log_warn "No raw for $SRR"; return 1; }
		
		# Trim
		trim_galore --cores 2 --paired "$raw1" "$raw2" --output_dir "$trim_dir"
		local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
		local tg_r2="$trim_dir/${SRR}_2_val_2.fq"
		[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz" "${tg_r2}.gz"
		
		trimmomatic PE -threads 2 "$tg_r1" "$tg_r2" \
			"${tg_r1}.tmp" /dev/null "${tg_r2}.tmp" /dev/null HEADCROP:${HEADCROP_BASES}
		mv "${tg_r1}.tmp" "$tg_r1"
		mv "${tg_r2}.tmp" "$tg_r2"
		
		verify_trimming_and_cleanup "$SRR" "$tg_r1" "$tg_r2" "$raw1" "$raw2"
	}
	export -f _parallel_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${JOBS:-2}" _parallel_worker {}
	gzip_trimmed_fastq_files
}

download_and_trim_srrs_wget_parallel() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS HEADCROP_BASES DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING
	export -f timestamp log log_info log_warn log_error
	
	_wget_worker() {
		local SRR="$1"
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		local trimmed1="$trim_dir/${SRR}_1_val_1.fq"
		[[ -f "$trimmed1" || -f "${trimmed1}.gz" ]] && return 0
		
		# Get ENA links
		local ena_links=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2)
		[[ -z "$ena_links" ]] && { log_warn "No ENA links for $SRR"; return 1; }
		
		wget -q -c -P "$raw_dir" "ftp://$(echo "$ena_links" | cut -f1 -d';')" "ftp://$(echo "$ena_links" | cut -f2 -d';')"
		
		trim_galore --cores 1 --paired "$raw_dir/${SRR}_1.fastq.gz" "$raw_dir/${SRR}_2.fastq.gz" --output_dir "$trim_dir"
		
		local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
		[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz" "$trim_dir/${SRR}_2_val_2.fq.gz"
		
		trimmomatic PE -threads 1 "$tg_r1" "$trim_dir/${SRR}_2_val_2.fq" \
			"${tg_r1}.tmp" /dev/null "$trim_dir/${SRR}_2.tmp" /dev/null HEADCROP:${HEADCROP_BASES}
		mv "${tg_r1}.tmp" "$tg_r1"
		mv "$trim_dir/${SRR}_2.tmp" "$trim_dir/${SRR}_2_val_2.fq"
		
		[[ "$DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING" == "TRUE" ]] && rm -f "$raw_dir"/*.fastq.gz
	}
	export -f _wget_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${JOBS:-2}" _wget_worker {}
	gzip_trimmed_fastq_files
}

download_kingfisher_and_trim_srrs() {
	local SRR_LIST=("$@")
	[[ ${#SRR_LIST[@]} -eq 0 ]] && SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	
	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS JOBS HEADCROP_BASES DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING
	export -f timestamp log log_info log_warn log_error run_with_space_time_log
	
	_kingfisher_worker() {
		local SRR="$1"
		local raw_dir="$RAW_DIR_ROOT/$SRR"
		local trim_dir="$TRIM_DIR_ROOT/$SRR"
		mkdir -p "$raw_dir" "$trim_dir"
		
		local trimmed1="$trim_dir/${SRR}_1_val_1.fq"
		[[ -f "$trimmed1" || -f "${trimmed1}.gz" ]] && return 0
		
		if [[ ! -f "$raw_dir/${SRR}_1.fastq.gz" ]]; then
			local kf_threads=$((THREADS / JOBS))
			[[ $kf_threads -lt 1 ]] && kf_threads=1
			kingfisher get --run-identifiers "$SRR" --output-directory "$raw_dir" \
				--download-threads "$kf_threads" --extraction-threads "$kf_threads" --gzip --check-md5sums || return 1
		fi
		
		local trim_cores=$((THREADS / JOBS))
		[[ $trim_cores -lt 1 ]] && trim_cores=1
		trim_galore --cores "$trim_cores" --paired "$raw_dir/${SRR}_1.fastq.gz" "$raw_dir/${SRR}_2.fastq.gz" --output_dir "$trim_dir"
		
		local tg_r1="$trim_dir/${SRR}_1_val_1.fq"
		[[ -f "${tg_r1}.gz" ]] && gunzip "${tg_r1}.gz" "$trim_dir/${SRR}_2_val_2.fq.gz"
		
		trimmomatic PE -threads "$trim_cores" "$tg_r1" "$trim_dir/${SRR}_2_val_2.fq" \
			"${tg_r1}.tmp" /dev/null "$trim_dir/${SRR}_2.tmp" /dev/null HEADCROP:${HEADCROP_BASES}
		mv "${tg_r1}.tmp" "$tg_r1"
		mv "$trim_dir/${SRR}_2.tmp" "$trim_dir/${SRR}_2_val_2.fq"
		
		[[ "$DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING" == "TRUE" ]] && rm -f "$raw_dir"/*.fastq.gz
	}
	export -f _kingfisher_worker
	
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${JOBS:-2}" _kingfisher_worker {}
	gzip_trimmed_fastq_files
}
