#!/bin/bash

# ==============================================================================
# PIPELINE UTILITIES
# ==============================================================================
# Core pipeline functions for data processing and analysis
# ==============================================================================

set -euo pipefail

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# Raw and Processed Data Directories
RAW_DIR_ROOT="1_RAW_SRR"
TRIM_DIR_ROOT="2_TRIMMED_SRR"
FASTQC_ROOT="3_FastQC"

# Method-specific output directories
HISAT2_REF_GUIDED_ROOT="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/5_stringtie_WD/a_Method_1_RAW_RESULTs"

# Method 2: HISAT2 De Novo
HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/5_stringtie_WD/a_Method_2_RAW_RESULTs"

# Method 3: Trinity De Novo
TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/5_stringtie_WD/a_Method_3_RAW_RESULTs"

# Method 4: Salmon SAF Quantification
SALMON_SAF_ROOT="4_POST_PROCESSING/4d_Method_4_Salmon_Saf_Quantification"

# Method 2: HISAT2 De Novo Assembly
HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/5_stringtie_WD/a_Method_2_RAW_RESULTs"

# Method 3: Trinity De Novo Assembly
TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/5_stringtie_WD/a_Method_3_RAW_RESULTs"

# Method 4: Salmon SAF Quantification
SALMON_SAF_ROOT="4_POST_PROCESSING/4d_Method_4_Salmon_Saf_Quantification"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/quant"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/matrices"

# Method 5: Bowtie2 RSEM Quantification
BOWTIE2_RSEM_ROOT="4_POST_PROCESSING/4e_Method_5_Bowtie2_Quantification"
RSEM_INDEX_ROOT="$BOWTIE2_RSEM_ROOT/index"
RSEM_QUANT_ROOT="$BOWTIE2_RSEM_ROOT/quant"
RSEM_MATRIX_ROOT="$BOWTIE2_RSEM_ROOT/matrices"

# Create required directories
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$TRINITY_DE_NOVO_ROOT" "$STRINGTIE_TRINITY_DE_NOVO_ROOT" \
	"$SALMON_SAF_ROOT" "$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
	"$BOWTIE2_RSEM_ROOT" "$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT"

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

show_pipeline_configuration() {
	# Display which pipelines are enabled for this run
	log_info "=== PIPELINE CONFIGURATION ==="
	log_info "Data Processing:"
	log_info "  Download & Trim SRR: $([ "$RUN_DOWNLOAD_and_TRIM_SRR" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Quality Control: $([ "$RUN_QUALITY_CONTROL" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	
	log_info "Analysis Methods:"
	log_info "  Method 1 - HISAT2 Ref-Guided: $([ "$RUN_METHOD_1_HISAT2_REF_GUIDED" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 2 - HISAT2 De Novo: $([ "$RUN_METHOD_2_HISAT2_DE_NOVO" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 3 - Trinity De Novo: $([ "$RUN_METHOD_3_TRINITY_DE_NOVO" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 4 - Salmon SAF: $([ "$RUN_METHOD_4_SALMON_SAF" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 5 - Bowtie2 + RSEM: $([ "$RUN_METHOD_5_BOWTIE2_RSEM" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	
	log_info "Validation & Comparison:"
	log_info "  Method Comparison: $([ "$RUN_METHOD_COMPARISON" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	
	# Count enabled methods
	local enabled_methods=0
	[[ "$RUN_METHOD_1_HISAT2_REF_GUIDED" = "TRUE" ]] && ((enabled_methods++))
	[[ "$RUN_METHOD_2_HISAT2_DE_NOVO" = "TRUE" ]] && ((enabled_methods++))
	[[ "$RUN_METHOD_3_TRINITY_DE_NOVO" = "TRUE" ]] && ((enabled_methods++))
	[[ "$RUN_METHOD_4_SALMON_SAF" = "TRUE" ]] && ((enabled_methods++))
	[[ "$RUN_METHOD_5_BOWTIE2_RSEM" = "TRUE" ]] && ((enabled_methods++))
	
	log_info "Total Methods Enabled: $enabled_methods"
	
	if [[ $enabled_methods -eq 0 ]]; then
		log_warn "WARNING: No analysis methods are enabled! Please enable at least one method."
	elif [[ $enabled_methods -gt 1 && "$RUN_METHOD_COMPARISON" = "TRUE" ]]; then
		log_info "EXCELLENT: Multiple methods enabled with comparison - excellent for validation!"
	fi
	
	log_info "============================="
}

mamba_install() {
	# Install required bioinformatics tools and dependencies
	log_info "Installing prerequisites via mamba..."

	# Check if mamba is available
	if ! command -v mamba >/dev/null 2>&1; then
		log_error "Mamba not found."
		return 1
	fi
	
	# Install packages
	mamba install -c conda-forge -c bioconda \
		aria2 parallel-fastq-dump sra-tools \
		hisat2 stringtie samtools bowtie2 rsem salmon trinity trim-galore \
		fastqc multiqc \
		parallel -y
	
	log_info "Prerequisites installation completed."
}

# ==============================================================================
# DATA DOWNLOAD AND PREPROCESSING FUNCTIONS
# ==============================================================================

gzip_trimmed_fastq_files() {
	# Compress all trimmed FASTQ files in the trimming directory using GNU parallel
	log_info "Compressing all trimmed FASTQ files in $TRIM_DIR_ROOT using GNU parallel..."
	find "$TRIM_DIR_ROOT" -type f -name "*.fq" -print0 | \
		parallel -0 -j "$JOBS" gzip {}
	log_info "Parallel compression of trimmed FASTQ files completed."
}

run_quality_control() {
	# Run quality control analysis on trimmed reads
	local SRR="$1"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	
	if [[ ! -d "$TrimGalore_DIR" ]]; then
		log_warn "Trimmed directory not found for $SRR. Skipping QC."
		return 1
	fi
	
	log_info "Running quality control for $SRR..."
	
	# Create QC output directory
	mkdir -p "$FASTQC_ROOT/$SRR"
	
	# FastQC on trimmed reads (skip if HTML report already exists)
	if command -v fastqc >/dev/null 2>&1; then
		outdir="$FASTQC_ROOT/${SRR}_trimmed"
		mkdir -p "$outdir"
		# If any FastQC HTML is present, assume QC was done and skip
		if compgen -G "$outdir/*_fastqc.html" >/dev/null; then
			log_info "FastQC HTML already exists for $SRR in $outdir. Skipping FastQC."
		else
			run_with_space_time_log \
				fastqc -t "${THREADS:-1}" -o "$outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || \
					log_warn "FastQC failed for $SRR"
		fi
	else
		log_warn "FastQC not available. Skipping read quality assessment."
	fi
	
	# MultiQC summary (if available)
	if command -v multiqc >/dev/null 2>&1; then
		run_with_space_time_log \
			multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || \
				log_warn "MultiQC failed. Individual FastQC reports still available."
	fi
}

download_and_trim_srrs() {
    # Download and trim RNA-seq data for each SRR sample
    local SRR_LIST=("$@")

    # Default to global project list if none provided
    if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
        SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
    fi

    for SRR in "${SRR_LIST[@]}"; do
        local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
        local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
        local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
        local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
        local raw1 raw2

        log_info "Working on $SRR..."
        mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

        # Check if trimmed files already exist
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Trimmed files for $SRR already exist. Skipping download and trimming."
            log_info "--------------------------------------------------"
            continue
        fi

        # Check for existing raw FASTQ files
        if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
            raw1="$raw_files_DIR/${SRR}_1.fastq"
            raw2="$raw_files_DIR/${SRR}_2.fastq"
            log_info "Raw FASTQ files for $SRR already exist. Skipping download."
        elif [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
            raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
            raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
            log_info "Compressed raw FASTQ files for $SRR already exist. Skipping download."
        else
            # Download with prefetch + fasterq-dump
            log_info "Downloading $SRR..."
            run_with_space_time_log prefetch "$SRR" --output-directory "$raw_files_DIR"
            
            # Important: include .sra file path explicitly
            run_with_space_time_log fasterq-dump --split-files --threads "${THREADS}" \
                "$raw_files_DIR/$SRR/$SRR.sra" -O "$raw_files_DIR"
			# Compress FASTQ files to save space
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
				log_info "Compressing FASTQ files for $SRR..."
				run_with_space_time_log gzip "$raw_files_DIR/${SRR}_1.fastq" "$raw_files_DIR/${SRR}_2.fastq"
			fi

            # Set raw file paths after download
            if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
                raw1="$raw_files_DIR/${SRR}_1.fastq"
                raw2="$raw_files_DIR/${SRR}_2.fastq"
            elif [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
                raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
                raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
            else
                log_warn "Raw FASTQ for $SRR not found after download; skipping."
                log_info "--------------------------------------------------"
                continue
            fi
        fi

        # Trim reads using Trim Galore
        log_info "Trimming $SRR..."
        run_with_space_time_log trim_galore --cores "${THREADS}" \
            --paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

        log_info "Done working on $SRR."
        log_info "--------------------------------------------------"

        # Verify trimming success and delete raw files if successful
        local trimming_success=false
        local trimmed_files_exist=false
        local trimmed_files_not_empty=false

        # Check if trimmed files exist (either uncompressed or compressed)
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            trimmed_files_exist=true
            
            # Check if files have content (not empty)
            if { [[ -s "$trimmed1" && -s "$trimmed2" ]] || [[ -s "${trimmed1}.gz" && -s "${trimmed2}.gz" ]]; }; then
                trimmed_files_not_empty=true
            fi
        fi

        # Verify overall trimming success
        if [[ "$trimmed_files_exist" == true && "$trimmed_files_not_empty" == true ]]; then
            trimming_success=true
            log_info "Trimming successfully completed for $SRR"
            
            # Delete raw files only if trimming was successful
            log_info "Deleting raw FASTQ files for $SRR after successful trimming..."
            rm -f "$raw1" "$raw2"
            log_info "Raw files deleted for $SRR"
        else
            log_warn "Trimming may have failed for $SRR - keeping raw files for potential reprocessing"
            if [[ "$trimmed_files_exist" == false ]]; then
                log_warn "Trimmed output files not found"
            elif [[ "$trimmed_files_not_empty" == false ]]; then
                log_warn "Trimmed output files are empty"
            fi
        fi
    done
	gzip_trimmed_fastq_files
	log_info "All SRR samples processed."
}

download_kingfisher_and_trim_srrs() {
	# Download and trim RNA-seq data using kingfisher with GNU parallel
	local SRR_LIST=("$@")

	# Default to global project list if none provided
	if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
		SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	fi

	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE JOBS
	export -f timestamp log log_info log_warn log_error run_with_space_time_log

	# Define per-SRR worker function for kingfisher
	_process_single_srr_kingfisher() {
		local SRR="$1"
		local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
		local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		local raw1 raw2

		log_info "Working on $SRR with kingfisher..."
		mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

		# Skip if trimmed files already exist
		if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
			log_info "Trimmed files for $SRR already exist. Skipping."
			return 0
		fi

		# Check for existing raw FASTQ files
		if [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
			raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
			log_info "Raw FASTQ files for $SRR already exist. Skipping download."
		else
			# Download with kingfisher (using fewer threads per job since we're running in parallel)
			log_info "Downloading $SRR with kingfisher..."
			local kingfisher_threads=$((THREADS / JOBS))
			# Add safety check to prevent division by zero
			[[ $kingfisher_threads -lt 1 ]] && kingfisher_threads=1
			
			run_with_space_time_log kingfisher get \
				--run-identifiers "$SRR" \
				--output-directory "$raw_files_DIR" \
				--download-threads "$kingfisher_threads" \
				--extraction-threads "$kingfisher_threads" \
				--gzip \
				--check-md5sums || {
				log_warn "Kingfisher download failed for $SRR"; return 1;
			}
			
			# Verify kingfisher output
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
				raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
				raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
				log_info "Kingfisher download successful for $SRR"
			else
				log_warn "Raw FASTQ for $SRR not found after kingfisher download; skipping."
				return 1
			fi
		fi

		# Trim reads using fewer cores per job
		log_info "Trimming $SRR..."
		local trim_cores=$((THREADS / JOBS))
		# Add safety check to prevent division by zero
		[[ $trim_cores -lt 1 ]] && trim_cores=1
		
		run_with_space_time_log trim_galore --cores "$trim_cores" \
			--paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

		# Verify trimming success and cleanup
		local trimming_success=false
		local trimmed_files_exist=false
		local trimmed_files_not_empty=false

		# Check if trimmed files exist and have content
		if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
			trimmed_files_exist=true
			
			if { [[ -s "$trimmed1" && -s "$trimmed2" ]] || [[ -s "${trimmed1}.gz" && -s "${trimmed2}.gz" ]]; }; then
				trimmed_files_not_empty=true
			fi
		fi

		# Process based on verification results
		if [[ "$trimmed_files_exist" == true && "$trimmed_files_not_empty" == true ]]; then
			trimming_success=true
			log_info "Trimming successfully completed for $SRR"
			log_info "Cleaning up raw FASTQs for $SRR..."
			rm -f "$raw1" "$raw2"
			log_info "Raw files deleted"
		else
			log_warn "Trimming verification failed for $SRR - keeping raw files"
			[[ "$trimmed_files_exist" == false ]] && log_warn "Trimmed files missing"
			[[ "$trimmed_files_not_empty" == false ]] && log_warn "Trimmed files are empty"
		fi

		log_info "Done with $SRR."
		log_info "--------------------------------------------------"
	}

	export -f _process_single_srr_kingfisher

	# Run jobs in parallel using GNU parallel
	log_info "Starting parallel kingfisher download and trimming..."
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${JOBS}" _process_single_srr_kingfisher {}
	log_info "All parallel kingfisher jobs completed."

	gzip_trimmed_fastq_files
	log_info "All SRR samples processed with kingfisher."
}

download_and_trim_srrs_parallel() {
	# Parallel Download and Trimming of SRR Samples

	local SRR_LIST=("$@")
	if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
		SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	fi

	export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE
	export -f timestamp log log_info log_warn log_error run_with_space_time_log

	# Define per-SRR worker function
	_process_single_srr() {
		local SRR="$1"
		local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
		local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		local raw1 raw2

		log_info "Working on $SRR..."
		mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

		# Skip if trimming already done
		if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
			log_info "Trimmed files for $SRR already exist. Skipping."
			log_info "--------------------------------------------------"
			return 0
		fi

		# Check for existing raw FASTQs
		if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq"
			raw2="$raw_files_DIR/${SRR}_2.fastq"
			log_info "Raw FASTQ files for $SRR already exist. Skipping download."
		elif [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
			raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
			log_info "Compressed FASTQ files for $SRR already exist. Skipping download."
		else
			# Download
			log_info "Downloading $SRR..."
			prefetch "$SRR" --output-directory "$raw_files_DIR" || {
				log_warn "Prefetch failed for $SRR"; return 1;
			}

			# Convert SRA → FASTQ
			fasterq-dump --split-files --threads "${THREADS}" \
				"$raw_files_DIR/$SRR/$SRR.sra" -O "$raw_files_DIR" || {
				log_warn "fasterq-dump failed for $SRR"; return 1;
			}

			# Compress FASTQ files to save space
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
				log_info "Compressing FASTQ files for $SRR..."
				gzip "$raw_files_DIR/${SRR}_1.fastq" "$raw_files_DIR/${SRR}_2.fastq"
			fi

			# Set raw file paths after download and compression
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
				raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
				raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
			elif [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
				raw1="$raw_files_DIR/${SRR}_1.fastq"
				raw2="$raw_files_DIR/${SRR}_2.fastq"
			else
				log_warn "Raw FASTQ for $SRR not found after download; skipping."
				log_info "--------------------------------------------------"
				return 1
			fi
		fi

		# Trim reads
		log_info "Trimming $SRR..."
		run_with_space_time_log trim_galore --cores "${THREADS}" \
			--paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

		# Verify trimming success and cleanup
		local trimming_success=false
		local trimmed_files_exist=false
		local trimmed_files_not_empty=false

		# Check existence and content of trimmed files
		if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
			trimmed_files_exist=true
			
			# Verify file sizes
			if { [[ -s "$trimmed1" && -s "$trimmed2" ]] || [[ -s "${trimmed1}.gz" && -s "${trimmed2}.gz" ]]; }; then
				trimmed_files_not_empty=true
			fi
		fi

		# Process based on verification results
		if [[ "$trimmed_files_exist" == true && "$trimmed_files_not_empty" == true ]]; then
			trimming_success=true
			log_info "Trimming successfully completed for $SRR"
			log_info "Cleaning up raw FASTQs for $SRR..."
			rm -f "$raw1" "$raw2"
			log_info "Raw files deleted"
		else
			log_warn "Trimming verification failed for $SRR - keeping raw files"
			if [[ "$trimmed_files_exist" == false ]]; then
				log_warn "Trimmed files missing"
			elif [[ "$trimmed_files_not_empty" == false ]]; then
				log_warn "Trimmed files are empty"
			fi
		fi

		log_info "Done with $SRR."
		log_info "--------------------------------------------------"
	}

	export -f _process_single_srr

	# Run jobs in parallel
	log_info "Starting parallel download and trimming..."
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr {}
	log_info "All parallel SRR jobs completed."

	gzip_trimmed_fastq_files
}

download_and_trim_srrs_wget_parallel() {
    # Parallel ENA FASTQ Download (wget) + Trim Galore

    local SRR_LIST=("$@")
    if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
        SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
    fi

    export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE
    export -f timestamp log log_info log_warn log_error run_with_space_time_log

    # Worker function for one SRR
    _process_single_srr_wget() {
        local SRR="$1"
        local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
        local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
        local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
        local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"

        log_info "PROCESSING: Working on $SRR..."
        mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

        # Skip if trimming already completed
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "SKIPPED: Trimmed files exist for $SRR - skipping."
            return 0
        fi

        # Fetch ENA FASTQ URLs dynamically
        local ena_links fq1 fq2
        ena_links=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" \
                     | tail -n +2)

        if [[ -z "$ena_links" ]]; then
            log_warn "WARNING: No ENA FASTQ links found for $SRR. Skipping."
            return 1
        fi

        fq1="ftp://$(echo "$ena_links" | cut -f1 -d';')"
        fq2="ftp://$(echo "$ena_links" | cut -f2 -d';')"

        log_info "DOWNLOADING: FASTQs for $SRR using wget..."
        wget -q -c -P "$raw_files_DIR" "$fq1" "$fq2"

        # Verify successful download
        if [[ ! -s "$raw_files_DIR/${SRR}_1.fastq.gz" || ! -s "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
            log_warn "FAILED: FASTQ download failed for $SRR. Check ENA mirrors."
            return 1
        fi

        # Trim Galore (low RAM, single-thread)
        log_info "TRIMMING: Adapters for $SRR..."
        run_with_space_time_log trim_galore --cores 1 \
            --paired "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz" \
            --output_dir "$TrimGalore_DIR"

        # Cleanup if trimming successful
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "CLEANUP: Removing raw FASTQs for $SRR..."
            rm -f "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz"
        fi

        log_info "COMPLETED: Finished $SRR."
        log_info "--------------------------------------------------"
    }

    export -f _process_single_srr_wget

    # Run jobs in parallel
    log_info "STARTING: Parallel ENA downloads and trimming..."
    printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr_wget {}
    log_info "COMPLETED: All parallel wget + trimming jobs complete."
}
