#!/bin/bash
# ==============================================================================
# METHOD 3: STAR SPLICE-AWARE ALIGNMENT PIPELINE
# ==============================================================================
# STAR splice-aware alignment with Salmon quantification and tximport for DESeq2
# Uses 2-pass mode for better junction detection
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M3_STAR_SOURCED:-}" == "true" ]] && return 0
export M3_STAR_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# STAR CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# STAR-specific settings (GPU-aware defaults)
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-$(get_star_genome_load 2>/dev/null || echo NoSharedMemory)}"
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"
STAR_STRAND_SPECIFIC="${STAR_STRAND_SPECIFIC:-intronMotif}"

# GTF annotation file for splice junction detection (required for --sjdbOverhang)
STAR_GTF_FILE="${STAR_GTF_FILE:-0_INPUT_FASTAs/gtf/reference/Eggplant_V4.1_function_IPR_final_formatted_v3.gtf}"

# ==============================================================================
# STAR TISSUE-SPECIFIC PIPELINE
# ==============================================================================

star_tissue_specific_pipeline() {
	local fasta="" rnaseq_list=() metadata_file=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--METADATA) metadata_file="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done
	
	[[ -z "$fasta" || ! -f "$fasta" ]] && { log_error "Valid FASTA required"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	
	declare -A sample_metadata
	if [[ -n "$metadata_file" && -f "$metadata_file" ]]; then
		load_sample_metadata "$metadata_file" sample_metadata || {
			log_warn "Metadata load failed - running pooled alignment"
			star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			return $?
		}
	else
		log_warn "No metadata - running pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	# Group samples by tissue
	declare -A tissue_samples
	for SRR in "${rnaseq_list[@]}"; do
		local tissue="${sample_metadata[${SRR}_condition]:-unknown}"
		tissue_samples["$tissue"]+="$SRR "
	done
	
	local tissue_count=${#tissue_samples[@]}
	log_info "[STAR TISSUE] Found $tissue_count tissue types"
	
	if [[ $tissue_count -lt 2 ]]; then
		log_warn "Single tissue detected - using pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	log_info "[STAR TISSUE] Running tissue-specific alignments for better isoform detection"
	
	# Run STAR per tissue
	for tissue in "${!tissue_samples[@]}"; do
		local tissue_srrs=(${tissue_samples[$tissue]})
		log_step "STAR alignment for tissue: $tissue (${#tissue_srrs[@]} samples)"
		
		star_alignment_pipeline \
			--FASTA "$fasta" \
			--RNASEQ_LIST "${tissue_srrs[@]}" \
			--TISSUE_TAG "$tissue"
	done
}

# ==============================================================================
# MAIN STAR ALIGNMENT PIPELINE
# ==============================================================================

star_alignment_pipeline() {
	local fasta="" rnaseq_list=() tissue_tag=""
	
	# Check for GNU parallel
	if ! command -v parallel >/dev/null 2>&1; then
		log_warn "[PERFORMANCE] GNU parallel not found - Salmon quantification will run sequentially"
	fi
	
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--TISSUE_TAG) tissue_tag="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done
	
	[[ -z "$fasta" ]] && { log_error "No FASTA file specified. Use --FASTA <fasta_file>."; return 1; }
	[[ ! -f "$fasta" ]] && { log_error "FASTA file '$fasta' not found."; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	[[ ${#rnaseq_list[@]} -eq 0 ]] && { log_error "No RNA-seq samples provided."; return 1; }

	local fasta_base fasta_tag star_index_dir star_genome_dir
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	
	# Convert STAR roots to absolute paths to avoid issues when running from different directories
	local abs_star_index_root abs_star_align_root
	abs_star_index_root="$(cd "$PROJECT_ROOT" 2>/dev/null && mkdir -p "$STAR_INDEX_ROOT" && cd "$STAR_INDEX_ROOT" && pwd)"
	abs_star_align_root="$(cd "$PROJECT_ROOT" 2>/dev/null && mkdir -p "$STAR_ALIGN_ROOT" && cd "$STAR_ALIGN_ROOT" && pwd)"
	
	# Fallback to original paths if absolute conversion fails
	[[ -z "$abs_star_index_root" ]] && abs_star_index_root="$STAR_INDEX_ROOT"
	[[ -z "$abs_star_align_root" ]] && abs_star_align_root="$STAR_ALIGN_ROOT"
	
	# Set directories based on tissue tag (using absolute paths)
	if [[ -n "$tissue_tag" ]]; then
		star_index_dir="${abs_star_index_root}/${fasta_tag}_${tissue_tag}_star_index"
		star_genome_dir="${abs_star_align_root}/${fasta_tag}_${tissue_tag}_alignments"
		log_info "[STAR] Tissue-specific alignment for: $tissue_tag (${#rnaseq_list[@]} samples)"
	else
		star_index_dir="${abs_star_index_root}/${fasta_tag}_star_index"
		star_genome_dir="${abs_star_align_root}/${fasta_tag}_alignments"
		log_info "[STAR] Pooled alignment: ${#rnaseq_list[@]} samples"
	fi
	
	# Log resolved paths for debugging
	log_info "[STAR] Index directory (absolute): $star_index_dir"
	log_info "[STAR] Output directory (absolute): $star_genome_dir"
	
	local star_overhang=$((STAR_READ_LENGTH - 1))
	log_info "[STAR] CPU allocation: Total=$THREADS threads"

	# STEP 1: BUILD STAR GENOME INDEX
	log_step "STAR genome index generation for $fasta_tag"
	
	if [[ -f "$star_index_dir/SAindex" ]]; then
		log_info "[STAR INDEX] Index for $fasta_tag already exists. Skipping."
	else
		log_info "[STAR INDEX] Building STAR genome index from transcriptome..."
		mkdir -p "$star_index_dir"
		
		log_file_size "$fasta" "Input transcriptome FASTA for STAR indexing"
		
		# Validate GTF file exists
		if [[ ! -f "$STAR_GTF_FILE" ]]; then
			log_error "GTF annotation file not found: $STAR_GTF_FILE"
			log_error "Set STAR_GTF_FILE to a valid GTF file path"
			return 1
		fi
		log_info "[STAR INDEX] Using GTF annotation: $STAR_GTF_FILE"
		
		run_with_space_time_log --input "$fasta" --output "$star_index_dir" \
			STAR --runMode genomeGenerate \
				--genomeDir "$star_index_dir" \
				--genomeFastaFiles "$fasta" \
				--sjdbGTFfile "$STAR_GTF_FILE" \
				--sjdbOverhang "$star_overhang" \
				--genomeSAindexNbases 12 \
				--runThreadN "$THREADS"
		
		[[ ! -f "$star_index_dir/SAindex" ]] && { log_error "STAR index generation failed"; return 1; }
		log_info "[STAR INDEX] Index built successfully: $star_index_dir"
	fi
	
	# STEP 2: ALIGN READS WITH STAR (2-PASS MODE)
	mkdir -p "$star_genome_dir"
	
	# Verify directory exists and is writable
	if [[ ! -d "$star_genome_dir" ]]; then
		log_error "Failed to create STAR output directory: $star_genome_dir"
		return 1
	fi
	if [[ ! -w "$star_genome_dir" ]]; then
		log_error "STAR output directory is not writable: $star_genome_dir"
		return 1
	fi
	log_info "[STAR] Output directory verified: $star_genome_dir"
	
	log_step "STAR splice-aware alignment for $fasta_tag samples"
	
	for SRR in "${rnaseq_list[@]}"; do
		local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"
		
		[[ -f "$bam_output" ]] && { log_info "[STAR] Alignment for $SRR already exists. Skipping."; continue; }
		
		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found; skipping."; continue; }
		
		log_info "[STAR] Aligning: $SRR"
		
		local star_reads=""
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			star_reads="$trimmed1 $trimmed2"
			log_info "[STAR] Processing paired-end reads for $SRR"
		else
			star_reads="$trimmed1"
			log_info "[STAR] Processing single-end reads for $SRR"
		fi
		
		# Create dedicated temp directory for STAR inside the output directory
		# Using absolute path within the alignment output to avoid permission/conflict issues
		local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
		rm -rf "$star_tmp_dir"  # STAR requires the tmp dir to not exist
		
		# Clean up any stale STAR temporary files from previous failed runs
		rm -rf "${star_genome_dir}/${SRR}_STARtmp" 2>/dev/null || true
		rm -rf "${star_genome_dir}/${SRR}_"*.tmp 2>/dev/null || true
		
		# CRITICAL: Ensure output directory exists and is writable RIGHT BEFORE STAR runs
		# This fixes "could not create output file" errors on HPC/server environments
		mkdir -p "$star_genome_dir"
		if [[ ! -d "$star_genome_dir" ]]; then
			log_error "[STAR] FATAL: Cannot create output directory: $star_genome_dir"
			return 1
		fi
		
		# Test write permissions by creating a test file
		local test_file="${star_genome_dir}/.write_test_$$"
		if ! touch "$test_file" 2>/dev/null; then
			log_error "[STAR] FATAL: Output directory not writable: $star_genome_dir"
			return 1
		fi
		rm -f "$test_file"
		
		log_info "[STAR] Output path verified writable: $star_genome_dir/${SRR}_*"
		
		run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$star_genome_dir" \
			STAR --runMode alignReads \
				--genomeDir "$star_index_dir" \
				--readFilesIn $star_reads \
				--readFilesCommand zcat \
				--outFileNamePrefix "$star_genome_dir/${SRR}_" \
				--outTmpDir "$star_tmp_dir" \
				--outSAMtype BAM SortedByCoordinate \
				--outSAMstrandField "$STAR_STRAND_SPECIFIC" \
				--outSAMunmapped Within \
				--outFilterType BySJout \
				--outFilterMultimapNmax 20 \
				--alignSJoverhangMin 8 \
				--alignSJDBoverhangMin 1 \
				--outFilterMismatchNmax 999 \
				--outFilterMismatchNoverReadLmax 0.04 \
				--alignIntronMin 20 \
				--alignIntronMax 1000000 \
				--alignMatesGapMax 1000000 \
				--limitBAMsortRAM 31000000000 \
				--twopassMode Basic \
				--runThreadN "$THREADS"
		
		[[ ! -f "$bam_output" ]] && { log_error "STAR alignment failed for $SRR"; return 1; }
		
		# Clean up temp directory after successful alignment
		rm -rf "$star_tmp_dir" 2>/dev/null || true
		
		log_info "[STAR] Successfully aligned: $SRR"
	done
	
	log_info "[STAR] All samples aligned successfully"
	
	# STEP 3: SALMON QUANTIFICATION
	local tag_suffix=""
	[[ -n "${tissue_tag:-}" ]] && tag_suffix="_${tissue_tag}"
	local salmon_idx="${abs_star_align_root}/${fasta_tag}${tag_suffix}_salmon_index"
	local quant_root="${abs_star_align_root}/${fasta_tag}${tag_suffix}_salmon_quant"
	mkdir -p "$quant_root"
	
	# Build Salmon index
	if [[ -f "$salmon_idx/versionInfo.json" ]]; then
		log_info "[SALMON INDEX] STAR Salmon index exists - skipping build"
	else
		log_step "Building Salmon index for transcriptome"
		run_with_space_time_log --input "$fasta" --output "$salmon_idx" \
			salmon index -t "$fasta" -i "$salmon_idx" -k 31 --threads "$THREADS"
	fi
	
	# Quantify samples
	log_info "[SALMON] Starting quantification for ${#rnaseq_list[@]} samples"
	
	for SRR in "${rnaseq_list[@]}"; do
		local quant_dir="$quant_root/$SRR"
		
		[[ -f "$quant_dir/quant.sf" ]] && { log_info "[SALMON] Quantification for $SRR exists - skipping"; continue; }
		
		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" ]] && { log_warn "No trimmed reads for $SRR - skipping"; continue; }
		
		log_info "[SALMON] Quantifying: $SRR"
		mkdir -p "$quant_dir"
		
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
				--validateMappings -l A -1 "$trimmed1" -2 "$trimmed2"
		else
			run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
				--validateMappings -l A -r "$trimmed1"
		fi
		
		[[ -f "$quant_dir/quant.sf" ]] && log_info "[SALMON] Successfully quantified: $SRR"
	done

	# STEP 4: PREPARE TXIMPORT FILES FOR DESEQ2
	log_step "Preparing tximport input for DESeq2 (STAR + Salmon)"
	
	local method3_root="${STAR_ALIGN_ROOT%/*}"
	local matrix_dir="$method3_root/6_matrices_from_stringtie"
	mkdir -p "$matrix_dir"
	
	local sample_metadata="$matrix_dir/sample_info.tsv"
	local tx2gene_file="$matrix_dir/tx2gene_${fasta_tag}${tag_suffix}.tsv"
	
	# Create tx2gene mapping
	if [[ ! -f "$tx2gene_file" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping"
		awk '/^>/{tx=$1; gsub(/^>/, "", tx); gene=tx; gsub(/\.[0-9]+$/, "", gene); print tx "\t" gene}' "$fasta" > "$tx2gene_file"
		log_info "[TXIMPORT] Created tx2gene mapping: $(wc -l < "$tx2gene_file") transcripts"
	fi
	
	# Verify quantifications
	local quant_count=0
	for SRR in "${rnaseq_list[@]}"; do
		[[ -f "$quant_root/$SRR/quant.sf" ]] && ((quant_count++))
	done
	
	[[ $quant_count -lt 2 ]] && { log_error "Insufficient quantifications: $quant_count (need ≥2)"; return 1; }
	
	# Create sample metadata
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" rnaseq_list[@]
	
	# Generate tximport R script
	local tximport_script="$matrix_dir/run_tximport_star_salmon.R"
	if [[ ! -f "$tximport_script" ]]; then
		generate_tximport_star_script "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" "$tximport_script"
	fi
	
	# Run tximport if R is available
	if command -v Rscript >/dev/null 2>&1; then
		log_step "Running tximport to import Salmon quantifications"
		if Rscript "$tximport_script" "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" 2>&1 | tee -a "$LOG_FILE"; then
			log_info "[TXIMPORT] Successfully imported counts for DESeq2"
			local count_matrix="$matrix_dir/gene_counts_tximport.tsv"
			[[ -f "$count_matrix" ]] && validate_count_matrix "$count_matrix" "gene" 2
		else
			log_warn "[TXIMPORT] tximport failed - check R dependencies"
		fi
	else
		log_warn "[TXIMPORT] Rscript not found - run manually: Rscript $tximport_script"
	fi
	
	log_step "STAR alignment pipeline completed for $fasta_tag"
	log_info "STAR alignments: $star_genome_dir"
	log_info "Salmon quantifications: $quant_root"
	log_info "Tximport outputs: $matrix_dir"
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Run tximport for STAR+Salmon using external R helper
# Usage: run_tximport_star <quant_dir> <metadata_file> <tx2gene_file> [output_dir]
run_tximport_star() {
	local quant_dir="$1"
	local metadata_file="$2"
	local tx2gene_file="$3"
	local output_dir="${4:-$(dirname "$metadata_file")}"
	local helper_script="$SCRIPT_DIR/helpers/tximport_star_helper.R"
	
	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi
	
	log_info "[TXIMPORT] Running STAR+Salmon import..."
	Rscript "$helper_script" "$quant_dir" "$metadata_file" "$tx2gene_file" "$output_dir"
}

# Generate tximport script (legacy - copies helper)
generate_tximport_star_script() {
	local quant_root="$1"
	local sample_metadata="$2"
	local tx2gene_file="$3"
	local matrix_dir="$4"
	local output_script="$5"
	local helper_script="$SCRIPT_DIR/helpers/tximport_star_helper.R"
	
	if [[ -f "$helper_script" ]]; then
		cp "$helper_script" "$output_script"
		chmod +x "$output_script"
		log_info "[TXIMPORT] Copied helper to: $output_script"
	else
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi
}
