#!/bin/bash

# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE: TISSUE
# ==============================================================================
# Description: RNA-seq analysis pipeline using multiple alignment approaches
# Author: Mark Cyril R. Mercado
# Version: v10
# Date: October 2025
# 
# Pipeline Methods:
# - Method 1: HISAT2 Reference Guided 
# - Method 2: HISAT2 De Novo Assembly
# - Method 3: Trinity De Novo Assembly
# - Method 4: Salmon SAF Quantification
# - Method 5: Bowtie2 + RSEM Quantification
#
# SCRIPT ORGANIZATION:
# 1. Configuration 
# 2. Pipeline Configuration Examples  
# 3. Input Files and Data Sources
# 4. Directory Structure and Output Paths  
# 5. Cleanup Options and Testing
# 6. Logging System and Utility Functions
# 7. Pipeline Functions:
#    	- Data Download and Preprocessing Functions
#		- Pipelines Functions	 
# 8. Main Execution Functions
# 9. Script Execution
# 10. Post-processing Options
# ==============================================================================

# Allow pipeline to continue on errors (individual methods handle their own failures)
set -uo pipefail
source "modules/logging_utils.sh"
source "modules/pipeline_utils.sh"

# ------------------------------------------------------------------------------
# Global safe defaults and guards
# ------------------------------------------------------------------------------
# Provide safe defaults for commonly referenced environment variables to avoid
# unbound variable errors under 'set -u' and ensure reasonable behavior.
THREADS="${THREADS:-4}"
keep_bam_global="${keep_bam_global:-n}"

# Ensure SRR arrays exist (some pipelines fall back to these lists)
if ! declare -p SRR_COMBINED_LIST >/dev/null 2>&1; then
	declare -a SRR_COMBINED_LIST=()
fi
if ! declare -p SRR_LIST_PRJNA328564 >/dev/null 2>&1; then
	declare -a SRR_LIST_PRJNA328564=()
fi

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================
# Standardized folder structure across all methods in 4_POST_PROCESSING/
# Pattern: 4_{alignment_tool}_WD, 5_{quant_tool}_WD, 6_matrices, 7_Heatmap_Outputs
# Common: a_gene_groups_input_list, b_modules_for_{method}, logs

# Raw and Processed Data Directories
RAW_DIR_ROOT="1_RAW_SRR"
TRIM_DIR_ROOT="2_TRIMMED_SRR"
FASTQC_ROOT="3_FastQC"

# Method 1: HISAT2 Reference Guided
HISAT2_REF_GUIDED_ROOT="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/5_stringtie_WD"
HISAT2_REF_GUIDED_MATRIX_ROOT="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided/6_matrices"

# Method 2: HISAT2 De Novo Assembly
HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/5_stringtie_WD"
HISAT2_DE_NOVO_MATRIX_ROOT="4_POST_PROCESSING/4b_Method_2_HISAT2_De_Novo/6_matrices"

# Method 3: Trinity De Novo Assembly
TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4_POST_PROCESSING/4c_Method_3_Trinity_De_Novo/5_stringtie_WD"

# Method 4: Salmon SAF Quantification
SALMON_SAF_ROOT="4_POST_PROCESSING/4d_Method_4_Salmon_Saf_Quantification"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/4_Salmon_WD/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/5_Salmon_Quant_WD"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/6_matrices_from_Salmon"
SALMON_MATRIX_ROOT="${SALMON_MATRIX_ROOT:-$SALMON_SAF_MATRIX_ROOT}"

# Method 5: Bowtie2 + RSEM Quantification
BOWTIE2_RSEM_ROOT="4_POST_PROCESSING/4e_Method_5_Bowtie2_Quantification"
RSEM_INDEX_ROOT="$BOWTIE2_RSEM_ROOT/4_Bowtie2_WD/index"
RSEM_QUANT_ROOT="$BOWTIE2_RSEM_ROOT/5_RSEM_Quant_WD"
RSEM_MATRIX_ROOT="$BOWTIE2_RSEM_ROOT/6_matrices_from_RSEM"

# Create required directories
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$TRINITY_DE_NOVO_ROOT" "$STRINGTIE_TRINITY_DE_NOVO_ROOT" \
	"$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
	"$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT"

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS FOR METADATA AND VALIDATION
# ------------------------------------------------------------------------------

# Generate tximport R script for proper RSEM/Salmon import into DESeq2
generate_tximport_script() {
	local method="$1"  # "rsem" or "salmon"
	local quant_dir="$2"
	local output_script="$3"
	local metadata_file="$4"
	
	# Path to the template R script
	local template_script="modules/tximport.R"
	
	# Check if template exists
	if [[ ! -f "$template_script" ]]; then
		log_error "[TXIMPORT] Template script not found: $template_script"
		log_error "[TXIMPORT] Please ensure tximport.R exists in the modules directory"
		return 1
	fi
	
	# Copy template to output location
	cp "$template_script" "$output_script"
	
	# Replace placeholders with actual values using sed
	sed -i "s|METHOD_PLACEHOLDER|$method|g" "$output_script"
	sed -i "s|QUANT_DIR_PLACEHOLDER|$quant_dir|g" "$output_script"
	sed -i "s|METADATA_FILE_PLACEHOLDER|$metadata_file|g" "$output_script"
	
	chmod +x "$output_script"
	log_info "[TXIMPORT] Generated tximport R script: $output_script"
	log_info "[TXIMPORT] Run this script for proper import into DESeq2: Rscript $output_script"
}

# Load sample metadata from external file for proper experimental design
load_sample_metadata() {
	local metadata_file="${1:-sample_conditions.txt}"
	local -n metadata_array=$2  # Use nameref to return associative array
	
	if [[ ! -f "$metadata_file" ]]; then
		log_warn "Metadata file not found: $metadata_file - Using default (all samples as 'treatment')"
		log_warn "⚠️  Differential expression analysis disabled without proper metadata"
		return 1
	fi
	
	# Validate metadata file structure
	local line_count=$(wc -l < "$metadata_file")
	if [[ $line_count -lt 2 ]]; then
		log_error "Metadata file must contain at least 2 samples"
		return 1
	fi
	
	# Parse metadata into associative array
	while IFS=$'\t' read -r srr condition batch; do
		[[ -z "$srr" || "$srr" =~ ^# ]] && continue  # Skip empty/comment lines
		metadata_array["${srr}_condition"]="$condition"
		metadata_array["${srr}_batch"]="$batch"
	done < "$metadata_file"
	
	# Validate that we have at least 2 different conditions
	local unique_conditions=$(printf '%s\n' "${metadata_array[@]}" | grep -v "_batch" | sort -u | wc -l)
	if [[ $unique_conditions -lt 2 ]]; then
		log_error "Need ≥2 experimental conditions for DE analysis (found: $unique_conditions)"
		return 1
	fi
	
	log_info "Loaded metadata: $(( ${#metadata_array[@]} / 2 )) samples, $unique_conditions conditions"
	return 0
}

# Create sample metadata CSV for DESeq2 with proper condition assignment
create_sample_metadata() {
	local metadata_file="$1"
	local -a sample_list=("${!2}")
	local metadata_source="${3:-sample_conditions.txt}"
	local output_format="${4:-csv}"  # csv or tsv
	
	# Try to load external metadata
	declare -A sample_metadata
	local has_external_metadata=false
	if load_sample_metadata "$metadata_source" sample_metadata 2>/dev/null; then
		has_external_metadata=true
		log_info "Using experimental conditions from: $metadata_source"
	fi
	
	# Determine delimiter based on output format
	local delim=","
	if [[ "$output_format" == "tsv" || "$metadata_file" == *.tsv ]]; then
		delim=$'\t'
	fi
	
	# Create metadata file with appropriate format
	echo -e "sample${delim}condition${delim}batch" > "$metadata_file"
	for SRR in "${sample_list[@]}"; do
		if [[ "$has_external_metadata" == "true" ]]; then
			local condition="${sample_metadata[${SRR}_condition]:-unknown}"
			local batch="${sample_metadata[${SRR}_batch]:-1}"
			echo -e "$SRR${delim}$condition${delim}$batch" >> "$metadata_file"
		else
			# Default fallback - assign placeholder conditions
			# Try to infer from SRR ID patterns, otherwise use generic names
			local condition="condition_${SRR}"
			echo -e "$SRR${delim}$condition${delim}1" >> "$metadata_file"
		fi
	done
	
	if [[ "$has_external_metadata" == "false" ]]; then
		log_warn "⚠️  Sample conditions need manual specification for DE analysis"
		log_warn "📝 Edit: $metadata_file or create: $metadata_source"
		log_warn "   Format (tab-separated): SRR_ID<tab>tissue/organ<tab>batch"
		log_warn "   Example: SRR3884597<tab>root<tab>1"
		log_warn "⚠️  DE requires ≥2 replicates/condition, ≥3 recommended"
		
		# Generate example metadata template
		local metadata_dir=$(dirname "$metadata_file")
		local template_file="$metadata_dir/sample_info_template.txt"
		if [[ ! -f "$template_file" ]]; then
			{
				echo "# SAMPLE METADATA TEMPLATE FOR TISSUE-SPECIFIC EXPRESSION PROFILING"
				echo "# Copy to 'sample_conditions.txt' in modules/ directory and edit condition/batch columns"
				echo "#"
				echo "# Format: SRR_ID<TAB>condition<TAB>batch"
				echo "# condition: Tissue/organ name (e.g., root, leaf, flower, stem)"
				echo "# batch: Experimental batch number (use 1 if all samples from same batch)"
				echo "#"
				for SRR in "${sample_list[@]}"; do
					echo "$SRR	tissue_name_here	1"
				done
			} > "$template_file"
			log_info "Generated template: $template_file"
		fi
	fi
}

# Validate count matrix quality and structure
validate_count_matrix() {
	local matrix="$1"
	local matrix_type="${2:-gene}"  # gene or transcript
	local min_samples="${3:-2}"
	
	[[ ! -f "$matrix" ]] && { log_error "Matrix not found: $matrix"; return 1; }
	
	log_step "Validating $matrix_type matrix: $(basename $matrix)"
	
	# Check number of samples
	local header=$(head -n1 "$matrix")
	local num_samples
	if [[ "$matrix" == *.csv ]]; then
		num_samples=$(echo "$header" | tr ',' '\n' | tail -n +2 | wc -l)
	else
		num_samples=$(echo "$header" | tr '\t' '\n' | tail -n +2 | wc -l)
	fi
	
	# Check sufficient samples
	if [[ $num_samples -lt $min_samples ]]; then
		log_error "Insufficient samples: $num_samples (need ≥$min_samples)"
		return 1
	fi
	
	# Check for duplicate IDs
	local total_rows=$(tail -n +2 "$matrix" | wc -l)
	local unique_ids
	if [[ "$matrix" == *.csv ]]; then
		unique_ids=$(tail -n +2 "$matrix" | cut -d',' -f1 | sort -u | wc -l)
	else
		unique_ids=$(tail -n +2 "$matrix" | cut -f1 | sort -u | wc -l)
	fi
	
	if [[ $total_rows -ne $unique_ids ]]; then
		log_error "Duplicate ${matrix_type} IDs detected!"
		log_error "Total rows: $total_rows, Unique IDs: $unique_ids"
		return 1
	fi
	
	# Check for zero-count genes/transcripts
	local zero_count
	if [[ "$matrix" == *.csv ]]; then
		zero_count=$(awk -F',' 'NR>1 {s=0; for(i=2;i<=NF;i++) s+=$i; if(s==0) count++} END {print count+0}' "$matrix")
	else
		zero_count=$(awk 'NR>1 {s=0; for(i=2;i<=NF;i++) s+=$i; if(s==0) count++} END {print count+0}' "$matrix")
	fi
	
	local zero_pct=$((total_rows > 0 ? zero_count * 100 / total_rows : 0))
	
	log_info "[VALIDATION] Passed - Samples: $num_samples, ${matrix_type^}s: $unique_ids, Zero: $zero_count ($zero_pct%)"
	
	[[ $zero_pct -gt 50 ]] && log_warn "⚠️  >50% ${matrix_type}s have zero counts - check mapping/depth"
	
	return 0
}

# Improved read length detection with fallback handling
detect_read_length() {
	local fastq="$1"
	local default_length="${2:-150}"
	
	[[ ! -f "$fastq" ]] && { echo "$default_length"; return 1; }
	
	# Determine decompression command
	local decompress_cmd="cat"
	case "$fastq" in
		*.gz) decompress_cmd="zcat" ;;
		*.bz2) decompress_cmd="bzcat" ;;
	esac
	
	# Sample first 1000 reads to get average length
	local avg_length=$($decompress_cmd "$fastq" 2>/dev/null | \
		awk 'NR%4==2 {sum+=length($0); count++} count==1000 {print int(sum/count); exit}')
	
	# Validate detected length
	if [[ -z "$avg_length" || $avg_length -lt 50 || $avg_length -gt 300 ]]; then
		log_warn "Could not detect valid read length from $fastq (got: ${avg_length:-N/A})"
		log_warn "Using default: $default_length bp"
		echo "$default_length"
		return 1
	fi
	
	echo "$avg_length"
	return 0
}

# ------------------------------------------------------------------------------
# MAIN PIPELINES FOR ALIGNMENT, ASSEMBLY, AND QUANTIFICATION
# ------------------------------------------------------------------------------

run_quality_control() {
	# Run quality control analysis on trimmed reads
	local SRR="$1"
	local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
	
	[[ ! -d "$TrimGalore_DIR" ]] && { log_warn "Trimmed dir not found for $SRR - skipping QC"; return 1; }
	
	log_step "QC: $SRR"
	
	# FastQC on trimmed reads
	if command -v fastqc >/dev/null 2>&1; then
		local outdir="$FASTQC_ROOT/${SRR}_trimmed"
		mkdir -p "$outdir"
		if compgen -G "$outdir/*_fastqc.html" >/dev/null; then
			log_info "[FASTQC] Report exists - skipping"
		else
			run_with_space_time_log \
				fastqc -t "${THREADS:-1}" -o "$outdir" \
					"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || \
					log_warn "[FASTQC] Failed for $SRR"
		fi
	else
		log_warn "[FASTQC] Not available - skipping read quality assessment"
	fi
	
	# MultiQC summary
	command -v multiqc >/dev/null 2>&1 && \
		run_with_space_time_log multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || \
		log_warn "[MULTIQC] Failed - individual reports available"
}

# Combined pipeline: Build HISAT2 reference-guided index, align reads, assemble, merge, and quantify transcripts
hisat2_ref_guided_pipeline() {
	# HISAT2 Reference Guided Pipeline using reference GTF and genome
	local fasta="" gtf="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)
				fasta="$2"; shift 2;;
			--GTF)
				gtf="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1")
					shift
				done
				;;
			*)
				log_error "Unknown option: $1"
				return 1;;
		esac
	done
	
	if [[ -z "$fasta" ]]; then
		log_error "No FASTA file specified. Use --FASTA <genome_fasta>."
		return 1
	fi
	if [[ ! -f "$fasta" ]]; then
		log_error "FASTA file '$fasta' not found."
		return 1
	fi
	if [[ -z "$gtf" ]]; then
		log_error "No GTF file specified. Use --GTF <annotation_gtf>."
		return 1
	fi
	if [[ ! -f "$gtf" ]]; then
		log_error "GTF file '$gtf' not found."
		return 1
	fi
	if [[ ${#rnaseq_list[@]} -eq 0 ]]; then
		rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	fi

	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_ref_guided"

	# BUILD HISAT2 REFERENCE-GUIDED INDEX
	mkdir -p "$HISAT2_REF_GUIDED_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "[INDEX] Ref-Guided index exists - skipping build"
	else
		log_step "Building HISAT2 Ref-Guided index: $fasta_base"
		
		# Extract splice sites and exons from GTF
		local splice_sites="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_splice_sites.txt"
		local exons="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_exons.txt"
		
		# Extract and validate
		hisat2_extract_splice_sites.py "$gtf" > "$splice_sites" || \
			{ log_error "Failed to extract splice sites - check GTF format/transcript features"; return 1; }
		
		hisat2_extract_exons.py "$gtf" > "$exons" || \
			{ log_error "Failed to extract exons - check GTF format/exon features"; return 1; }
		
		[[ ! -s "$splice_sites" ]] && { log_error "Empty splice sites file - GTF/FASTA mismatch"; return 1; }
		[[ ! -s "$exons" ]] && { log_error "Empty exons file - GTF/FASTA mismatch"; return 1; }
		
		log_info "[INDEX] Extracted $(wc -l < "$splice_sites") splice sites, $(wc -l < "$exons") exons"
		
		# Build index with splice sites and exons
		log_file_size "$fasta" "Input FASTA for HISAT2 index - $fasta_tag"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_REF_GUIDED_INDEX_DIR" \
			hisat2-build \
				-p "${THREADS}" \
				--ss "$splice_sites" \
				--exon "$exons" \
				"$fasta" \
				"$index_prefix"
		log_file_size "$HISAT2_REF_GUIDED_INDEX_DIR" "HISAT2 index output - $fasta_tag"
	fi

	# ALIGNMENT, SORTING, AND STRINGTIE ASSEMBLY for each SRR
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="" trimmed2=""
		mkdir -p "$HISAT2_DIR"
		
		# Find trimmed FASTQ files (paired-end or single-end)
		if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
		elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1 && compgen -G "$TrimGalore_DIR/${SRR}*val_2.*" >/dev/null 2>&1; then
			local files1=( "$TrimGalore_DIR"/${SRR}*val_1.* )
			local files2=( "$TrimGalore_DIR"/${SRR}*val_2.* )
			trimmed1="${files1[0]}"
			trimmed2="${files2[0]}"
		elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq" ]]; then
			# Single-end reads (uncompressed)
			trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq"
			trimmed2=""
		elif [[ -f "$TrimGalore_DIR/${SRR}_trimmed.fq.gz" ]]; then
			# Single-end reads (compressed)
			trimmed1="$TrimGalore_DIR/${SRR}_trimmed.fq.gz"
			trimmed2=""
		elif compgen -G "$TrimGalore_DIR/${SRR}*trimmed.fq*" >/dev/null 2>&1; then
			# Generic single-end pattern
			local files=( "$TrimGalore_DIR"/${SRR}*trimmed.fq* )
			trimmed1="${files[0]}"
			trimmed2=""
		fi
		
		if [[ -z "$trimmed1" ]]; then
			log_warn "Trimmed FASTQ not found for $SRR - skipping"
			continue
		fi

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
		local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped.sam"
		
		# Align and sort if BAM doesn't exist
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "[ALIGN] BAM exists for $SRR/$fasta_tag - skipping"
		else
			log_step "Aligning: $SRR -> $fasta_tag (HISAT2 Ref-Guided)"
			# Check if paired-end or single-end reads
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				# Paired-end alignment
				log_input_output_size "$trimmed1" "$HISAT2_DIR" "HISAT2 alignment input - $SRR"
				run_with_space_time_log --input "$TrimGalore_DIR" --output "$HISAT2_DIR" \
					hisat2 -p "${THREADS}" --dta \
						-x "$index_prefix" \
						-1 "$trimmed1" \
						-2 "$trimmed2" \
						-S "$sam"
			else
				# Single-end alignment
				run_with_space_time_log \
					hisat2 -p "${THREADS}" --dta \
						-x "$index_prefix" \
						-U "$trimmed1" \
						-S "$sam"
			fi
			
			log_info "[SAMTOOLS] Converting to sorted BAM..."
			log_file_size "$sam" "SAM file before sorting - $SRR"
			run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			log_file_size "$bam" "Sorted BAM file - $SRR"
			run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
			log_info "[CLEANUP] Deleting SAM file to save space"
			rm -f "$sam"
		fi

		# StringTie assembly with reference GTF
		local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		local out_gene_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv"
		local ballgown_dir="$out_dir/ballgown"
		mkdir -p "$out_dir" "$ballgown_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "[STRINGTIE] Assembly exists for $SRR/$fasta_tag - skipping"
		else
			log_step "Assembling transcripts: $SRR -> $fasta_tag (ref GTF)"
			log_input_output_size "$bam" "$out_dir" "StringTie assembly - $SRR"
			run_with_space_time_log --input "$bam" --output "$out_dir" \
				stringtie -p "$THREADS" "$bam" \
					-G "$gtf" \
					-o "$out_gtf" \
					-A "$out_gene_abundances_tsv" \
					-B \
					-C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf"
			log_file_size "$out_gtf" "StringTie output GTF - $SRR"
		fi
	done
	
	# Create merged GTF for all samples
	log_step "Creating merged GTF file for reference-guided assembly"
	local merge_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/merged"
	local merged_gtf="$merge_dir/${fasta_tag}_ref_guided_merged.gtf"
	local gtf_list="$merge_dir/gtf_list.txt"
	mkdir -p "$merge_dir"
	
	# Create list of GTF files for merging
	true > "$gtf_list"
	for SRR in "${rnaseq_list[@]}"; do
		local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		if [[ -f "$out_gtf" ]]; then
			echo "$out_gtf" >> "$gtf_list"
		fi
	done
	
	if [[ -f "$merged_gtf" ]]; then
		log_info "[STRINGTIE MERGE] Merged GTF for $fasta_tag (Ref-Guided) already exists. Skipping merge."
	else
		log_info "[STRINGTIE MERGE] Merging GTF files for $fasta_tag (Ref-Guided)..."
		run_with_space_time_log \
			stringtie --merge \
				-p "$THREADS" \
				-G "$gtf" \
				-o "$merged_gtf" \
				"$gtf_list"
	fi
	
	# Re-estimate abundances with merged GTF for better quantification
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
		local final_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/final"
		local final_gtf="$final_dir/${SRR}_${fasta_tag}_ref_guided_final.gtf"
		local final_abundances="$final_dir/${SRR}_${fasta_tag}_ref_guided_final_abundances.tsv"
		local ballgown_dir="$final_dir/ballgown"
		
		mkdir -p "$final_dir" "$ballgown_dir"
		
		# Skip if BAM was deleted and final quantification already exists
		if [[ ! -f "$bam" && -f "$final_gtf" ]]; then
			log_info "[STRINGTIE QUANT] Final quantification for $SRR (Ref-Guided) already exists. Skipping re-estimation."
			continue
		fi
		
		if [[ -f "$bam" ]]; then
			log_step "Re-estimating abundances for $SRR with merged GTF (Ref-Guided)"
			run_with_space_time_log \
				stringtie \
					-p "$THREADS" \
					-e -B \
					-G "$merged_gtf" \
					-A "$final_abundances" \
					-o "$final_gtf" \
					"$bam"
			
			# Cleanup BAM files after final quantification if specified
			if [[ "$keep_bam_global" != "y" ]]; then
				log_info "[CLEANUP] Deleting BAM file for $SRR (Ref-Guided) after final quantification."
				rm -f "$bam" "${bam}.bai"
			fi
		else
			log_warn "BAM file not found for $SRR. Cannot re-estimate abundances."
		fi
	done
	
	# PREPARE COUNT MATRICES FOR DESEQ2 USING prepDE.py
	log_step "Preparing count matrices for DESeq2 using prepDE.py"
	local deseq2_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/deseq2_input"
	local prepde_sample_list="$deseq2_dir/sample_list.txt"
	local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
	local transcript_count_matrix="$deseq2_dir/transcript_count_matrix.csv"
	
	mkdir -p "$deseq2_dir"
	
	# Create sample list file for prepDE.py
	true > "$prepde_sample_list"
	local samples_found=0
	for SRR in "${rnaseq_list[@]}"; do
		local final_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/final"
		local final_gtf="$final_dir/${SRR}_${fasta_tag}_ref_guided_final.gtf"
		
		if [[ -f "$final_gtf" ]]; then
			echo "$SRR $final_gtf" >> "$prepde_sample_list"
			((samples_found++))
		else
			log_warn "Final GTF not found for $SRR. Skipping from count matrix preparation."
		fi
	done
	
	# Validate we have enough samples
	if [[ $samples_found -lt 2 ]]; then
		log_error "Insufficient samples for count matrix generation (found: $samples_found, need: ≥2)"
		return 1
	fi
	log_info "[PREPDE] Found $samples_found samples for count matrix generation"
	
	# Check if count matrices already exist
	if [[ -f "$gene_count_matrix" && -f "$transcript_count_matrix" ]]; then
		log_info "[PREPDE] Count matrices for $fasta_tag (Ref-Guided) already exist. Skipping prepDE.py."
	else
		# Auto-detect read length from first available trimmed FASTQ file using improved function
		local read_length=150  # Default fallback
		local read_length_detected=false
		for SRR in "${rnaseq_list[@]}"; do
			local trim_dir="$TRIM_DIR_ROOT/$SRR"
			# Try multiple patterns for both paired-end and single-end trimmed files
			local r1=$(find "$trim_dir" -type f \( \
				-name "*_1_val_1.fq.gz" -o -name "*_1_val_1.fq" \
				-o -name "*_val_1.fq.gz" -o -name "*_val_1.fq" \
				-o -name "*_1.fq.gz" -o -name "*_1.fq" \
				-o -name "*_trimmed.fq.gz" -o -name "*_trimmed.fq" \
			\) 2>/dev/null | head -n1)
			
			if [[ -f "$r1" ]]; then
				# Attempt to detect read length
				read_length=$(detect_read_length "$r1" 150)
				if [[ $? -eq 0 && "$read_length" -gt 0 ]]; then
					read_length_detected=true
					log_info "[PREPDE] Auto-detected read length: $read_length bp (from $SRR: $(basename "$r1"))"
					break
				else
					log_warn "[PREPDE] Could not detect valid read length from $(basename "$r1"), trying next sample..."
				fi
			fi
		done
		
		if [[ "$read_length_detected" == "false" ]]; then
			log_warn "[PREPDE] Could not auto-detect read length from any sample. Using default: $read_length bp"
			log_warn "[PREPDE] For better accuracy, verify read length matches your data"
			log_warn "[PREPDE] You can manually specify read length by modifying the pipeline or prepDE.py command"
		fi
		
		# Run prepDE.py to generate count matrices
		if command -v prepDE.py >/dev/null 2>&1; then
			log_info "[PREPDE] Running prepDE.py to generate count matrices (read length: $read_length bp)..."
			run_with_space_time_log \
				prepDE.py \
					-i "$prepde_sample_list" \
					-g "$gene_count_matrix" \
					-t "$transcript_count_matrix" \
					-l "$read_length"
		elif command -v python >/dev/null 2>&1 && python -c "import prepDE" 2>/dev/null; then
			log_info "[PREPDE] Running prepDE.py via python to generate count matrices (read length: $read_length bp)..."
			run_with_space_time_log \
				python -m prepDE \
					-i "$prepde_sample_list" \
					-g "$gene_count_matrix" \
					-t "$transcript_count_matrix" \
					-l "$read_length"
		else
			log_error "prepDE.py is required but not found in PATH"
			log_error "Install from: https://github.com/gpertea/stringtie/blob/master/prepDE.py"
			log_error "Or download: wget https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py && chmod +x prepDE.py"
			log_error "Then add to PATH or place in working directory"
			return 1
		fi
		
		# Verify prepDE.py generated output files
		if [[ ! -f "$gene_count_matrix" || ! -f "$transcript_count_matrix" ]]; then
			log_error "prepDE.py failed to generate count matrices"
			log_error "Check sample_list file: $prepde_sample_list"
			return 1
		fi
	fi
	
	# Create sample metadata file for DESeq2 using improved function
	local sample_metadata="$deseq2_dir/sample_metadata.csv"
	if [[ ! -f "$sample_metadata" ]]; then
		create_sample_metadata "$sample_metadata" rnaseq_list[@]
	fi
	
	# Validate count matrix quality using validation function
	if [[ -f "$gene_count_matrix" ]]; then
		validate_count_matrix "$gene_count_matrix" "gene" 2
	fi
	
	if [[ -f "$transcript_count_matrix" ]]; then
		validate_count_matrix "$transcript_count_matrix" "transcript" 2
	fi
	
	# POST-PROCESSING - GENERATE VISUALIZATIONS
	log_step "Post-processing: Generating visualizations (heatmaps and bar graphs)"
	
	local method1_dir="4_POST_PROCESSING/4a_Method_1_HISAT2_Ref_Guided"
	local heatmap_basic_script="$method1_dir/b_modules_for_Method_1/b_make_heatmap_of_matrices.R"
	local heatmap_cv_script="$method1_dir/b_modules_for_Method_1/c_make_heatmap_with_CV_of_matrices.R"
	local bargraph_script="$method1_dir/b_modules_for_Method_1/d_make_BarGraph_of_matrices.R"
	local gene_groups_dir="$method1_dir/a_gene_groups_input_list"
	
	mkdir -p "$method1_dir/7_Heatmap_Outputs/I_Basic_Heatmap"
	mkdir -p "$method1_dir/7_Heatmap_Outputs/II_Heatmap_with_CV"
	mkdir -p "$method1_dir/7_Heatmap_Outputs/III_Bar_Graphs"
	
	if [[ -f "$heatmap_basic_script" && -f "$sample_metadata" && -f "$gene_count_matrix" ]]; then
		log_info "[VISUALIZATION] Checking for gene group files..."
		
		local gene_group_files=()
		if [[ -d "$gene_groups_dir" ]]; then
			while IFS= read -r -d '' file; do
				gene_group_files+=("$file")
			done < <(find "$gene_groups_dir" -name "*.txt" -type f -print0 2>/dev/null)
		fi
		
		if [[ ${#gene_group_files[@]} -eq 0 ]]; then
			log_warn "[VISUALIZATION] No gene group files found in $gene_groups_dir"
		else
			log_info "[VISUALIZATION] Found ${#gene_group_files[@]} gene group file(s)"
			
			for gene_group_file in "${gene_group_files[@]}"; do
				local group_name=$(basename "$gene_group_file" .txt)
				log_info "[VISUALIZATION] Processing: $group_name"
				
				[[ -f "$heatmap_basic_script" ]] && \
					Rscript --vanilla "$heatmap_basic_script" -d "$method1_dir" -f "$fasta_tag" -c "counts" -g "$gene_group_file" 2>&1 | tee -a "$LOG_FILE"
				
				[[ -f "$heatmap_cv_script" ]] && \
					Rscript --vanilla "$heatmap_cv_script" -d "$method1_dir" -f "$fasta_tag" -c "counts" -g "$gene_group_file" 2>&1 | tee -a "$LOG_FILE"
				
				[[ -f "$bargraph_script" ]] && \
					Rscript --vanilla "$bargraph_script" -d "$method1_dir" -f "$fasta_tag" -c "counts" -g "$gene_group_file" 2>&1 | tee -a "$LOG_FILE"
			done
		fi
	fi
	
	log_step "HISAT2 reference-guided pipeline completed for $fasta_tag"
	log_info "Merged GTF: $merged_gtf"
	log_info "Final quantifications in: $STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/*/final/"
	log_info "DESeq2 input files:"
	log_info "  - Gene count matrix: $gene_count_matrix"
	log_info "  - Transcript count matrix: $transcript_count_matrix" 
	log_info "  - Sample metadata: $sample_metadata"
}

# Combined pipeline: Build HISAT2 index, align reads, assemble, merge, and quantify transcripts
hisat2_de_novo_pipeline() {
	local fasta="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)
				fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1")
					shift
				done
				;;
			*)
				log_error "Unknown option: $1"
				return 1;;
		esac
	done
	
	if [[ -z "$fasta" ]]; then
		log_error "No FASTA file specified. Use --FASTA <fasta_file>."
		return 1
	fi
	if [[ ! -f "$fasta" ]]; then
		log_error "FASTA file '$fasta' not found."
		return 1
	fi
	if [[ ${#rnaseq_list[@]} -eq 0 ]]; then
		rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	fi

	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"

	# BUILD HISAT2 INDEX
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "[HISAT2 INDEX] De novo index for $fasta_base already exists. Skipping build."
	else
		log_step "Building HISAT2 de novo index from $fasta"
		log_file_size "$fasta" "Input FASTA for HISAT2 de novo index - $fasta_tag"
		run_with_space_time_log --input "$fasta" --output "$HISAT2_DE_NOVO_INDEX_DIR" hisat2-build -p "${THREADS}" "$fasta" "$index_prefix"
		log_file_size "$HISAT2_DE_NOVO_INDEX_DIR" "HISAT2 de novo index output - $fasta_tag"
	fi

	# ALIGNMENT, SORTING, AND STRINGTIE ASSEMBLY for each SRR
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="" trimmed2=""
		mkdir -p "$HISAT2_DIR"
		
		# Find trimmed FASTQ files
		if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
		elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1 && compgen -G "$TrimGalore_DIR/${SRR}*val_2.*" >/dev/null 2>&1; then
			local files1=( "$TrimGalore_DIR"/${SRR}*val_1.* )
			local files2=( "$TrimGalore_DIR"/${SRR}*val_2.* )
			trimmed1="${files1[0]}"
			trimmed2="${files2[0]}"
		fi
		
		if [[ -z "$trimmed1" || -z "$trimmed2" ]]; then
			log_warn "Trimmed FASTQ for $SRR not found in $TrimGalore_DIR; skipping."
			continue
		fi

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped.sam"
		
		# Align and sort if BAM doesn't exist
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "[HISAT2 ALIGN] BAM file for $SRR and $fasta_tag (de novo) already exists. Skipping alignment."
		else
			log_step "Aligning $fasta_tag with $SRR using HISAT2 De Novo Approach"
			# Check if paired-end or single-end reads
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				# Paired-end alignment
				log_input_output_size "$trimmed1" "$HISAT2_DIR" "HISAT2 de novo alignment - $SRR"
				run_with_space_time_log --input "$TrimGalore_DIR" --output "$HISAT2_DIR" \
					hisat2 -p "${THREADS}" --dta \
						-x "$index_prefix" \
						-1 "$trimmed1" \
						-2 "$trimmed2" \
						-S "$sam"
			else
				# Single-end alignment
				run_with_space_time_log \
					hisat2 -p "${THREADS}" --dta \
						-x "$index_prefix" \
						-U "$trimmed1" \
						-S "$sam"
			fi
			
			log_info "[SAMTOOLS] Converting SAM to sorted BAM for $fasta_tag with $SRR (de novo)..."
			log_file_size "$sam" "SAM file before sorting (de novo) - $SRR"
			run_with_space_time_log --input "$sam" --output "$bam" samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			log_file_size "$bam" "Sorted BAM file (de novo) - $SRR"
			run_with_space_time_log samtools index -@ "${THREADS}" "$bam"
			log_info "[CLEANUP] Deleting SAM file to save disk space"
			rm -f "$sam"
			log_info "[HISAT2 ALIGN] Done aligning $fasta_tag with $SRR (de novo)."
		fi

		# StringTie assembly
		local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
		local out_gene_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"
		mkdir -p "$out_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "[STRINGTIE] De novo assembly for $fasta_tag/$SRR already exists. Skipping."
		else
			log_step "Assembling transcripts for $fasta_tag with $SRR (de novo)"
			log_input_output_size "$bam" "$out_dir" "StringTie de novo assembly - $SRR"
			run_with_space_time_log --input "$bam" --output "$out_dir" \
				stringtie -p "$THREADS" "$bam" \
					-o "$out_gtf" \
					-A "$out_gene_abundances_tsv"
			log_file_size "$out_gtf" "StringTie de novo output GTF - $SRR"
		fi
		
		# Cleanup BAM and SAM files after quantification (unless explicitly kept)
		if [[ "$keep_bam_global" != "y" ]]; then
			log_info "[CLEANUP] Deleting BAM/SAM files to save disk space"
			rm -f "$bam" "${bam}.bai" "$sam"
		else
			log_info "[KEEP_BAM] Preserving BAM file for further analysis"
		fi
		log_info "[STRINGTIE] Done processing $fasta_tag with $SRR (de novo)."
		log_info "--------------------------------------------------"
	done
}

# Trinity de novo with tissue-specific assembly option
trinity_tissue_specific_pipeline() {
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
			log_warn "Metadata load failed - running pooled assembly"
			trinity_de_novo_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			return $?
		}
	else
		log_warn "No metadata - running pooled assembly"
		trinity_de_novo_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	# Group samples by tissue
	declare -A tissue_samples
	for SRR in "${rnaseq_list[@]}"; do
		local tissue="${sample_metadata[${SRR}_condition]:-unknown}"
		tissue_samples["$tissue"]+="$SRR "
	done
	
	local tissue_count=${#tissue_samples[@]}
	log_info "[TRINITY TISSUE] Found $tissue_count tissue types"
	
	if [[ $tissue_count -lt 2 ]]; then
		log_warn "Single tissue detected - using pooled assembly"
		trinity_de_novo_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	log_info "[TRINITY TISSUE] Running tissue-specific assemblies for better isoform detection"
	
	# Run Trinity per tissue
	for tissue in "${!tissue_samples[@]}"; do
		local tissue_srrs=(${tissue_samples[$tissue]})
		log_step "Trinity assembly for tissue: $tissue (${#tissue_srrs[@]} samples)"
		
		trinity_de_novo_alignment_pipeline \
			--FASTA "$fasta" \
			--RNASEQ_LIST "${tissue_srrs[@]}" \
			--TISSUE_TAG "$tissue"
	done
}

# Trinity de novo assembly with Salmon quantification and tximport for DESeq2
# 
# STATISTICAL CONSIDERATIONS FOR TISSUE-SPECIFIC EXPRESSION PROFILING:
# 1. Trinity pools all samples into ONE assembly - may miss tissue-specific isoforms
# 2. In-silico normalization DISABLED to preserve biological variation
# 3. Salmon quasi-mapping for fast, accurate quantification
# 4. Tximport properly handles transcript-level uncertainty for gene-level inference
# 5. Strand-specific protocols supported via TRINITY_STRAND_SPECIFIC env var
#
# WORKFLOW:
# 1. Trinity de novo assembly from trimmed reads
# 2. Salmon indexing and quantification (pseudo-alignment)
# 3. Tximport for proper statistical aggregation to gene level
# 4. DESeq2-ready outputs with tx2gene mapping
#
# OPTIMIZATIONS:
# - Salmon index built once and reused
# - Memory configurable via TRINITY_MAX_MEMORY (default: 50G)
# - Min contig length filter reduces spurious transcripts
# - Full cleanup enabled to save disk space
#
trinity_de_novo_alignment_pipeline() {
	local fasta="" rnaseq_list=() tissue_tag=""
	
	# Check for GNU parallel and warn if not available
	if ! command -v parallel >/dev/null 2>&1; then
		log_warn "[PERFORMANCE] GNU parallel not found - Salmon quantification will run sequentially"
		log_warn "[PERFORMANCE] Install with: sudo apt-get install parallel (Ubuntu/Debian)"
		log_warn "[PERFORMANCE] Or: brew install parallel (macOS)"
		log_warn "[PERFORMANCE] Parallel processing can speed up quantification by 2-4x"
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
	
	if [[ -z "$fasta" ]]; then
		log_error "No FASTA file specified. Use --FASTA <fasta_file>."
		return 1
	fi
	if [[ ! -f "$fasta" ]]; then
		log_error "FASTA file '$fasta' not found."
		return 1
	fi
	if [[ ${#rnaseq_list[@]} -eq 0 ]]; then
		rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	fi
	if [[ ${#rnaseq_list[@]} -eq 0 ]]; then
		log_error "No RNA-seq samples provided for Trinity pipeline. Cannot proceed."
		return 1
	fi

	local fasta_base fasta_tag trinity_out_dir trinity_fasta
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	
	# Append tissue tag if provided
	if [[ -n "$tissue_tag" ]]; then
		trinity_out_dir="$TRINITY_DE_NOVO_ROOT/${fasta_tag}_${tissue_tag}_trinity_out_dir"
		log_info "[TRINITY] Tissue-specific assembly for: $tissue_tag (${#rnaseq_list[@]} samples)"
	else
		trinity_out_dir="$TRINITY_DE_NOVO_ROOT/${fasta_tag}_trinity_out_dir"
		log_info "[TRINITY] Pooled assembly: ${#rnaseq_list[@]} samples"
		log_warn "⚠️  Pooled assembly may miss tissue-specific isoforms"
		log_warn "   Use trinity_tissue_specific_pipeline with --METADATA for better detection"
	fi
	
	trinity_fasta="$trinity_out_dir/Trinity.fasta"
	
	# Trinity parameters for tissue-specific expression profiling
	local trinity_max_memory="${TRINITY_MAX_MEMORY:-50G}"
	local trinity_min_contig_length="${TRINITY_MIN_CONTIG_LENGTH:-200}"
	local trinity_min_kmer_cov="${TRINITY_MIN_KMER_COV:-2}"
	# For strand-specific libraries (set based on your protocol):
	# --SS_lib_type RF (dUTP/NSR/NNSR), FR (Ligation), F/R (single-end)
	local trinity_strand="${TRINITY_STRAND_SPECIFIC:-}"  # Leave empty for unstranded
	
	# Trinity parallelization settings
	local trinity_bfly_cpu=$((THREADS / 4))
	[[ $trinity_bfly_cpu -lt 1 ]] && trinity_bfly_cpu=1
	log_info "[TRINITY] CPU allocation: Total=$THREADS, Butterfly=$trinity_bfly_cpu"

	# STEP 1: DE NOVO TRANSCRIPTOME ASSEMBLY WITH TRINITY
	log_step "Trinity de novo transcriptome assembly for $fasta_tag"
	
	if [[ -f "$trinity_fasta" || -f "${trinity_fasta}.gz" ]]; then
		log_info "[TRINITY] Assembly for $fasta_tag already exists. Skipping."
		# Set trinity_fasta to the existing file
		[[ -f "${trinity_fasta}.gz" ]] && trinity_fasta="${trinity_fasta}.gz"
	else
		log_info "[TRINITY] Starting de novo assembly for $fasta_tag..."
		mkdir -p "$trinity_out_dir"
		
		log_info "[TRINITY] Collecting trimmed FASTQ files from: $TRIM_DIR_ROOT"
		
		# Collect all trimmed FASTQ files for Trinity
		local left_reads=() right_reads=() single_reads=()
		local paired_count=0 single_count=0
		
		log_info "[TRINITY] Processing ${#rnaseq_list[@]} SRR samples for Trinity input..."
		
		for SRR in "${rnaseq_list[@]}"; do
			local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
			local trimmed1="" trimmed2=""
			
			# Find trimmed FASTQ files
			if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
				trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
				trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
			elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
				trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
				trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
			elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1 && compgen -G "$TrimGalore_DIR/${SRR}*val_2.*" >/dev/null 2>&1; then
				local files1=( "$TrimGalore_DIR"/${SRR}*val_1.* )
				local files2=( "$TrimGalore_DIR"/${SRR}*val_2.* )
				trimmed1="${files1[0]}"
				trimmed2="${files2[0]}"
			elif compgen -G "$TrimGalore_DIR/${SRR}*trimmed.*" >/dev/null 2>&1; then
				# Single-end trimmed files (common pattern)
				local files=( "$TrimGalore_DIR"/${SRR}*trimmed.* )
				trimmed1="${files[0]}"
			elif compgen -G "$TrimGalore_DIR/${SRR}*.fq*" >/dev/null 2>&1; then
				# Generic single-end pattern
				local files=( "$TrimGalore_DIR"/${SRR}*.fq* )
				trimmed1="${files[0]}"
			fi
			
			# Determine if paired-end or single-end
			if [[ -n "$trimmed1" && -n "$trimmed2" && -f "$trimmed1" && -f "$trimmed2" ]]; then
				# Paired-end reads
				left_reads+=("$trimmed1")
				right_reads+=("$trimmed2")
				paired_count=$((paired_count + 1))
				log_info "[TRINITY] Added paired-end reads for $SRR: $trimmed1, $trimmed2"
			elif [[ -n "$trimmed1" && -f "$trimmed1" ]]; then
				# Single-end reads
				single_reads+=("$trimmed1")
				single_count=$((single_count + 1))
				log_info "[TRINITY] Added single-end reads for $SRR: $trimmed1"
			else
				log_warn "Trimmed FASTQ for $SRR not found; skipping from Trinity input."
			fi
		done
		
		log_info "[TRINITY] Read collection completed: Paired=$paired_count, Single=$single_count"
		log_info "Left reads array size: ${#left_reads[@]}"
		log_info "Right reads array size: ${#right_reads[@]}"
		log_info "Single reads array size: ${#single_reads[@]}"
		
		if [[ ${#left_reads[@]} -eq 0 && ${#single_reads[@]} -eq 0 ]]; then
			log_error "No trimmed FASTQ files found for Trinity assembly."
			log_error "Searched in: $TRIM_DIR_ROOT/<SRR>/"
			return 1
		fi
		
		log_info "Found $paired_count paired-end samples and $single_count single-end samples"
		
		# Run Trinity based on read type
		if [[ ${#left_reads[@]} -gt 0 && ${#right_reads[@]} -gt 0 ]]; then
			# Mixed or paired-end reads
			local left_files=$(IFS=,; echo "${left_reads[*]}")
			local right_files=$(IFS=,; echo "${right_reads[*]}")
			
			if [[ ${#single_reads[@]} -gt 0 ]]; then
				# Mixed reads: both paired and single
				local single_files=$(IFS=,; echo "${single_reads[*]}")
				log_info "[TRINITY] Running Trinity with mixed paired/single-end reads"
				
				# Build Trinity command with optional strand parameter
				local trinity_cmd=(Trinity --seqType fq 
					--left "$left_files" 
					--right "$right_files" 
					--single "$single_files"
					--CPU "$THREADS"
					--bflyCalculateCPU
					--bflyCPU "$trinity_bfly_cpu"
					--max_memory "$trinity_max_memory"
					--output "$trinity_out_dir"
					--min_contig_length "$trinity_min_contig_length"
					--min_kmer_cov "$trinity_min_kmer_cov"
					--full_cleanup
					--no_version_check)
				
				# Add strand-specific parameter if specified
				[[ -n "$trinity_strand" ]] && trinity_cmd+=(--SS_lib_type "$trinity_strand")
				
				# Note: --normalize_reads removed to preserve biological variation for tissue comparison
				log_warn "⚠️  In-silico normalization disabled - preserves tissue-specific expression differences"
				
				log_file_size "$TRIM_DIR_ROOT" "Input trimmed reads directory for Trinity"
				run_with_space_time_log --input "$TRIM_DIR_ROOT" --output "$trinity_out_dir" "${trinity_cmd[@]}"
				log_file_size "$trinity_out_dir" "Trinity de novo assembly output"
			else
				# Only paired-end reads
				log_info "[TRINITY] Running Trinity with paired-end reads only"
				
				local trinity_cmd=(Trinity --seqType fq 
					--left "$left_files" 
					--right "$right_files"
					--CPU "$THREADS"
					--bflyCalculateCPU
					--bflyCPU "$trinity_bfly_cpu"
					--max_memory "$trinity_max_memory"
					--output "$trinity_out_dir"
					--min_contig_length "$trinity_min_contig_length"
					--min_kmer_cov "$trinity_min_kmer_cov"
					--full_cleanup
					--no_version_check)
				
				[[ -n "$trinity_strand" ]] && trinity_cmd+=(--SS_lib_type "$trinity_strand")
				
				log_file_size "$TRIM_DIR_ROOT" "Input trimmed reads directory for Trinity"
				run_with_space_time_log --input "$TRIM_DIR_ROOT" --output "$trinity_out_dir" "${trinity_cmd[@]}"
				log_file_size "$trinity_out_dir" "Trinity de novo assembly output"
			fi
		elif [[ ${#single_reads[@]} -gt 0 ]]; then
			# Only single-end reads
			local single_files=$(IFS=,; echo "${single_reads[*]}")
			log_info "[TRINITY] Running Trinity with single-end reads only"
			
			local trinity_cmd=(Trinity --seqType fq 
				--single "$single_files"
				--CPU "$THREADS"
				--bflyCalculateCPU
				--bflyCPU "$trinity_bfly_cpu"
				--max_memory "$trinity_max_memory"
				--output "$trinity_out_dir"
				--min_contig_length "$trinity_min_contig_length"
				--min_kmer_cov "$trinity_min_kmer_cov"
				--full_cleanup
				--no_version_check)
			
			[[ -n "$trinity_strand" ]] && trinity_cmd+=(--SS_lib_type "$trinity_strand")
			
			run_with_space_time_log "${trinity_cmd[@]}"
		else
			log_error "No valid reads found for Trinity assembly."
			return 1
		fi
		
		# Check if Trinity assembly was successful
		if [[ ! -f "$trinity_fasta" ]]; then
			log_error "Trinity assembly failed - Trinity.fasta not found at: $trinity_fasta"
			log_error "Check Trinity log files in: $trinity_out_dir"
			return 1
		fi
		
		log_info "[TRINITY] Assembly completed successfully: $trinity_fasta"
		
		# Generate assembly statistics and quality metrics
		log_info "[TRINITY] Generating assembly statistics and quality metrics..."
		local stats_file="$trinity_out_dir/assembly_stats.txt"
		local quality_report="$trinity_out_dir/assembly_quality_report.txt"
		
		# Run TrinityStats.pl if available
		if command -v TrinityStats.pl >/dev/null 2>&1; then
			TrinityStats.pl "$trinity_fasta" > "$stats_file" 2>&1 || \
				log_warn "Failed to generate detailed Trinity statistics"
			
			# Check assembly quality metrics
			if [[ -f "$stats_file" ]]; then
				local n50=$(grep -oP 'N50:\s*\K[0-9]+' "$stats_file" 2>/dev/null || echo "0")
				local median_contig=$(grep -oP 'Median contig length:\s*\K[0-9]+' "$stats_file" 2>/dev/null || echo "0")
				local total_assembled=$(grep -oP 'Total assembled bases:\s*\K[0-9]+' "$stats_file" 2>/dev/null || echo "0")
				
				# Quality assessment
				{
					echo "═══════════════════════════════════════════════════════════"
					echo "TRINITY ASSEMBLY QUALITY REPORT"
					echo "Assembly: $fasta_tag"
					echo "Date: $(date)"
					echo "═══════════════════════════════════════════════════════════"
					echo ""
					
					if [[ ${n50:-0} -lt 500 ]]; then
						echo "⚠️  WARNING: Low N50 ($n50 bp)"
						echo "   - Assembly may be highly fragmented"
						echo "   - Consider increasing read depth or adjusting assembly parameters"
						echo "   - Recommended N50: >1000 bp for most transcriptomes"
					elif [[ ${n50:-0} -lt 1000 ]]; then
						echo "✓ Acceptable N50: $n50 bp"
						echo "   - Assembly quality is moderate"
					else
						echo "✓✓ Good N50: $n50 bp"
						echo "   - Assembly quality is high"
					fi
					echo ""
					
					if [[ ${median_contig:-0} -lt 300 ]]; then
						echo "⚠️  WARNING: Low median contig length ($median_contig bp)"
						echo "   - Many short contigs detected"
					else
						echo "✓ Median contig length: $median_contig bp"
					fi
					echo ""
					
					if [[ ${total_assembled:-0} -lt 10000000 ]]; then
						echo "⚠️  WARNING: Low total assembled bases ($total_assembled bp)"
						echo "   - Transcriptome may be incomplete"
						echo "   - Expected: 10-100 Mbp for most eukaryotic transcriptomes"
					else
						echo "✓ Total assembled bases: $total_assembled bp"
					fi
					echo ""
					
					echo "See full statistics in: $stats_file"
					echo "═══════════════════════════════════════════════════════════"
				} > "$quality_report"
				
				cat "$quality_report"
				log_info "[TRINITY] Assembly N50: $n50 bp"
			fi
		fi
		
		# Basic statistics using standard tools
		local num_transcripts=$(grep -c "^>" "$trinity_fasta" || echo "0")
		local total_length=$(grep -v "^>" "$trinity_fasta" | tr -d '\n' | wc -c)
		local avg_length=0
		[[ $num_transcripts -gt 0 ]] && avg_length=$((total_length / num_transcripts))
		
		log_info "[TRINITY] Assembly: $num_transcripts transcripts, avg $avg_length bp"
		
		# Quality validation
		[[ $num_transcripts -lt 1000 ]] && \
			log_warn "⚠️  Low transcript count ($num_transcripts) - check read depth/quality"
		
		# Save basic stats
		{
			echo "Trinity Assembly Statistics for $fasta_tag"
			echo "==========================================="
			echo "Assembly file: $trinity_fasta"
			echo "Number of transcripts: $num_transcripts"
			echo "Total sequence length: $total_length bp"
			echo "Average transcript length: $avg_length bp"
			echo "Generated: $(date)"
		} > "$trinity_out_dir/assembly_summary.txt"
		
		# Clean up Trinity intermediate files to save disk space
		log_info "[TRINITY CLEANUP] Removing large intermediate assembly files"
		log_info "[TRINITY CLEANUP] Removing chrysalis directories..."
		find "$trinity_out_dir" -name "chrysalis" -type d -exec rm -rf {} + 2>/dev/null || true
		log_info "[TRINITY CLEANUP] Removing jellyfish kmer files..."
		find "$trinity_out_dir" -name "jellyfish.kmers*" -delete 2>/dev/null || true
		log_info "[TRINITY CLEANUP] Removing inchworm kmer files..."
		find "$trinity_out_dir" -name "inchworm.kmer*" -delete 2>/dev/null || true
		log_info "[TRINITY CLEANUP] Removing recursive_trinity.cmds files..."
		find "$trinity_out_dir" -name "recursive_trinity.cmds*" -delete 2>/dev/null || true
		log_info "[TRINITY CLEANUP] Trinity intermediate file cleanup complete"
	fi

	# Use Trinity assembly directly (keep uncompressed during pipeline for efficiency)
	local trinity_ref="$trinity_fasta"
	
	# Handle if assembly is already compressed from previous run
	if [[ "$trinity_fasta" == *.gz ]]; then
		trinity_ref="${trinity_fasta%.gz}"
		if [[ ! -f "$trinity_ref" ]]; then
			log_info "[TRINITY] Decompressing assembly for downstream analysis..."
			if ! gunzip -k "$trinity_fasta"; then
				log_error "Failed to decompress Trinity assembly"
				return 1
			fi
		fi
	fi

	# STEP 2: BUILD SALMON INDEX AND QUANTIFY TRANSCRIPTS
	# Use tissue-aware tags to avoid collisions when running per tissue
	local tag_suffix=""
	[[ -n "${tissue_tag:-}" ]] && tag_suffix="_${tissue_tag}"
	local salmon_idx="$TRINITY_DE_NOVO_ROOT/${fasta_tag}${tag_suffix}_salmon_index"
	local quant_root="$TRINITY_DE_NOVO_ROOT/${fasta_tag}${tag_suffix}_salmon_quant"
	mkdir -p "$quant_root"
	
	# Build Salmon index once
	if [[ -f "$salmon_idx/versionInfo.json" ]]; then
		log_info "[SALMON INDEX] Trinity Salmon index exists - skipping build"
	else
		log_step "Building Salmon index for Trinity assembly"
		log_file_size "$trinity_ref" "Trinity assembly for Salmon indexing"
		run_with_space_time_log --input "$trinity_ref" --output "$salmon_idx" salmon index -t "$trinity_ref" -i "$salmon_idx" -k 31 --threads "$THREADS"
		log_file_size "$salmon_idx" "Salmon index for Trinity assembly"
	fi
	
	# Quantify each sample with Salmon using GNU parallel
	log_info "[SALMON] Starting parallel quantification for ${#rnaseq_list[@]} samples"
	
	# Export required variables for parallel execution
	export salmon_idx quant_root TRIM_DIR_ROOT THREADS trinity_strand LOG_FILE
	export -f log_info log_warn log_error log_step run_with_space_time_log
	
	# Create parallel quantification function
	quantify_sample() {
		set +e  # Disable exit on error for this function
		local SRR="$1"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="" trimmed2=""
		
		# Find trimmed FASTQ files
		if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
		elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1 && compgen -G "$TrimGalore_DIR/${SRR}*val_2.*" >/dev/null 2>&1; then
			local files1=( "$TrimGalore_DIR"/${SRR}*val_1.* )
			local files2=( "$TrimGalore_DIR"/${SRR}*val_2.* )
			trimmed1="${files1[0]}"
			trimmed2="${files2[0]}"
		elif compgen -G "$TrimGalore_DIR/${SRR}*trimmed.*" >/dev/null 2>&1; then
			local files=( "$TrimGalore_DIR"/${SRR}*trimmed.* )
			trimmed1="${files[0]}"
		elif compgen -G "$TrimGalore_DIR/${SRR}*.fq*" >/dev/null 2>&1; then
			local files=( "$TrimGalore_DIR"/${SRR}*.fq* )
			trimmed1="${files[0]}"
		fi
		
		if [[ -z "$trimmed1" ]]; then
			echo "[WARN] No trimmed reads for $SRR - skipping"
			set -e
			return 0
		fi
		
		local quant_dir="$quant_root/$SRR"
		if [[ -f "$quant_dir/quant.sf" ]]; then
			echo "[SALMON] Quantification for $SRR exists - skipping"
			set -e
			return 0
		fi
		
		echo "[SALMON] Quantifying: $SRR"
		mkdir -p "$quant_dir"
		
		# Calculate threads per sample (divide total threads by number of parallel jobs)
		local threads_per_sample=$((THREADS / 2))
		[[ $threads_per_sample -lt 1 ]] && threads_per_sample=1
		
		# Build Salmon command
		local salmon_cmd="salmon quant -i $salmon_idx -o $quant_dir -p $threads_per_sample --validateMappings"
		
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			salmon_cmd="$salmon_cmd -l A -1 $trimmed1 -2 $trimmed2"
		else
			salmon_cmd="$salmon_cmd -l A -r $trimmed1"
		fi
		
		# Add strand-specific library type if specified
		if [[ -n "$trinity_strand" ]]; then
			case "$trinity_strand" in
				RF) salmon_cmd="${salmon_cmd// -l A/ -l ISR}" ;;
				FR) salmon_cmd="${salmon_cmd// -l A/ -l ISF}" ;;
				F) salmon_cmd="${salmon_cmd// -l A/ -l SF}" ;;
				R) salmon_cmd="${salmon_cmd// -l A/ -l SR}" ;;
			esac
		fi
		
		# Run Salmon and capture exit code
		eval $salmon_cmd 2>&1
		local salmon_exit=$?
		
		# Validate output
		if [[ -f "$quant_dir/quant.sf" && $salmon_exit -eq 0 ]]; then
			local num_transcripts=$(tail -n +2 "$quant_dir/quant.sf" | wc -l)
			local expressed=$(awk 'NR>1 && $5>0' "$quant_dir/quant.sf" | wc -l)
			local total_reads=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_dir/quant.sf")
			echo "[SALMON] $SRR: $expressed/$num_transcripts expressed, $total_reads total counts"
			[[ $total_reads -lt 500000 ]] && echo "[WARN] Low counts: $total_reads (expected >1M)"
			set -e
			return 0
		else
			echo "[ERROR] Salmon quantification failed for $SRR (exit code: $salmon_exit)"
			set -e
			return 1
		fi
	}
	export -f quantify_sample
	
	# Run parallel quantification (2 samples at a time to balance I/O and compute)
	if command -v parallel >/dev/null 2>&1; then
		printf "%s\n" "${rnaseq_list[@]}" | parallel -j 2 --line-buffer --halt soon,fail=1 quantify_sample {} 2>&1 | tee -a "$LOG_FILE" || {
			local parallel_exit=$?
			log_error "[SALMON] Parallel quantification failed with exit code: $parallel_exit"
			return 1
		}
		log_info "[SALMON] Parallel quantification completed"
	else
		log_warn "[SALMON] GNU parallel not found, falling back to sequential processing"
		local failed_samples=0
		for SRR in "${rnaseq_list[@]}"; do
			if ! quantify_sample "$SRR" 2>&1 | tee -a "$LOG_FILE"; then
				log_error "[SALMON] Quantification failed for $SRR"
				((failed_samples++))
			fi
		done
		if [[ $failed_samples -gt 0 ]]; then
			log_error "[SALMON] $failed_samples sample(s) failed during quantification"
			return 1
		fi
	fi


	# STEP 3: PREPARE TXIMPORT FILES FOR DESEQ2
	log_step "Preparing tximport input for DESeq2 (Trinity + Salmon)"
	
	local method3_root="${TRINITY_DE_NOVO_ROOT%/*}"
	local matrix_dir="$method3_root/6_matrices_from_stringtie"
	mkdir -p "$matrix_dir"
	
	local sample_metadata="$matrix_dir/sample_info.tsv"
	# Include tissue tag in tx2gene filename to avoid overwriting when running multiple tissues
	local tx2gene_file="$matrix_dir/tx2gene_${fasta_tag}${tag_suffix}.tsv"
	
	# Create tx2gene mapping from Trinity FASTA (TRINITY_*: gene == transcript ID without final _i[0-9]+)
	# Validate tissue_tag is set properly to avoid filename issues
	[[ -z "${tissue_tag+x}" ]] && tissue_tag=""
	if [[ ! -f "$tx2gene_file" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping"
		awk '/^>/{
			tx=$1; gsub(/^>/, "", tx);
			gene=tx; if (match(gene, /^(.+)_i[0-9]+$/, arr)) { gene=arr[1]; }
			print tx "\t" gene
		}' "$trinity_ref" > "$tx2gene_file"
		log_info "[TXIMPORT] Created tx2gene mapping: $(wc -l < "$tx2gene_file") transcripts"
	fi
	
	# Verify Salmon quantifications
	local quant_count=0
	for SRR in "${rnaseq_list[@]}"; do
		[[ -f "$quant_root/$SRR/quant.sf" ]] && ((quant_count++))
	done
	
	[[ $quant_count -lt 2 ]] && { log_error "Insufficient quantifications: $quant_count (need ≥2)"; return 1; }
	log_info "[TXIMPORT] Found $quant_count samples with Salmon quantifications"
	
	# Create sample metadata file
	if [[ ! -f "$sample_metadata" ]]; then
		create_sample_metadata "$sample_metadata" rnaseq_list[@]
	fi
	
	# Generate tximport R script
	local tximport_script="$matrix_dir/run_tximport_trinity_salmon.R"
	if [[ ! -f "$tximport_script" ]]; then
		cat > "$tximport_script" << 'RSCRIPT'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
quant_dir <- if(length(args) > 0) args[1] else "QUANT_DIR_PLACEHOLDER"
metadata_file <- if(length(args) > 1) args[2] else "METADATA_FILE_PLACEHOLDER"
tx2gene_file <- if(length(args) > 2) args[3] else "TX2GENE_FILE_PLACEHOLDER"
output_dir <- if(length(args) > 3) args[4] else dirname(metadata_file)

cat("Trinity + Salmon tximport for DESeq2\n")
cat("Quant dir:", quant_dir, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Tx2gene:", tx2gene_file, "\n")

# Load metadata
coldata <- read.delim(metadata_file, header=TRUE, stringsAsFactors=FALSE)
samples <- coldata$sample
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples
all(file.exists(files)) || stop("Missing quant.sf files")

# Load tx2gene mapping
tx2gene <- read.delim(tx2gene_file, header=FALSE, col.names=c("TXNAME", "GENEID"))

# Import with tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)

# Save outputs
saveRDS(txi, file.path(output_dir, "tximport_trinity_salmon.rds"))
saveRDS(dds, file.path(output_dir, "deseq2_dataset_trinity.rds"))

# Export count matrices
write.table(txi$counts, file.path(output_dir, "gene_counts_tximport.tsv"), 
            sep="\t", quote=FALSE, col.names=NA)
write.table(txi$abundance, file.path(output_dir, "gene_tpm_tximport.tsv"), 
            sep="\t", quote=FALSE, col.names=NA)

cat("\nTximport completed successfully!\n")
cat("Outputs saved to:", output_dir, "\n")
cat("- tximport_trinity_salmon.rds\n")
cat("- deseq2_dataset_trinity.rds\n")
cat("- gene_counts_tximport.tsv\n")
cat("- gene_tpm_tximport.tsv\n")
RSCRIPT
		
		# Replace placeholders
		sed -i "s|QUANT_DIR_PLACEHOLDER|$quant_root|g" "$tximport_script"
		sed -i "s|METADATA_FILE_PLACEHOLDER|$sample_metadata|g" "$tximport_script"
		sed -i "s|TX2GENE_FILE_PLACEHOLDER|$tx2gene_file|g" "$tximport_script"
		
		chmod +x "$tximport_script"
		log_info "[TXIMPORT] Generated R script: $tximport_script"
	fi
	
	# Run tximport if R is available
	if command -v Rscript >/dev/null 2>&1; then
		log_step "Running tximport to import Salmon quantifications"
		if Rscript "$tximport_script" "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" 2>&1 | tee -a "$LOG_FILE"; then
			log_info "[TXIMPORT] ✓ Successfully imported counts for DESeq2"
			
			# Validate output
			local count_matrix="$matrix_dir/gene_counts_tximport.tsv"
			[[ -f "$count_matrix" ]] && validate_count_matrix "$count_matrix" "gene" 2
		else
			log_warn "[TXIMPORT] ✗ tximport failed - check R dependencies"
			log_warn "Install: BiocManager::install(c('tximport', 'DESeq2'))"
		fi
	else
		log_warn "[TXIMPORT] Rscript not found - run manually: Rscript $tximport_script"
	fi
	
	# STEP 4: POST-PROCESSING - GENERATE VISUALIZATIONS
	log_step "Post-processing: Generating visualizations"
	
	local gene_groups_dir="$method3_root/a_gene_groups_input_list"
	mkdir -p "$method3_root/7_Heatmap_Outputs/"{I_Basic_Heatmap,II_Heatmap_with_CV,III_Bar_Graphs}
	
	# Check if count matrix exists from tximport
	local count_matrix="$matrix_dir/gene_counts_tximport.tsv"
	if [[ ! -f "$count_matrix" ]]; then
		log_warn "[VISUALIZATION] Count matrix not found - run tximport first"
		log_info "[VISUALIZATION] Execute: Rscript $matrix_dir/run_tximport_trinity_salmon.R"
	elif [[ ! -d "$gene_groups_dir" ]] || [[ -z "$(ls -A "$gene_groups_dir"/*.txt 2>/dev/null)" ]]; then
		log_warn "[VISUALIZATION] No gene group files in: $gene_groups_dir"
		log_info "[VISUALIZATION] Create gene group files (one gene ID per line) to generate plots"
	else
		log_info "[VISUALIZATION] Visualizations require custom DESeq2 analysis"
		log_info "[VISUALIZATION] Use the count matrix: $count_matrix"
		log_info "[VISUALIZATION] Use the sample info: $sample_metadata"
		log_info "[VISUALIZATION] Use gene groups from: $gene_groups_dir"
		log_info ""
		log_info "[VISUALIZATION] Run DESeq2 normalization and generate heatmaps in R:"
		log_info "  library(DESeq2); library(ComplexHeatmap)"
		log_info "  dds <- readRDS('$matrix_dir/deseq2_dataset_trinity.rds')"
		log_info "  dds <- DESeq(dds)"
		log_info "  # Then create custom visualizations with normalized counts"
	fi
	
	# FINAL SUMMARY
	log_step "Trinity de novo pipeline completed for $fasta_tag"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	log_info "📁 OUTPUT STRUCTURE:"
	log_info "  Trinity assembly: $trinity_fasta"
	log_info "  Salmon index: $salmon_idx"
	log_info "  Quantifications: $quant_root"
	log_info "  Tximport outputs: $matrix_dir"
	log_info ""
	log_info "📊 DESEQ2 INPUT FILES:"
	log_info "  1. Sample metadata: $sample_metadata"
	log_info "  2. Tx2gene mapping: $tx2gene_file"
	log_info "  3. Tximport RDS: $matrix_dir/tximport_trinity_salmon.rds"
	log_info "  4. DESeq2 dataset: $matrix_dir/deseq2_dataset_trinity.rds"
	log_info "  5. Count matrix: $matrix_dir/gene_counts_tximport.tsv"
	log_info ""
	log_info "🔬 NEXT STEPS FOR TISSUE-SPECIFIC EXPRESSION PROFILING:"
	log_info "  1. Edit sample metadata with tissue information: $sample_metadata"
	log_info "  2. Create gene group files: $gene_groups_dir (one gene ID per line)"
	log_info "  3. Run tximport if not done: Rscript $matrix_dir/run_tximport_trinity_salmon.R"
	log_info "  4. Perform differential expression in R:"
	log_info "     dds <- readRDS('$matrix_dir/deseq2_dataset_trinity.rds')"
	log_info "     dds <- DESeq(dds)"
	log_info "     results(dds, contrast=c('condition', 'treatment', 'control'))"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	
	# Compress Trinity assembly to save space
	if [[ -f "$trinity_ref" && "$trinity_ref" != *.gz ]]; then
		log_info "[CLEANUP] Compressing Trinity assembly..."
		gzip -f "$trinity_ref" && log_info "[CLEANUP] Compressed: ${trinity_ref}.gz" || \
			log_warn "Compression failed but pipeline completed"
	fi
}

# Quantify expression using decoy-aware Salmon (Selective Alignment)
# Fast and accurate 
salmon_saf_pipeline() {
    local fasta="" genome="" rnaseq_list=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --FASTA) fasta="$2"; shift 2;;
            --GENOME) genome="$2"; shift 2;;
            --RNASEQ_LIST)
                shift
                while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
            *) log_error "Unknown arg: $1"; return 1;;
        esac
    done

    [[ -z "$fasta" || -z "$genome" ]] && { log_error "Usage: --FASTA genes.fa --GENOME genome.fa"; return 1; }
    [[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")

    local tag="$(basename "${fasta%.*}")"
    local work="tmp_${tag}_gentrome"
    local idx_dir="$SALMON_INDEX_ROOT/${tag}_decoySAF"
    local quant_root="$SALMON_QUANT_ROOT/$tag"
    local matrix_dir="$SALMON_MATRIX_ROOT/$tag"

    mkdir -p "$SALMON_INDEX_ROOT" "$quant_root" "$matrix_dir" "$work"

    # --- Build decoy-aware index ---
    if [[ -f "$idx_dir/versionInfo.json" ]]; then
        log_info "[SALMON INDEX] Decoy index already exists. Skipping."
    else
        log_step "Building decoy-aware Salmon index for $tag"
        awk '/^>/{print substr($0,2); next}{next}' "$genome" > "$work/decoys.txt"
        cat "$fasta" "$genome" > "$work/gentrome.fa"
        log_file_size "$work/gentrome.fa" "Gentrome FASTA for Salmon - $tag"
        log_file_size "$work/decoys.txt" "Decoy list for Salmon - $tag"
        run_with_space_time_log --input "$work" --output "$idx_dir" salmon index \
            -t "$work/gentrome.fa" \
            -d "$work/decoys.txt" \
            -i "$idx_dir" \
            -k 31 -p "$THREADS"
        log_file_size "$idx_dir" "Salmon index output - $tag"
        log_info "[CLEANUP] Removing temporary gentrome work directory"
        rm -rf "$work"
    fi

    # --- Quantification per SRR ---
    for SRR in "${rnaseq_list[@]}"; do
        local tdir="$TRIM_DIR_ROOT/$SRR"
        local r1=$(ls "$tdir"/${SRR}*val_1.* 2>/dev/null | head -n1)
        local r2=$(ls "$tdir"/${SRR}*val_2.* 2>/dev/null | head -n1)
        
        # Check for single-end reads if paired-end not found
        if [[ -z "$r1" ]]; then
            r1=$(ls "$tdir"/${SRR}*trimmed.* 2>/dev/null | head -n1)
        fi
        
        [[ -z "$r1" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

        local out_dir="$quant_root/$SRR"
        mkdir -p "$out_dir"
        if [[ -f "$out_dir/quant.sf" ]]; then
            log_info "[SALMON QUANT] Quantification for $SRR already exists. Skipping."
            continue
        fi

		log_step "Quantifying expression for $SRR with Salmon"
		# Determine if we have paired-end or single-end reads
		if [[ -n "$r2" && -f "$r2" ]]; then
			# Paired-end reads
			log_info "[SALMON QUANT] Using paired-end reads for $SRR"
			log_input_output_size "$r1" "$out_dir" "Salmon quantification input - $SRR"
			run_with_space_time_log --input "$tdir" --output "$out_dir" salmon quant \
				-i "$idx_dir" -l A \
				-1 "$r1" -2 "$r2" \
				-p "$THREADS" \
				--validateMappings \
				--seqBias --gcBias --posBias \
				--numBootstraps 30 \
				-o "$out_dir"
			log_file_size "$out_dir/quant.sf" "Salmon quantification output - $SRR"
		else
			# Single-end reads
			log_info "[SALMON QUANT] Using single-end reads for $SRR"
			run_with_space_time_log salmon quant \
				-i "$idx_dir" -l A \
				-r "$r1" \
				-p "$THREADS" \
				--validateMappings \
				--seqBias --gcBias --posBias \
				--numBootstraps 30 \
				-o "$out_dir"
		fi
    done

    # --- Merge matrices for Salmon ---
    log_step "Generating gene and transcript matrices (Salmon)"
    
    # Check if gene_trans_map exists, create if needed
    local gene_trans_map="${fasta}.gene_trans_map"
    if [[ ! -f "$gene_trans_map" ]]; then
        log_info "[SALMON MATRIX] Creating gene-transcript mapping file..."
        
        # Detect if this is a Trinity assembly by checking header format
        local is_trinity=false
        if grep -q "^>TRINITY_" "$fasta" 2>/dev/null; then
            is_trinity=true
            log_info "[SALMON MATRIX] Detected Trinity assembly format"
        fi
        
        # Extract gene info from FASTA headers with format-specific handling
        if [[ "$is_trinity" == "true" ]]; then
            # Trinity format: >TRINITY_DN1000_c0_g1_i1 [additional_info]
            # Gene ID: TRINITY_DN1000_c0_g1
            # Transcript ID: TRINITY_DN1000_c0_g1_i1
            awk '/^>/ {
                # Extract full transcript ID
                trans = $1
                gsub(/^>/, "", trans)
                
                # Extract gene ID by removing isoform suffix (_i[0-9]+)
                gene = trans
                if (match(gene, /^(.+)_i[0-9]+$/, arr)) {
                    gene = arr[1]
                }
                
                print gene "\t" trans
            }' "$fasta" > "$gene_trans_map"
        else
            # Generic format handling for non-Trinity assemblies
            grep "^>" "$fasta" | sed 's/^>//' | awk '{
                trans=$1
                # Try to extract gene ID from different header formats
                if (match($0, /gene=([^ ]+)/, arr)) {
                    gene=arr[1]
                } else if (match($0, /gene_id[=:]([^ ]+)/, arr)) {
                    gene=arr[1]
                } else if (match(trans, /^([^|]+)\|/, arr)) {
                    gene=arr[1]  # Format: gene|transcript
                } else {
                    gene=trans  # Fallback: use transcript as gene
                }
                print gene "\t" trans
            }' > "$gene_trans_map"
        fi
        
        local unique_genes=$(cut -f1 "$gene_trans_map" | sort -u | wc -l)
        local total_transcripts=$(wc -l < "$gene_trans_map")
        log_info "[SALMON MATRIX] Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
        
        # Validate mapping
        if [[ $total_transcripts -eq 0 ]]; then
            log_error "Failed to create gene-transcript mapping!"
            return 1
        fi
    fi
    
    # Generate count matrices using Trinity's abundance_estimates_to_matrix.pl
    if command -v abundance_estimates_to_matrix.pl >/dev/null 2>&1; then
        log_info "[SALMON MATRIX] Running abundance_estimates_to_matrix.pl for Salmon quantifications..."
        run_with_space_time_log abundance_estimates_to_matrix.pl \
            --est_method salmon \
            --gene_trans_map "$gene_trans_map" \
            --out_prefix "$matrix_dir/genes" \
            --name_sample_by_basedir "$quant_root"/*/quant.sf
    else
        log_warn "abundance_estimates_to_matrix.pl not found. Creating manual count matrix..."
        
        # Manual matrix creation from Salmon quant.sf files
        local temp_gene_ids="$matrix_dir/temp_gene_ids.txt"
        local temp_counts="$matrix_dir/temp_counts.txt"
        
        # Get gene IDs from first successful quantification
        local first_sample=""
        for SRR in "${rnaseq_list[@]}"; do
            if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
                first_sample="$SRR"
                awk 'NR>1 {print $1}' "$quant_root/$SRR/quant.sf" > "$temp_gene_ids"
                break
            fi
        done
        
        if [[ -n "$first_sample" ]]; then
            # Extract NumReads (column 5) for each sample
            for SRR in "${rnaseq_list[@]}"; do
                if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
                    awk 'NR>1 {print int($5 + 0.5)}' "$quant_root/$SRR/quant.sf" > "$matrix_dir/${SRR}_counts.tmp"
                    
                    # Validate count quality
                    local total_counts=$(awk '{sum+=$1} END {print sum}' "$matrix_dir/${SRR}_counts.tmp")
                    if [[ $total_counts -lt 100000 ]]; then
                        log_warn "Low count total for $SRR: $total_counts (expected >1M for typical RNA-seq)"
                    fi
                else
                    log_warn "Quantification not found for $SRR. Using zeros."
                    local num_genes=$(wc -l < "$temp_gene_ids")
                    yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
                fi
            done
            
            # Combine all count files
            paste "$temp_gene_ids" "$matrix_dir"/*_counts.tmp > "$temp_counts"
            
            # Create count matrix
            echo -n "gene_id" > "$matrix_dir/genes.counts.matrix"
            for SRR in "${rnaseq_list[@]}"; do
                echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"
            done
            echo "" >> "$matrix_dir/genes.counts.matrix"
            cat "$temp_counts" >> "$matrix_dir/genes.counts.matrix"
            
            # Cleanup temporary files created during matrix generation
            log_info "[CLEANUP] Removing temporary matrix generation files"
            rm -f "$temp_gene_ids" "$temp_counts" "$matrix_dir"/*_counts.tmp
        fi
    fi
    
    # --- Prepare DESeq2-compatible outputs ---
    log_step "Preparing DESeq2-compatible count matrix for Salmon pipeline"
    local deseq2_dir="$matrix_dir/deseq2_input"
    local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
    local sample_metadata="$deseq2_dir/sample_metadata.csv"
    mkdir -p "$deseq2_dir"
    
    # Verify we have quantifications
    local quant_count=0
    for SRR in "${rnaseq_list[@]}"; do
        [[ -f "$quant_root/$SRR/quant.sf" ]] && ((quant_count++))
    done
    
    if [[ $quant_count -lt 2 ]]; then
        log_error "Insufficient Salmon quantifications (found: $quant_count, need: ≥2)"
        return 1
    fi
    log_info "[SALMON] Found $quant_count samples with successful quantifications"
    
    # Convert tab-delimited matrix to CSV for DESeq2
    if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
        if [[ ! -f "$gene_count_matrix" ]]; then
            log_info "[SALMON MATRIX] Converting count matrix to CSV format for DESeq2..."
            # Convert tabs to commas and rename first column
            sed 's/\t/,/g' "$matrix_dir/genes.counts.matrix" | \
                sed '1s/gene_id/Gene_ID/' > "$gene_count_matrix"
        fi
    else
        log_warn "Count matrix not found. Cannot create DESeq2 input."
    fi
    
    # Create sample metadata file using improved function
    if [[ ! -f "$sample_metadata" ]]; then
        create_sample_metadata "$sample_metadata" rnaseq_list[@]
    fi
    
    # Generate tximport R script for Salmon (RECOMMENDED for publication)
    local tximport_script="$deseq2_dir/run_tximport_salmon.R"
    if [[ ! -f "$tximport_script" ]]; then
        generate_tximport_script "salmon" "$quant_root" "$tximport_script" "$sample_metadata"
        log_info "[TXIMPORT] Script generated: $tximport_script (recommended for DESeq2)"
    fi
    
    # Create TPM matrix for exploratory analysis
    if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
        local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
        if [[ ! -f "$tpm_matrix" ]]; then
            log_info "[SALMON MATRIX] Converting TPM matrix to CSV format..."
            sed 's/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" | \
                sed '1s/gene_id/Gene_ID/' > "$tpm_matrix"
        fi
    fi
    
    # Create aggregated summary statistics
    local summary_file="$deseq2_dir/salmon_summary.txt"
    if [[ ! -f "$summary_file" ]]; then
        log_info "[SALMON SUMMARY] Creating Salmon quantification summary..."
        {
            echo "==================================================================="
            echo "Salmon SAF Quantification Summary for $tag"
            echo "==================================================================="
            echo "Date: $(date)"
            echo ""
            echo "Samples processed: ${#rnaseq_list[@]}"
            echo "Quantification method: Salmon Selective Alignment with decoy-aware indexing"
            echo ""
            echo "Per-sample statistics:"
            echo "-------------------------------------------------------------------"
            
            for SRR in "${rnaseq_list[@]}"; do
                if [[ -f "$quant_root/$SRR/quant.sf" ]]; then
                    local total_transcripts=$(awk 'NR>1' "$quant_root/$SRR/quant.sf" | wc -l)
                    local expressed_transcripts=$(awk 'NR>1 && $5>0' "$quant_root/$SRR/quant.sf" | wc -l)
                    local total_reads=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_root/$SRR/quant.sf")
                    
                    echo "$SRR:"
                    echo "  Total transcripts: $total_transcripts"
                    echo "  Expressed transcripts (NumReads > 0): $expressed_transcripts"
                    echo "  Total mapped reads: $total_reads"
                    
                    # Extract mapping rate if available (using POSIX-compatible awk instead of grep -P)
                    if [[ -f "$quant_root/$SRR/aux_info/meta_info.json" ]]; then
                        # Extract percent_mapped without requiring grep -P
                        local mapping_rate
                        mapping_rate=$(awk -F: '/"percent_mapped"/ {gsub(/[ ,"]/, "", $2); print $2}' "$quant_root/$SRR/aux_info/meta_info.json" 2>/dev/null || echo "N/A")
                        echo "  Mapping rate: ${mapping_rate}%"
                    fi
                fi
            done
            
            echo ""
            echo "Output files:"
            echo "  - Count matrix: $gene_count_matrix"
            echo "  - Sample metadata: $sample_metadata"
            if [[ -f "$deseq2_dir/gene_tpm_matrix.csv" ]]; then
                echo "  - TPM matrix: $deseq2_dir/gene_tpm_matrix.csv"
            fi
            echo "==================================================================="
        } > "$summary_file"
        
        log_info "[SALMON SUMMARY] Summary saved to: $summary_file"
    fi
    
    # Validate count matrix quality using validation function
    if [[ -f "$gene_count_matrix" ]]; then
        validate_count_matrix "$gene_count_matrix" "gene" 2
    fi
    
    log_step "COMPLETED: Salmon-SAF pipeline for $tag"
    log_info "Count matrices: $matrix_dir/"
    log_info "DESeq2 input files:"
    log_info "  - Gene count matrix: $gene_count_matrix"
    log_info "  - Sample metadata: $sample_metadata"
    if [[ -f "$deseq2_dir/gene_tpm_matrix.csv" ]]; then
        log_info "  - TPM matrix: $deseq2_dir/gene_tpm_matrix.csv"
    fi
    log_info "  - Summary: $summary_file"
    
    # Apply normalization to expression matrices
    normalize_expression_data "$matrix_dir" "salmon"
}

# Quantify expression using Bowtie2 + RSEM
# Reviewer Preferred Method
bowtie2_rsem_pipeline() {
    local fasta="" rnaseq_list=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --FASTA) fasta="$2"; shift 2;;
            --RNASEQ_LIST)
                shift
                while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
            *) log_error "Unknown arg: $1"; return 1;;
        esac
    done

    [[ -z "$fasta" ]] && { log_error "Usage: --FASTA genes.fa"; return 1; }
    [[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	
	# Convert line endings if dos2unix is available
	if command -v dos2unix >/dev/null 2>&1; then
		dos2unix "$fasta" 2>/dev/null || true
	fi

    local tag="$(basename "${fasta%.*}")"
    local rsem_idx="$RSEM_INDEX_ROOT/${tag}_rsem"
    local quant_root="$RSEM_QUANT_ROOT/$tag"
    local matrix_dir="$RSEM_MATRIX_ROOT/$tag"

    mkdir -p "$RSEM_INDEX_ROOT" "$quant_root" "$matrix_dir"

    # --- Prepare reference ---
    if [[ -f "${rsem_idx}.grp" ]]; then
        log_info "[RSEM INDEX] RSEM reference already exists. Skipping."
    else
        log_step "Building RSEM reference for $tag"
        log_file_size "$fasta" "Input FASTA for RSEM index - $tag"
        run_with_space_time_log --input "$fasta" --output "$RSEM_INDEX_ROOT" rsem-prepare-reference --bowtie2 "$fasta" "$rsem_idx"
        log_file_size "$RSEM_INDEX_ROOT" "RSEM index output - $tag"
    fi

    # --- Quantify each SRR ---
    for SRR in "${rnaseq_list[@]}"; do
        local tdir="$TRIM_DIR_ROOT/$SRR"
        local r1=$(ls "$tdir"/${SRR}*val_1.* 2>/dev/null | head -n1)
        local r2=$(ls "$tdir"/${SRR}*val_2.* 2>/dev/null | head -n1)
        
        # Check for single-end reads if paired-end not found
        if [[ -z "$r1" ]]; then
            r1=$(ls "$tdir"/${SRR}*trimmed.* 2>/dev/null | head -n1)
        fi
        
        [[ -z "$r1" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

        local out_dir="$quant_root/$SRR"
        mkdir -p "$out_dir"
        if [[ -f "$out_dir/${SRR}.genes.results" ]]; then
            log_info "[RSEM QUANT] RSEM results for $SRR already exist. Skipping."
            continue
        fi

        log_step "Running Bowtie2 + RSEM for $SRR"
		# Determine if we have paired-end or single-end reads
		if [[ -n "$r2" && -f "$r2" ]]; then
			# Paired-end reads
			log_info "[RSEM QUANT] Using paired-end reads for $SRR"
			log_input_output_size "$r1" "$out_dir" "RSEM quantification input - $SRR"
			run_with_space_time_log --input "$tdir" --output "$out_dir" \
				rsem-calculate-expression \
					--paired-end \
					--bowtie2 \
					--num-threads "$THREADS" \
					"$r1" "$r2" "$rsem_idx" "$out_dir/$SRR"
			log_file_size "$out_dir/${SRR}.genes.results" "RSEM gene results - $SRR"
			
			# Cleanup BAM files after RSEM quantification (unless explicitly kept)
			if [[ "$keep_bam_global" != "y" ]]; then
				log_info "[CLEANUP] Deleting RSEM BAM files to save disk space"
				rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
					  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
			else
				log_info "[KEEP_BAM] Preserving RSEM BAM files for further analysis"
			fi
		else
			# Single-end reads
			log_info "[RSEM QUANT] Using single-end reads for $SRR"
			run_with_space_time_log \
				rsem-calculate-expression \
					--bowtie2 \
					--num-threads "$THREADS" \
					"$r1" "$rsem_idx" "$out_dir/$SRR"
			
			# Cleanup BAM files after RSEM quantification (unless explicitly kept)
			if [[ "$keep_bam_global" != "y" ]]; then
				log_info "[CLEANUP] Deleting RSEM BAM files to save disk space"
				rm -f "$out_dir/${SRR}.transcript.bam" "$out_dir/${SRR}.genome.bam" \
					  "$out_dir/${SRR}.transcript.sorted.bam" "$out_dir/${SRR}.transcript.sorted.bam.bai"
			else
				log_info "[KEEP_BAM] Preserving RSEM BAM files for further analysis"
			fi
		fi
    done

    # --- Merge matrices for RSEM ---
    log_step "Generating gene and transcript matrices (RSEM)"
    
    # Check if gene_trans_map exists, create if needed
    local gene_trans_map="${fasta}.gene_trans_map"
    if [[ ! -f "$gene_trans_map" ]]; then
        log_info "[RSEM MATRIX] Creating gene-transcript mapping file..."
        
        # Detect if this is a Trinity assembly by checking header format
        local is_trinity=false
        if grep -q "^>TRINITY_" "$fasta" 2>/dev/null; then
            is_trinity=true
            log_info "[RSEM MATRIX] Detected Trinity assembly format"
        fi
        
        # Extract gene info from FASTA headers with format-specific handling
        if [[ "$is_trinity" == "true" ]]; then
            # Trinity format: >TRINITY_DN1000_c0_g1_i1 [additional_info]
            # Gene ID: TRINITY_DN1000_c0_g1
            # Transcript ID: TRINITY_DN1000_c0_g1_i1
            awk '/^>/ {
                # Extract full transcript ID
                trans = $1
                gsub(/^>/, "", trans)
                
                # Extract gene ID by removing isoform suffix (_i[0-9]+)
                gene = trans
                if (match(gene, /^(.+)_i[0-9]+$/, arr)) {
                    gene = arr[1]
                }
                
                print gene "\t" trans
            }' "$fasta" > "$gene_trans_map"
        else
            # Generic format handling for non-Trinity assemblies
            grep "^>" "$fasta" | sed 's/^>//' | awk '{
                trans=$1
                # Try to extract gene ID from different header formats
                if (match($0, /gene=([^ ]+)/, arr)) {
                    gene=arr[1]
                } else if (match($0, /gene_id[=:]([^ ]+)/, arr)) {
                    gene=arr[1]
                } else if (match(trans, /^([^|]+)\|/, arr)) {
                    gene=arr[1]  # Format: gene|transcript
                } else {
                    gene=trans  # Fallback: use transcript as gene
                }
                print gene "\t" trans
            }' > "$gene_trans_map"
        fi
        
        local unique_genes=$(cut -f1 "$gene_trans_map" | sort -u | wc -l)
        local total_transcripts=$(wc -l < "$gene_trans_map")
        log_info "[RSEM MATRIX] Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
        
        # Validate mapping
        if [[ $total_transcripts -eq 0 ]]; then
            log_error "Failed to create gene-transcript mapping!"
            return 1
        fi
    fi
    
    run_with_space_time_log \
		abundance_estimates_to_matrix.pl \
			--est_method RSEM \
			--gene_trans_map "$gene_trans_map" \
			--out_prefix "$matrix_dir/genes" \
			--name_sample_by_basedir "$quant_root"/*/*.genes.results || {
		log_warn "abundance_estimates_to_matrix.pl failed. Creating manual count matrix..."
		
		# Manual matrix creation from RSEM .genes.results files
		local temp_gene_ids="$matrix_dir/temp_gene_ids.txt"
		local temp_counts="$matrix_dir/temp_counts.txt"
		local temp_tpm="$matrix_dir/temp_tpm.txt"
		local temp_fpkm="$matrix_dir/temp_fpkm.txt"
		
		# Get gene IDs from first successful quantification
		local first_sample=""
		for SRR in "${rnaseq_list[@]}"; do
			if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
				first_sample="$SRR"
				awk 'NR>1 {print $1}' "$quant_root/$SRR/${SRR}.genes.results" > "$temp_gene_ids"
				break
			fi
		done
		
		if [[ -n "$first_sample" ]]; then
			# Extract expected_count (column 5), TPM (column 6), and FPKM (column 7) for each sample
			for SRR in "${rnaseq_list[@]}"; do
				if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
					# Expected counts (rounded to integers)
					awk 'NR>1 {print int($5 + 0.5)}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_counts.tmp"
					# TPM values
					awk 'NR>1 {print $6}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_tpm.tmp"
					# FPKM values
					awk 'NR>1 {print $7}' "$quant_root/$SRR/${SRR}.genes.results" > "$matrix_dir/${SRR}_fpkm.tmp"
				else
					log_warn "RSEM results not found for $SRR. Using zeros."
					local num_genes=$(wc -l < "$temp_gene_ids")
					yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_counts.tmp"
					yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_tpm.tmp"
					yes 0 | head -n "$num_genes" > "$matrix_dir/${SRR}_fpkm.tmp"
				fi
			done
			
			# Combine all count files
			paste "$temp_gene_ids" "$matrix_dir"/*_counts.tmp > "$temp_counts"
			paste "$temp_gene_ids" "$matrix_dir"/*_tpm.tmp > "$temp_tpm"
			paste "$temp_gene_ids" "$matrix_dir"/*_fpkm.tmp > "$temp_fpkm"
			
			# Create count matrix
			echo -n "gene_id" > "$matrix_dir/genes.counts.matrix"
			for SRR in "${rnaseq_list[@]}"; do
				echo -ne "\t$SRR" >> "$matrix_dir/genes.counts.matrix"
			done
			echo "" >> "$matrix_dir/genes.counts.matrix"
			cat "$temp_counts" >> "$matrix_dir/genes.counts.matrix"
			
			# Create TPM matrix
			echo -n "gene_id" > "$matrix_dir/genes.TPM.not_cross_norm"
			for SRR in "${rnaseq_list[@]}"; do
				echo -ne "\t$SRR" >> "$matrix_dir/genes.TPM.not_cross_norm"
			done
			echo "" >> "$matrix_dir/genes.TPM.not_cross_norm"
			cat "$temp_tpm" >> "$matrix_dir/genes.TPM.not_cross_norm"
			
			# Create FPKM matrix
			echo -n "gene_id" > "$matrix_dir/genes.FPKM.not_cross_norm"
			for SRR in "${rnaseq_list[@]}"; do
				echo -ne "\t$SRR" >> "$matrix_dir/genes.FPKM.not_cross_norm"
			done
			echo "" >> "$matrix_dir/genes.FPKM.not_cross_norm"
			cat "$temp_fpkm" >> "$matrix_dir/genes.FPKM.not_cross_norm"
			
			# Cleanup temporary files created during matrix generation
			log_info "[CLEANUP] Removing temporary RSEM matrix generation files"
			rm -f "$temp_gene_ids" "$temp_counts" "$temp_tpm" "$temp_fpkm" \
				  "$matrix_dir"/*_counts.tmp "$matrix_dir"/*_tpm.tmp "$matrix_dir"/*_fpkm.tmp
		fi
	}
    
    # --- Prepare DESeq2-compatible outputs ---
	log_step "Preparing DESeq2-compatible count matrix for RSEM pipeline"
	# Concise note: prefer tximport for DESeq2; a script will be generated below
	log_info "[NOTE] RSEM reports expected counts; for DESeq2, prefer tximport (script will be generated)."
    local deseq2_dir="$matrix_dir/deseq2_input"
    local gene_count_matrix="$deseq2_dir/gene_count_matrix.csv"
    local sample_metadata="$deseq2_dir/sample_metadata.csv"
    mkdir -p "$deseq2_dir"
    
    # Verify we have quantifications
    local quant_count=0
    for SRR in "${rnaseq_list[@]}"; do
        [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]] && ((quant_count++))
    done
    
    if [[ $quant_count -lt 2 ]]; then
        log_error "Insufficient RSEM quantifications (found: $quant_count, need: ≥2)"
        return 1
    fi
    log_info "[RSEM] Found $quant_count samples with successful quantifications"
    
    # Convert tab-delimited matrix to CSV for DESeq2
    if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
        if [[ ! -f "$gene_count_matrix" ]]; then
            log_info "[RSEM MATRIX] Converting count matrix to CSV format for DESeq2..."
            # Convert tabs to commas and rename first column
            sed 's/\t/,/g' "$matrix_dir/genes.counts.matrix" | \
                sed '1s/gene_id/Gene_ID/' > "$gene_count_matrix"
        fi
    else
        log_warn "Count matrix not found. Cannot create DESeq2 input."
    fi
    
    # Create sample metadata file using improved function
    if [[ ! -f "$sample_metadata" ]]; then
        create_sample_metadata "$sample_metadata" rnaseq_list[@]
    fi
    
    # Generate tximport R script for RSEM (RECOMMENDED for publication)
    local tximport_script="$deseq2_dir/run_tximport_rsem.R"
    if [[ ! -f "$tximport_script" ]]; then
        generate_tximport_script "rsem" "$quant_root" "$tximport_script" "$sample_metadata"
        log_info "[TXIMPORT] Script generated: $tximport_script (recommended for DESeq2)"
    fi
    
    # Create TPM and FPKM matrices for exploratory analysis
    if [[ -f "$matrix_dir/genes.TPM.not_cross_norm" ]]; then
        local tpm_matrix="$deseq2_dir/gene_tpm_matrix.csv"
        if [[ ! -f "$tpm_matrix" ]]; then
            log_info "[RSEM MATRIX] Converting TPM matrix to CSV format..."
            sed 's/\t/,/g' "$matrix_dir/genes.TPM.not_cross_norm" | \
                sed '1s/gene_id/Gene_ID/' > "$tpm_matrix"
        fi
    fi
    
    if [[ -f "$matrix_dir/genes.FPKM.not_cross_norm" ]]; then
        local fpkm_matrix="$deseq2_dir/gene_fpkm_matrix.csv"
        if [[ ! -f "$fpkm_matrix" ]]; then
            log_info "[RSEM MATRIX] Converting FPKM matrix to CSV format..."
            sed 's/\t/,/g' "$matrix_dir/genes.FPKM.not_cross_norm" | \
                sed '1s/gene_id/Gene_ID/' > "$fpkm_matrix"
        fi
    fi
    
    # Create aggregated summary statistics
    local summary_file="$deseq2_dir/rsem_summary.txt"
    if [[ ! -f "$summary_file" ]]; then
        log_info "[RSEM SUMMARY] Creating RSEM quantification summary..."
        {
            echo "==================================================================="
            echo "RSEM Quantification Summary for $tag"
            echo "==================================================================="
            echo "Date: $(date)"
            echo ""
            echo "Samples processed: ${#rnaseq_list[@]}"
            echo "Quantification method: RSEM with Bowtie2 alignment"
            echo ""
            echo "Per-sample statistics:"
            echo "-------------------------------------------------------------------"
            
            for SRR in "${rnaseq_list[@]}"; do
                if [[ -f "$quant_root/$SRR/${SRR}.genes.results" ]]; then
                    local total_genes=$(awk 'NR>1' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
                    local expressed_genes=$(awk 'NR>1 && $5>0' "$quant_root/$SRR/${SRR}.genes.results" | wc -l)
                    local total_counts=$(awk 'NR>1 {sum+=$5} END {print int(sum)}' "$quant_root/$SRR/${SRR}.genes.results")
                    
                    echo "$SRR:"
                    echo "  Total genes: $total_genes"
                    echo "  Expressed genes (count > 0): $expressed_genes"
                    echo "  Total expected counts: $total_counts"
                fi
            done
            
            echo ""
            echo "Output files:"
            echo "  - Count matrix: $gene_count_matrix"
            echo "  - Sample metadata: $sample_metadata"
            if [[ -f "$deseq2_dir/gene_tpm_matrix.csv" ]]; then
                echo "  - TPM matrix: $deseq2_dir/gene_tpm_matrix.csv"
            fi
            if [[ -f "$deseq2_dir/gene_fpkm_matrix.csv" ]]; then
                echo "  - FPKM matrix: $deseq2_dir/gene_fpkm_matrix.csv"
            fi
            echo "==================================================================="
        } > "$summary_file"
        
        log_info "[RSEM SUMMARY] Summary saved to: $summary_file"
    fi

    # Validate count matrix quality using validation function
    if [[ -f "$gene_count_matrix" ]]; then
        validate_count_matrix "$gene_count_matrix" "gene" 2
    fi

    log_step "COMPLETED: Bowtie2-RSEM pipeline for $tag"
    log_info "Count matrices: $matrix_dir/"
    log_info "DESeq2 input files:"
    log_info "  - Gene count matrix: $gene_count_matrix"
    log_info "  - Sample metadata: $sample_metadata"
    if [[ -f "$deseq2_dir/gene_tpm_matrix.csv" ]]; then
        log_info "  - TPM matrix: $deseq2_dir/gene_tpm_matrix.csv"
    fi
    if [[ -f "$deseq2_dir/gene_fpkm_matrix.csv" ]]; then
        log_info "  - FPKM matrix: $deseq2_dir/gene_fpkm_matrix.csv"
    fi
    log_info "  - Summary: $summary_file"
    
    # Apply normalization to expression matrices
    normalize_expression_data "$matrix_dir" "RSEM"
}

# ------------------------------------------------------------------------------
# METHOD COMPARISON AND VALIDATION FUNCTIONS
# ------------------------------------------------------------------------------

# IMPORTANT: Sample Metadata Configuration
# =========================================
# All pipeline methods create a default sample_metadata.csv file with:
#   - All samples assigned to "treatment" condition
#   - All samples in batch "1"
# 
# THIS MUST BE CUSTOMIZED for differential expression analysis!
# 
# Option 1: Create a sample_conditions.txt file with format:
#   SRR_ID condition batch
#   SRR3884597 control 1
#   SRR3884653 treatment 1
#
# Option 2: Manually edit the generated sample_metadata.csv files:
#   - Located in: */deseq2_input/sample_metadata.csv
#   - Assign appropriate experimental conditions and batch information
#
# Without proper condition assignment, DESeq2 cannot perform differential
# expression analysis!
# =========================================

normalize_expression_data() {
	# Apply normalization to expression matrices
	local matrix_dir="$1"
	local method="$2"
	
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		log_step "Applying TMM normalization for $method"
		if command -v normalize_matrix.pl >/dev/null 2>&1; then
			run_with_space_time_log normalize_matrix.pl "$matrix_dir/genes.counts.matrix" \
				--est_method "$method" \
				--out_prefix "$matrix_dir/genes.TMM" || \
				log_warn "TMM normalization failed for $method"
		else
			log_warn "normalize_matrix.pl not found. Skipping TMM normalization."
		fi
	fi
}

# Cross-method validation and comparison
compare_methods_summary() {
	local fasta_tag="$1"
	
	log_step "Cross-Method Validation Summary for $fasta_tag"
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
	
	# Method 1: HISAT2 Reference-Guided
	local m1_matrix="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
	if [[ -f "$m1_matrix" ]]; then
		local m1_genes=$(tail -n +2 "$m1_matrix" | wc -l)
		local m1_samples=$(head -n1 "$m1_matrix" | tr ',' '\n' | tail -n +2 | wc -l)
		log_info "✅ Method 1 (HISAT2 Ref-Guided): $m1_genes genes, $m1_samples samples"
		log_info "   Status: BEST for publication (true read counts via prepDE.py)"
	fi
	
	# Method 3: Trinity De Novo
	local m3_matrix="$STRINGTIE_TRINITY_DE_NOVO_ROOT/$fasta_tag/deseq2_input/gene_count_matrix.csv"
	local m3_tximport="$STRINGTIE_TRINITY_DE_NOVO_ROOT/$fasta_tag/deseq2_input/run_tximport_trinity.R"
	if [[ -f "$m3_matrix" ]]; then
		local m3_genes=$(tail -n +2 "$m3_matrix" | wc -l)
		local m3_samples=$(head -n1 "$m3_matrix" | tr ',' '\n' | tail -n +2 | wc -l)
		log_info "⚠️  Method 3 (Trinity): $m3_genes genes, $m3_samples samples"
		log_info "   Status: Uses RSEM expected_count (statistical estimates)"
		if [[ -f "$m3_tximport" ]]; then
			log_info "   ✅ tximport script available: $m3_tximport"
		fi
	fi
	
	# Method 4: Salmon SAF
	local m4_matrix="$SALMON_SAF_ROOT/$fasta_tag/matrices/deseq2_input/gene_count_matrix.csv"
	local m4_tximport="$SALMON_SAF_ROOT/$fasta_tag/matrices/deseq2_input/run_tximport_salmon.R"
	if [[ -f "$m4_matrix" ]]; then
		local m4_genes=$(tail -n +2 "$m4_matrix" | wc -l)
		local m4_samples=$(head -n1 "$m4_matrix" | tr ',' '\n' | tail -n +2 | wc -l)
		log_info "✅ Method 4 (Salmon SAF): $m4_genes genes, $m4_samples samples"
		log_info "   Status: EXCELLENT for publication (fast, accurate, modern)"
		if [[ -f "$m4_tximport" ]]; then
			log_info "   ✅ tximport script available: $m4_tximport"
		fi
	fi
	
	# Method 5: Bowtie2 + RSEM
	local m5_matrix="$BOWTIE2_RSEM_ROOT/$fasta_tag/matrices/deseq2_input/gene_count_matrix.csv"
	local m5_tximport="$BOWTIE2_RSEM_ROOT/$fasta_tag/matrices/deseq2_input/run_tximport_rsem.R"
	if [[ -f "$m5_matrix" ]]; then
		local m5_genes=$(tail -n +2 "$m5_matrix" | wc -l)
		local m5_samples=$(head -n1 "$m5_matrix" | tr ',' '\n' | tail -n +2 | wc -l)
		log_info "⚠️  Method 5 (Bowtie2+RSEM): $m5_genes genes, $m5_samples samples"
		log_info "   Status: Uses RSEM expected_count (statistical estimates)"
		if [[ -f "$m5_tximport" ]]; then
			log_info "   ✅ tximport script available: $m5_tximport"
		fi
	fi
	
	log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}

