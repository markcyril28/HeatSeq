#!/bin/bash

# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE - METHOD 2: HISAT2 DE NOVO
# ==============================================================================
# Description: RNA-seq analysis pipeline using HISAT2 de novo assembly approach
# Author: Mark Cyril R. Mercado
# Version: v10
# Date: October 2025
# 
# Pipeline Methods:
# - Method 1: HISAT2 Reference Guided (not in this script)
# - Method 2: HISAT2 De Novo Assembly (this script)
# - Method 3: Trinity De Novo Assembly
# - Method 4: Salmon SAF Quantification
# - Method 5: Bowtie2 + RSEM Quantification
#
# SCRIPT ORGANIZATION:
# 1. Configuration and Runtime Switches (All 5 Methods + QC Options)
# 2. Pipeline Configuration Examples  
# 3. Input Files and Data Sources
# 4. RNA-seq Data Sources (SRR Lists)
# 5. Directory Structure and Output Paths  
# 6. Cleanup Options and Testing
# 7. Logging System and Utility Functions
# 8. Pipeline Functions:
#    	- Data Download and Preprocessing Functions
#		- Pipelines Functions	 
# 9. Main Execution Functions
# 10. Script Execution
# 11. Post-processing Options
# ==============================================================================

set -euo pipefail


# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# Raw and Processed Data Directories
RAW_DIR_ROOT="1_RAW_SRR"                # Raw SRR download directory
TRIM_DIR_ROOT="2_TRIMMED_SRR"           # Trimmed reads directory
FASTQC_ROOT="3_FastQC"                  # FastQC quality control reports

# Method 1: HISAT2 Reference Guided 
HISAT2_REF_GUIDED_ROOT="4_OUTPUTS/4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="4_OUTPUTS/4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="4_OUTPUTS/4a_Method_1_HISAT2_Ref_Guided/5_stringtie_WD/a_Method_1_RAW_RESULTs"

# Method 2: HISAT2 De Novo Assembly (main method in this script)
HISAT2_DE_NOVO_ROOT="4_OUTPUTS/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="4_OUTPUTS/4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4_OUTPUTS/4b_Method_2_HISAT2_De_Novo/5_stringtie_WD/a_Method_2_RAW_RESULTs"

# Method 3: Trinity De Novo Assembly
TRINITY_DE_NOVO_ROOT="4_OUTPUTS/4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4_OUTPUTS/4c_Method_3_Trinity_De_Novo/5_stringtie_WD/a_Method_3_RAW_RESULTs"

# Method 4: Salmon SAF Quantification
SALMON_SAF_ROOT="4_OUTPUTS/4d_Method_4_Salmon_Saf_Quantification"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/quant"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/matrices"

# Method 5: Bowtie2 + RSEM Quantification
BOWTIE2_RSEM_ROOT="4_OUTPUTS/4e_Method_5_Bowtie2_Quantification"
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

# ------------------------------------------------------------------------------
# MAIN PIPELINES FOR ALIGNMENT, ASSEMBLY, AND QUANTIFICATION
# ------------------------------------------------------------------------------

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
		log_info "HISAT2 reference-guided index already exists for $fasta_base. Skipping build."
	else
		log_info "Building HISAT2 reference-guided index from $fasta with GTF $gtf..."
		
		# Extract splice sites and exons from GTF
		local splice_sites="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_splice_sites.txt"
		local exons="$HISAT2_REF_GUIDED_INDEX_DIR/${fasta_tag}_exons.txt"
		
		run_with_time_to_log hisat2_extract_splice_sites.py "$gtf" > "$splice_sites"
		run_with_time_to_log hisat2_extract_exons.py "$gtf" > "$exons"
		
		# Build index with splice sites and exons
		run_with_time_to_log hisat2-build \
			-p "${THREADS}" \
			--ss "$splice_sites" \
			--exon "$exons" \
			"$fasta" \
			"$index_prefix"
	fi

	# ALIGNMENT, SORTING, AND STRINGTIE ASSEMBLY for each SRR
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
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

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped_sorted.bam"
		local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_ref_guided_mapped.sam"
		
		# Align and sort if BAM doesn't exist
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "BAM for $SRR and $fasta_tag (ref-guided) already exists. Skipping alignment."
		else
			log_info "Aligning $fasta_tag with $SRR using reference-guided approach..."
			# Check if paired-end or single-end reads
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				# Paired-end alignment
				run_with_time_to_log \
					hisat2 -p "${THREADS}" --dta \
						$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--rna-strandness $RNA_STRAND_PROTOCOL") \
						-x "$index_prefix" \
						-1 "$trimmed1" \
						-2 "$trimmed2" \
						-S "$sam"
			else
				# Single-end alignment
				run_with_time_to_log \
					hisat2 -p "${THREADS}" --dta \
						$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--rna-strandness $RNA_STRAND_PROTOCOL") \
						-x "$index_prefix" \
						-U "$trimmed1" \
						-S "$sam"
			fi
			
			log_info "Converting SAM to sorted BAM for $fasta_tag with $SRR (ref-guided)..."
			run_with_time_to_log samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_time_to_log samtools index -@ "${THREADS}" "$bam"
			rm -f "$sam"
			log_info "Done aligning $fasta_tag with $SRR (ref-guided)."
		fi

		# StringTie assembly with reference GTF
		local out_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		local out_gene_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_ref_guided_gene_abundances.tsv"
		local ballgown_dir="$out_dir/ballgown"
		mkdir -p "$out_dir" "$ballgown_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "Reference-guided assembly exists for $fasta_tag/$SRR. Skipping."
		else
			log_info "Assembling transcripts for $fasta_tag with $SRR using reference GTF..."
			run_with_time_to_log \
				stringtie -p "$THREADS" "$bam" \
					-G "$gtf" \
					-o "$out_gtf" \
					-A "$out_gene_abundances_tsv" \
					-B \
					-e \
					-C "$out_dir/${SRR}_${fasta_tag}_ref_guided_cov_refs.gtf"
		fi
		
		# Cleanup BAM files if specified
		if [[ "$keep_bam_global" != "y" ]]; then
			log_info "Deleting the BAM file for $SRR (ref-guided)."
			rm -f "$bam" "${bam}.bai"
		fi
		
		log_info "Done processing $fasta_tag with $SRR (ref-guided)."
		log_info "--------------------------------------------------"
	done
	
	# Create merged GTF for all samples
	log_info "Creating merged GTF file for reference-guided assembly..."
	local merge_dir="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/merged"
	local merged_gtf="$merge_dir/${fasta_tag}_ref_guided_merged.gtf"
	local gtf_list="$merge_dir/gtf_list.txt"
	mkdir -p "$merge_dir"
	
	# Create list of GTF files for merging
	> "$gtf_list"
	for SRR in "${rnaseq_list[@]}"; do
		local out_gtf="$STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/$SRR/${SRR}_${fasta_tag}_ref_guided_stringtie_assembled.gtf"
		if [[ -f "$out_gtf" ]]; then
			echo "$out_gtf" >> "$gtf_list"
		fi
	done
	
	if [[ -f "$merged_gtf" ]]; then
		log_info "Merged GTF already exists for $fasta_tag (ref-guided). Skipping merge."
	else
		log_info "Merging GTF files for $fasta_tag (ref-guided)..."
		run_with_time_to_log stringtie --merge \
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
		
		mkdir -p "$final_dir"
		
		# Skip if BAM was deleted and final quantification already exists
		if [[ ! -f "$bam" && -f "$final_gtf" ]]; then
			log_info "Final quantification already exists for $SRR (ref-guided). Skipping re-estimation."
			continue
		fi
		
		if [[ -f "$bam" ]]; then
			log_info "Re-estimating abundances for $SRR with merged GTF (ref-guided)..."
			run_with_time_to_log stringtie \
				-p "$THREADS" \
				-e -B \
				-G "$merged_gtf" \
				-A "$final_abundances" \
				-o "$final_gtf" \
				"$bam"
		fi
	done
	
	log_info "HISAT2 reference-guided pipeline completed for $fasta_tag."
	log_info "Merged GTF: $merged_gtf"
	log_info "Final quantifications in: $STRINGTIE_HISAT2_REF_GUIDED_ROOT/$fasta_tag/*/final/"
}

hisat2_de_novo_pipeline() {
	# Combined pipeline: Build HISAT2 index, align reads, assemble, merge, and quantify transcripts
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
		rnaseq_list=("${SRR_LIST_PRJNA328564[@]}")
	fi

	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"

	# BUILD HISAT2 INDEX
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "HISAT2 index already exists for $fasta_base. Skipping build."
	else
		log_info "Building HISAT2 index from $fasta ..."
		run_with_time_to_log hisat2-build -p "${THREADS}" "$fasta" "$index_prefix"
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
			log_info "BAM for $SRR and $fasta_tag already exists. Skipping alignment."
		else
			log_info "Aligning $fasta_tag with $SRR..."
			# Check if paired-end or single-end reads
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				# Paired-end alignment
				run_with_time_to_log \
					hisat2 -p "${THREADS}" --dta \
						$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--rna-strandness $RNA_STRAND_PROTOCOL") \
						-x "$index_prefix" \
						-1 "$trimmed1" \
						-2 "$trimmed2" \
						-S "$sam"
			else
				# Single-end alignment
				run_with_time_to_log \
					hisat2 -p "${THREADS}" --dta \
						$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--rna-strandness $RNA_STRAND_PROTOCOL") \
						-x "$index_prefix" \
						-U "$trimmed1" \
						-S "$sam"
			fi
			
			log_info "Converting SAM to sorted BAM for $fasta_tag with $SRR..."
			run_with_time_to_log samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_time_to_log samtools index -@ "${THREADS}" "$bam"
			log_info "Deleting the SAM file."
			rm -f "$sam"
			log_info "Done aligning $fasta_tag with $SRR."
		fi

		# StringTie assembly
		local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
		local out_gene_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo.tsv"
		mkdir -p "$out_dir"
		
		if [[ -f "$out_gtf" ]]; then
			log_info "Assembly exists for $fasta_tag/$SRR. Skipping."
		else
			log_info "Assembling transcripts for $fasta_tag with $SRR ..."
			run_with_time_to_log \
				stringtie -p "$THREADS" "$bam" \
					-o "$out_gtf" \
					-A "$out_gene_abundances_tsv"
		fi
		
		log_info "Deleting the SAM and BAM file."
		rm -f "$bam" "${bam}.bai"
		log_info "Done processing $fasta_tag with $SRR."
		log_info "--------------------------------------------------"
	done
}

# Trinity de novo alignment, stringtie quantification that can an input to DeSeq2 pipeline
trinity_de_novo_alignment_pipeline() {
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
		rnaseq_list=("${SRR_LIST_PRJNA328564[@]}")
	fi

	local fasta_base fasta_tag trinity_out_dir trinity_fasta
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	trinity_out_dir="$TRINITY_DE_NOVO_ROOT/${fasta_tag}_trinity_out_dir"
	trinity_fasta="$trinity_out_dir/Trinity.fasta"

	# STEP 1: DE NOVO TRANSCRIPTOME ASSEMBLY WITH TRINITY
	log_step "Trinity de novo transcriptome assembly for $fasta_tag"
	
	if [[ -f "$trinity_fasta" ]]; then
		log_info "Trinity assembly already exists for $fasta_tag. Skipping assembly."
	else
		log_info "Starting Trinity de novo assembly for $fasta_tag..."
		mkdir -p "$trinity_out_dir"
		
		# Collect all trimmed FASTQ files for Trinity
		local left_reads=() right_reads=() single_reads=()
		local paired_count=0 single_count=0
		
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
				((paired_count++))
				log_info "Added paired-end reads for $SRR: $trimmed1, $trimmed2"
			elif [[ -n "$trimmed1" && -f "$trimmed1" ]]; then
				# Single-end reads
				single_reads+=("$trimmed1")
				((single_count++))
				log_info "Added single-end reads for $SRR: $trimmed1"
			else
				log_warn "Trimmed FASTQ for $SRR not found; skipping from Trinity input."
			fi
		done
		
		if [[ ${#left_reads[@]} -eq 0 && ${#single_reads[@]} -eq 0 ]]; then
			log_error "No trimmed FASTQ files found for Trinity assembly."
			return 1
		fi
		
		log_info "Found $paired_count paired-end samples and $single_count single-end samples"
		
		# Run Trinity based on read type
		if [[ ${#left_reads[@]} -gt 0 && ${#right_reads[@]} -gt 0 ]]; then
			# Only paired-end reads
			local left_files=$(IFS=,; echo "${left_reads[*]}")
			local right_files=$(IFS=,; echo "${right_reads[*]}")
			log_info "Running Trinity with paired-end reads only..."
			run_with_time_to_log Trinity \
				--seqType fq \
				--left "$left_files" \
				--right "$right_files" \
				--CPU "$THREADS" \
				--max_memory 20G \
				--output "$trinity_out_dir" \
				--normalize_reads \
				--full_cleanup \
				--no_version_check
		elif [[ ${#single_reads[@]} -gt 0 ]]; then
			# Only single-end reads
			local single_files=$(IFS=,; echo "${single_reads[*]}")
			log_info "Running Trinity with single-end reads only..."
			run_with_time_to_log Trinity \
				--seqType fq \
				--single "$single_files" \
				--CPU "$THREADS" \
				--max_memory 20G \
				--output "$trinity_out_dir" \
				--normalize_reads \
				--full_cleanup \
				--no_version_check
		else
			log_error "No valid reads found for Trinity assembly."
			return 1
		fi
		
		# Optional: Compress Trinity assembly to save space
		if [[ -f "$trinity_fasta" ]]; then
			log_info "Compressing Trinity assembly to save disk space..."
			gzip -f "$trinity_fasta"
			trinity_fasta="${trinity_fasta}.gz"
		fi
		
		# Optional: Remove Trinity intermediate directories if only the assembly is needed
		log_info "Cleaning up Trinity intermediate files..."
		find "$trinity_out_dir" -name "chrysalis" -type d -exec rm -rf {} + 2>/dev/null || true
		find "$trinity_out_dir" -name "jellyfish.kmers*" -delete 2>/dev/null || true
	fi

	# STEP 2: ALIGN READS TO TRINITY ASSEMBLY AND RUN STRINGTIE FOR QUANTIFICATION
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/${fasta_tag}_trinity/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="" trimmed2=""
		mkdir -p "$HISAT2_DIR"
		
		# Find trimmed FASTQ files (same logic as above)
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
			# Single-end trimmed files
			local files=( "$TrimGalore_DIR"/${SRR}*trimmed.* )
			trimmed1="${files[0]}"
		elif compgen -G "$TrimGalore_DIR/${SRR}*.fq*" >/dev/null 2>&1; then
			# Generic single-end pattern
			local files=( "$TrimGalore_DIR"/${SRR}*.fq* )
			trimmed1="${files[0]}"
		fi
		
		if [[ -z "$trimmed1" ]]; then
			log_warn "Trimmed FASTQ for $SRR not found in $TrimGalore_DIR; skipping."
			continue
		fi

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trinity_mapped_sorted.bam"
		local abundance_dir="$HISAT2_DIR/trinity_abundance_${SRR}_${fasta_tag}"
		
		# Use Trinity's native alignment and quantification workflow
		if [[ -f "$abundance_dir/RSEM.genes.results" ]]; then
			log_info "Trinity abundance estimation for $SRR and $fasta_tag already exists. Skipping."
		else
			log_info "Running Trinity align_and_estimate_abundance for $SRR to Trinity assembly of $fasta_tag..."
			
			# Trinity's recommended workflow for quantification
			mkdir -p "$abundance_dir"
			# Determine if we have paired-end or single-end reads
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				# Paired-end reads
				log_info "Using paired-end reads for $SRR"
				run_with_time_to_log align_and_estimate_abundance.pl \
					--transcripts "$trinity_fasta" \
					--seqType fq \
					--left "$trimmed1" \
					--right "$trimmed2" \
					--est_method RSEM \
					--aln_method bowtie2 \
					--trinity_mode \
					--prep_reference \
					--thread_count "$THREADS" \
					--output_dir "$abundance_dir"
			else
				# Single-end reads
				log_info "Using single-end reads for $SRR"
				run_with_time_to_log align_and_estimate_abundance.pl \
					--transcripts "$trinity_fasta" \
					--seqType fq \
					--single "$trimmed1" \
					--est_method RSEM \
					--aln_method bowtie2 \
					--trinity_mode \
					--prep_reference \
					--thread_count "$THREADS" \
					--output_dir "$abundance_dir"
			fi
			
			# Create BAM file link for StringTie compatibility (if needed)
			if [[ -f "$abundance_dir/bowtie2.bam" ]]; then
				ln -sf "$abundance_dir/bowtie2.bam" "$bam"
				ln -sf "$abundance_dir/bowtie2.bam.bai" "${bam}.bai" 2>/dev/null || samtools index "$bam"
			fi
		fi

		# StringTie quantification for DeSeq2 compatibility (using Trinity abundance if BAM exists)
		local out_dir="$STRINGTIE_TRINITY_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_trinity_stringtie_quantified.gtf"
		local out_gene_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_trinity_gene_abundances.tsv"
		local out_transcript_abundances_tsv="$out_dir/${SRR}_${fasta_tag}_trinity_transcript_abundances.tsv"
		mkdir -p "$out_dir"
		
		if [[ -f "$out_gtf" && -f "$out_gene_abundances_tsv" ]]; then
			log_info "StringTie quantification exists for Trinity assembly of $fasta_tag/$SRR. Skipping."
		else
			log_info "Quantifying transcripts with StringTie for Trinity assembly of $fasta_tag with $SRR..."
			
			# First pass: assemble transcripts without reference
			run_with_time_to_log \
				stringtie -p "$THREADS" "$bam" \
					-o "$out_gtf" \
					-A "$out_gene_abundances_tsv" \
					-C "$out_dir/${SRR}_${fasta_tag}_trinity_cov_refs.gtf" \
					-B \
					-e \
					-G "$out_gtf"
					
			# Create transcript abundance file for DeSeq2
			run_with_time_to_log \
				stringtie -p "$THREADS" "$bam" \
					-e -B \
					-G "$out_gtf" \
					-A "$out_transcript_abundances_tsv" \
					-o "$out_dir/${SRR}_${fasta_tag}_trinity_final.gtf"
		fi
		
		# Cleanup BAM files if specified
		if [[ "$keep_bam_global" != "y" ]]; then
			log_info "Deleting the BAM file for $SRR."
			rm -f "$bam" "${bam}.bai"
		fi
		
		log_info "Done processing Trinity assembly of $fasta_tag with $SRR."
		log_info "--------------------------------------------------"
	done
	
	# STEP 3: PREPARE DESEQ2 INPUT FILES
	log_info "Preparing count matrices for DeSeq2..."
	local deseq2_dir="$STRINGTIE_TRINITY_DE_NOVO_ROOT/$fasta_tag/deseq2_input"
	mkdir -p "$deseq2_dir"
	
	# Create sample table for DESeq2
	local sample_table="$deseq2_dir/sample_table.csv"
	echo "sample,condition,path" > "$sample_table"
	
	for SRR in "${rnaseq_list[@]}"; do
		local out_dir="$STRINGTIE_TRINITY_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local gene_count_file="$out_dir/${SRR}_${fasta_tag}_trinity_gene_abundances.tsv"
		
		if [[ -f "$gene_count_file" ]]; then
			# Extract condition from SRR (you may need to customize this)
			local condition="sample"  # Default condition, customize based on your metadata
			echo "$SRR,$condition,$gene_count_file" >> "$sample_table"
		fi
	done
	
	log_info "Trinity de novo pipeline completed for $fasta_tag."
	log_info "DeSeq2 input files prepared in: $deseq2_dir"
	log_info "Sample table: $sample_table"
}

salmon_saf_pipeline() {
    # Quantify expression using decoy-aware Salmon (Selective Alignment)
	# Fast and accurate 
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
        log_info "Salmon decoy index already exists. Skipping."
    else
        log_info "Building decoy-aware Salmon index for $tag..."
        awk '/^>/{print substr($0,2); next}{next}' "$genome" > "$work/decoys.txt"
        cat "$fasta" "$genome" > "$work/gentrome.fa"
        run_with_time_to_log salmon index \
            -t "$work/gentrome.fa" \
            -d "$work/decoys.txt" \
            -i "$idx_dir" \
            -k 31 -p "$THREADS"
        rm -rf "$work"
    fi

    # --- Quantification per SRR ---
    for SRR in "${rnaseq_list[@]}"; do
        local tdir="$TRIM_DIR_ROOT/$SRR"
        local r1=$(ls "$tdir"/${SRR}*val_1.* 2>/dev/null | head -n1)
        local r2=$(ls "$tdir"/${SRR}*val_2.* 2>/dev/null | head -n1)
        [[ -z "$r1" || -z "$r2" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

        local out_dir="$quant_root/$SRR"
        mkdir -p "$out_dir"
        if [[ -f "$out_dir/quant.sf" ]]; then
            log_info "Quantification already exists for $SRR. Skipping."
            continue
        fi

        log_info "Quantifying expression for $SRR..."
		# Determine if we have paired-end or single-end reads
		if [[ -n "$r2" && -f "$r2" ]]; then
			# Paired-end reads
			log_info "Using paired-end reads for $SRR"
			run_with_time_to_log salmon quant \
			-i "$idx_dir" -l A \
			-1 "$r1" -2 "$r2" \
			-p "$THREADS" \
			--validateMappings \
			--seqBias --gcBias --posBias \
			--rangeFactorizationBins 4 \
			--numBootstraps 100 \
			--numGibbsSamples 20 \
			--thinningFactor 16 \
			$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--libType $([ "$RNA_STRAND_PROTOCOL" = "RF" ] && echo "ISR" || echo "ISF")") \
			-o "$out_dir"
		else
			# Single-end reads
			log_info "Using single-end reads for $SRR"
			run_with_time_to_log salmon quant \
			-i "$idx_dir" -l A \
			-r "$r1" \
			-p "$THREADS" \
			--validateMappings \
			--seqBias --gcBias --posBias \
			--rangeFactorizationBins 4 \
			--numBootstraps 100 \
			--numGibbsSamples 20 \
			--thinningFactor 16 \
			$([ "$RNA_STRAND_PROTOCOL" != "unstranded" ] && echo "--libType $([ "$RNA_STRAND_PROTOCOL" = "RF" ] && echo "SR" || echo "SF")") \
			-o "$out_dir"
		fi
    done

    # --- Merge matrices ---
    log_info "Generating gene and transcript matrices..."
    
    # Check if gene_trans_map exists, create if needed
    local gene_trans_map="${fasta}.gene_trans_map"
    if [[ ! -f "$gene_trans_map" ]]; then
        log_info "Creating gene-transcript mapping file..."
        grep "^>" "$fasta" | sed 's/^>//' | awk '{print $1 "\t" $1}' > "$gene_trans_map"
    fi
    
    run_with_time_to_log abundance_estimates_to_matrix.pl \
        --est_method salmon \
        --gene_trans_map "$gene_trans_map" \
        --out_prefix "$matrix_dir/genes" \
        --name_sample_by_basedir "$quant_root"/*/quant.sf

    log_info "COMPLETED: Salmon-SAF pipeline for $tag"
    log_info "Outputs: $matrix_dir/"
}

bowtie2_rsem_pipeline() {
    # Quantify expression using Bowtie2 + RSEM
	# Reviewer Preferred Method
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

    local tag="$(basename "${fasta%.*}")"
    local rsem_idx="$RSEM_INDEX_ROOT/${tag}_rsem"
    local quant_root="$RSEM_QUANT_ROOT/$tag"
    local matrix_dir="$RSEM_MATRIX_ROOT/$tag"

    mkdir -p "$RSEM_INDEX_ROOT" "$quant_root" "$matrix_dir"

    # --- Prepare reference ---
    if [[ -f "${rsem_idx}.grp" ]]; then
        log_info "RSEM reference already exists. Skipping."
    else
        log_info "Building RSEM reference for $tag..."
        run_with_time_to_log rsem-prepare-reference --bowtie2 "$fasta" "$rsem_idx"
    fi

    # --- Quantify each SRR ---
    for SRR in "${rnaseq_list[@]}"; do
        local tdir="$TRIM_DIR_ROOT/$SRR"
        local r1=$(ls "$tdir"/${SRR}*val_1.* 2>/dev/null | head -n1)
        local r2=$(ls "$tdir"/${SRR}*val_2.* 2>/dev/null | head -n1)
        [[ -z "$r1" || -z "$r2" ]] && { log_warn "Missing trimmed reads for $SRR. Skipping."; continue; }

        local out_dir="$quant_root/$SRR"
        mkdir -p "$out_dir"
        if [[ -f "$out_dir/${SRR}.genes.results" ]]; then
            log_info "RSEM results already exist for $SRR. Skipping."
            continue
        fi

        log_info "Running Bowtie2 + RSEM for $SRR..."
		# Determine if we have paired-end or single-end reads
		if [[ -n "$r2" && -f "$r2" ]]; then
			# Paired-end reads
			log_info "Using paired-end reads for $SRR"
			run_with_time_to_log rsem-calculate-expression \
			--paired-end \
			--bowtie2 \
			--num-threads "$THREADS" \
			"$r1" "$r2" "$rsem_idx" "$out_dir/$SRR"
		else
			# Single-end reads
			log_info "Using single-end reads for $SRR"
			run_with_time_to_log rsem-calculate-expression \
			--bowtie2 \
			--num-threads "$THREADS" \
			"$r1" "$rsem_idx" "$out_dir/$SRR"
		fi
    done

    # --- Merge matrices ---
    log_info "Generating gene and transcript matrices..."
    
    # Check if gene_trans_map exists, create if needed
    local gene_trans_map="${fasta}.gene_trans_map"
    if [[ ! -f "$gene_trans_map" ]]; then
        log_info "Creating gene-transcript mapping file..."
        grep "^>" "$fasta" | sed 's/^>//' | awk '{print $1 "\t" $1}' > "$gene_trans_map"
    fi
    
    run_with_time_to_log abundance_estimates_to_matrix.pl \
        --est_method RSEM \
        --gene_trans_map "$gene_trans_map" \
        --out_prefix "$matrix_dir/genes" \
        --name_sample_by_basedir "$quant_root"/*/*.genes.results

    log_info "COMPLETED: Bowtie2-RSEM pipeline for $tag"
    log_info "Outputs: $matrix_dir/"
}

# ------------------------------------------------------------------------------
# METHOD COMPARISON AND VALIDATION FUNCTIONS
# ------------------------------------------------------------------------------

normalize_expression_data() {
	# Apply normalization to expression matrices
	local matrix_dir="$1"
	local method="$2"
	
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		log_info "Applying TMM normalization for $method..."
		if command -v normalize_matrix.pl >/dev/null 2>&1; then
			run_with_time_to_log normalize_matrix.pl "$matrix_dir/genes.counts.matrix" \
				--est_method "$method" \
				--out_prefix "$matrix_dir/genes.TMM" || \
				log_warn "TMM normalization failed for $method"
		else
			log_warn "normalize_matrix.pl not found. Skipping TMM normalization."
		fi
	fi
}
