#!/bin/bash

set -euo pipefail

# Chmod and Converts the file to Unix line endings
#chmod +x ./*.sh   # Ensure the run script is executable
#dos2unix ./*.sh   # Convert Windows line endings to Unix
# ============================================================================== 
# CONFIGURATION
# ============================================================================== 
#read -p "How many CPU threads to use? [default: 4]: " user_threads
THREADS="${user_threads:-64}"  # Number of threads to use for parallel operations

# ==============================================================================
# INPUT FILES
# ==============================================================================
#All_SmelGIF_GTF_FILE="0_INPUT_FASTAs/All_SmelDMP_Head_Gene_Name_v4.gtf"
#Eggplant_V4_1_transcripts_function_FASTA_FILE="0_INPUT_FASTAs/Eggplant_V4_1_transcripts_function.fa"

ALL_FASTA_FILES=(
	# List of FASTA files to process
	"0_INPUT_FASTAs/TEST.fasta"
	#"0_INPUT_FASTAs/SmelGIF_with_Cell_Cycle_Control_genes.fasta"
	#"0_INPUT_FASTAs/SmelDMP_CDS_Control_Best.fasta"
	#"0_INPUT_FASTAs/SmelGIF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/SmelGRF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/SmelGRF-GIF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/Control_Genes_Puta.fasta"
	#"0_INPUT_FASTAs/SmelGRF_with_Cell_Cycle_Control_genes.fasta"
)

SRR_LIST_PRJNA328564=(
	# List of SRR sample IDs to process
	#SRR3884664 # Fruits Calyx Stage 2
	SRR3884653 # Fruits Flesh Stage 2
	SRR3884631 # Fruits 6 cm
	#SRR3884677 # Cotyledons
	#SRR3884679 # Pistils
	SRR3884597 # Flowers
	SRR3884686 # Buds 0.7 cm
	SRR3884687 # Buds, Opened Buds
	SRR3884689 # Leaves
	SRR3884690 # Stems
	#SRR3884685 # Radicles
	SRR3884675 # Roots
)

SRR_LIST_SAMN28540077=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540077&o=acc_s%3Aa&s=SRR20722234,SRR20722233,SRR20722232,SRR20722230,SRR20722225,SRR20722226,SRR20722227,SRR20722228,SRR20722229
	SRR2072232	# mature_fruits
	SRR20722226	# young_fruits
	SRR20722234	# flowers

	#SRR20722228	# sepals
	SRR20722230	# mature_leaves
	#SRR20722225	# young_leaves
	SRR20722227	# stems
	#SRR20722233	# leaf_buds
	SRR20722229	# roots
)

SRR_LIST_SAMN28540068=(
	#Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540068&o=acc_s%3Aa
	# Solanum virginianum. wild eggplant or Thai green eggplant, but it is not the same species as the common eggplant (Solanum melongena).
)

OTHER_SRR_LIST=(
	# Possible Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP390977&o=acc_s%3Aa
	SRR34564302	# Fruits (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR34564302&display=metadata)
	SRR34848077 # Leaves (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR34848077&display=metadata)
		# Cotyledons (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884677&display=metadata)
	SRR3479277 # Pistil (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3479277&display=metadata)
	SRR3884597 # Flowers (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884597&display=metadata)
		SRR3884686 # Buds 0.7 cm
	 # Buds, Opened Buds
	 # Leaves
	 # Stems
	 # Radicles
	SRR20722229 # Roots (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR20722229&display=metadata)
	# Add other SRR IDs here if needed
)

SRR_COMBINED_LIST=(
	"${SRR_LIST_PRJNA328564[@]}"
	"${SRR_LIST_SAMN28540077[@]}"
	"${OTHER_SRR_LIST[@]}"
)

for SRR in "${SRR_COMBINED_LIST[@]}"; do
	echo "$SRR"
done

#: << 'OFF'
RAW_DIR_ROOT="1_RAW_SRR"
TRIM_DIR_ROOT="2_TRIMMED_SRR"
FASTQC_ROOT="3_FastQC"

HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/4_HISAT2_RESULTS"
HISAT2_DE_NOVO_INDEX_DIR="4b_Method_2_HISAT2_De_Novo/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/5_stringtie/a_Method_2_Results"

mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" "$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT"

# SWTICHES. 
RUN_DOWNLOAD_SRR=TRUE
RUN_TRIMMING=TRUE
RUN_HISAT2_INDEX_ALIGN_SORT=TRUE
RUN_STRINGTIE_ASSEMBLE_MERGE_QUANTIFY=TRUE
RUN_CLEANUP_BAM=TRUE

#HISAT2_DE_NOVO_ROOT="02_HISAT2/"
#HISAT2_DE_NOVO_INDEX_DIR="02_HISAT2/index"
# ==============================================================================
# CLEANUP OPTIONS and Testing Essentials
#rm -rf $TRIM_DIR_ROOT
#rm -rf $HISAT2_DE_NOVO_ROOT      # Remove previous HISAT2 results
#rm -rf $HISAT2_DE_NOVO_INDEX_DIR # Remove previous HISAT2 index
#rm -rf $STRINGTIE_HISAT2_DE_NOVO_ROOT   # Remove previous StringTie results

# ==============================================================================
# LOGGING
# ==============================================================================
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_full_log.log}"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; }
log_error() { log ERROR "$@"; }
log_step() { log INFO "=============== $* ==============="; }

setup_logging() {
	# Set up logging and output redirection
	# Prompt for BAM cleanup option after logging setup
	#read -p "Do you want to keep the BAM files after StringTie assembly? (y/n) [default: y]: " keep_bam_global
	keep_bam_global="${keep_bam_global:-n}"

	mkdir -p "$LOG_DIR"
	#rm -f "$LOG_DIR"/*.log
	#echo "Choose log output:"
	#echo "1) Console and file"
	#echo "2) File only"
	#read -p "Enter 1 or 2 [default: 1]: " log_choice
	log_choice="${log_choice:-1}"
	if [[ "$log_choice" == "2" ]]; then
		exec >"$LOG_FILE" 2>&1
	else
		exec > >(tee -a "$LOG_FILE") 2>&1
	fi
	log_info "Logging to: $LOG_FILE"


}

trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
trap 'log_info "Script finished. See log: $LOG_FILE"' EXIT

run_with_time_to_log() {
	# Run a command and log resource usage (tracks time and memory)
	/usr/bin/time -v "$@" >> "$LOG_FILE" 2>&1
}

# ==============================================================================
# FUNCTIONS
# ==============================================================================
download_srrs() {
	# Download RNA-seq data for each SRR sample using fasterq-dump
	local SRR_LIST=("$@")
	if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
		SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	fi
	for SRR in "${SRR_LIST[@]}"; do
		local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
		log_info "Working on $SRR."
		mkdir -p "$raw_files_DIR"
		if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]] \
		   || [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
			log_info "Raw fastq files for $SRR already exist. Skipping download."
		else
			fasterq-dump --split-files "$SRR" -O "$raw_files_DIR"
		fi
		log_info "Done working on $SRR."
		log_info "--------------------------------------------------"
	done
}

trim_reads() {
	# Trim reads for each sample using Trim Galore
	local SRR_LIST=("$@")
	if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
		SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
	fi
	for SRR in "${SRR_LIST[@]}"; do
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
		local raw1 raw2 trimmed1 trimmed2
		mkdir -p "$TrimGalore_DIR"
		if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq"
			raw2="$raw_files_DIR/${SRR}_2.fastq"
		elif [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
			raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
		else
			log_warn "Raw FASTQ for $SRR not found in $raw_files_DIR; skipping."
			continue
		fi
		trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
		trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		if { [[ -f "$trimmed1" && -f "$trimmed2" ]]; } || { [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
			log_info "Trimmed files for $SRR already exist. Skipping trimming."
			continue
		fi
		log_info "Trimming $SRR..."
		run_with_time_to_log trim_galore --cores "${THREADS}" --paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"
		log_info "Done trimming $SRR."
	done
}

hisat2_index_align_sort() {
	# Build HISAT2 index from FASTA and align trimmed reads for each sample
	local fasta="" rnaseq_list=()
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)      fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1")
					shift
				done
				;;
			*)         echo "Unknown option: $1" >&2; return 1;;
		esac
	done
	if [[ -z "$fasta" ]]; then
		log_error "No FASTA file specified. Use --FASTA <fasta_file>."
		return 1
	fi
	if [[ ! -f "$fasta" ]]; then
		log_error "Fasta file '$fasta' not found."
		return 1
	fi
	if [[ ${#rnaseq_list[@]} -eq 0 ]]; then
		rnaseq_list=("${SRR_LIST_PRJNA328564[@]}")
	fi
	mkdir -p "$HISAT2_DE_NOVO_INDEX_DIR"
	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_DE_NOVO_INDEX_DIR/${fasta_tag}_index"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "HISAT2 index already exists at $HISAT2_DE_NOVO_INDEX_DIR for $fasta_base. Skipping build."
	else
		log_info "Building HISAT2 index at $HISAT2_DE_NOVO_INDEX_DIR from $fasta ..."
		run_with_time_to_log hisat2-build -p "${THREADS}" "$fasta" "$index_prefix"
	fi
	export HISAT2_DE_NOVO_INDEX_DIR
	export HISAT2_INDEX_PREFIX="$index_prefix"
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
		local trimmed1="" trimmed2=""
		mkdir -p "$HISAT2_DIR"
		# Check for uncompressed trimmed FASTQ files
		if [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"
		# Check for compressed trimmed FASTQ files
		elif [[ -f "$TrimGalore_DIR/${SRR}_1_val_1.fq.gz" && -f "$TrimGalore_DIR/${SRR}_2_val_2.fq.gz" ]]; then
			trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq.gz"
			trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq.gz"
		# Fallback: match any file pattern for trimmed FASTQ
		elif compgen -G "$TrimGalore_DIR/${SRR}*val_1.*" >/dev/null 2>&1 && compgen -G "$TrimGalore_DIR/${SRR}*val_2.*" >/dev/null 2>&1; then
			local _a=( $TrimGalore_DIR/${SRR}*val_1.* ); trimmed1="${_a[0]}"
			local _b=( $TrimGalore_DIR/${SRR}*val_2.* ); trimmed2="${_b[0]}"
		fi
		if [[ -z "$trimmed1" || -z "$trimmed2" ]]; then
			log_warn "Trimmed FASTQ for $SRR not found in $TrimGalore_DIR; skipping alignment."
			continue
		fi
		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped.sam"
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "BAM for $SRR and $fasta_tag already exists. Skipping alignment."
			continue
		fi
		log_info "Aligning $fasta_tag with $SRR..."
		run_with_time_to_log \
			hisat2 -p "${THREADS}" --dta \
				-x "$index_prefix" \
				-1 "$trimmed1" \
				-2 "$trimmed2" \
				-S "$sam"
		log_info "Converting SAM to sorted BAM for $fasta_tag with $SRR..."
		run_with_time_to_log \
			samtools sort -@ "${THREADS}" -o "$bam" "$sam"
		run_with_time_to_log \
			samtools index -@ "${THREADS}" "$bam"
		rm -f "$sam"
		log_info "Done aligning $fasta_tag with $SRR."
	done
}

stringtie_assemble() {
	# Assemble transcripts for each sample using StringTie
	local fasta rnaseq_list=()

	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--RNASEQ_LIST) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
			*) shift;;
		esac
	done

	local fasta_tag="${fasta##*/}"; fasta_tag="${fasta_tag%.*}"
	for SRR in "${rnaseq_list[@]}"; do
		local out_dir="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local bam="$HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local out_gtf="$out_dir/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"

		mkdir -p "$out_dir"
		[[ -f "$out_gtf" ]] && log_info "Assembly exists for $fasta_tag/$SRR. Skipping." && continue

		log_info "Assembling transcripts for $fasta_tag with $SRR ..."
		run_with_time_to_log \
			stringtie -p "$THREADS" "$bam" \
				-o "$out_gtf" \
				-A "$out_dir/${SRR}_${fasta_tag}_gene_abundances_de_novo_v1.tsv"
	done
}

stringtie_merge() {
	# Merge all StringTie assemblies for a given FASTA into a single GTF
	local fasta
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--RNASEQ_LIST) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do shift; done;;
			*) shift;;
		esac
	done

	local fasta_tag="${fasta##*/}"; fasta_tag="${fasta_tag%.*}"
	local MERGELIST="$STRINGTIE_HISAT2_DE_NOVO_ROOT/mergelist_${fasta_tag}.txt"
	mkdir -p "$STRINGTIE_HISAT2_DE_NOVO_ROOT"
	find "$(realpath "$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag")" -type f \
		-name "*${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf" > "$MERGELIST"

	local merged_gtf="$STRINGTIE_HISAT2_DE_NOVO_ROOT/merged_transcripts_de_novo_${fasta_tag}.gtf"
	[[ -s "$merged_gtf" ]] && log_info "Merged GTF for $fasta_tag already exists. Skipping." && return

	log_info "Merging StringTie assemblies for $fasta_tag into $merged_gtf ..."
	run_with_time_to_log \
		stringtie --merge -p "$THREADS" -o "$merged_gtf" "$MERGELIST"
}

stringtie_quantify() {
	# Quantify gene expression for each sample using merged GTF
	local fasta rnaseq_list=()

	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--RNASEQ_LIST) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do rnaseq_list+=("$1"); shift; done;;
			*) shift;;
		esac
	done

	local fasta_tag="${fasta##*/}"; fasta_tag="${fasta_tag%.*}"
	local merged_gtf="$STRINGTIE_HISAT2_DE_NOVO_ROOT/merged_transcripts_de_novo_${fasta_tag}.gtf"
	[[ ! -s "$merged_gtf" ]] && log_error "Merged GTF not found at $merged_gtf. Run stringtie_merge first." && return 1

	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local STRINGTIE_DIR="$STRINGTIE_HISAT2_DE_NOVO_ROOT/$fasta_tag/$SRR"
		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local out_gtf="$STRINGTIE_DIR/${SRR}_${fasta_tag}_expression_estimates_de_novo_v2.gtf"

		mkdir -p "$STRINGTIE_DIR"
		[[ -f "$out_gtf" ]] && log_info "Quantification exists for $fasta_tag/$SRR. Skipping." && continue

		log_info "Quantifying expression for $fasta_tag with $SRR using merged GTF ..."
		run_with_time_to_log \
			stringtie -p "$THREADS" -e -B "$bam" \
				-G "$merged_gtf" \
				-o "$out_gtf" \
				-A "$STRINGTIE_DIR/${SRR}_${fasta_tag}_gene_abundances_de_novo_v2.tsv" \
				-C "$STRINGTIE_DIR/${SRR}_${fasta_tag}_transcripts_with_coverage_de_novo.gtf"
	done
}

cleanup_bam_files() {
	# Optionally clean up BAM files after processing
	# Accepts --FASTA argument for per-fasta cleanup
	local fasta=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)
				fasta="$2"; shift 2;;
			*)
				shift;;
		esac
	done
	local keep_bam="$keep_bam_global"
	if [[ ! "$keep_bam" =~ ^[yYnN]$ ]]; then
		log_error "Invalid input for BAM cleanup: $keep_bam. Please enter 'y' or 'n'."
		exit 1
	fi
	keep_bam="${keep_bam,,}"

	local fasta_base fasta_tag
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	local bam_dir="$HISAT2_DE_NOVO_ROOT/$fasta_tag"

	if [[ "$keep_bam" == "n" ]]; then
		log_info "Deleting BAM files as per user request..."
		run_with_time_to_log find "$bam_dir" -type f -name "*.bam" -exec rm -f {} +
		run_with_time_to_log find "$bam_dir" -type f -name "*.bam.bai" -exec rm -f {} +
		log_info "BAM files deleted."
	else
		log_info "Keeping BAM files as per user request."
	fi
}

# ==============================================================================
# RUN / ENTRYPOINT
# ==============================================================================
run_all() {
	# Main pipeline entrypoint: runs all steps for each FASTA and RNA-seq list
	# Steps: logging, tool check, download, trim, align, assemble, merge, quantify, cleanup
		local fasta=""
		local rnaseq_list=()
		# Parse arguments
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
					shift;;
			esac
		done

		local start_time end_time elapsed formatted_elapsed
		start_time=$(date +%s)
		setup_logging
		log_step "Script started at: $(date -d @$start_time)"

		if [[ $RUN_DOWNLOAD_SRR == "TRUE" ]]; then
			log_step "STEP 00: Download SRR files"
			download_srrs "${rnaseq_list[@]}"
		fi	

		if [[ $RUN_TRIMMING == "TRUE" ]]; then
			log_step "STEP 01: Trimming"
			trim_reads "${rnaseq_list[@]}"
		fi

		if [[ $RUN_HISAT2_INDEX_ALIGN_SORT == "TRUE" ]]; then
			log_step "STEP 02: Build HISAT2 index and align"
			hisat2_index_align_sort --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		fi

		if [[ $RUN_STRINGTIE_ASSEMBLE_MERGE_QUANTIFY == "TRUE" ]]; then
			log_step "STEP 03: StringTie assembly, merge, and quantification"
			stringtie_assemble --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			stringtie_merge --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			stringtie_quantify --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		fi

		if [[ $RUN_CLEANUP_BAM == "TRUE" ]]; then
			log_step "STEP 04: Optional BAM cleanup"
			cleanup_bam_files --FASTA "$fasta"
		else
			log_info "Skipping BAM cleanup as per configuration."
		fi
		
		end_time=$(date +%s)
		log_step "Final timing"
		log_info "Script ended at: $(date -d @$end_time)"
		elapsed=$((end_time - start_time))
		formatted_elapsed=$(date -u -d @${elapsed} +%H:%M:%S)
		log_info "Elapsed time: $formatted_elapsed"	
}

# Run the pipeline for each FASTA input and SRR list.
for fasta_input in "${ALL_FASTA_FILES[@]}"; do
	# Run the pipeline for each FASTA file and all SRR samples
	#run_all --FASTA "$fasta_input" --RNASEQ_LIST "${SRR_LIST_PRJNA328564[@]}"
	run_all --FASTA "$fasta_input" --RNASEQ_LIST "${SRR_COMBINED_LIST[@]}"
done

# Zip all the content of this folder: STRINGTIE_HISAT2_DE_NOVO_ROOT="03_stringtie/"
#
#tar -czvf 03_stringtie__$(date +%Y%m%d_%H%M%S).tar.gz 03_stringtie/  # Archive all StringTie results for sharing or backup

#OFF