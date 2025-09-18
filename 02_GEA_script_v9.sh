#!/bin/bash

set -euo pipefail

chmod +x ./02_GEA_run.sh

# ============================================================================== 
# CONFIGURATION
# ============================================================================== 
#read -p "How many CPU threads to use? [default: 4]: " user_threads
THREADS="${user_threads:-64}"

# ==============================================================================
# INPUT FILES
# ==============================================================================
#All_SmelGIF_GTF_FILE="00_INPUTS/All_SmelDMP_Head_Gene_Name_v4.gtf"
#Eggplant_V4_1_transcripts_function_FASTA_FILE="00_INPUTS/Eggplant_V4_1_transcripts_function.fa"

ALL_FASTA_FILES=(
	"00_INPUTS/TEST.fa"
	"00_INPUTS/All_Control.fa"
	"00_INPUTS/SmelGRF.fasta"
	"00_INPUTS/SmelGIF.fasta"
)

SRR_LIST_PRJNA328564=(
	#SRR3884631 # Fruits Â¯ 6 cm
	#SRR3884664 # Fruits Calyx Stage 2
	#SRR3884653 # Fruits Flesh Stage 2
	#SRR3884677 # Cotyledons
	#SRR3884679 # Pistils /
	#SRR3884597 # Flowers
	SRR3884686 # Buds Â¯ 0\,7 cm /
	#SRR3884687 # Buds, Opened Buds
	#SRR3884689 # Leaves
	#SRR3884690 # Stems
	#SRR3884685 # Radicles
	#SRR3884675 # Roots
)

RAW_DIR_ROOT="00_Raw_Files_and_FastQC/PRJNA328564"
TRIM_DIR_ROOT="01_Trimmed_TrimGalore_Ver/PRJNA328564"
HISAT2_ROOT="02_HISAT2/TrimGalore_Ver"
STRINGTIE_ROOT="03_stringtie/TrimGalore_Ver"
HISAT2_INDEX_DIR="02_HISAT2/index"

# Testing Essentials
#rm -rf $TRIM_DIR_ROOT	
#rm -rf $HISAT2_ROOT
#rm -rf $HISAT2_INDEX_DIR
rm -rf $STRINGTIE_ROOT


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
log_step() { log INFO "===== $* ====="; }

setup_logging() {
	# Set up logging and output redirection
	# Prompt for BAM cleanup option after logging setup
	#read -p "Do you want to keep the BAM files after StringTie assembly? (y/n) [default: y]: " keep_bam_global
	keep_bam_global="${keep_bam_global:-y}"

	mkdir -p "$LOG_DIR"
	#rm -f "$LOG_DIR"/*.log
	echo "Choose log output:"
	echo "1) Console and file"
	echo "2) File only"
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
	# Run a command and log resource usage
    /usr/bin/time -v "$@" >> "$LOG_FILE" 2>&1
}

# ==============================================================================
# TOOLING
# ==============================================================================
require_tools() {
	# Ensure required tools and conda environment are available
	log_info "Ensuring conda environment 'GEA_ENV' and required tools..."
	local env_name="GEA_env"
	local required_cmds=(fasterq-dump trim_galore hisat2 samtools stringtie)
	if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
		log_info "Creating environment '$env_name' ..."
		conda create -y -n "$env_name" -c conda-forge -c bioconda \
			sra-tools trim-galore hisat2 samtools stringtie
	fi
	declare -A pkg_for_cmd=(
		[fasterq-dump]="sra-tools"
		[trim_galore]="trim-galore"
		[hisat2]="hisat2"
		[samtools]="samtools"
		[stringtie]="stringtie"
	)
	local missing_pkgs=()
	for cmd in "${required_cmds[@]}"; do
		if [[ ! -x "${CONDA_PREFIX:-}/bin/$cmd" ]]; then
			missing_pkgs+=( "${pkg_for_cmd[$cmd]}" )
		fi
	done
	if (( ${#missing_pkgs[@]} > 0 )); then
		local uniq=() seen=""
		for p in "${missing_pkgs[@]}"; do
			if [[ " $seen " != *" $p "* ]]; then
				uniq+=( "$p" ); seen+=" $p"
			fi
		done
		log_info "Installing missing packages in '$env_name': ${uniq[*]}"
		conda install -y -n "$env_name" -c conda-forge -c bioconda "${uniq[@]}"
	fi
	local ok=1
	for cmd in "${required_cmds[@]}"; do
		if [[ -x "${CONDA_PREFIX:-}/bin/$cmd" ]]; then
			continue
		elif command -v "$cmd" >/dev/null 2>&1; then
			log_warn "'$cmd' is not in '$env_name' but found in PATH; consider installing it in the env."
		else
			log_error "Required command '$cmd' not found even after installation."
			ok=0
		fi
	done
	if [[ $ok -ne 1 ]]; then
		exit 1
	fi
	log_info "Using conda env: ${CONDA_DEFAULT_ENV:-unknown} (${CONDA_PREFIX:-})"
}

# ==============================================================================
# FUNCTIONS
# ==============================================================================
download_srrs() {
	# Download RNA-seq data for each SRR sample
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
		log_info "-------------------------"
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
	# Build HISAT2 index and align reads for each sample
	# Build HISAT2 index and align SRR reads
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
	mkdir -p "$HISAT2_INDEX_DIR"
	local fasta_base fasta_tag index_prefix
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	index_prefix="$HISAT2_INDEX_DIR/${fasta_tag}_index"
	if ls "${index_prefix}".*.ht2 >/dev/null 2>&1; then
		log_info "HISAT2 index already exists at $HISAT2_INDEX_DIR for $fasta_base. Skipping build."
	else
		log_info "Building HISAT2 index at $HISAT2_INDEX_DIR from $fasta ..."
		run_with_time_to_log hisat2-build -p "${THREADS}" "$fasta" "$index_prefix"
	fi
	export HISAT2_INDEX_DIR
	export HISAT2_INDEX_PREFIX="$index_prefix"
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_ROOT/$fasta_tag/$SRR"
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
		if [[ -f "$bam" && -f "${bam}.bai" && -f "${sam}" ]]; then
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
		run_with_time_to_log \
			samtools sort -@ "${THREADS}" -o "$bam" "$sam"
		run_with_time_to_log \
			samtools index -@ "${THREADS}" "$bam"
		rm -f "$sam"
		log_info "Done aligning $fasta_tag with $SRR."
	done
}

stringtie_assemble() {
	# Run StringTie assembly for each sample
	local fasta=""
	local rnaseq_list=()
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
	local fasta_base fasta_tag
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	for SRR in "${rnaseq_list[@]}"; do
		local stringtie_chr_DIR="$STRINGTIE_ROOT/$fasta_tag/$SRR"
		local HISAT2_DIR="$HISAT2_ROOT/$fasta_tag/$SRR"
		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		local out_gtf="$stringtie_chr_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf"
		mkdir -p "$stringtie_chr_DIR"
		if [[ -f "$out_gtf" ]]; then
			log_info "StringTie assembly exists for $fasta_tag/$SRR. Skipping."
			continue
		fi
		run_with_time_to_log \
			stringtie \
				-p "$THREADS" "$bam" \
				-o "$out_gtf" \
				-A "$stringtie_chr_DIR/${SRR}_${fasta_tag}_gene_abundances_de_novo_v1.tsv"
	done
}

stringtie_merge() {
	# Merge StringTie assemblies for each FASTA
	local fasta=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)
				fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do shift; done;;
			*)
				shift;;
		esac
	done
	local fasta_base fasta_tag
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	local MERGELIST="$STRINGTIE_ROOT/mergelist_${fasta_tag}.txt"
	mkdir -p "$STRINGTIE_ROOT"
	find "$(realpath "$STRINGTIE_ROOT/$fasta_tag")" -type f -name "*${fasta_tag}_trimmed_mapped_sorted_stringtie_assembled_de_novo.gtf" > "$MERGELIST"
	local merged_gtf="$STRINGTIE_ROOT/merged_transcripts_de_novo_${fasta_tag}.gtf"
	if [[ -s "$merged_gtf" ]]; then
		log_info "Merged transcripts GTF for $fasta_tag already exists. Skipping merge."
		return
	fi
	log_info "Merging transcripts listed in $MERGELIST for $fasta_tag ..."
	run_with_time_to_log \
		stringtie --merge -p "$THREADS" -o "$merged_gtf" "$MERGELIST"
}

stringtie_quantify() {
	# Quantify expression using merged GTF for each sample
	local fasta=""
	local rnaseq_list=()
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
	local fasta_base fasta_tag
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	local merged_gtf="$STRINGTIE_ROOT/merged_transcripts_de_novo_${fasta_tag}.gtf"
	if [[ ! -s "$merged_gtf" ]]; then
		log_error "Merged GTF not found at $merged_gtf for $fasta_tag. Run stringtie_merge first."
		return
	fi
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_ROOT/$fasta_tag/$SRR"
		local STRINGTIE_DIR="$STRINGTIE_ROOT/$fasta_tag/$SRR"
		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trimmed_mapped_sorted.bam"
		mkdir -p "$STRINGTIE_DIR"
		local out_gtf="$STRINGTIE_DIR/${SRR}_${fasta_tag}_expression_estimates_de_novo_v2.gtf"
		if [[ -f "$out_gtf" ]]; then
			log_info "StringTie quantify exists for $fasta_tag/$SRR. Skipping."
			continue
		fi
		run_with_time_to_log stringtie -p "$THREADS" -e -B "$bam" \
			-G "$merged_gtf" \
			-o "$out_gtf" \
			-A "$STRINGTIE_DIR/${SRR}_${fasta_tag}_gene_abundances_de_novo_v2.tsv" \
			-C "$STRINGTIE_DIR/${SRR}_${fasta_tag}_transcripts_with_coverage_de_novo.gtf"
	done
}

cleanup_bam_files() {
	# Optionally clean up BAM files after processing
	# Use the global keep_bam_global variable set after logging setup
	local keep_bam="$keep_bam_global"
	if [[ ! "$keep_bam" =~ ^[yYnN]$ ]]; then
		log_error "Invalid input for BAM cleanup: $keep_bam. Please enter 'y' or 'n'."
		exit 1
	fi
	keep_bam="${keep_bam,,}"

	if [[ "$keep_bam" == "n" ]]; then
		log_info "Deleting BAM files as per user request..."
		run_with_time_to_log find "$HISAT2_ROOT" -type f -name "*.bam" -exec rm -f {} +
		run_with_time_to_log find "$HISAT2_ROOT" -type f -name "*.bam.bai" -exec rm -f {} +
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
		require_tools

		log_step "STEP 00: Download SRR files"
		download_srrs "${rnaseq_list[@]}"

		log_step "STEP 01: Trimming"
		trim_reads "${rnaseq_list[@]}"

		log_step "STEP 02: Build HISAT2 index and align"
		hisat2_index_align_sort --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"

		log_step "STEP 03: StringTie assembly, merge, and quantification"
		stringtie_assemble --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		stringtie_merge --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		stringtie_quantify --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"

		cleanup_bam_files

		end_time=$(date +%s)
		log_step "Final timing"
		log_info "Script ended at: $(date -d @$end_time)"
		elapsed=$((end_time - start_time))
		formatted_elapsed=$(date -u -d @${elapsed} +%H:%M:%S)
		log_info "Elapsed time: $formatted_elapsed"
}

# Run the pipeline for each FASTA input and SRR list.
for fasta_input in "${ALL_FASTA_FILES[@]}"; do
	run_all --FASTA "$fasta_input" --RNASEQ_LIST "${SRR_LIST_PRJNA328564[@]}"
done

# Zip all the content of this folder: STRINGTIE_ROOT="03_stringtie/TrimGalore_Ver"
#tar -czvf 03_stringtie_TrimGalore_Ver_$(date +%Y%m%d_%H%M%S).tar.gz 03_stringtie/TrimGalore_Ver
