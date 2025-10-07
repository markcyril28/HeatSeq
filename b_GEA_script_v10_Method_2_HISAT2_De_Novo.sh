#!/bin/bash

set -euo pipefail
#conda install -c conda-forge aria2 -y
#conda install -c bioconda -c conda-forge parallel-fastq-dump -y

# ============================================================================== 
# CONFIGURATION and SWITCHES
# ============================================================================== 
THREADS=8  # Number of threads to use for parallel operations
JOBS=3      # Number of parallel jobs for GNU Parallel 
RUN_DOWNLOAD_and_TRIM_SRR=TRUE
RUN_HISAT2_INDEX_ALIGN_SORT_STRINGTIE=FALSE

# ==============================================================================
# INPUT FILES
# ==============================================================================

#All_SmelGIF_GTF_FILE="0_INPUT_FASTAs/All_SmelDMP_Head_Gene_Name_v4.gtf"
#Eggplant_V4_1_transcripts_function_FASTA_FILE="0_INPUT_FASTAs/Eggplant_V4_1_transcripts_function.fa"

ALL_FASTA_FILES=(
	# List of FASTA files to process
	"0_INPUT_FASTAs/All_Smel_Genes.fasta"
	#"0_INPUT_FASTAs/Eggplant_V4.1_transcripts.function.fa"
	#"0_INPUT_FASTAs/TEST.fasta"
	#"0_INPUT_FASTAs/SmelGIF_with_Cell_Cycle_Control_genes.fasta"
	#"0_INPUT_FASTAs/SmelDMP_CDS_Control_Best.fasta"
	#"0_INPUT_FASTAs/SmelGIF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/SmelGRF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/SmelGRF-GIF_with_Best_Control_Cyclo.fasta"
	#"0_INPUT_FASTAs/Control_Genes_Puta.fasta"
	#"0_INPUT_FASTAs/SmelGRF_with_Cell_Cycle_Control_genes.fasta"
)

# SRA Run Selector Tool link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa

SRR_LIST_PRJNA328564=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa
	#SRR3884653 	# Fruits Flesh Stage 2 #SRR3884664 # Fruits Calyx Stage 2
	SRR3884631 	# Fruits 6 cm #SRR3884677 # Cotyledons #SRR3884679 # Pistils
	SRR3884597 	# Flowers
	SRR3884686 	# Buds 0.7 cm #SRR3884687 	# Buds, Opened Buds
	SRR3884689 	# Leaves
	SRR3884690 	# Stems
	SRR3884675 	# Roots #SRR3884685 # Radicles
)

SRR_LIST_SAMN28540077=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540077&o=acc_s%3Aa&s=SRR20722234,SRR20722233,SRR20722232,SRR20722230,SRR20722225,SRR20722226,SRR20722227,SRR20722228,SRR20722229
	SRR2072232	# mature_fruits #SRR20722226	# young_fruits
	SRR20722234	# flowers #SRR20722228	# sepals
	SRR21010466 # Buds, Nonparthenocarpy ID: PRJNA865018
	SRR20722230	# mature_leaves #SRR20722233	# leaf_buds
	SRR20722227	# stems
	SRR20722229	# roots
)

SRR_LIST_SAMN28540068=(
	#Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540068&o=acc_s%3Aa
	SRR20722387 # mature_fruits
	SRR23909863 # Fully Develop (FD) Flower ID: PRJNA941250
	SRR20722297 # flower_buds #SRR20722385 # sepals
	SRR20722386 # mature_leaves #SRR20722383 # young_leaves #SRR20722296 # leaf_buds
	SRR20722384 # stems
	SRR31755282 # Roots (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP552204&o=acc_s%3Aa)
)

# A Good Dataset for SmelDMP GEA: 
# 	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=%20%20PRJNA865018&o=acc_s%3Aa PRJNA865018
#	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA941250&o=acc_s%3Aa PRJNA941250 # Buds, Opened Buds
OTHER_SRR_LIST=(
	# Possible Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP390977&o=acc_s%3Aa
	SRR34564302	# Fruits (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR34564302&display=metadata)
	SRR34848077 # Leaves (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR34848077&display=metadata)
	PRJNA613773 # Leaf (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA613773&o=acc_s%3Aa)
		# Cotyledons (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884677&display=metadata)
	SRR3479277 # Pistil (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3479277&display=metadata)
	SRR3884597 # Flowers (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884597&display=metadata)

	PRJNA341784 # Flower buds lang (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA341784&o=acc_s%3Aa)
	PRJNA477924 # Leaf and Root (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA477924&o=acc_s%3Aa)
	 # Leaves
	 # Stems
	 # Radicles
	
	# Add other SRR IDs here if needed
)

SRR_COMBINED_LIST=(
	"${SRR_LIST_PRJNA328564[@]}"
	"${SRR_LIST_SAMN28540077[@]}"
	"${SRR_LIST_SAMN28540068[@]}"
	#"${OTHER_SRR_LIST[@]}"
)

#:<< 'OFF'
RAW_DIR_ROOT="1_RAW_SRR"
TRIM_DIR_ROOT="2_TRIMMED_SRR"
FASTQC_ROOT="3_FastQC"

HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/5_stringtie_WD/a_Method_2_RAW_RESULTs"

TRINITY_DE_NOVO_ROOT="4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4c_Method_3_Trinity_De_Novo/5_stringtie_WD/a_Method_3_RAW_RESULTs"

mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$TRINITY_DE_NOVO_ROOT" "$STRINGTIE_TRINITY_DE_NOVO_ROOT"

# ==============================================================================
# CLEANUP OPTIONS and Testing Essentials
rm -rf "$RAW_DIR_ROOT"               	# Remove previous raw SRR files
#rm -rf $FASTQC_ROOT               		 # Remove previous FastQC results
#rm -rf $HISAT2_DE_NOVO_ROOT      		 # Remove previous HISAT2 results
#rm -rf $HISAT2_DE_NOVO_INDEX_DIR 		 # Remove previous HISAT2 index
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

        # --------------------------------------------------
        # Check if trimmed files already exist
        # --------------------------------------------------
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Trimmed files for $SRR already exist. Skipping download and trimming."
            log_info "--------------------------------------------------"
            continue
        fi

        # --------------------------------------------------
        # Check for existing raw FASTQ files
        # --------------------------------------------------
        if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
            raw1="$raw_files_DIR/${SRR}_1.fastq"
            raw2="$raw_files_DIR/${SRR}_2.fastq"
            log_info "Raw FASTQ files for $SRR already exist. Skipping download."
        elif [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
            raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
            raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
            log_info "Compressed raw FASTQ files for $SRR already exist. Skipping download."
        else
            # --------------------------------------------------
            # Download with prefetch + fasterq-dump
            # --------------------------------------------------
            log_info "Downloading $SRR..."
            prefetch "$SRR" --output-directory "$raw_files_DIR"
            
            # Important: include .sra file path explicitly
            fasterq-dump --split-files --threads "${THREADS}" \
                "$raw_files_DIR/$SRR/$SRR.sra" -O "$raw_files_DIR"

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

        # --------------------------------------------------
        # Trim reads using Trim Galore
        # --------------------------------------------------
        log_info "Trimming $SRR..."
        run_with_time_to_log trim_galore --cores "${THREADS}" \
            --paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

        log_info "Done working on $SRR."
        log_info "--------------------------------------------------"

        # --------------------------------------------------
        # Delete raw FASTQ files after successful trimming
        # --------------------------------------------------
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Deleting raw FASTQ files for $SRR after successful trimming..."
            rm -f "$raw1" "$raw2"
            log_info "Raw files deleted for $SRR."
        fi
    done
}

download_and_trim_srrs_parallel() {
    # ============================================
    # Parallel Download and Trimming of SRR Samples
    # ============================================

    local SRR_LIST=("$@")
    if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
        SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
    fi

    export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE
    export -f timestamp log log_info log_warn log_error run_with_time_to_log

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

            raw1="$raw_files_DIR/${SRR}_1.fastq"
            raw2="$raw_files_DIR/${SRR}_2.fastq"
        fi

        # Trim reads
        log_info "Trimming $SRR..."
        run_with_time_to_log trim_galore --cores "${THREADS}" \
            --paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

        # Clean up raw files
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Trimming successful; removing raw FASTQs for $SRR..."
            rm -f "$raw1" "$raw2"
        fi

        log_info "Done with $SRR."
        log_info "--------------------------------------------------"
    }

    export -f _process_single_srr

    # ============================================
    # Run jobs in parallel
    # ============================================
    log_info "Starting parallel download and trimming..."
    printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr {}
    log_info "All parallel SRR jobs completed."
}

download_and_trim_srrs_parallel_fastqdump() {
    # =========================================================
    # Parallel SRA Download and Trimming using parallel-fastq-dump
    # =========================================================

    local SRR_LIST=("$@")
    if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
        SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
    fi

    export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE
    export -f timestamp log log_info log_warn log_error run_with_time_to_log

    # ---------------------------------------------------------
    # Inner worker function for one SRR
    # ---------------------------------------------------------
    _process_single_srr_fastqdump() {
        local SRR="$1"
        local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
        local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
        local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
        local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"

        log_info "Working on $SRR..."
        mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

        # Skip if trimmed files already exist
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Trimmed files for $SRR already exist. Skipping."
            log_info "--------------------------------------------------"
            return 0
        fi

        # ---------------------------------------------------------
        # Download and extract FASTQs using parallel-fastq-dump
        # ---------------------------------------------------------
        log_info "Downloading and converting $SRR with parallel-fastq-dump..."
        parallel-fastq-dump --sra-id "$SRR" \
            --threads "${THREADS}" \
            --split-files \
            --outdir "$raw_files_DIR" \
            --gzip \
            --tmpdir "$raw_files_DIR/tmp" || {
                log_warn "parallel-fastq-dump failed for $SRR"
                return 1
            }

        # ---------------------------------------------------------
        # Verify FASTQ presence
        # ---------------------------------------------------------
        if [[ ! -f "$raw_files_DIR/${SRR}_1.fastq.gz" || ! -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
            log_warn "FASTQs not found after parallel-fastq-dump for $SRR"
            return 1
        fi

        # ---------------------------------------------------------
        # Trim reads with Trim Galore
        # ---------------------------------------------------------
        log_info "Trimming $SRR..."
        run_with_time_to_log trim_galore --cores "${THREADS}" \
            --paired "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz" \
            --output_dir "$TrimGalore_DIR"

        log_info "Done working on $SRR."
        log_info "--------------------------------------------------"

        # ---------------------------------------------------------
        # Cleanup raw files if trimming successful
        # ---------------------------------------------------------
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "Trimming successful. Deleting raw FASTQs for $SRR..."
            rm -f "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz"
        fi
    }

    export -f _process_single_srr_fastqdump

    # ---------------------------------------------------------
    # Batch execution across SRRs in parallel
    # ---------------------------------------------------------
    log_info "Starting batch processing with parallel-fastq-dump..."
    printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr_fastqdump {}
    log_info "All SRRs processed with parallel-fastq-dump."
}

download_and_trim_srrs_wget_parallel() {
    # ============================================
    # Parallel ENA FASTQ Download (wget) + Trim Galore
    # ============================================

    local SRR_LIST=("$@")
    if [[ ${#SRR_LIST[@]} -eq 0 ]]; then
        SRR_LIST=("${SRR_LIST_PRJNA328564[@]}")
    fi

    export RAW_DIR_ROOT TRIM_DIR_ROOT THREADS LOG_FILE
    export -f timestamp log log_info log_warn log_error run_with_time_to_log

    # --------------------------------------------
    # Worker function for one SRR
    # --------------------------------------------
    _process_single_srr_wget() {
        local SRR="$1"
        local raw_files_DIR="$RAW_DIR_ROOT/$SRR"
        local TrimGalore_DIR="$TRIM_DIR_ROOT/$SRR"
        local trimmed1="$TrimGalore_DIR/${SRR}_1_val_1.fq"
        local trimmed2="$TrimGalore_DIR/${SRR}_2_val_2.fq"

        log_info "▶ Working on $SRR..."
        mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

        # Skip if trimming already completed
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "✔ Trimmed files exist for $SRR — skipping."
            return 0
        fi

        # ----------------------------------------
        # Fetch ENA FASTQ URLs dynamically
        # ----------------------------------------
        local ena_links fq1 fq2
        ena_links=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" \
                     | tail -n +2)

        if [[ -z "$ena_links" ]]; then
            log_warn "⚠ No ENA FASTQ links found for $SRR. Skipping."
            return 1
        fi

        fq1="ftp://$(echo "$ena_links" | cut -f1 -d';')"
        fq2="ftp://$(echo "$ena_links" | cut -f2 -d';')"

        log_info "⬇ Downloading FASTQs for $SRR using wget..."
        wget -q -c -P "$raw_files_DIR" "$fq1" "$fq2"

        # Verify successful download
        if [[ ! -s "$raw_files_DIR/${SRR}_1.fastq.gz" || ! -s "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
            log_warn "❌ FASTQ download failed for $SRR. Check ENA mirrors."
            return 1
        fi

        # ----------------------------------------
        # Trim Galore (low RAM, single-thread)
        # ----------------------------------------
        log_info "✂ Trimming adapters for $SRR..."
        run_with_time_to_log trim_galore --cores 1 \
            --paired "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz" \
            --output_dir "$TrimGalore_DIR"

        # ----------------------------------------
        # Cleanup if trimming successful
        # ----------------------------------------
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "🧹 Removing raw FASTQs for $SRR..."
            rm -f "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz"
        fi

        log_info "✅ Finished $SRR."
        log_info "--------------------------------------------------"
    }

    export -f _process_single_srr_wget

    # --------------------------------------------
    # Run jobs in parallel
    # --------------------------------------------
    log_info "🚀 Starting parallel ENA downloads and trimming..."
    printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr_wget {}
    log_info "🎯 All parallel wget + trimming jobs complete."
}

hisat2_de_novo_index_align_sort_stringtie_pipeline() {
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
			run_with_time_to_log \
				hisat2 -p "${THREADS}" --dta \
					-x "$index_prefix" \
					-1 "$trimmed1" \
					-2 "$trimmed2" \
					-S "$sam"
			
			log_info "Converting SAM to sorted BAM for $fasta_tag with $SRR..."
			run_with_time_to_log samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_time_to_log samtools index -@ "${THREADS}" "$bam"
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
		
		log_info "Deleting the BAM file."
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
		local left_reads=() right_reads=()
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
			fi
			
			if [[ -n "$trimmed1" && -n "$trimmed2" ]]; then
				left_reads+=("$trimmed1")
				right_reads+=("$trimmed2")
			else
				log_warn "Trimmed FASTQ for $SRR not found; skipping from Trinity input."
			fi
		done
		
		if [[ ${#left_reads[@]} -eq 0 ]]; then
			log_error "No trimmed FASTQ files found for Trinity assembly."
			return 1
		fi
		
		# Run Trinity
		local left_files=$(IFS=,; echo "${left_reads[*]}")
		local right_files=$(IFS=,; echo "${right_reads[*]}")
		
		log_info "Running Trinity with ${#left_reads[@]} sample pairs..."
		run_with_time_to_log Trinity \
			--seqType fq \
			--left "$left_files" \
			--right "$right_files" \
			--CPU "$THREADS" \
			--max_memory 20G \
			--output "$trinity_out_dir" \
			--normalize_reads \
			--full_cleanup
	fi

	# STEP 2: ALIGN READS TO TRINITY ASSEMBLY AND RUN STRINGTIE FOR QUANTIFICATION
	for SRR in "${rnaseq_list[@]}"; do
		local HISAT2_DIR="$HISAT2_DE_NOVO_ROOT/${fasta_tag}_trinity/$SRR"
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

		local bam="$HISAT2_DIR/${SRR}_${fasta_tag}_trinity_mapped_sorted.bam"
		
		# Align reads directly to Trinity assembly using bowtie2 (Trinity's preferred aligner)
		if [[ -f "$bam" && -f "${bam}.bai" ]]; then
			log_info "BAM for $SRR and Trinity assembly of $fasta_tag already exists. Skipping alignment."
		else
			log_info "Building bowtie2 index and aligning $SRR to Trinity assembly of $fasta_tag..."
			
			# Build bowtie2 index for Trinity assembly
			local trinity_index="$HISAT2_DIR/${fasta_tag}_trinity_index"
			if [[ ! -f "${trinity_index}.1.bt2" ]]; then
				run_with_time_to_log bowtie2-build "$trinity_fasta" "$trinity_index"
			fi
			
			# Align with bowtie2
			local sam="$HISAT2_DIR/${SRR}_${fasta_tag}_trinity_mapped.sam"
			run_with_time_to_log \
				bowtie2 -p "${THREADS}" --local --no-unal \
					-x "$trinity_index" \
					-1 "$trimmed1" \
					-2 "$trimmed2" \
					-S "$sam"
			
			log_info "Converting SAM to sorted BAM for Trinity assembly of $fasta_tag with $SRR..."
			run_with_time_to_log samtools sort -@ "${THREADS}" -o "$bam" "$sam"
			run_with_time_to_log samtools index -@ "${THREADS}" "$bam"
			rm -f "$sam"
		fi

		# StringTie quantification for DeSeq2 compatibility
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

# ==============================================================================
# RUN / ENTRYPOINT
# ==============================================================================
run_all() {
	# Main pipeline entrypoint: runs all steps for each FASTA and RNA-seq list
	# Steps: Logging, Download and  Trim, HISAT2 alignment to Stringtie, and Cleanup. 
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

	log_info "SRR samples to process:"
	for SRR in "${rnaseq_list[@]}"; do
		log_info "$SRR"
	done

	if [[ $RUN_DOWNLOAD_and_TRIM_SRR == "TRUE" ]]; then
		log_step "STEP 01: Download and trim RNA-seq data"
		#download_and_trim_srrs "${rnaseq_list[@]}"
		download_and_trim_srrs_parallel "${rnaseq_list[@]}"
		#download_and_trim_srrs_parallel_fastqdump "${rnaseq_list[@]}"
		#download_and_trim_srrs_wget_parallel "${rnaseq_list[@]}"
	fi

	if [[ $RUN_HISAT2_INDEX_ALIGN_SORT_STRINGTIE == "TRUE" ]]; then
		log_step "STEP 02: Build HISAT2 index and align"
		hisat2_de_novo_index_align_sort_stringtie_pipeline \
			--FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
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
	run_all \
		--FASTA "$fasta_input" \
		--RNASEQ_LIST "${SRR_COMBINED_LIST[@]}"
done

# Zip all the content of this folder: STRINGTIE_HISAT2_DE_NOVO_ROOT="03_stringtie/"

#tar -czvf 03_stringtie__$(date +%Y%m%d_%H%M%S).tar.gz 03_stringtie/  # Archive all StringTie results for sharing or backup

#OFF
