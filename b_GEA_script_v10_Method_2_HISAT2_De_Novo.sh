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
# CONFIGURATION AND RUNTIME SWITCHES
# ============================================================================== 

# Runtime Configuration
THREADS=32                               # Number of threads to use for parallel operations
JOBS=4                                  # Number of parallel jobs for GNU Parallel 

# RNA-seq Library Configuration
RNA_STRAND_PROTOCOL="RF"                # RNA-seq strand protocol: "RF" (dUTP), "FR" (ligation), or "unstranded"
                                       # RF = first read reverse complement, second read forward (most common)
                                       # FR = first read forward, second read reverse complement  
                                       # unstranded = no strand specificity

# Pipeline Control Switches
RUN_MAMBA_INSTALLATION=FALSE
RUN_DOWNLOAD_and_TRIM_SRR=TRUE
RUN_GZIP_TRIMMED_FILES=TRUE
RUN_HEATMAP_WRAPPER_for_HISAT2_DE_NOVO=FALSE
RUN_ZIP_RESULTS=FALSE

# GEA Methods 
RUN_METHOD_2_HISAT2_DE_NOVO=FALSE	
RUN_METHOD_1_HISAT2_REF_GUIDED=FALSE
RUN_METHOD_3_TRINITY_DE_NOVO=FALSE
RUN_METHOD_4_SALMON_SAF=FALSE
RUN_METHOD_5_BOWTIE2_RSEM=FALSE

# Quality Control and Analysis Options
RUN_QUALITY_CONTROL=FALSE
RUN_METHOD_COMPARISON=FALSE

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

# Legacy GTF and FASTA references (commented out)
#All_SmelGIF_GTF_FILE="0_INPUT_FASTAs/All_SmelDMP_Head_Gene_Name_v4.gtf"
#Eggplant_V4_1_transcripts_function_FASTA_FILE="0_INPUT_FASTAs/Eggplant_V4.1_transcripts_function.fa"

# FASTA Files for Analysis
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

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS)
# ==============================================================================

# SRA Run Selector Tool link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa

SRR_LIST_PRJNA328564=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa
	# Developmental stages arranged from early to late
	SRR3884685	# Radicles (earliest - germination)
	SRR3884677	# Cotyledons (seed leaves)
	SRR3884675	# Roots (root development)
	SRR3884690	# Stems (vegetative growth)
	SRR3884689	# Leaves (vegetative growth)
	SRR3884684	# Senescent_leaves (leaf aging)
	SRR3884686	# Buds_0.7cm (flower bud initiation)
	SRR3884687	# Opened_Buds (flower development)
	SRR3884597	# Flowers (anthesis)
	SRR3884679	# Pistils (female reproductive parts)
	SRR3884608	# Fruits_1cm (early fruit development)
	SRR3884620	# Fruits_Stage_1 (early fruit stage)
	SRR3884631	# Fruits_6cm (fruit enlargement)
	SRR3884642	# Fruits_Skin_Stage_2 (mid fruit development)
	SRR3884653	# Fruits_Flesh_Stage_2 (mid fruit development)
	SRR3884664	# Fruits_Calyx_Stage_2 (mid fruit development)
	SRR3884680	# Fruits_Skin_Stage_3 (late fruit development)
	SRR3884681	# Fruits_Flesh_Stage_3 (late fruit development)
	SRR3884678	# Fruits_peduncle (fruit attachment)
)

SRR_LIST_SAMN28540077=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540077&o=acc_s%3Aa&s=SRR20722234,SRR20722233,SRR20722232,SRR20722230,SRR20722225,SRR20722226,SRR20722227,SRR20722228,SRR20722229
	SRR2072232	# mature_fruits #SRR20722226	# young_fruits
	SRR20722234	# flowers #SRR20722228	# sepals
	SRR21010466 # Buds, Nonparthenocarpy ID: PRJNA865018 
	SRR20722233	# leaf_buds #SRR20722230	# mature_leaves
	SRR20722227	# stems
	SRR20722229	# roots
)

SRR_LIST_SAMN28540068=(
	#Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540068&o=acc_s%3Aa
	SRR20722387 # mature_fruits
	SRR23909863 # Fully Develop (FD) Flower ID: PRJNA941250
	SRR20722297 # flower_buds #SRR20722385 # sepals
	SRR20722296 # leaf_buds #SRR20722386 # mature_leaves #SRR20722383 # young_leaves
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
)

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# Raw and Processed Data Directories
RAW_DIR_ROOT="1_RAW_SRR"                # Raw SRR download directory
TRIM_DIR_ROOT="2_TRIMMED_SRR"           # Trimmed reads directory
FASTQC_ROOT="3_FastQC"                  # FastQC quality control reports

# Method 1: HISAT2 Reference Guided 
HISAT2_REF_GUIDED_ROOT="4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="4a_Method_1_HISAT2_Ref_Guided/4_HISAT2_WD/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="4a_Method_1_HISAT2_Ref_Guided/5_stringtie_WD/a_Method_1_RAW_RESULTs"

# Method 2: HISAT2 De Novo Assembly (main method in this script)
HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="4b_Method_2_HISAT2_De_Novo/4_HISAT2_WD/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="4b_Method_2_HISAT2_De_Novo/5_stringtie_WD/a_Method_2_RAW_RESULTs"

# Method 3: Trinity De Novo Assembly
TRINITY_DE_NOVO_ROOT="4c_Method_3_Trinity_De_Novo/4_Trinity_WD"
STRINGTIE_TRINITY_DE_NOVO_ROOT="4c_Method_3_Trinity_De_Novo/5_stringtie_WD/a_Method_3_RAW_RESULTs"

# Method 4: Salmon SAF Quantification
SALMON_SAF_ROOT="4d_Method_4_Salmon_Saf_Quantification"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/quant"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/matrices"

# Method 5: Bowtie2 + RSEM Quantification
BOWTIE2_RSEM_ROOT="4e_Method_5_Bowtie2_Quantification"
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
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================

rm -rf "$RAW_DIR_ROOT"                   # Remove previous raw SRR files
#rm -rf "$FASTQC_ROOT"                   # Remove previous FastQC results
#rm -rf "$HISAT2_DE_NOVO_ROOT"           # Remove previous HISAT2 results
#rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"      # Remove previous HISAT2 index
#rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT" # Remove previous StringTie results

# ==============================================================================
# LOGGING SYSTEM, ENVIRONMENT, and PROGRAM
# ==============================================================================

# Logging Configuration
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_full_log.log}"

# Logging Functions
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
	log_choice="${log_choice:-1}"
	if [[ "$log_choice" == "2" ]]; then
		exec >"$LOG_FILE" 2>&1
	else
		exec > >(tee -a "$LOG_FILE") 2>&1
	fi
	log_info "Logging to: $LOG_FILE"
}

# Error handling and cleanup traps
trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
trap 'log_info "Script finished. See log: $LOG_FILE"' EXIT

# Utility function for command timing and logging
run_with_time_to_log() {
	# Run a command and log resource usage (tracks time and memory)
	/usr/bin/time -v "$@" >> "$LOG_FILE" 2>&1
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
	
	# FastQC on trimmed reads
	if command -v fastqc >/dev/null 2>&1; then
		fastqc -t "$THREADS" -o "$FASTQC_ROOT/$SRR" \
			"$TrimGalore_DIR"/${SRR}*val*.fq* 2>/dev/null || \
			log_warn "FastQC failed for $SRR"
	else
		log_warn "FastQC not available. Skipping read quality assessment."
	fi
	
	# MultiQC summary (if available)
	if command -v multiqc >/dev/null 2>&1; then
		multiqc "$FASTQC_ROOT" -o "$FASTQC_ROOT/summary" --force 2>/dev/null || \
			log_warn "MultiQC failed. Individual FastQC reports still available."
	fi
}


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

# Prerequisites (install if needed):
#conda install -c conda-forge aria2 -y
#conda install -c bioconda -c conda-forge parallel-fastq-dump -y

# Prerequisites installation - install required bioinformatics tools
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
# PIPELINE FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# DATA DOWNLOAD AND PREPROCESSING FUNCTIONS
# ------------------------------------------------------------------------------
gzip_trimmed_fastq_files() {
	# Compress all trimmed FASTQ files in the trimming directory using GNU parallel
	log_info "Compressing all trimmed FASTQ files in $TRIM_DIR_ROOT using GNU parallel..."
	find "$TRIM_DIR_ROOT" -type f -name "*.fq" -print0 | \
		parallel -0 -j "$JOBS" gzip {}
	log_info "Parallel compression of trimmed FASTQ files completed."
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
			# Compress FASTQ files to save space
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq" && -f "$raw_files_DIR/${SRR}_2.fastq" ]]; then
				log_info "Compressing FASTQ files for $SRR..."
				gzip "$raw_files_DIR/${SRR}_1.fastq" "$raw_files_DIR/${SRR}_2.fastq"
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

        # --------------------------------------------------
        # Trim reads using Trim Galore
        # --------------------------------------------------
        log_info "Trimming $SRR..."
        run_with_time_to_log trim_galore --cores "${THREADS}" \
            --paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

        log_info "Done working on $SRR."
        log_info "--------------------------------------------------"

        # --------------------------------------------------
        # Verify trimming success and delete raw files if successful
        # --------------------------------------------------
        # Check trimmed output files and their sizes
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
	# Download and trim RNA-seq data using kingfisher
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

		log_info "Working on $SRR with kingfisher..."
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
		# Check if trimmed files already exist - skip if so
		if [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; then
			log_info "Trimmed files for $SRR already exist. Skipping download and trimming."
			log_info "--------------------------------------------------"
			continue
		fi
		
		# Check for existing raw FASTQ files
		if [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
			raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
			raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
			log_info "Raw FASTQ files for $SRR already exist. Skipping download."
		else
			# --------------------------------------------------
			# Download with kingfisher (tries multiple sources automatically)
			# --------------------------------------------------
			log_info "Downloading $SRR with kingfisher..."
			run_with_time_to_log kingfisher get \
				--run-identifiers "$SRR" \
				--output-directory "$raw_files_DIR" \
				--download-threads "$THREADS" \
				--extraction-threads "$THREADS" \
				--gzip \
				--check-md5sums
			
			# Kingfisher outputs files as SRR_1.fastq.gz and SRR_2.fastq.gz
			if [[ -f "$raw_files_DIR/${SRR}_1.fastq.gz" && -f "$raw_files_DIR/${SRR}_2.fastq.gz" ]]; then
				raw1="$raw_files_DIR/${SRR}_1.fastq.gz"
				raw2="$raw_files_DIR/${SRR}_2.fastq.gz"
				log_info "Kingfisher download successful for $SRR"
			else
				log_warn "Raw FASTQ for $SRR not found after kingfisher download; skipping."
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
		# Verify trimming success and delete raw files if successful
		# --------------------------------------------------
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
	log_info "All SRR samples processed with kingfisher."
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
		run_with_time_to_log trim_galore --cores "${THREADS}" \
			--paired "$raw1" "$raw2" --output_dir "$TrimGalore_DIR"

		# --------------------------------------------------
		# Verify trimming success and cleanup
		# --------------------------------------------------
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

	# ============================================
	# Run jobs in parallel
	# ============================================
	log_info "Starting parallel download and trimming..."
	printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr {}
	log_info "All parallel SRR jobs completed."

	gzip_trimmed_fastq_files
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
        # Verify trimming success and cleanup raw files
        # ---------------------------------------------------------
        local trimming_success=false
        local trimmed_files_exist=false
        local trimmed_files_not_empty=false
        local trimmed_files_complete=false

        # Check for existence of trimmed files
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            trimmed_files_exist=true
            
            # Verify file sizes are non-zero
            if { [[ -s "$trimmed1" && -s "$trimmed2" ]] || [[ -s "${trimmed1}.gz" && -s "${trimmed2}.gz" ]]; }; then
                trimmed_files_not_empty=true
            fi

            # Additional check for complete files (no .part files left)
            if ! compgen -G "${TrimGalore_DIR}/${SRR}*.part" > /dev/null; then
                trimmed_files_complete=true
            fi
        fi

        # Final verification and cleanup
        if [[ "$trimmed_files_exist" == true && "$trimmed_files_not_empty" == true && "$trimmed_files_complete" == true ]]; then
            trimming_success=true
            log_info "Trimming successfully verified for $SRR"
            log_info "Cleaning up raw FASTQs..."
            rm -f "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz"
            log_info "Raw files deleted"
        else
            log_warn "Trimming verification failed for $SRR - keeping raw files"
            [[ "$trimmed_files_exist" == false ]] && log_warn "Trimmed files missing"
            [[ "$trimmed_files_not_empty" == false ]] && log_warn "Trimmed files are empty"
            [[ "$trimmed_files_complete" == false ]] && log_warn "Trimming process may be incomplete (.part files found)"
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

        log_info "PROCESSING: Working on $SRR..."
        mkdir -p "$raw_files_DIR" "$TrimGalore_DIR"

        # Skip if trimming already completed
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "SKIPPED: Trimmed files exist for $SRR - skipping."
            return 0
        fi

        # ----------------------------------------
        # Fetch ENA FASTQ URLs dynamically
        # ----------------------------------------
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

        # ----------------------------------------
        # Trim Galore (low RAM, single-thread)
        # ----------------------------------------
        log_info "TRIMMING: Adapters for $SRR..."
        run_with_time_to_log trim_galore --cores 1 \
            --paired "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz" \
            --output_dir "$TrimGalore_DIR"

        # ----------------------------------------
        # Cleanup if trimming successful
        # ----------------------------------------
        if { [[ -f "$trimmed1" && -f "$trimmed2" ]] || [[ -f "${trimmed1}.gz" && -f "${trimmed2}.gz" ]]; }; then
            log_info "CLEANUP: Removing raw FASTQs for $SRR..."
            rm -f "$raw_files_DIR/${SRR}_1.fastq.gz" "$raw_files_DIR/${SRR}_2.fastq.gz"
        fi

        log_info "COMPLETED: Finished $SRR."
        log_info "--------------------------------------------------"
    }

    export -f _process_single_srr_wget

    # --------------------------------------------
    # Run jobs in parallel
    # --------------------------------------------
    log_info "STARTING: Parallel ENA downloads and trimming..."
    printf "%s\n" "${SRR_LIST[@]}" | parallel -j "${PARALLEL_JOBS:-$JOBS}" _process_single_srr_wget {}
    log_info "COMPLETED: All parallel wget + trimming jobs complete."
}

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

# ==============================================================================
# MAIN EXECUTION FUNCTIONS
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
	
	# Show pipeline configuration
	#show_pipeline_configuration

	log_info "SRR samples to process:"
	for SRR in "${rnaseq_list[@]}"; do
		log_info "$SRR"
	done

	if [[ $RUN_DOWNLOAD_and_TRIM_SRR == "TRUE" ]]; then
		log_step "STEP 01: Download and trim RNA-seq data"
		#download_and_trim_srrs "${rnaseq_list[@]}"
		#download_and_trim_srrs_parallel "${rnaseq_list[@]}"
		#download_and_trim_srrs_parallel_fastqdump "${rnaseq_list[@]}"
		#download_and_trim_srrs_wget_parallel "${rnaseq_list[@]}"
		download_kingfisher_and_trim_srrs "${rnaseq_list[@]}"
		
		if [[ $RUN_QUALITY_CONTROL == "TRUE" ]]; then
			log_step "STEP 01b: Quality control analysis"
			for SRR in "${rnaseq_list[@]}"; do
				run_quality_control "$SRR"
			done
		fi
	fi


	# Method 1: HISAT2 Reference-Guided Pipeline
	if [[ $RUN_METHOD_1_HISAT2_REF_GUIDED == "TRUE" ]]; then
		log_step "STEP 02a: HISAT2 Reference-Guided Pipeline"
		# Note: Requires GTF file - add --GTF parameter when calling
		log_warn "HISAT2 Reference-Guided pipeline requires GTF file parameter"
		# hisat2_ref_guided_pipeline --FASTA "$fasta" --GTF "$gtf_file" --RNASEQ_LIST "${rnaseq_list[@]}"
	fi

	# Method 2: HISAT2 De Novo Pipeline (Main method)
	if [[ $RUN_METHOD_2_HISAT2_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 02b: HISAT2 De Novo Pipeline"
		hisat2_de_novo_pipeline \
			--FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
	fi

	# Method 3: Trinity De Novo Pipeline
	if [[ $RUN_METHOD_3_TRINITY_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 03: Trinity De Novo Assembly and Quantification"
		trinity_de_novo_alignment_pipeline \
			--FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
	fi

	# Method 4: Salmon SAF Quantification
	if [[ $RUN_METHOD_4_SALMON_SAF == "TRUE" ]]; then
		log_step "STEP 04: Salmon SAF Quantification"
		# Note: Requires genome file for decoy-aware indexing
		local genome_file="${fasta%.*}_genome.fa"  # Assumes genome file exists
		if [[ -f "$genome_file" ]]; then
			salmon_saf_pipeline \
				--FASTA "$fasta" --GENOME "$genome_file" --RNASEQ_LIST "${rnaseq_list[@]}"
		else
			log_warn "Genome file '$genome_file' not found. Skipping Salmon SAF pipeline."
			log_warn "Please provide genome file for decoy-aware Salmon quantification."
		fi
	fi

	# Method 5: Bowtie2 + RSEM Quantification
	if [[ $RUN_METHOD_5_BOWTIE2_RSEM == "TRUE" ]]; then
		log_step "STEP 05: Bowtie2 + RSEM Quantification"
		bowtie2_rsem_pipeline \
			--FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
	fi

	end_time=$(date +%s)
	log_step "Final timing"
	log_info "Script ended at: $(date -d @$end_time)"
	elapsed=$((end_time - start_time))
	formatted_elapsed=$(date -u -d @${elapsed} +%H:%M:%S)
	log_info "Elapsed time: $formatted_elapsed"
}

# ==============================================================================
# SCRIPT EXECUTION
# ==============================================================================

# Call the function if installation is enabled
if [[ $RUN_MAMBA_INSTALLATION == "TRUE" ]]; then
	mamba_install
fi

if [[ $RUN_GZIP_TRIMMED_FILES == "TRUE" ]]; then
	log_step "Gzipping trimmed FASTQ files to save space"
	gzip_trimmed_fastq_files
fi

# Execute the pipeline for each FASTA input file
for fasta_input in "${ALL_FASTA_FILES[@]}"; do
	# Run the complete pipeline for each FASTA file with all SRR samples
	run_all \
		--FASTA "$fasta_input" \
		--RNASEQ_LIST "${SRR_COMBINED_LIST[@]}"
done

# ==============================================================================
# POST-PROCESSING: HEATMAP WRAPPER EXECUTION for HISAT2 DE NOVO
# ==============================================================================

if [[ $RUN_HEATMAP_WRAPPER_for_HISAT2_DE_NOVO == "TRUE" ]]; then
	log_step "Heatmap Wrapper post-processing enabled"
	# Execute the Heatmap Wrapper script for post-processing
	if [[ -f "4b_Method_2_HISAT2_De_Novo/0_Heatmap_Wrapper.sh" ]]; then
		log_step "Executing Heatmap Wrapper post-processing script"
		chmod +x 4b_Method_2_HISAT2_De_Novo/*.sh
		chmod +x 4b_Method_2_HISAT2_De_Novo/0_Heatmap_Wrapper.sh

		if bash 4b_Method_2_HISAT2_De_Novo/0_Heatmap_Wrapper.sh 2>&1; then
			log_info "Heatmap Wrapper completed successfully"
		else
			log_error "Heatmap Wrapper failed with exit code $?"
		fi
	else
		log_warn "Heatmap Wrapper script not found - skipping"
	fi
fi


# ==============================================================================
# POST-PROCESSING OPTIONS (COMMENTED OUT)
# ==============================================================================

if [[ $RUN_ZIP_RESULTS == "TRUE" ]]; then
	# Optional: Archive StringTie results for sharing or backup
	#tar -czvf "stringtie_results_$(date +%Y%m%d_%H%M%S).tar.gz" "$STRINGTIE_HISAT2_DE_NOVO_ROOT"
	#tar -czvf HISAT2_DE_NOVO_ROOT_HPC_$(date +%Y%m%d_%H%M%S).tar.gz $HISAT2_DE_NOVO_ROOT
	tar -czvf 4b_Method_2_HISAT2_De_Novo_$(date +%Y%m%d_%H%M%S).tar.gz 4b_Method_2_HISAT2_De_Novo/
fi

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
echo "END OF SCRIPT"