#!/bin/bash
# ==============================================================================
# CENTRALIZED CONFIGURATION
# ==============================================================================
# Single source of truth for all pipeline paths and settings
# Sourced by: pipeline_utils.sh, methods/methods.sh
# ==============================================================================

set -euo pipefail

# Guard against double-sourcing
[[ "${CONFIG_SOURCED:-}" == "true" ]] && return 0
export CONFIG_SOURCED="true"

# ==============================================================================
# GLOBAL RUNTIME SETTINGS
# ==============================================================================
THREADS="${THREADS:-4}"
JOBS="${JOBS:-2}"
keep_bam_global="${keep_bam_global:-n}"

# Bowtie2 alignment mode
BOWTIE2_MODE="${BOWTIE2_MODE:-very-sensitive-local}"

# ==============================================================================
# DIRECTORY STRUCTURE
# ==============================================================================
RAW_DIR_ROOT="1_RAW_SRR"
TRIM_DIR_ROOT="2_TRIMMED_SRR"
FASTQC_ROOT="3_FastQC"
POST_PROCESSING_ROOT="4_POST_PROCESSING"

# Method 1: HISAT2 Reference Guided
HISAT2_REF_GUIDED_ROOT="$POST_PROCESSING_ROOT/4a_M1_HISAT2_Ref_Guided/4_HISAT2_WD"
HISAT2_REF_GUIDED_INDEX_DIR="$HISAT2_REF_GUIDED_ROOT/index"
STRINGTIE_HISAT2_REF_GUIDED_ROOT="$POST_PROCESSING_ROOT/4a_M1_HISAT2_Ref_Guided/5_stringtie_WD"
HISAT2_REF_GUIDED_MATRIX_ROOT="$POST_PROCESSING_ROOT/4a_M1_HISAT2_Ref_Guided/6_matrices"

# Method 2: HISAT2 De Novo
HISAT2_DE_NOVO_ROOT="$POST_PROCESSING_ROOT/4b_M2_HISAT2_De_Novo/4_HISAT2_WD"
HISAT2_DE_NOVO_INDEX_DIR="$HISAT2_DE_NOVO_ROOT/index"
STRINGTIE_HISAT2_DE_NOVO_ROOT="$POST_PROCESSING_ROOT/4b_M2_HISAT2_De_Novo/5_stringtie_WD"
HISAT2_DE_NOVO_MATRIX_ROOT="$POST_PROCESSING_ROOT/4b_M2_HISAT2_De_Novo/6_matrices"

# Method 3: STAR Splice-Aware Alignment
STAR_ALIGN_ROOT="$POST_PROCESSING_ROOT/4c_M3_STAR_Alignment"
STAR_INDEX_ROOT="$STAR_ALIGN_ROOT/star_index"
STAR_GENOME_DIR="$STAR_ALIGN_ROOT/alignments"

# Method 4: Salmon SAF
SALMON_SAF_ROOT="$POST_PROCESSING_ROOT/4d_M4_Salmon_Saf_Quant"
SALMON_INDEX_ROOT="$SALMON_SAF_ROOT/4_Salmon_WD/index"
SALMON_QUANT_ROOT="$SALMON_SAF_ROOT/5_Salmon_Quant_WD"
SALMON_SAF_MATRIX_ROOT="$SALMON_SAF_ROOT/6_matrices_from_Salmon"
SALMON_MATRIX_ROOT="${SALMON_MATRIX_ROOT:-$SALMON_SAF_MATRIX_ROOT}"

# Method 5: Bowtie2 + RSEM
BOWTIE2_RSEM_ROOT="$POST_PROCESSING_ROOT/4e_M5_RSEM_Bowtie2_Quant"
RSEM_INDEX_ROOT="$BOWTIE2_RSEM_ROOT/4_Bowtie2_WD/index"
RSEM_QUANT_ROOT="$BOWTIE2_RSEM_ROOT/5_RSEM_Quant_WD"
RSEM_MATRIX_ROOT="$BOWTIE2_RSEM_ROOT/6_matrices_from_RSEM"

# ==============================================================================
# STAR CONFIGURATION
# ==============================================================================
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-NoSharedMemory}"
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"
STAR_STRAND_SPECIFIC="${STAR_STRAND_SPECIFIC:-intronMotif}"

# ==============================================================================
# SRR SAMPLE ARRAYS (initialize if not set)
# ==============================================================================
if ! declare -p SRR_COMBINED_LIST &>/dev/null; then
	declare -a SRR_COMBINED_LIST=()
fi
if ! declare -p SRR_LIST_PRJNA328564 &>/dev/null; then
	declare -a SRR_LIST_PRJNA328564=()
fi

# ==============================================================================
# INITIALIZE DIRECTORIES
# ==============================================================================
init_directories() {
	mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
		"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
		"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
		"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" "$STAR_GENOME_DIR" \
		"$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
		"$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT"
}

# ==============================================================================
# DISPLAY CONFIGURATION
# ==============================================================================
show_pipeline_configuration() {
	log_info "=== PIPELINE CONFIGURATION ==="
	log_info "Data Processing:"
	log_info "  Download SRR: $([ "${RUN_DOWNLOAD_SRR:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Trim SRR: $([ "${RUN_TRIM_SRR:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Quality Control: $([ "${RUN_QUALITY_CONTROL:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "Analysis Methods:"
	log_info "  Method 1 - HISAT2 Ref-Guided: $([ "${RUN_METHOD_1_HISAT2_REF_GUIDED:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 2 - HISAT2 De Novo: $([ "${RUN_METHOD_2_HISAT2_DE_NOVO:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 3 - STAR Alignment: $([ "${RUN_METHOD_3_STAR_ALIGNMENT:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 4 - Salmon SAF: $([ "${RUN_METHOD_4_SALMON_SAF:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "  Method 5 - Bowtie2 + RSEM: $([ "${RUN_METHOD_5_BOWTIE2_RSEM:-FALSE}" = "TRUE" ] && echo "[ENABLED]" || echo "[DISABLED]")"
	log_info "============================="
}
