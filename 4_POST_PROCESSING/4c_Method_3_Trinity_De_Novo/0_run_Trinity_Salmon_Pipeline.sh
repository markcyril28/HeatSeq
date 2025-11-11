#!/bin/bash
# ==============================================================================
# TRINITY + SALMON PIPELINE WRAPPER FOR METHOD 3
# ==============================================================================
# Description: Wrapper script to run Trinity de novo assembly with Salmon
#              quantification and tximport for DESeq2 analysis
# Author: Mark Cyril R. Mercado
# Version: v2
# Date: November 2025
#
# This script uses the updated Trinity pipeline with:
# - Salmon quantification (instead of RSEM)
# - Tximport for proper statistical handling
# - Direct DESeq2 integration
# ==============================================================================

set -euo pipefail

# Navigate to project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Source required modules
source "modules/logging_utils.sh"
source "modules/pipeline_utils.sh"
source "modules/methods.sh"

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# FASTA files to analyze
FASTA_FILES=(
	"0_INPUT_FASTAs/All_Smel_Genes.fasta"
	# Add more FASTA files here
)

# RNA-seq samples to process
SRR_LIST=(
	SRR3884686  # Buds_0.7cm
	SRR3884687  # Opened_Buds
	SRR3884597  # Flowers
	# Add more SRR IDs here
)

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

setup_logging

log_step "Starting Trinity + Salmon Pipeline (Method 3)"
log_info "Workspace: $PROJECT_ROOT"
log_info "FASTA files: ${#FASTA_FILES[@]}"
log_info "RNA-seq samples: ${#SRR_LIST[@]}"

for fasta in "${FASTA_FILES[@]}"; do
	if [[ ! -f "$fasta" ]]; then
		log_error "FASTA file not found: $fasta"
		continue
	fi
	
	log_step "Processing FASTA: $(basename "$fasta")"
	
	# Run Trinity + Salmon pipeline
	trinity_de_novo_alignment_pipeline \
		--FASTA "$fasta" \
		--RNASEQ_LIST "${SRR_LIST[@]}"
	
	log_info "✓ Completed pipeline for $(basename "$fasta")"
done

log_step "Trinity + Salmon Pipeline completed"
log_info "Next steps:"
log_info "1. Edit sample metadata: 4c_Method_3_Trinity_De_Novo/6_matrices_from_stringtie/sample_info.tsv"
log_info "2. Run tximport if not done: Rscript run_tximport_trinity_salmon.R"
log_info "3. Perform DESeq2 analysis in R"
log_info "4. Generate custom visualizations"
