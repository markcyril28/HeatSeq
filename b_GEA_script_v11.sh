#!/bin/bash

# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE
# ==============================================================================
# Description: RNA-seq analysis (Tissue/Organ-specific Expression Profiling) pipeline using 
# using various methods. 
# Author: Mark Cyril R. Mercado
# Version: v11
# Date: October 2025

# ==============================================================================

set -euo pipefail

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
#conda activate GEA_ENV
conda activate geaheat

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#source "modules/all_functions.sh"
source "modules/logging_utils.sh"
source "modules/pipeline_utils.sh"
source "modules/methods.sh"
bash init_setup.sh

# ============================================================================== 
# CONFIGURATION AND RUNTIME SWITCHES
# ============================================================================== 

# Runtime Configuration
THREADS=96                               # Number of threads to use for parallel operations
JOBS=4                                  # Number of parallel jobs for GNU Parallel 

# Export variables for function access
export THREADS JOBS

# Pipeline Control Switches, Quality Control and Analysis Options
# Pipeline Control Switches - Enable/disable pipeline stages
PIPELINE_STAGES=(
	#"MAMBA_INSTALLATION"
	#"DOWNLOAD_and_TRIM_SRR"
	#"GZIP_TRIMMED_FILES"
	#"QUALITY_CONTROL"
	
	## GEA Methods
	"METHOD_1_HISAT2_REF_GUIDED"
	"METHOD_2_HISAT2_DE_NOVO"
	"METHOD_4_SALMON_SAF"
	"METHOD_5_BOWTIE2_RSEM"
	#"METHOD_3_TRINITY_DE_NOVO"

	#"HEATMAP_WRAPPER"
	#"ZIP_RESULTS"
)

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

# Legacy GTF and FASTA references (commented out)
#All_SmelGIF_GTF_FILE="0_INPUT_FASTAs/All_SmelDMP_Head_Gene_Name_v4.gtf"
Eggplant_V4_1_transcripts_function_FASTA_FILE="0_INPUT_FASTAs/Eggplant_V4.1_transcripts.function.fa"
ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE="0_INPUT_FASTAs/All_Smel_Genes_Full_Name_reformatted.gtf"
decoy="0_INPUT_FASTAs/TEST.fasta"
gtf_file="${ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE}"

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
	SRR3884686	# Buds_0.7cm (flower bud initiation) [MAIN INTEREST]
	SRR3884687	# Opened_Buds (flower development) 	 [MAIN INTEREST]
	SRR3884597	# Flowers (anthesis)/				 [MAIN INTEREST]	
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

SRR_LIST_PRJNA865018=(
# A Good Dataset for SmelDMP GEA: 
# 	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=%20%20PRJNA865018&o=acc_s%3Aa PRJNA865018
#	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA941250&o=acc_s%3Aa PRJNA941250 # Buds, Opened Buds
	SRR21010454	# buds_3
	SRR21010456	# buds_2
	SRR21010466	# buds_1
	SRR21010458	# flowers_3
	SRR21010460	# flowers_2
	SRR21010462	# flowers_1
	SRR21010450	# fruits_2
	SRR21010452	# fruits_1
	SRR21010464	# fruits_3
)

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
	#"${SRR_LIST_PRJNA328564[@]}"	# Main Dataset for GEA. 
	"${SRR_LIST_SAMN28540077[@]}"	# Chinese Dataset for replicability. 
	"${SRR_LIST_SAMN28540068[@]}"	# Chinese Dataset for replicability. 
	#"${SRR_LIST_PRJNA865018[@]}"	# Good Dataset for SmelDMP GEA.
)


# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# Create required directories (including log directories)
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$TRINITY_DE_NOVO_ROOT" "$STRINGTIE_TRINITY_DE_NOVO_ROOT" \
	"$SALMON_SAF_ROOT" "$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
	"$BOWTIE2_RSEM_ROOT" "$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT" \
	"logs/log_files" "logs/time_files"

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================

#rm -rf "$RAW_DIR_ROOT"                   # Remove previous raw SRR files
#rm -rf "$FASTQC_ROOT"                   # Remove previous FastQC results
#rm -rf "$HISAT2_DE_NOVO_ROOT"           # Remove previous HISAT2 results
#rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"      # Remove previous HISAT2 index
#rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT" # Remove previous StringTie results

# ==============================================================================
# MAIN EXECUTION FUNCTIONS
# ==============================================================================

# Convert array to boolean flags for backward compatibility
RUN_MAMBA_INSTALLATION=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^MAMBA_INSTALLATION$" && echo "TRUE" || echo "FALSE")
RUN_DOWNLOAD_and_TRIM_SRR=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^DOWNLOAD_and_TRIM_SRR$" && echo "TRUE" || echo "FALSE")
RUN_GZIP_TRIMMED_FILES=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^GZIP_TRIMMED_FILES$" && echo "TRUE" || echo "FALSE")
RUN_QUALITY_CONTROL=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^QUALITY_CONTROL$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_1_HISAT2_REF_GUIDED=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_1_HISAT2_REF_GUIDED$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_2_HISAT2_DE_NOVO=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_2_HISAT2_DE_NOVO$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_3_TRINITY_DE_NOVO=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_3_TRINITY_DE_NOVO$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_4_SALMON_SAF=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_4_SALMON_SAF$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_5_BOWTIE2_RSEM=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_5_BOWTIE2_RSEM$" && echo "TRUE" || echo "FALSE")
RUN_HEATMAP_WRAPPER=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^HEATMAP_WRAPPER$" && echo "TRUE" || echo "FALSE")
RUN_ZIP_RESULTS=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^ZIP_RESULTS$" && echo "TRUE" || echo "FALSE")



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
		download_and_trim_srrs_parallel "${rnaseq_list[@]}"
		#download_and_trim_srrs_parallel_fastqdump "${rnaseq_list[@]}"
		#download_and_trim_srrs_wget_parallel "${rnaseq_list[@]}"
		#download_kingfisher_and_trim_srrs "${rnaseq_list[@]}"
		
		if [[ $RUN_QUALITY_CONTROL == "TRUE" ]]; then
			log_step "STEP 01b: Quality Control analysis"
			for SRR in "${rnaseq_list[@]}"; do
				run_quality_control "$SRR"
			done
		fi
	fi

	if [[ $RUN_QUALITY_CONTROL == "TRUE" ]]; then
			log_step "STEP 01b: Quality Control analysis"
			for SRR in "${rnaseq_list[@]}"; do
				run_quality_control "$SRR"
			done
	fi


	# Method 1: HISAT2 Reference-Guided Pipeline
	if [[ $RUN_METHOD_1_HISAT2_REF_GUIDED == "TRUE" ]]; then
		log_step "STEP 02a: HISAT2 Reference-Guided Pipeline"
		log_warn "HISAT2 Reference-Guided pipeline requires GTF file parameter"
		if hisat2_ref_guided_pipeline --FASTA "$fasta" --GTF "$gtf_file" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 1 completed successfully"
		else
			log_error "Method 1 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 2: HISAT2 De Novo Pipeline (Main method)
	if [[ $RUN_METHOD_2_HISAT2_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 02b: HISAT2 De Novo Pipeline"
		if hisat2_de_novo_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 2 completed successfully"
		else
			log_error "Method 2 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 3: Trinity De Novo Pipeline
	if [[ $RUN_METHOD_3_TRINITY_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 03: Trinity De Novo Assembly and Quantification"
		if trinity_de_novo_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 3 completed successfully"
		else
			log_error "Method 3 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 4: Salmon SAF Quantification
	if [[ $RUN_METHOD_4_SALMON_SAF == "TRUE" ]]; then
		log_step "STEP 04: Salmon SAF Quantification"
		local genome_file="$decoy"
		if [[ -f "$genome_file" ]]; then
			if salmon_saf_pipeline --FASTA "$fasta" --GENOME "$genome_file" --RNASEQ_LIST "${rnaseq_list[@]}"; then
				log_info "Method 4 completed successfully"
			else
				log_error "Method 4 failed (exit code: $?) - continuing with remaining methods"
			fi
		else
			log_warn "Genome file '$genome_file' not found. Skipping Salmon SAF pipeline."
			log_warn "Please provide genome file for decoy-aware Salmon quantification."
		fi
	fi

	# Method 5: Bowtie2 + RSEM Quantification
	if [[ $RUN_METHOD_5_BOWTIE2_RSEM == "TRUE" ]]; then
		log_step "STEP 05: Bowtie2 + RSEM Quantification"
		if bowtie2_rsem_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 5 completed successfully"
		else
			log_error "Method 5 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Generate cross-method validation summary
	local fasta_base="$(basename "$fasta")"
	local fasta_tag="${fasta_base%.*}"
	compare_methods_summary "$fasta_tag"

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

if [[ $RUN_HEATMAP_WRAPPER == "TRUE" ]]; then
	log_step "Heatmap Wrapper post-processing enabled"
	
	# Navigate to post-processing directory
	if [[ ! -d "4_POST_PROCESSING" ]]; then
		log_error "Directory '4_POST_PROCESSING' not found"
	else
		cd 4_POST_PROCESSING || {
			log_error "Failed to change to 4_POST_PROCESSING directory"
			exit 1
		}
		
		# Execute the Heatmap Wrapper script for post-processing
		if [[ -f "run_all_post_processing.sh" ]]; then
			log_step "Executing Heatmap Wrapper post-processing script"
			chmod +x ./*.sh
			chmod +x run_all_post_processing.sh
			
			if bash "run_all_post_processing.sh" 2>&1; then
				log_info "Heatmap Wrapper completed successfully"
			else
				local exit_code=$?
				log_error "Heatmap Wrapper failed with exit code $exit_code"
			fi
		else
			log_warn "Heatmap Wrapper script 'run_all_post_processing.sh' not found - skipping"
		fi
		
		# Return to original directory
		cd - > /dev/null || log_warn "Failed to return to previous directory"
	fi
fi

# ==============================================================================
# POST-PROCESSING OPTIONS (COMMENTED OUT)
# ==============================================================================

if [[ $RUN_ZIP_RESULTS == "TRUE" ]]; then
	# Optional: Archive StringTie results for sharing or backup
	#tar -czvf "stringtie_results_$(date +%Y%m%d_%H%M%S).tar.gz" "$STRINGTIE_HISAT2_DE_NOVO_ROOT"
	#tar -czvf HISAT2_DE_NOVO_ROOT_HPC_$(date +%Y%m%d_%H%M%S).tar.gz $HISAT2_DE_NOVO_ROOT
	#tar -czvf 4b_Method_2_HISAT2_De_Novo_$(date +%Y%m%d_%H%M%S).tar.gz 4b_Method_2_HISAT2_De_Novo/
	log_step "Creating compressed archive for folders: 4_POST_PROCESSING, 4_OUTPUTS, and logs"
	TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
	tar -czf "CMSC244_${TIMESTAMP}.tar.gz" 4_POST_PROCESSING 4_OUTPUTS logs
	log_info "Archive created: CMSC244_${TIMESTAMP}.tar.gz"
fi

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
echo "END OF SCRIPT"
