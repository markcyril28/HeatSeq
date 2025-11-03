# ==============================================================================
# RSEM TO DESEQ2 PREPARATION SCRIPT (METHOD 5)
# ==============================================================================
# Description: This R script processes RSEM quantification files (`.genes.results`)
#              generated from Bowtie2 + RSEM Quantification (Method 5) and prepares
#              them for differential expression analysis using DESeq2.
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: October 2025
#
# Input: RSEM `.genes.results` files for each sample.
# Output: A combined gene count matrix and a sample information (colData) table,
#         both in TSV format, suitable for direct import into DESeq2.
#
# Required packages: dplyr, tibble, readr
# ==============================================================================

# Load required libraries
supressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory for the project
BASE_DIR <- getwd()

# Directory where RSEM quantification results for Method 5 are stored
RSEM_QUANT_ROOT <- file.path(BASE_DIR, "4e_Method_5_Bowtie2_Quantification", "quant")

# Output directory for DESeq2 input files
DESEQ2_OUTPUT_DIR <- file.path(BASE_DIR, "4e_Method_5_Bowtie2_Quantification", "DESeq2_input")

# List of FASTA tags used in the analysis (from ALL_FASTA_FILES in bash script)
FASTA_TAGS <- c("TEST") # Example: c("All_Smel_Genes", "TEST")

# List of SRR IDs used in the analysis (from SRR_COMBINED_LIST in bash script)
SRR_IDS <- c("SRR3884686", "SRR3884687", "SRR3884597")

# Create output directory if it doesn't exist
dir.create(DESEQ2_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to read and process a single RSEM results file
read_rsem_results <- function(file_path, srr_id) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Read the TSV file
  data <- read_tsv(file_path, show_col_types = FALSE)
  
  # RSEM .genes.results typically has columns like 'gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM'.
  # For DESeq2, we need integer counts, so 'expected_count' is the most appropriate.
  data_filtered <- data %>%
    dplyr::select(gene_id, expected_count) %>%
    dplyr::rename(!!srr_id := expected_count) # Rename expected_count column to SRR ID
  
  return(data_filtered)
}

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("RSEM TO DESEQ2 PREPARATION (METHOD 5)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

all_count_matrices <- list()

for (fasta_tag in FASTA_TAGS) {
  cat(paste0("Processing FASTA tag: ", fasta_tag, "\n"))
  
  # Initialize a list to hold data frames for current fasta_tag
  current_fasta_matrices <- list()
  
  for (srr_id in SRR_IDS) {
    file_name <- paste0(srr_id, ".genes.results")
    file_path <- file.path(RSEM_QUANT_ROOT, fasta_tag, srr_id, file_name)
    
    cat(paste0("  Reading: ", file_path, "\n"))
    srr_data <- read_rsem_results(file_path, srr_id)
    
    if (!is.null(srr_data)) {
      current_fasta_matrices[[srr_id]] <- srr_data
    }
  }
  
  # Merge all data frames for the current fasta_tag
  if (length(current_fasta_matrices) > 0) {
    # Use full_join to ensure all genes are included, filling missing values with 0
    merged_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), current_fasta_matrices)
    merged_matrix[is.na(merged_matrix)] <- 0
    
    # Set gene_id as row names and remove the column
    merged_matrix <- merged_matrix %>%
      column_to_rownames(var = "gene_id")
    
    # Ensure all counts are integers for DESeq2
    merged_matrix <- round(merged_matrix)
    
    all_count_matrices[[fasta_tag]] <- merged_matrix
    
    # Save the combined count matrix
    output_matrix_file <- file.path(DESEQ2_OUTPUT_DIR, paste0("gene_count_matrix_rsem_", fasta_tag, ".tsv"))
    write_tsv(as_tibble(merged_matrix, rownames = "gene_id"), output_matrix_file)
    cat(paste0("  Saved combined count matrix to: ", output_matrix_file, "\n"))
  } else {
    cat(paste0("  No data found for FASTA tag: ", fasta_tag, "\n"))
  }
}

# Create sample information (colData) table
cat("Creating sample information table...\n")

# This is a placeholder. Users will need to customize 'condition' based on their experimental design.
# For now, we'll create a dummy condition based on SRR ID order.

sample_data <- tibble(sample = SRR_IDS) %>%
  mutate(condition = factor(rep(c("cond1", "cond2"), length.out = length(SRR_IDS)))) # Example conditions

output_coldata_file <- file.path(DESEQ2_OUTPUT_DIR, "sample_info_rsem.tsv")
write_tsv(sample_data, output_coldata_file)
cat(paste0("Saved sample information table to: ", output_coldata_file, "\n"))

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("RSEM DESEQ2 PREPARATION COMPLETED.\n")
cat("Please review 'sample_info_rsem.tsv' and update 'condition' column based on your experimental design.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

