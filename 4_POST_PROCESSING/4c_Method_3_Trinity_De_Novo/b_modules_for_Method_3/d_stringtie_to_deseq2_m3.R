# ==============================================================================
# STRINGTIE TO DESEQ2 PREPARATION SCRIPT FOR TRINITY DE NOVO (METHOD 3)
# ==============================================================================
# Description: This R script processes StringTie gene abundance files generated
#              from Trinity de novo assembly (Method 3) and prepares them for
#              differential expression analysis using DESeq2.
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: October 2025
#
# Input: StringTie gene abundance TSV files (e.g., SRRXXXXXX_FASTA_TAG_trinity_gene_abundances.tsv)
# Output: A combined gene count matrix and a sample information (colData) table,
#         both in TSV format, suitable for direct import into DESeq2.
#
# Required packages: dplyr, tibble, readr
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# ==============================================================================
# COMMAND-LINE ARGUMENT PARSING
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript d_stringtie_to_deseq2_m3.R <base_dir> <fasta_tags_comma_separated> <srr_ids_comma_separated>", call. = FALSE)
}

BASE_DIR <- args[1]
FASTA_TAGS <- strsplit(args[2], ",")[[1]]
SRR_IDS <- strsplit(args[3], ",")[[1]]

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Directory where StringTie results for Method 3 are stored
STRINGTIE_M3_ROOT <- file.path(BASE_DIR, "5_stringtie_WD", "a_Method_3_RAW_RESULTs")

# Output directory for DESeq2 input files
# Path sync with modules/methods.sh: 6_matrices_from_stringtie
DESEQ2_OUTPUT_DIR <- file.path(BASE_DIR, "6_matrices_from_stringtie")


# Create output directory if it doesn't exist
dir.create(DESEQ2_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to read and process a single StringTie abundance file
read_stringtie_abundance <- function(file_path, srr_id) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Read the TSV file
  data <- read_tsv(file_path, show_col_types = FALSE)
  
  # Select relevant columns: gene_id and FPKM (or TPM, depending on what DESeq2 expects)
  # DESeq2 typically works with raw counts. StringTie's gene_abundances.tsv contains
  # 'FPKM' and 'TPM', but for DESeq2 we often need raw counts. If raw counts are not
  # directly available, FPKM/TPM can be used with appropriate DESeq2 settings or converted.
  # For simplicity, we'll use 'num_reads' if available, otherwise FPKM as a proxy.
  # NOTE: StringTie's gene_abundances.tsv usually contains 'TPM', 'FPKM', 'coverage', 'num_reads'.
  # 'num_reads' is the closest to raw counts. If not present, FPKM/TPM can be used but DESeq2
  # expects integer counts. A common workaround is to round FPKM/TPM or use a different tool
  # for raw count extraction (e.g., featureCounts on StringTie's BAMs).
  
  # Assuming 'num_reads' column is available and represents raw counts
  if ("num_reads" %in% colnames(data)) {
    data_filtered <- data %>%
      dplyr::select(gene_id, num_reads) %>%
      dplyr::rename(!!srr_id := num_reads) # Rename num_reads column to SRR ID
  } else if ("FPKM" %in% colnames(data)) {
    warning(paste("num_reads not found for", srr_id, ". Using FPKM as proxy for counts. Consider using featureCounts for raw counts."))
    data_filtered <- data %>%
      dplyr::select(gene_id, FPKM) %>%
      dplyr::rename(!!srr_id := FPKM) %>% # Rename FPKM column to SRR ID
      dplyr::mutate(!!srr_id := round(!!sym(srr_id))) # Round FPKM to nearest integer for DESeq2
  } else {
    stop(paste("Neither 'num_reads' nor 'FPKM' found in", file_path, ". Cannot proceed."))
  }
  
  return(data_filtered)
}

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("STRINGTIE TO DESEQ2 PREPARATION (METHOD 3)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

all_count_matrices <- list()

for (fasta_tag in FASTA_TAGS) {
  cat(paste0("Processing FASTA tag: ", fasta_tag, "\n"))
  
  # Initialize a list to hold data frames for current fasta_tag
  current_fasta_matrices <- list()
  
  for (srr_id in SRR_IDS) {
    file_name <- paste0(srr_id, "_", fasta_tag, "_trinity_gene_abundances.tsv")
    file_path <- file.path(STRINGTIE_M3_ROOT, fasta_tag, srr_id, file_name)
    
    cat(paste0("  Reading: ", file_path, "\n"))
    srr_data <- read_stringtie_abundance(file_path, srr_id)
    
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
    output_matrix_file <- file.path(DESEQ2_OUTPUT_DIR, paste0("gene_count_matrix_", fasta_tag, ".tsv"))
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

output_coldata_file <- file.path(DESEQ2_OUTPUT_DIR, "sample_info.tsv")
write_tsv(sample_data, output_coldata_file)
cat(paste0("Saved sample information table to: ", output_coldata_file, "\n"))

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("DESEQ2 PREPARATION COMPLETED.\n")
cat("Please review 'sample_info.tsv' and update 'condition' column based on your experimental design.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
