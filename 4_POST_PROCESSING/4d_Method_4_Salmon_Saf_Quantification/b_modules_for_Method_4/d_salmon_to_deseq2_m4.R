# ==============================================================================
# SALMON TO DESEQ2 PREPARATION SCRIPT (METHOD 4)
# ==============================================================================
# Description: This R script processes Salmon quantification files, imports them
#              using tximeta, and prepares a DESeqDataSet object for downstream
#              differential expression analysis.
#
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: October 2025
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tximeta)
  library(DESeq2)
  library(dplyr)
  library(readr)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# This script assumes the following directory structure:
# ./quant/<SRR_ID>/quant.sf

# Get SRR IDs from the directory names in 'quant'
srr_ids <- list.dirs("quant", full.names = FALSE, recursive = FALSE)
if (length(srr_ids) == 0) {
  stop("No sample directories found in 'quant'. Please ensure your Salmon output is organized as 'quant/<SRR_ID>/quant.sf'")
}

# Create a sample information table (colData)
coldata <- data.frame(
  sample = srr_ids,
  condition = factor(rep(c("cond1", "cond2"), length.out = length(srr_ids))), # Placeholder
  stringsAsFactors = FALSE
)

# Create file paths to quant.sf files
files <- file.path("quant", coldata$sample, "quant.sf")
names(files) <- coldata$sample
coldata$files <- files
coldata$names <- coldata$sample

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SALMON TO DESEQ2 PREPARATION (METHOD 4)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# --- 1. Import Salmon data with tximeta ---
cat("Importing Salmon quantification data with tximeta...\n")

# Use tximeta to import data and summarize to gene level
# This requires a transcript-to-gene mapping. If not available, tximeta will
# attempt to find one automatically. For this example, we assume it can.
se <- tximeta(coldata)
gse <- summarizeToGene(se)

# --- 2. Create DESeqDataSet ---
cat("Creating DESeqDataSet object...\n")

dds <- DESeqDataSet(gse, design = ~ condition)

# --- 3. Save DESeq2 input files ---
cat("Saving DESeq2 input files...\n")

DESEQ2_INPUT_DIR <- "matrices"
dir.create(DESEQ2_INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Save count matrix
count_matrix <- counts(dds)
output_matrix_file <- file.path(DESEQ2_INPUT_DIR, "gene_count_matrix_salmon.tsv")
write.table(count_matrix, file = output_matrix_file, sep = "\t", quote = FALSE, col.names = NA)
cat(paste0("Saved count matrix to: ", output_matrix_file, "\n"))

# Save sample info
output_coldata_file <- file.path(DESEQ2_INPUT_DIR, "sample_info_salmon.tsv")
write_tsv(as_tibble(colData(dds)), output_coldata_file)
cat(paste0("Saved sample information to: ", output_coldata_file, "\n"))

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("DESEQ2 PREPARATION COMPLETED.\n")
cat("Please review 'sample_info_salmon.tsv' and update 'condition' column.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")