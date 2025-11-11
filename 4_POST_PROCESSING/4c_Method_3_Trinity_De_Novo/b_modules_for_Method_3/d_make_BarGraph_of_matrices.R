# ==============================================================================
# BAR GRAPH GENERATION SCRIPT (METHOD 3)
# ==============================================================================
# Description: This R script generates bar graphs for individual genes from
#              normalized expression data.
#
# Author: Mark Cyril R. Mercado
# Version: v2
# Date: October 2025
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(getopt)
})

# ==============================================================================
# COMMAND-LINE ARGUMENT PARSING
# ==============================================================================

spec <- matrix(c(
  'base_dir',   'd', 1, "character", "Base directory of the project",
  'fasta_tag',  'f', 1, "character", "FASTA tag for the analysis",
  'count_type', 'c', 1, "character", "Count type (e.g., 'counts')",
  'gene_group', 'g', 2, "character", "Name of the gene group file (optional)"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

BASE_DIR <- opt$base_dir
FASTA_TAG <- opt$fasta_tag
COUNT_TYPE <- opt$count_type
GENE_GROUP_FILE <- opt$gene_group

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DESEQ2_INPUT_DIR <- file.path(BASE_DIR, "7_DESeq2_input")
# Output directory in 7_Heatmap_Outputs structure
VISUALIZATION_OUTPUT_DIR <- file.path(BASE_DIR, "7_Heatmap_Outputs", "III_Bar_Graphs")
dir.create(VISUALIZATION_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

COUNT_MATRIX_FILE <- file.path(DESEQ2_INPUT_DIR, paste0("gene_count_matrix_", FASTA_TAG, ".tsv"))
SAMPLE_INFO_FILE <- file.path(DESEQ2_INPUT_DIR, "sample_info.tsv")

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("Loading data...\n")
count_data <- read_tsv(COUNT_MATRIX_FILE, show_col_types = FALSE) %>%
  column_to_rownames("gene_id")
col_data <- read_tsv(SAMPLE_INFO_FILE, show_col_types = FALSE) %>%
  column_to_rownames("sample")
if (!all(colnames(count_data) == rownames(col_data))) {
    col_data <- col_data[colnames(count_data), , drop = FALSE]
}

cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

if (!is.null(GENE_GROUP_FILE)) {
  cat(paste("Filtering by gene group:", GENE_GROUP_FILE, "\n"))
  gene_group <- read.table(GENE_GROUP_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
  normalized_counts <- normalized_counts[rownames(normalized_counts) %in% gene_group, ]
}

cat("Generating bar graphs...\n")
gene_group_name <- if (!is.null(GENE_GROUP_FILE)) tools::file_path_sans_ext(basename(GENE_GROUP_FILE)) else "All_Genes"
output_dir_group <- file.path(VISUALIZATION_OUTPUT_DIR, gene_group_name)
dir.create(output_dir_group, recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(normalized_counts)) {
  gene_name <- rownames(normalized_counts)[i]
  df <- data.frame(
    sample = colnames(normalized_counts),
    counts = normalized_counts[i, ]
  )
  
  p <- ggplot(df, aes(x = sample, y = counts, fill = sample)) +
    geom_bar(stat = "identity") +
    labs(title = gene_name, x = "Sample", y = "Normalized Counts") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  output_file <- file.path(output_dir_group, paste0(gene_name, "_bargraph.png"))
  ggsave(output_file, plot = p, width = 8, height = 6, units = "in")
}

cat(paste0("Bar graphs saved to: ", output_dir_group, "\n"))