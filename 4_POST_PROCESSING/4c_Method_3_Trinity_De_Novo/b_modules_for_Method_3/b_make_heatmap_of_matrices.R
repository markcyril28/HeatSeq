# ==============================================================================
# BASIC HEATMAP GENERATION SCRIPT (METHOD 3)
# ==============================================================================
# Description: This R script takes a gene count matrix and sample information,
#              performs DESeq2 analysis, and generates a basic heatmap of
#              normalized gene expression values.
#
# Author: Mark Cyril R. Mercado
# Version: v2
# Date: October 2025
#
# Input:
#   - Gene count matrix (TSV format)
#   - Sample information table (TSV format)
#   - Gene group file (optional)
#
# Output:
#   - Heatmap image (PNG format)
#   - Data table with normalized values (TSV format)
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(dplyr)
  library(tibble)
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

# Input directory for DESeq2 files
DESEQ2_INPUT_DIR <- file.path(BASE_DIR, "7_DESeq2_input")

# Output directory for visualizations
# Output directory in 7_Heatmap_Outputs structure
VISUALIZATION_OUTPUT_DIR <- file.path(BASE_DIR, "7_Heatmap_Outputs", "I_Basic_Heatmap")
dir.create(VISUALIZATION_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Input files
COUNT_MATRIX_FILE <- file.path(DESEQ2_INPUT_DIR, paste0("gene_count_matrix_", FASTA_TAG, ".tsv"))
SAMPLE_INFO_FILE <- file.path(DESEQ2_INPUT_DIR, "sample_info.tsv")

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("BASIC HEATMAP GENERATION (METHOD 3)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# --- 1. Load Data ---
cat("Loading count matrix and sample info...\n")

if (!file.exists(COUNT_MATRIX_FILE)) {
  stop(paste("Count matrix file not found:", COUNT_MATRIX_FILE))
}
if (!file.exists(SAMPLE_INFO_FILE)) {
  stop(paste("Sample info file not found:", SAMPLE_INFO_FILE))
}

count_data <- read_tsv(COUNT_MATRIX_FILE, show_col_types = FALSE) %>%
  column_to_rownames("gene_id")

col_data <- read_tsv(SAMPLE_INFO_FILE, show_col_types = FALSE) %>%
  column_to_rownames("sample")

# Ensure row and column names match
if (!all(colnames(count_data) %in% rownames(col_data))) {
  stop("Column names of count matrix do not match row names of sample info.")
}
if (!all(colnames(count_data) == rownames(col_data))) {
    col_data <- col_data[colnames(count_data), , drop = FALSE]
}

# --- 2. DESeq2 Analysis ---
cat("Running DESeq2 analysis...\n")

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
normalized_counts <- assay(vsd)

# --- 3. Filter by Gene Group (if provided) ---
if (!is.null(GENE_GROUP_FILE)) {
  cat(paste("Filtering by gene group:", GENE_GROUP_FILE, "\n"))
  gene_group <- read.table(GENE_GROUP_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
  normalized_counts <- normalized_counts[rownames(normalized_counts) %in% gene_group, ]
}

# --- 4. Heatmap Generation ---
cat("Generating heatmap...\n")

scaled_data <- t(scale(t(normalized_counts)))
color_palette <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

gene_group_name <- if (!is.null(GENE_GROUP_FILE)) tools::file_path_sans_ext(basename(GENE_GROUP_FILE)) else "All_Genes"
output_prefix <- paste0(gene_group_name, "_", COUNT_TYPE, "_geneName_Organ_from_", FASTA_TAG, "_vst_normalized")
output_png <- file.path(VISUALIZATION_OUTPUT_DIR, paste0(output_prefix, ".png"))
output_tsv <- file.path(VISUALIZATION_OUTPUT_DIR, paste0(output_prefix, ".tsv"))

png(output_png, width = 8, height = 12, units = "in", res = 300)
Heatmap(scaled_data,
        name = "Z-score",
        col = color_palette,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        cluster_rows = TRUE,
        cluster_columns = TRUE)
dev.off()

cat(paste0("Saved heatmap to: ", output_png, "\n"))

# --- 5. Save Data ---
cat("Saving normalized data...\n")
write_tsv(as_tibble(scaled_data, rownames = "gene_id"), output_tsv)
cat(paste0("Saved data to: ", output_tsv, "\n"))

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("BASIC HEATMAP GENERATION COMPLETED.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")