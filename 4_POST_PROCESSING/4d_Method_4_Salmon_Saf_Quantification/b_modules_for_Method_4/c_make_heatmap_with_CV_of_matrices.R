# ==============================================================================
# HEATMAP WITH CV GENERATION SCRIPT (METHOD 4)
# ==============================================================================
# Description: This R script generates a heatmap with coefficient of variation (CV)
#              for gene expression data after DESeq2 analysis.
#
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: October 2025
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
  'gene_group', 'g', 2, "character", "Name of the gene group file (optional)"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

GENE_GROUP_FILE <- opt$gene_group

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DESEQ2_INPUT_DIR <- "matrices"
VISUALIZATION_OUTPUT_DIR <- "9_Heatmap_with_Covariability_Visualizations"
dir.create(VISUALIZATION_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

COUNT_MATRIX_FILE <- file.path(DESEQ2_INPUT_DIR, "gene_count_matrix_salmon.tsv")
SAMPLE_INFO_FILE <- file.path(DESEQ2_INPUT_DIR, "sample_info_salmon.tsv")

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("Loading data...\n")
count_data <- as.matrix(read.table(COUNT_MATRIX_FILE, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
col_data <- read_tsv(SAMPLE_INFO_FILE, show_col_types = FALSE) %>%
  filter(sample %in% colnames(count_data)) %>%
  column_to_rownames("sample")

if (!all(colnames(count_data) == rownames(col_data))) {
    col_data <- col_data[colnames(count_data), , drop = FALSE]
}

cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = round(count_data), colData = col_data, design = ~ condition)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

if (!is.null(GENE_GROUP_FILE)) {
  cat(paste("Filtering by gene group:", GENE_GROUP_FILE, "\n"))
  gene_group <- read.table(GENE_GROUP_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
  normalized_counts <- normalized_counts[rownames(normalized_counts) %in% gene_group, ]
}

cat("Calculating CV and generating heatmap...\n")
cv_values <- apply(normalized_counts, 1, function(x) sd(x) / mean(x))
scaled_data <- t(scale(t(normalized_counts)))
color_palette <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ha <- rowAnnotation(CV = anno_barplot(cv_values))
gene_group_name <- if (!is.null(GENE_GROUP_FILE)) tools::file_path_sans_ext(basename(GENE_GROUP_FILE)) else "All_Genes"
output_prefix <- paste0(gene_group_name, "_salmon_cv_heatmap")
output_png <- file.path(VISUALIZATION_OUTPUT_DIR, paste0(output_prefix, ".png"))

png(output_png, width = 10, height = 12, units = "in", res = 300)
Heatmap(scaled_data, name = "Z-score", col = color_palette, cluster_rows = TRUE, cluster_columns = TRUE, right_annotation = ha)
dev.off()

cat(paste0("Saved CV heatmap to: ", output_png, "\n"))
