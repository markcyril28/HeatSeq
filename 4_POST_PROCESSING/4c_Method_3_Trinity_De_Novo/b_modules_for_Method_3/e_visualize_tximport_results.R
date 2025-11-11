#!/usr/bin/env Rscript
# ==============================================================================
# VISUALIZATION SCRIPT FOR TRINITY + SALMON TXIMPORT RESULTS
# ==============================================================================
# Description: Generate heatmaps and bar graphs from tximport results
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: November 2025
#
# Input: Tximport RDS file and gene group lists
# Output: Heatmaps and bar graphs in 7_Heatmap_Outputs/
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(getopt)
})

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# ==============================================================================

spec <- matrix(c(
  'base_dir',    'd', 1, "character", "Base directory (4c_Method_3_Trinity_De_Novo)",
  'fasta_tag',   'f', 1, "character", "FASTA tag",
  'gene_group',  'g', 1, "character", "Gene group file path",
  'output_type', 't', 2, "character", "Output type: heatmap, cv_heatmap, bargraph (default: all)"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (is.null(opt$base_dir) || is.null(opt$fasta_tag) || is.null(opt$gene_group)) {
  stop("Usage: Rscript e_visualize_tximport_results.R -d <base_dir> -f <fasta_tag> -g <gene_group_file>")
}

BASE_DIR <- opt$base_dir
FASTA_TAG <- opt$fasta_tag
GENE_GROUP_FILE <- opt$gene_group
OUTPUT_TYPE <- if (!is.null(opt$output_type)) opt$output_type else "all"

# ==============================================================================
# CONFIGURATION
# ==============================================================================

MATRIX_DIR <- file.path(BASE_DIR, "6_matrices_from_stringtie")
TXIMPORT_RDS <- file.path(MATRIX_DIR, "tximport_trinity_salmon.rds")
DESEQ2_RDS <- file.path(MATRIX_DIR, "deseq2_dataset_trinity.rds")
COUNT_MATRIX <- file.path(MATRIX_DIR, "gene_counts_tximport.tsv")
SAMPLE_INFO <- file.path(MATRIX_DIR, "sample_info.tsv")

OUTPUT_HEATMAP_DIR <- file.path(BASE_DIR, "7_Heatmap_Outputs", "I_Basic_Heatmap")
OUTPUT_CV_DIR <- file.path(BASE_DIR, "7_Heatmap_Outputs", "II_Heatmap_with_CV")
OUTPUT_BARGRAPH_DIR <- file.path(BASE_DIR, "7_Heatmap_Outputs", "III_Bar_Graphs")

dir.create(OUTPUT_HEATMAP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_CV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_BARGRAPH_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading tximport results...\n")

if (!file.exists(DESEQ2_RDS)) {
  stop("DESeq2 dataset not found. Run tximport first: ", DESEQ2_RDS)
}

dds <- readRDS(DESEQ2_RDS)
dds <- DESeq(dds, quiet = TRUE)

# Load gene groups
gene_list <- read.table(GENE_GROUP_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
gene_group_name <- tools::file_path_sans_ext(basename(GENE_GROUP_FILE))

cat("Gene group:", gene_group_name, "\n")
cat("Genes to visualize:", length(gene_list), "\n")

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Filter for genes of interest
genes_in_data <- intersect(gene_list, rownames(normalized_counts))
if (length(genes_in_data) == 0) {
  stop("No genes from gene group found in count matrix")
}

cat("Genes found in data:", length(genes_in_data), "\n")

filtered_counts <- normalized_counts[genes_in_data, , drop = FALSE]

# ==============================================================================
# BASIC HEATMAP
# ==============================================================================

if (OUTPUT_TYPE %in% c("all", "heatmap")) {
  cat("Generating basic heatmap...\n")
  
  # Log2 transform
  log_counts <- log2(filtered_counts + 1)
  
  # Z-score normalization
  mat_scaled <- t(scale(t(log_counts)))
  
  output_file <- file.path(OUTPUT_HEATMAP_DIR, 
                           paste0(gene_group_name, "_", FASTA_TAG, "_heatmap.png"))
  
  png(output_file, width = 1200, height = 800, res = 120)
  
  ht <- Heatmap(mat_scaled,
                name = "Z-score",
                col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = paste("Gene Expression:", gene_group_name),
                row_names_gp = gpar(fontsize = 8))
  
  draw(ht)
  dev.off()
  
  cat("Saved heatmap:", output_file, "\n")
}

# ==============================================================================
# CV HEATMAP
# ==============================================================================

if (OUTPUT_TYPE %in% c("all", "cv_heatmap")) {
  cat("Generating CV heatmap...\n")
  
  # Calculate coefficient of variation
  cv_data <- apply(filtered_counts, 1, function(x) {
    sd(x) / mean(x) * 100
  })
  
  cv_df <- data.frame(
    Gene = names(cv_data),
    CV = cv_data,
    Mean_Expression = rowMeans(filtered_counts)
  )
  
  # Add CV to matrix
  mat_with_cv <- cbind(filtered_counts, CV = cv_data)
  mat_log <- log2(mat_with_cv + 1)
  
  output_file <- file.path(OUTPUT_CV_DIR,
                           paste0(gene_group_name, "_", FASTA_TAG, "_cv_heatmap.png"))
  
  png(output_file, width = 1400, height = 800, res = 120)
  
  ht <- Heatmap(mat_log[, -ncol(mat_log)],
                name = "log2(count+1)",
                col = colorRamp2(c(0, 5, 10), c("white", "yellow", "red")),
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = paste("Gene Expression with CV:", gene_group_name),
                row_names_gp = gpar(fontsize = 8)) +
       Heatmap(mat_log[, ncol(mat_log), drop = FALSE],
               name = "CV %",
               col = colorRamp2(c(0, 50, 100), c("green", "yellow", "red")),
               width = unit(1, "cm"),
               show_column_names = TRUE)
  
  draw(ht)
  dev.off()
  
  cat("Saved CV heatmap:", output_file, "\n")
}

# ==============================================================================
# BAR GRAPHS
# ==============================================================================

if (OUTPUT_TYPE %in% c("all", "bargraph")) {
  cat("Generating bar graphs...\n")
  
  for (gene_id in genes_in_data) {
    gene_data <- data.frame(
      Sample = colnames(filtered_counts),
      Expression = as.numeric(filtered_counts[gene_id, ])
    )
    
    output_file <- file.path(OUTPUT_BARGRAPH_DIR,
                             paste0(gene_id, "_", FASTA_TAG, "_bargraph.png"))
    
    p <- ggplot(gene_data, aes(x = Sample, y = Expression)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Expression of", gene_id),
           x = "Sample",
           y = "Normalized Count") +
      geom_hline(yintercept = mean(gene_data$Expression), 
                 linetype = "dashed", color = "red")
    
    ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  }
  
  cat("Saved", length(genes_in_data), "bar graphs to:", OUTPUT_BARGRAPH_DIR, "\n")
}

cat("\nVisualization completed successfully!\n")
