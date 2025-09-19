#!/usr/bin/env Rscript

# ===============================================
# Libraries
# ===============================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
})

# ===============================================
# Configuration
# ===============================================
count_dir <- "/mnt/e/HeatSeq/count_matrices"   # <-- output of bash script
out_dir   <- "/mnt/e/HeatSeq/heatmaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================
# Helper function
# ===============================================
process_group <- function(file) {
  group_name <- gsub("_counts.tsv$", "", basename(file))
  message("Processing group: ", group_name)

  # Read count matrix
  counts <- read.delim(file, row.names = 1, check.names = FALSE)

  # Remove rows with all zeros
  counts <- counts[rowSums(counts) > 0, ]

  # Build dummy metadata (all samples same condition)
  # You should edit this with real experimental design
  samples <- colnames(counts)
  condition <- rep("ConditionA", length(samples))
  coldata <- data.frame(row.names = samples, condition)

  # DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData   = coldata,
                                design    = ~ condition)

  # Variance-stabilizing transform
  vsd <- vst(dds, blind = TRUE)
  mat <- assay(vsd)

  # Select top 50 most variable genes
  topVarGenes <- head(order(rowVars(mat), decreasing = TRUE), 50)

  # Heatmap
  pdf(file.path(out_dir, paste0(group_name, "_heatmap.pdf")),
      width = 8, height = 10)
  pheatmap(mat[topVarGenes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = coldata,
           show_rownames = FALSE,
           fontsize = 10,
           main = paste("Top variable genes -", group_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  dev.off()

  message("Saved heatmap: ", file.path(out_dir, paste0(group_name, "_heatmap.pdf")))
}

# ===============================================
# Main
# ===============================================
files <- list.files(count_dir, pattern = "_counts.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No *_counts.tsv files found in ", count_dir)
}

lapply(files, process_group)

message("All heatmaps saved in: ", out_dir)
