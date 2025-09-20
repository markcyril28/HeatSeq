# Load required libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

library(DESeq2)
library(pheatmap)

# Load count matrix from prepDE.py
count_matrix <- read.csv("04_DESeq2/transcript_count_matrix.csv", row.names = 1)

# Remove non-numeric columns if present
count_matrix <- count_matrix[sapply(count_matrix, is.numeric)]

# Round counts to integers
count_matrix <- round(count_matrix)

# Create dummy metadata assuming all samples are in one condition (no DE analysis)
sample_ids <- colnames(count_matrix)
col_data <- data.frame(row.names = sample_ids, condition = rep("sample", length(sample_ids)))

# Construct DESeq2 dataset (no differential expression design)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~1
)

# Filter low count genes (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Apply variance stabilizing transformation (recommended for small datasets)
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract transformed expression matrix
vst_matrix <- assay(vst_data)

# Select top 50 most variable genes for heatmap
top_var_genes <- head(order(rowVars(vst_matrix), decreasing = TRUE), 50)
vst_subset <- vst_matrix[top_var_genes, ]

# Export heatmap as JPEG
jpeg("04_DESeq2/heatmap_top50.jpeg", width = 1200, height = 1000, res = 150)
pheatmap(
  vst_subset,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Most Variable Genes"
)
dev.off()

