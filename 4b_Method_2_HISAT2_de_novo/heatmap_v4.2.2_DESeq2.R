# Load required libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(DESeq2)
library(pheatmap)

# Load count matrix generated from prepDE.py
count_matrix <- read.csv("04_DESeq2/transcript_count_matrix.csv", row.names = 1)

# Ensure only numeric columns are retained
count_matrix <- count_matrix[sapply(count_matrix, is.numeric)]

# Round counts to nearest integers
count_matrix <- round(count_matrix)

# Construct dummy sample metadata (e.g., single condition)
sample_ids <- colnames(count_matrix)
col_data <- data.frame(row.names = sample_ids, condition = rep("sample", length(sample_ids)))

# Create DESeq2 dataset without a design formula (transformation only)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~1
)

# Filter out low-count transcripts (e.g., total counts < 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Apply variance stabilizing transformation (VST)
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

# Select top 50 most variable genes by row variance
top_var_genes <- head(order(rowVars(vst_matrix), decreasing = TRUE), 50)
vst_subset <- vst_matrix[top_var_genes, ]

# Transpose matrix to have samples as rows and genes as columns
vst_transposed <- t(vst_subset)

# Exclude specific genes by their MSTRG IDs (before suffix stripping)
excluded_genes <- c("MSTRG.9", "MSTRG.12", "MSTRG.5", "MSTRG.8")
gene_prefixes <- sub("\\.\\d+$", "", colnames(vst_transposed))
vst_transposed <- vst_transposed[, !gene_prefixes %in% excluded_genes]

# Normalize each gene column to [0, 10]
normalize_0_10 <- function(x) {
  if (max(x) == min(x)) return(rep(5, length(x)))
  (x - min(x)) / (max(x) - min(x)) * 10
}
vst_normalized <- apply(vst_transposed, 2, normalize_0_10)
rownames(vst_normalized) <- rownames(vst_transposed)
colnames(vst_normalized) <- colnames(vst_transposed)

# Sort genes by total normalized expression
gene_order <- order(colSums(vst_normalized))
vst_normalized <- vst_normalized[, gene_order]

# Updated gene ID to reference name map (excluding removed entries)
gene_id_map <- c(
  "MSTRG.1" = "SmelDMP01.990",
  "MSTRG.2" = "SmelDMP01.730",
  "MSTRG.3" = "SmelDMP04",
  "MSTRG.4" = "SmelDMP10.200",
  "MSTRG.6" = "SmelDMP02",
  "MSTRG.7" = "SmelDMP12",
  "MSTRG.10" = "Smel_18s_rRNA06.900",
  "MSTRG.11" = "Smel_Cyclophilin01"
)

# Rename columns using the updated gene map
new_colnames <- sapply(colnames(vst_normalized), function(gene_id) {
  prefix <- sub("\\.\\d+$", "", gene_id)
  if (prefix %in% names(gene_id_map)) gene_id_map[prefix] else gene_id
})
colnames(vst_normalized) <- new_colnames

# Define sample order (by SRR accession)
sample_order <- c(
  "SRR3884631", 
  "SRR3884664", 
  "SRR3884653", 
  "SRR3884677",
  "SRR3884679", 
  "SRR3884597", 
  "SRR3884687",
  "SRR3884686", 
  "SRR3884689", 
  "SRR3884690", 
  "SRR3884685", 
  "SRR3884675"
)
vst_normalized <- vst_normalized[sample_order, ]

# Assign human-readable sample labels
sample_labels <- c(
  "SRR3884686" = "Buds ¯ 0.7 cm /",
  "SRR3884687" = "Buds, Opened Buds /",
  "SRR3884677" = "Cotyledons",
  "SRR3884597" = "Flowers /",
  "SRR3884631" = "Fruits ¯ 6 cm",
  "SRR3884664" = "Fruits Calyx Stage 2",
  "SRR3884653" = "Fruits Flesh Stage 2",
  "SRR3884689" = "Leaves",
  "SRR3884679" = "Pistils /",
  "SRR3884685" = "Radicles",
  "SRR3884675" = "Roots",
  "SRR3884690" = "Stems"
)
rownames(vst_normalized) <- sample_labels[rownames(vst_normalized)]

# Define violet/eggplant-themed color palette
violet_palette <- colorRampPalette(c("#f2e6ff", "#b266ff", "#4b0082"))(100)

# Export heatmap as JPEG
jpeg("04_DESeq2/heatmap_DESeq2_V2.jpeg", width = 1200, height = 1000, res = 150)
pheatmap(
  vst_normalized,
  color = violet_palette,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90,           # Rotate column labels
  scale = "none",           # Use fixed 0–10 scale
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "In Silico Gene Expression Analysis of SmelDMPs\n(DESeq2 VST normalized, then Min-Max [0–10] Normalized and Sorted by Expression)"
)
dev.off()
