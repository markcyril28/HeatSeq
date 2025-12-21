#!/usr/bin/env Rscript
# Bar Graph for Method 1 (HISAT2 Ref-Guided)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(getopt)
})

spec <- matrix(c(
  'base_dir',   'd', 1, "character", "Base directory",
  'fasta_tag',  'f', 1, "character", "FASTA tag",
  'count_type', 'c', 1, "character", "Count type",
  'gene_group', 'g', 2, "character", "Gene group file"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

BASE_DIR <- opt$base_dir
FASTA_TAG <- opt$fasta_tag
GENE_GROUP_FILE <- opt$gene_group

setwd(BASE_DIR)

# Path sync with modules/methods.sh: 5_stringtie_WD/<fasta_tag>/deseq2_input
matrix_path <- file.path("..", "5_stringtie_WD", FASTA_TAG, "deseq2_input", "gene_count_matrix.csv")
sample_info_path <- file.path("..", "5_stringtie_WD", FASTA_TAG, "deseq2_input", "sample_metadata.csv")

if (!file.exists(matrix_path)) stop("Matrix not found: ", matrix_path)
if (!file.exists(sample_info_path)) stop("Sample info not found: ", sample_info_path)

count_data <- read.csv(matrix_path, row.names=1, check.names=FALSE)
sample_info <- read.csv(sample_info_path, row.names=1)

if (!is.null(GENE_GROUP_FILE) && file.exists(GENE_GROUP_FILE)) {
  gene_ids <- read.table(GENE_GROUP_FILE, header=FALSE, stringsAsFactors=FALSE)[,1]
  count_data <- count_data[rownames(count_data) %in% gene_ids, ]
}

dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_info, design=~condition)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized=TRUE)

# Output directory in 7_Heatmap_Outputs structure
output_dir <- file.path("..", "7_Heatmap_Outputs", "III_Bar_Graphs")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

for (gene in rownames(norm_counts)) {
  df <- data.frame(Sample=colnames(norm_counts), Count=norm_counts[gene,])
  p <- ggplot(df, aes(x=Sample, y=Count)) + geom_bar(stat="identity", fill="steelblue") +
       theme_minimal() + labs(title=gene) + theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(out_dir, paste0(gene, "_", FASTA_TAG, "_bar.png")), p, width=8, height=6)
}

cat("Bar graphs saved in:", out_dir, "\n")
