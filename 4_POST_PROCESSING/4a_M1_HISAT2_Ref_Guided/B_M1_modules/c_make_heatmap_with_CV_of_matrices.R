#!/usr/bin/env Rscript
# CV Heatmap for Method 1 (HISAT2 Ref-Guided)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
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

cv_values <- apply(norm_counts, 1, function(x) (sd(x)/mean(x))*100)

# Output directory in 7_Heatmap_Outputs structure
output_dir <- file.path("..", "7_Heatmap_Outputs", "II_Heatmap_with_CV")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

group_name <- if (!is.null(GENE_GROUP_FILE)) tools::file_path_sans_ext(basename(GENE_GROUP_FILE)) else "all_genes"
out_file <- file.path(out_dir, paste0(group_name, "_", FASTA_TAG, "_CV_heatmap.png"))

png(out_file, width=900, height=600)
Heatmap(log2(norm_counts + 1), name="log2(count+1)",
        col=colorRamp2(c(0, 5, 10), c("blue", "white", "red")),
        row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=10),
        right_annotation=rowAnnotation(CV=anno_barplot(cv_values)))
dev.off()

cat("CV Heatmap saved:", out_file, "\n")
