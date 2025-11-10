#!/usr/bin/env Rscript

# ===============================================
# RSEM QUANTIFICATION TO MATRICES - METHOD 5
# ===============================================
# Processes RSEM quantification output directly
# RSEM provides properly normalized TPM and FPKM values
# NO DESeq2 needed - RSEM normalization is already correct

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ===============================================
# CONFIGURATION
# ===============================================

QUANT_DIR <- "3_quant"
MASTER_REFERENCE <- "All_Smel_Genes"
MATRICES_OUTPUT_DIR <- "4_matrices"
GENE_GROUPS_DIR <- "2_gene_groups"

# Sample information
SAMPLE_IDS <- c(
  "SRR3884675", "SRR3884690", "SRR3884689", "SRR3884684",
  "SRR3884686", "SRR3884687", "SRR3884597"
)

SAMPLE_LABELS <- c(
    # Roots
    "SRR3884675" = "Roots",      # PRJNA328564
    #"SRR20722229" = "Roots_2",     # SAMN28540077
    #"SRR31755282" = "Roots_3",     # SAMN28540068

    # Stems
    "SRR3884690" = "Stems",      # PRJNA328564
    #"SRR20722227" = "Stems_2",     # SAMN28540077
    #"SRR20722384" = "Stems_3",     # SAMN28540068

    # Leaves
    "SRR3884689" = "Leaves",     # PRJNA328564
    #"SRR20722230" = "Leaves_2",    # SAMN28540077
    #"SRR20722386" = "Leaves_3",    # SAMN28540068
    "SRR3884684" = "Senescent_leaves", # PRJNA328564

    # Buds
    "SRR3884686" = "Buds",       # PRJNA328564
    #"SRR21010466" = "Buds_2",      # SAMN28540077
    #"SRR20722297" = "Buds_3",      # SAMN28540068

    # Opened Buds
    "SRR3884687" = "Opened_Buds", # PRJNA328564

    # Flowers
    "SRR3884597" = "Flowers",    # PRJNA328564
    #"SRR20722234" = "Flowers_2",   # SAMN28540077
    #"SRR23909863" = "Flowers_3",   # SAMN28540068

    # Fruits
    "SRR3884631" = "Fruits",     # PRJNA328564
    #"SRR2072232" = "Fruits_2",     # SAMN28540077
    #"SRR20722387" = "Fruits_3",    # SAMN28540068
    "SRR3884608" = "Fruits_1cm",   # PRJNA328564
    "SRR3884620" = "Fruits_Stage_1", # PRJNA328564
    "SRR3884642" = "Fruits_Skin_Stage_2", # PRJNA328564
    "SRR3884653" = "Fruits_Flesh_Stage_2", # PRJNA328564
    "SRR3884664" = "Fruits_Calyx_Stage_2", # PRJNA328564
    "SRR3884680" = "Fruits_Skin_Stage_3", # PRJNA328564
    "SRR3884681" = "Fruits_Flesh_Stage_3", # PRJNA328564
    "SRR3884678" = "Fruits_peduncle", # PRJNA328564

    # Other organs
    "SRR3884685" = "Radicles",     # PRJNA328564
    "SRR3884677" = "Cotyledons",   # PRJNA328564
    "SRR3884679" = "Pistils"       # PRJNA328564
)


# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("RSEM QUANTIFICATION TO COUNT MATRICES - METHOD 5\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("RSEM provides properly normalized values:\n")
cat("  • TPM (Transcripts Per Million) - Normalized, cross-sample comparable\n")
cat("  • FPKM (Fragments Per Kilobase Million) - Normalized for length & depth\n")
cat("  • expected_count - Raw probabilistic counts\n\n")

cat("NO DESeq2 needed - RSEM normalization is statistically sound\n\n")

# ===============================================
# STEP 1: LOCATE RSEM OUTPUT FILES
# ===============================================

cat("Step 1: Locating RSEM output files...\n")

rsem_quant_dir <- file.path(QUANT_DIR, MASTER_REFERENCE)

# Build paths to RSEM .genes.results files
files <- file.path(rsem_quant_dir, SAMPLE_IDS, paste0(SAMPLE_IDS, ".genes.results"))
names(files) <- SAMPLE_IDS

# Check which files exist
files_exist <- file.exists(files)
if (sum(files_exist) == 0) {
  stop("ERROR: No RSEM quantification files found in ", rsem_quant_dir)
}

if (sum(files_exist) < length(files)) {
  missing <- SAMPLE_IDS[!files_exist]
  cat("WARNING: Missing quantification for samples:", paste(missing, collapse = ", "), "\n")
  files <- files[files_exist]
  SAMPLE_IDS <- SAMPLE_IDS[files_exist]
}

cat("Found", length(files), "RSEM quantification files\n")
cat("Samples:", paste(SAMPLE_IDS, collapse = ", "), "\n\n")

# ===============================================
# STEP 2: READ RSEM RESULTS
# ===============================================

cat("Step 2: Reading RSEM quantification results...\n")

# Function to read single RSEM results file
read_rsem_results <- function(file_path) {
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(data)
}

# Read all files
rsem_data_list <- list()
for (sample_id in names(files)) {
  rsem_data_list[[sample_id]] <- read_rsem_results(files[sample_id])
  cat("  Read:", sample_id, "-", nrow(rsem_data_list[[sample_id]]), "genes\n")
}

# Get gene IDs from first sample
gene_ids <- rsem_data_list[[1]]$gene_id

cat("\nTotal genes:", length(gene_ids), "\n\n")

# ===============================================
# STEP 3: BUILD COUNT MATRICES
# ===============================================

cat("Step 3: Building count matrices...\n")

# Initialize matrices
expected_count_matrix <- matrix(NA, nrow = length(gene_ids), ncol = length(SAMPLE_IDS))
tpm_matrix <- matrix(NA, nrow = length(gene_ids), ncol = length(SAMPLE_IDS))
fpkm_matrix <- matrix(NA, nrow = length(gene_ids), ncol = length(SAMPLE_IDS))

rownames(expected_count_matrix) <- gene_ids
rownames(tpm_matrix) <- gene_ids
rownames(fpkm_matrix) <- gene_ids

colnames(expected_count_matrix) <- SAMPLE_IDS
colnames(tpm_matrix) <- SAMPLE_IDS
colnames(fpkm_matrix) <- SAMPLE_IDS

# Fill matrices
for (i in seq_along(SAMPLE_IDS)) {
  sample_id <- SAMPLE_IDS[i]
  data <- rsem_data_list[[sample_id]]

  # Match gene IDs
  match_idx <- match(gene_ids, data$gene_id)

  expected_count_matrix[, i] <- data$expected_count[match_idx]
  tpm_matrix[, i] <- data$TPM[match_idx]
  fpkm_matrix[, i] <- data$FPKM[match_idx]
}

cat("Built 3 count matrices:\n")
cat("  • expected_count:", nrow(expected_count_matrix), "x", ncol(expected_count_matrix), "\n")
cat("  • TPM:", nrow(tpm_matrix), "x", ncol(tpm_matrix), "\n")
cat("  • FPKM:", nrow(fpkm_matrix), "x", ncol(fpkm_matrix), "\n\n")

# ===============================================
# HELPER FUNCTION: SAVE MATRIX
# ===============================================

save_matrix_file <- function(matrix_data, output_path, label_type = "SRR") {
  # Create dataframe
  df <- as.data.frame(matrix_data)

  # Relabel columns if needed
  if (label_type == "Organ") {
    colnames(df) <- SAMPLE_LABELS[colnames(df)]
  }

  # Add gene ID column
  df <- tibble::rownames_to_column(df, "GeneID")

  # Write to file
  write.table(df, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ===============================================
# STEP 4: SAVE FULL GENOME MATRICES
# ===============================================

cat("Step 4: Saving full genome count matrices...\n")

count_types <- c("expected_count", "TPM", "FPKM")
matrices_list <- list(
  expected_count = expected_count_matrix,
  TPM = tpm_matrix,
  FPKM = fpkm_matrix
)

for (count_type in count_types) {
  matrix_data <- matrices_list[[count_type]]

  # Save with SRR IDs
  output_file_srr <- file.path(output_dir,
    paste0(MASTER_REFERENCE, "_", count_type, "_geneID_from_", MASTER_REFERENCE, ".tsv"))
  save_matrix_file(matrix_data, output_file_srr, label_type = "SRR")
  cat("  Saved:", basename(output_file_srr), "\n")

  # Save with Organ labels
  output_file_organ <- file.path(output_dir,
    paste0(MASTER_REFERENCE, "_", count_type, "_geneName_from_", MASTER_REFERENCE, ".tsv"))
  save_matrix_file(matrix_data, output_file_organ, label_type = "Organ")
  cat("  Saved:", basename(output_file_organ), "\n")
}

cat("\n")

# ===============================================
# STEP 5: PROCESS GENE GROUPS
# ===============================================

cat("Step 5: Processing gene groups...\n")

# Find all gene group files
gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.tsv$", full.names = TRUE)

if (length(gene_group_files) == 0) {
  cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
} else {
  cat("Found", length(gene_group_files), "gene group files\n\n")

  for (gene_group_file in gene_group_files) {
    gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
    cat("  Processing:", gene_group_name, "\n")

    # Read gene list
    gene_list <- readLines(gene_group_file)
    gene_list <- gene_list[!grepl("^#", gene_list) & nzchar(gene_list)]
    gene_list <- trimws(gene_list)

    if (length(gene_list) == 0) {
      cat("    Skipping: No genes in list\n")
      next
    }

    # Match genes in data
    genes_in_data <- intersect(gene_list, gene_ids)

    if (length(genes_in_data) == 0) {
      cat("    Skipping: No matching genes found\n")
      next
    }

    cat("    Found", length(genes_in_data), "/", length(gene_list), "genes\n")

    # Create gene group output directory
    gene_group_dir <- file.path(output_dir, gene_group_name)
    dir.create(gene_group_dir, recursive = TRUE, showWarnings = FALSE)

    # Extract and save subsets for each count type
    for (count_type in count_types) {
      matrix_data <- matrices_list[[count_type]]
      subset_matrix <- matrix_data[genes_in_data, , drop = FALSE]

      # Save with SRR IDs
      subset_file_srr <- file.path(gene_group_dir,
        paste0(gene_group_name, "_", count_type, "_geneID_from_", MASTER_REFERENCE, ".tsv"))
      save_matrix_file(subset_matrix, subset_file_srr, label_type = "SRR")

      # Save with Organ labels
      subset_file_organ <- file.path(gene_group_dir,
        paste0(gene_group_name, "_", count_type, "_geneName_from_", MASTER_REFERENCE, ".tsv"))
      save_matrix_file(subset_matrix, subset_file_organ, label_type = "Organ")
    }

    cat("    Saved all count type matrices\n")
  }
}

# ===============================================
# STEP 6: GENERATE SUMMARY STATISTICS
# ===============================================

cat("\nStep 6: Generating summary statistics...\n")

summary_file <- file.path(output_dir, "rsem_quantification_summary.txt")

sink(summary_file)
cat("RSEM Quantification Summary - Method 5\n")
cat("========================================\n\n")

cat("Reference:", MASTER_REFERENCE, "\n")
cat("Total Genes:", length(gene_ids), "\n")
cat("Total Samples:", length(SAMPLE_IDS), "\n\n")

cat("Samples:\n")
for (sample_id in SAMPLE_IDS) {
  cat("  ", sample_id, "-", SAMPLE_LABELS[sample_id], "\n")
}
cat("\n")

cat("Count Types Generated:\n")
cat("  • expected_count - Raw probabilistic counts from RSEM EM algorithm\n")
cat("  • TPM - Transcripts Per Million (normalized, cross-sample comparable)\n")
cat("  • FPKM - Fragments Per Kilobase Million (normalized for length & depth)\n\n")

cat("TPM Summary Statistics:\n")
for (sample_id in SAMPLE_IDS) {
  tpm_values <- tpm_matrix[, sample_id]
  cat("  ", sample_id, ":\n")
  cat("    Mean TPM:", round(mean(tpm_values, na.rm = TRUE), 2), "\n")
  cat("    Median TPM:", round(median(tpm_values, na.rm = TRUE), 2), "\n")
  cat("    Genes with TPM > 1:", sum(tpm_values > 1, na.rm = TRUE), "\n")
}
cat("\n")

cat("Gene Groups Processed:", length(gene_group_files), "\n")
for (gene_group_file in gene_group_files) {
  gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
  gene_list <- readLines(gene_group_file)
  gene_list <- gene_list[!grepl("^#", gene_list) & nzchar(gene_list)]
  genes_in_data <- intersect(gene_list, gene_ids)
  cat("  ", gene_group_name, ":", length(genes_in_data), "/", length(gene_list), "genes found\n")
}

sink()

cat("Saved summary:", summary_file, "\n\n")

# ===============================================
# SUMMARY
# ===============================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("RSEM MATRIX GENERATION COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Summary:\n")
cat("  • Processed", length(SAMPLE_IDS), "samples\n")
cat("  • Extracted", length(gene_ids), "genes\n")
cat("  • Generated 3 count types: expected_count, TPM, FPKM\n")
cat("  • Processed", length(gene_group_files), "gene groups\n")
cat("  • Output directory:", output_dir, "\n\n")

cat("Normalization Approach:\n")
cat("  ✅ RSEM provides pre-normalized TPM and FPKM values\n")
cat("  ✅ TPM is normalized across samples (sum to 1 million)\n")
cat("  ✅ FPKM is normalized for gene length and sequencing depth\n")
cat("  ✅ NO DESeq2 needed - RSEM normalization is correct\n\n")

cat("Matrices ready for heatmap generation\n\n")

cat("Technical Note:\n")
cat("RSEM uses Expectation-Maximization algorithm to handle:\n")
cat("  • Multi-mapping reads\n")
cat("  • Isoform ambiguity\n")
cat("  • Probabilistic assignment\n")
cat("The TPM and FPKM values are already properly normalized.\n\n")
