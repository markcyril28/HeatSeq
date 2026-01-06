#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: RSEM TO DESEQ2 MATRICES - GENERIC
# ===============================================
# Processes RSEM quantification using tximport
# for statistically sound gene expression analysis

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# Source shared utilities (provides ensure_output_dir, convert_to_organ_labels, etc.)
SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))

# ===============================================
# RSEM-SPECIFIC HELPER
# ===============================================

save_count_matrix_rsem <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix, sample_labels) {
  # Save with SRR IDs
  output_file_srr <- file.path(output_dir, 
    paste0(base_name, "_expected_count_geneID_from_", master_ref, level_suffix, ".tsv"))
  
  counts_df <- as.data.frame(counts_matrix)
  counts_df <- tibble::rownames_to_column(counts_df, "GeneID")
  write.table(counts_df, output_file_srr, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", output_file_srr, "\n")
  
  # Save with Organ labels using shared function
  output_file_organ <- file.path(output_dir,
    paste0(base_name, "_expected_count_geneName_from_", master_ref, level_suffix, ".tsv"))
  
  counts_organ <- convert_to_organ_labels(counts_matrix)
  counts_organ_df <- as.data.frame(counts_organ)
  counts_organ_df <- tibble::rownames_to_column(counts_organ_df, "GeneID")
  write.table(counts_organ_df, output_file_organ, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", output_file_organ, "\n")
}

# ===============================================
# CONFIGURATION
# ===============================================

QUANT_DIR <- "5_RSEM_Quant_WD"
MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", "All_Smel_Genes")
MATRICES_OUTPUT_DIR <- "6_matrices"
# Use shared GENE_GROUPS_DIR from 0_shared_config.R (already sourced)

# Toggle to generate both gene-level and isoform-level matrices
GENERATE_GENE_LEVEL <- TRUE      # Use .genes.results (aggregated)
GENERATE_ISOFORM_LEVEL <- TRUE   # Use .isoforms.results (transcript-specific)

# Use SAMPLE_IDS and SAMPLE_LABELS from shared config (0_shared_config.R)
# Already sourced above - no duplication needed

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
ensure_output_dir(output_dir)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("TXIMPORT: RSEM QUANTIFICATION TO MATRICES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Configuration:\n")
cat("  • Master Reference:", MASTER_REFERENCE, "\n")
cat("  • Generate gene-level matrices:", GENERATE_GENE_LEVEL, "\n")
cat("  • Generate isoform-level matrices:", GENERATE_ISOFORM_LEVEL, "\n\n")

# Define processing levels
processing_levels <- list()
if (GENERATE_GENE_LEVEL) {
  processing_levels[["gene_level"]] <- list(
    file_type = ".genes.results",
    tx_in = FALSE,
    tx_out = FALSE,
    label = "Gene-Level",
    output_suffix = "_gene_level"
  )
}
if (GENERATE_ISOFORM_LEVEL) {
  processing_levels[["isoform_level"]] <- list(
    file_type = ".isoforms.results",
    tx_in = TRUE,
    tx_out = TRUE,
    label = "Isoform-Level",
    output_suffix = "_isoform_level"
  )
}

if (length(processing_levels) == 0) {
  stop("ERROR: At least one of GENERATE_GENE_LEVEL or GENERATE_ISOFORM_LEVEL must be TRUE")
}

# ===============================================
# PROCESS EACH LEVEL
# ===============================================

for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PROCESSING:", level_config$label, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # ===============================================
  # STEP 1: LOCATE RSEM OUTPUT FILES
  # ===============================================
  
  cat("Step 1: Locating RSEM", level_config$label, "output files...\n")
  
  rsem_quant_dir <- file.path(QUANT_DIR, MASTER_REFERENCE)
  
  # Build paths to RSEM files
  files <- file.path(rsem_quant_dir, SAMPLE_IDS, paste0(SAMPLE_IDS, level_config$file_type))
  names(files) <- SAMPLE_IDS
  
  # Check which files exist
  files_exist <- file.exists(files)
  if (sum(files_exist) == 0) {
    cat("ERROR: No", level_config$label, "quantification files found in", rsem_quant_dir, "\n")
    cat("Skipping this level...\n\n")
    next
  }
  
  current_sample_ids <- SAMPLE_IDS
  if (sum(files_exist) < length(files)) {
    missing <- SAMPLE_IDS[!files_exist]
    cat("WARNING: Missing quantification for samples:", paste(missing, collapse = ", "), "\n")
    files <- files[files_exist]
    current_sample_ids <- SAMPLE_IDS[files_exist]
  }
  
  cat("Found", length(files), "RSEM quantification files\n")
  cat("Samples:", paste(current_sample_ids, collapse = ", "), "\n\n")
  
  # ===============================================
  # STEP 2: IMPORT WITH TXIMPORT
  # ===============================================
  
  cat("Step 2: Importing RSEM data with tximport...\n")
  
  # tximport for RSEM
  txi <- tximport(files, type = "rsem", txIn = level_config$tx_in, txOut = level_config$tx_out)
  
  entity_type <- if (level_config$tx_out) "transcripts" else "genes"
  cat("Successfully imported data for", ncol(txi$counts), "samples\n")
  cat("Total", entity_type, ":", nrow(txi$counts), "\n\n")

  # ===============================================
  # STEP 3: CREATE SAMPLE METADATA
  # ===============================================
  
  cat("Step 3: Creating sample metadata...\n")
  
  # Get conditions from SAMPLE_LABELS, fallback to sample ID if not found
  conditions <- SAMPLE_LABELS[current_sample_ids]
  missing_labels <- is.na(conditions)
  if (any(missing_labels)) {
    cat("WARNING: No labels found for samples:", paste(current_sample_ids[missing_labels], collapse = ", "), "\n")
    cat("Using sample IDs as condition labels for these samples\n")
    conditions[missing_labels] <- current_sample_ids[missing_labels]
  }
  
  sample_data <- data.frame(
    SampleID = current_sample_ids,
    Condition = conditions,
    row.names = current_sample_ids,
    stringsAsFactors = FALSE
  )
  
  cat("Sample metadata:\n")
  print(sample_data)
  cat("\n")
  
  # ===============================================
  # STEP 4: CHECK FOR REPLICATES
  # ===============================================
  
  cat("Step 4: Checking for biological replicates...\n")
  
  condition_counts <- table(sample_data$Condition)
  has_replicates <- any(condition_counts > 1)
  
  cat("Samples per condition:\n")
  print(condition_counts)
  cat("\n")
  
  if (!has_replicates) {
    cat("WARNING: No biological replicates detected!\n")
    cat("Using TPM normalization for visualization purposes\n\n")
  }
  
  # ===============================================
  # STEP 5: NORMALIZATION
  # ===============================================
  
  if (has_replicates) {
    cat("Step 5: Running DESeq2 normalization...\n")
    
    dds <- DESeqDataSetFromTximport(
      txi = txi,
      colData = sample_data,
      design = ~ Condition
    )
    
    dds <- DESeq(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    cat("DESeq2 normalization complete\n")
    
  } else {
    cat("Step 5: Computing TPM values...\n")
    normalized_counts <- txi$abundance
    cat("TPM normalization complete (visualization only)\n")
  }
  
  rownames(normalized_counts) <- rownames(normalized_counts)
  cat("Matrix dimensions:", nrow(normalized_counts), entity_type, "x", ncol(normalized_counts), "samples\n\n")
  
  # ===============================================
  # STEP 6: SAVE FULL NORMALIZED MATRIX
  # ===============================================
  
  cat("Step 6: Saving normalized count matrix...\n")
  
  level_output_dir <- file.path(output_dir, level_name)
  ensure_output_dir(level_output_dir)
  
  save_count_matrix(normalized_counts, level_output_dir, MASTER_REFERENCE, 
                   MASTER_REFERENCE, level_config$output_suffix, SAMPLE_LABELS)
  cat("\n")
  
  # ===============================================
  # STEP 7: PROCESS GENE GROUPS
  # ===============================================
  
  cat("Step 7: Processing gene groups...\n")
  
  gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)
  
  if (length(gene_group_files) == 0) {
    cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group files\n\n")
  
    for (gene_group_file in gene_group_files) {
      gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
      cat("  Processing:", gene_group_name, "\n")
  
      # Read gene list from CSV (first column is Gene_ID)
      tryCatch({
        if (grepl("\\.csv$", gene_group_file, ignore.case = TRUE)) {
          gene_df <- read.csv(gene_group_file, stringsAsFactors = FALSE, header = TRUE)
          if ("Gene_ID" %in% colnames(gene_df)) {
            gene_list <- gene_df$Gene_ID
          } else {
            gene_list <- gene_df[[1]]  # Fallback to first column
          }
        } else {
          gene_list <- suppressWarnings(readLines(gene_group_file))
          gene_list <- gene_list[!grepl("^#|^Gene_ID", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
        }
        gene_list <- trimws(gene_list)
      }, error = function(e) {
        cat("    Error reading file:", e$message, "\n")
        gene_list <- character(0)
      })
  
      if (length(gene_list) == 0) {
        cat("    Skipping: No genes in list\n")
        next
      }
  
      genes_in_data <- intersect(gene_list, rownames(normalized_counts))
  
      if (length(genes_in_data) == 0) {
        cat("    Skipping: No matching", entity_type, "found\n")
        next
      }
  
      cat("    Found", length(genes_in_data), "/", length(gene_list), entity_type, "\n")
  
      gene_group_dir <- file.path(level_output_dir, gene_group_name)
      ensure_output_dir(gene_group_dir)
  
      subset_counts <- normalized_counts[genes_in_data, , drop = FALSE]
      save_count_matrix(subset_counts, gene_group_dir, gene_group_name,
                       MASTER_REFERENCE, level_config$output_suffix, SAMPLE_LABELS)
      cat("    Saved subset matrices\n")
    }
  }
  
  # ===============================================
  # LEVEL SUMMARY
  # ===============================================
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(level_config$label, "PROCESSING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

# ===============================================
# FINAL SUMMARY
# ===============================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL PROCESSING COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Generated matrices for", length(processing_levels), "level(s)\n")
cat("Output directory:", output_dir, "\n")
cat("Matrices ready for heatmap generation\n\n")
