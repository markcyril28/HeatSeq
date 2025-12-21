#!/usr/bin/env Rscript

# ===============================================
# TXIMPORT: RSEM TO DESEQ2 MATRICES - METHOD 5
# ===============================================
# Processes RSEM quantification using tximport
# for statistically sound gene expression analysis

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# ===============================================
# LOAD SHARED CONFIGURATION
# ===============================================

source("B_M5_modules/0_shared_config.R")

# ===============================================
# HELPER FUNCTIONS
# ===============================================

# Convert column names to organ labels
convert_to_organ_labels <- function(counts_matrix) {
  colnames_organ <- SAMPLE_LABELS[colnames(counts_matrix)]
  colnames_organ[is.na(colnames_organ)] <- colnames(counts_matrix)[is.na(colnames_organ)]
  result <- counts_matrix
  colnames(result) <- colnames_organ
  result
}

# Save matrix with both SRR and Organ labels
save_count_matrix <- function(counts_matrix, output_dir, base_name, master_ref, level_suffix) {
  # Save with SRR IDs
  output_file_srr <- file.path(output_dir, 
    paste0(base_name, "_expected_count_geneID_from_", master_ref, level_suffix, ".tsv"))
  
  counts_df <- as.data.frame(counts_matrix)
  counts_df <- tibble::rownames_to_column(counts_df, "GeneID")
  write.table(counts_df, output_file_srr, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved:", output_file_srr, "\n")
  
  # Save with Organ labels
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
MATRICES_OUTPUT_DIR <- "6_matrices"

# Gene groups directory - points to subdirectory based on MASTER_REFERENCE
GENE_GROUPS_DIR <- paste0("A_GeneGroups_InputList/for_", MASTER_REFERENCE)

# Toggle to generate both gene-level and isoform-level matrices
GENERATE_GENE_LEVEL <- TRUE      # Use .genes.results (aggregated)
GENERATE_ISOFORM_LEVEL <- TRUE   # Use .isoforms.results (transcript-specific)

# Read SRR filter from environment variable (set by wrapper script)
srr_filter_env <- Sys.getenv("SRR_FILTER_LIST", unset = "")
if (nzchar(srr_filter_env)) {
  requested_samples <- trimws(strsplit(srr_filter_env, " ")[[1]])
  # Filter to only include requested samples that exist in SAMPLE_IDS
  SAMPLE_IDS <- SAMPLE_IDS[SAMPLE_IDS %in% requested_samples]
  cat("SRR Filter applied:", length(SAMPLE_IDS), "samples selected\n")
  cat("Selected samples:", paste(SAMPLE_IDS, collapse = ", "), "\n\n")
} else {
  cat("Processing all available samples (no filter)\n\n")
}

# Create output directories
output_dir <- file.path(MATRICES_OUTPUT_DIR, MASTER_REFERENCE)
ensure_output_dir(output_dir)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("TXIMPORT: RSEM QUANTIFICATION TO MATRICES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Configuration:\n")
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
  
  # Create sample metadata data frame
  sample_data <- data.frame(
    SampleID = current_sample_ids,
    Condition = SAMPLE_LABELS[current_sample_ids],
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
  
  # Count samples per condition
  condition_counts <- table(sample_data$Condition)
  has_replicates <- any(condition_counts > 1)
  
  cat("Samples per condition:\n")
  print(condition_counts)
  cat("\n")
  
  if (!has_replicates) {
    cat("WARNING: No biological replicates detected!\n")
    cat("Cannot use DESeq2 normalization without replicates.\n")
    cat("Falling back to TPM normalization for visualization purposes.\n\n")
  }
  
  # ===============================================
  # STEP 5: NORMALIZATION
  # ===============================================
  
  if (has_replicates) {
    cat("Step 5: Running DESeq2 normalization (replicates available)...\n")
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromTximport(
      txi = txi,
      colData = sample_data,
      design = ~ Condition
    )
    
    # Run DESeq2 for normalization
    dds <- DESeq(dds)
    
    # Extract normalized counts
    normalized_counts <- counts(dds, normalized = TRUE)
    
    cat("DESeq2 normalization complete\n")
    
  } else {
    cat("Step 5: Computing TPM values (no replicates available)...\n")
    
    # Calculate TPM from tximport abundance
    # tximport already provides abundance in TPM units for RSEM
    normalized_counts <- txi$abundance
    
    cat("TPM normalization complete\n")
    cat("Note: Using TPM (Transcripts Per Million) for visualization\n")
    cat("      This is appropriate for heatmaps and expression profiles\n")
    cat("      but NOT suitable for differential expression analysis\n")
  }
  
  # Keep row names as-is (preserve transcript IDs with suffixes for isoform level)
  rownames(normalized_counts) <- rownames(normalized_counts)
  
  cat("Normalized count matrix dimensions:", nrow(normalized_counts), entity_type, "x", ncol(normalized_counts), "samples\n\n")
  
  # ===============================================
  # STEP 6: SAVE FULL NORMALIZED MATRIX
  # ===============================================
  
  cat("Step 6: Saving normalized count matrix...\n")
  
  # Create level-specific output directory
  level_output_dir <- file.path(output_dir, level_name)
  ensure_output_dir(level_output_dir)
  
  # Save matrices using helper function
  save_count_matrix(normalized_counts, level_output_dir, MASTER_REFERENCE, 
                   MASTER_REFERENCE, level_config$output_suffix)
  cat("\n")
  
  # ===============================================
  # STEP 7: PROCESS GENE GROUPS
  # ===============================================
  
  cat("Step 7: Processing gene groups...\n")
  
  # Find all gene group files (.txt or .tsv)
  gene_group_files <- list.files(GENE_GROUPS_DIR, pattern = "\\.(txt|tsv)$", full.names = TRUE)
  
  if (length(gene_group_files) == 0) {
    cat("No gene group files found in", GENE_GROUPS_DIR, "\n")
  } else {
    cat("Found", length(gene_group_files), "gene group files\n\n")
  
    for (gene_group_file in gene_group_files) {
      gene_group_name <- tools::file_path_sans_ext(basename(gene_group_file))
      cat("  Processing:", gene_group_name, "\n")
  
      # Read gene list (suppress incomplete final line warning)
      gene_list <- suppressWarnings(readLines(gene_group_file))
      gene_list <- gene_list[!grepl("^#", gene_list) & nzchar(gene_list)]
      gene_list <- trimws(gene_list)
  
      if (length(gene_list) == 0) {
        cat("    Skipping: No genes in list\n")
        next
      }
  
      # Match genes in normalized counts
      genes_in_data <- intersect(gene_list, rownames(normalized_counts))
  
      if (length(genes_in_data) == 0) {
        cat("    Skipping: No matching", entity_type, "found\n")
        next
      }
  
      cat("    Found", length(genes_in_data), "/", length(gene_list), entity_type, "\n")
  
      # Create gene group output directory
      gene_group_dir <- file.path(level_output_dir, gene_group_name)
      ensure_output_dir(gene_group_dir)
  
      # Extract subset and save
      subset_counts <- normalized_counts[genes_in_data, , drop = FALSE]
      save_count_matrix(subset_counts, gene_group_dir, gene_group_name,
                       MASTER_REFERENCE, level_config$output_suffix)
      cat("    Saved subset matrices\n")
    }
  }
  
  # ===============================================
  # STEP 8: SAVE DESEQ2 OBJECT FOR FURTHER ANALYSIS (IF AVAILABLE)
  # ===============================================
  
  if (has_replicates && exists("dds")) {
    cat("\nStep 8: Saving DESeq2 object for downstream analysis...\n")
    
    deseq2_output_dir <- file.path(level_output_dir, "deseq2_objects")
    ensure_output_dir(deseq2_output_dir)
    
    dds_file <- file.path(deseq2_output_dir, paste0("dds_tximport", level_config$output_suffix, ".rds"))
    saveRDS(dds, file = dds_file)
    
    cat("Saved DESeq2 object:", dds_file, "\n")
    cat("This object can be used for differential expression analysis\n\n")
  } else {
    cat("\nStep 8: Skipping DESeq2 object save (no replicates available)\n")
    cat("Note: Differential expression analysis requires biological replicates\n\n")
  }
  
  # ===============================================
  # LEVEL SUMMARY
  # ===============================================
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(level_config$label, "PROCESSING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Summary:\n")
  cat("  • Imported", ncol(txi$counts), "samples with tximport\n")
  cat("  • Processed", nrow(normalized_counts), entity_type, "\n")
  
  if (has_replicates) {
    cat("  • Applied DESeq2 normalization\n")
  } else {
    cat("  • Applied TPM normalization (no replicates available)\n")
  }
  
  cat("  • Saved normalized count matrices\n")
  cat("  • Processed", length(gene_group_files), "gene groups\n")
  cat("  • Output directory:", level_output_dir, "\n\n")

} # End of processing level loop

# ===============================================
# FINAL SUMMARY
# ===============================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ALL PROCESSING COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Generated matrices for", length(processing_levels), "level(s):\n")
for (level_name in names(processing_levels)) {
  level_config <- processing_levels[[level_name]]
  cat("  •", level_config$label, "-", level_config$file_type, "\n")
}

cat("\nOutput directory:", output_dir, "\n")
cat("Matrices ready for heatmap generation\n\n")

cat("Important Notes:\n")
cat("  • Gene-level: Better for overall gene expression patterns\n")
cat("  • Isoform-level: Required when gene lists use transcript IDs (e.g., SmelGRF_08.140)\n")
cat("  • TPM values are suitable for visualization but NOT statistical analysis\n")
cat("  • For differential expression, biological replicates are required\n\n")
