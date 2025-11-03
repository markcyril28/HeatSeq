# ==============================================================================
# SALMON SAF TO DESEQ2 PREPARATION SCRIPT (METHOD 4)
# ==============================================================================
# Description: This R script processes Salmon quantification files (`quant.sf`)
#              generated from Salmon SAF Quantification (Method 4), aggregates
#              transcript-level counts to gene-level, and prepares them for
#              differential expression analysis using DESeq2.
# Author: Mark Cyril R. Mercado
# Version: v1
# Date: October 2025
#
# Input: Salmon `quant.sf` files for each sample.
# Output: A combined gene count matrix and a sample information (colData) table,
#         both in TSV format, suitable for direct import into DESeq2.
#
# Required packages: dplyr, tibble, readr, tximport, AnnotationDbi, org.Mm.eg.db (example)
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(tximport)
  # For gene ID mapping, you might need an annotation package like org.Mm.eg.db
  # library(AnnotationDbi)
  # library(org.Mm.eg.db) # Example for mouse. Change as needed.
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base directory for the project
BASE_DIR <- getwd()

# Directory where Salmon quantification results for Method 4 are stored
SALMON_QUANT_ROOT <- file.path(BASE_DIR, "4d_Method_4_Salmon_Saf_Quantification", "quant")

# Output directory for DESeq2 input files
DESEQ2_OUTPUT_DIR <- file.path(BASE_DIR, "4d_Method_4_Salmon_Saf_Quantification", "DESeq2_input")

# List of FASTA tags used in the analysis (from ALL_FASTA_FILES in bash script)
FASTA_TAGS <- c("TEST") # Example: c("All_Smel_Genes", "TEST")

# List of SRR IDs used in the analysis (from SRR_COMBINED_LIST in bash script)
SRR_IDS <- c("SRR3884686", "SRR3884687", "SRR3884597")

# Create output directory if it doesn't exist
dir.create(DESEQ2_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to create a transcript-to-gene mapping table
# This function needs to be customized based on your transcriptome FASTA header format.
# Example: If transcript IDs are like 'TRINITY_DN1234_c0_g1_i1', and gene IDs are 'TRINITY_DN1234_c0_g1'
# Or if your FASTA headers contain both transcript and gene IDs.
# For Salmon, the 'Name' column in quant.sf is the transcript ID.
# You need a way to map these transcript IDs to gene IDs.
# This often comes from the FASTA file used to build the Salmon index, or a separate annotation file.

# Placeholder for tx2gene creation. CUSTOMIZE THIS PART.
# This example assumes transcript IDs are in the format 'geneID|transcriptID'
# If your transcript IDs are just gene IDs, then tx2gene is not strictly needed for gene-level counts.
# If your transcript IDs are complex, you'll need to parse them or use an external mapping file.

# For the 'TEST.fasta' example, let's assume transcript IDs are directly gene IDs for simplicity
# If your actual data has transcript IDs that need mapping, you'll need to adjust this.
# A common scenario for de novo assemblies is that Trinity outputs transcript IDs like
# 'TRINITY_DNXXXX_cY_gZ_iA' where 'TRINITY_DNXXXX_cY_gZ' is the gene ID.

# For now, we'll create a dummy tx2gene that maps transcript IDs to themselves (gene-level counts)
# if the 'Name' column in quant.sf is already gene-level. If it's transcript-level, you need a proper mapping.

# Example of how to create tx2gene if transcript IDs are like 'geneID_transcriptID'
# tx2gene_df <- data.frame(
#   TXNAME = c("geneA_t1", "geneA_t2", "geneB_t1"),
#   GENEID = c("geneA", "geneA", "geneB"),
#   stringsAsFactors = FALSE
# )

# For the purpose of this script, we'll assume a simple case or that tximport can handle it
# if the 'Name' column in quant.sf is already gene-level or easily parsed.

# If you have a separate file mapping transcript IDs to gene IDs, load it here.
# For example: tx2gene_df <- read_tsv("path/to/tx2gene_map.tsv")

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SALMON SAF TO DESEQ2 PREPARATION (METHOD 4)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

for (fasta_tag in FASTA_TAGS) {
  cat(paste0("Processing FASTA tag: ", fasta_tag, "\n"))
  
  # Construct file paths for tximport
  files <- file.path(SALMON_QUANT_ROOT, fasta_tag, SRR_IDS, "quant.sf")
  names(files) <- SRR_IDS
  
  # Check if all quant.sf files exist
  if (!all(file.exists(files))) {
    missing_files <- files[!file.exists(files)]
    warning(paste("Missing Salmon quant.sf files:", paste(missing_files, collapse = ", ")))
    cat(paste0("  Skipping FASTA tag: ", fasta_tag, " due to missing files.\n"))
    next
  }
  
  # Placeholder for tx2gene. This needs to be properly defined based on your data.
  # For now, we'll try without tx2gene, which means tximport will try to guess
  # or you'll get transcript-level counts. For gene-level, a proper tx2gene is crucial.
  # If your 'Name' column in quant.sf is already gene IDs, you might not need tx2gene.
  
  # Example of a simple tx2gene if transcript IDs are like 'geneID|transcriptID'
  # This part needs to be adapted to your specific FASTA header format.
  # For 'TEST.fasta', if it contains gene sequences, then Salmon quantifies genes directly.
  # If it contains transcript sequences, you need a mapping.
  
  # For now, let's assume the 'Name' column in quant.sf contains gene IDs directly
  # or that we want transcript-level counts if a proper tx2gene is not provided.
  
  # If you have a tx2gene data frame, uncomment and use it:
  # txi <- tximport(files, type = "salmon", tx2gene = tx2gene_df)
  
  # Attempt tximport without tx2gene (will give transcript-level counts)
  # or if 'Name' column in quant.sf is already gene IDs, it will work for gene-level.
  cat("  Importing Salmon quantification data with tximport...\n")
  txi <- tryCatch({
    tximport(files, type = "salmon", txOut = FALSE, countsFromAbundance = "scaledTPM")
  }, error = function(e) {
    message("Error during tximport: ", e$message)
    message("Attempting with txOut=TRUE to get transcript-level data. You may need to provide a tx2gene mapping.")
    tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM")
  })
  
  # Extract gene-level counts (if txOut=FALSE was successful, or aggregate manually if txOut=TRUE)
  if (is.list(txi) && "counts" %in% names(txi)) {
    gene_counts <- as.data.frame(txi$counts)
    gene_counts <- round(gene_counts) # Ensure integer counts for DESeq2
    
    # Set gene IDs as row names
    gene_counts <- gene_counts %>%
      rownames_to_column(var = "gene_id") %>%
      column_to_rownames(var = "gene_id")
    
    # Save the combined gene count matrix
    output_matrix_file <- file.path(DESEQ2_OUTPUT_DIR, paste0("gene_count_matrix_salmon_", fasta_tag, ".tsv"))
    write_tsv(as_tibble(gene_counts, rownames = "gene_id"), output_matrix_file)
    cat(paste0("  Saved combined gene count matrix to: ", output_matrix_file, "\n"))
  } else {
    cat(paste0("  Could not extract gene counts for FASTA tag: ", fasta_tag, ". Check tximport output.\n"))
  }
}

# Create sample information (colData) table
cat("Creating sample information table...\n")

# This is a placeholder. Users will need to customize 'condition' based on their experimental design.
# For now, we'll create a dummy condition based on SRR ID order.

sample_data <- tibble(sample = SRR_IDS) %>%
  mutate(condition = factor(rep(c("cond1", "cond2"), length.out = length(SRR_IDS)))) # Example conditions

output_coldata_file <- file.path(DESEQ2_OUTPUT_DIR, "sample_info_salmon.tsv")
write_tsv(sample_data, output_coldata_file)
cat(paste0("Saved sample information table to: ", output_coldata_file, "\n"))

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SALMON SAF DESEQ2 PREPARATION COMPLETED.\n")
cat("Please review 'sample_info_salmon.tsv' and update 'condition' column based on your experimental design.\n")
cat("If gene-level counts are incorrect, ensure a proper tx2gene mapping is provided for tximport.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
