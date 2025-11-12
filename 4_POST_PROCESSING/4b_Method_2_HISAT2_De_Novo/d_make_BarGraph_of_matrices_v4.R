# ===============================================
# BAR GRAPH GENERATION FOR GENE EXPRESSION
# ===============================================
#
# This script generates bar graphs showing individual gene expression across developmental stages.
# Output is saved to "8_Bar_Graph_Visualizations" directory.
#
# Important: Normalization handling for StringTie metrics
# - Coverage: RAW read depth (NOT library-size normalized) - REQUIRES CPM normalization
# - FPKM/TPM: Already library-size and gene-length normalized - only needs log2
#
# Features:
# - Individual bar plots for each gene
# - Expression values across developmental stages (organs)
# - Color-coded by expression level using custom color scale
# - Professional layout matching reference design
# - Saves both PNG bar graphs and TSV data files
# - Multiple normalization schemes - Each gene generates TWO bar graphs:
#   * sorted_by_organ: Expression values in original developmental stage order (x-axis)
#   * sorted_by_expression: Expression values sorted by intensity (lowest to highest from left to right)
#   * Color-coded bars by expression level using custom color scale
#   * Professional layout matching reference design
#   * Y-axis labeled according to normalization scheme
# - Output structure: 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_organ/{gene_id}_{normalization_scheme}_bargraph.png
#                   : 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_expression/{gene_id}_{normalization_scheme}_bargraph.png
#
# Normalization Schemes:
# - raw: No transformation, original count values
#   * Use when: Inspecting raw count distributions or QC
# - raw_normalized: CPM normalization without log transformation
#   * Use when: Need to compare samples with different sequencing depths on linear scale
# - count_type_normalized: Count-type specific (Coverage: CPM+log2, FPKM/TPM: log2 only)
#   * Uses consistent pseudocount of 1 across all count types for standardization
#   * Note: Both coverage (post-CPM) and FPKM/TPM use log2(x+1) transformation
#   * Use when: RECOMMENDED DEFAULT for most visualizations and analyses
# - zscore: Count-type normalization + Z-score standardization per gene
#   * Requires >1 gene and >1 sample; warns and falls back to count-type norm otherwise
#   * Use when: Comparing expression patterns across different genes (removes gene-specific baseline)
# - zscore_scaled_to_ten: Z-score normalized then scaled to [0,10] range PER GENE
#   * Genes with constant expression (zero variance) will appear as 0
#   * Each gene independently scaled to use full [0,10] range
#   * Use when: Creating heatmaps or visualizations where all genes should use same color scale
# - cpm: CPM normalization + log2 (for coverage data only, warns if used with FPKM/TPM)
#   * Automatically redirects to count-type normalization for FPKM/TPM inputs
#   * Use when: Working exclusively with raw coverage counts and need library-size normalization
#
# Required packages: ggplot2, dplyr, RColorBrewer, grid, scales
# ===============================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(grid)
  library(scales)
  library(gridExtra)
})

# ===============================================
# CONFIGURATION
# ===============================================

# Input and output directories
BASE_DIR <- getwd()
MATRICES_DIR <- file.path("5_stringtie_WD", "b_Method_2_COUNT_MATRICES")
BARGRAPH_OUT_DIR <- "8_Bar_Graph_Visualizations"

# File naming configuration
MASTER_REFERENCE <- "All_Smel_Genes"
MASTER_REFERENCE_SUFFIX <- paste0("_from_", MASTER_REFERENCE)

# Gene groups (Boilerplates)
FASTA_GROUPS <- c(
  # Control Gene Groups
  "Best_Cell_Cycle_Associated_Control_Genes",
  "Best_Control_Genes",
  
  # Individual Gene Groups
  "SmelDMPs",
  #"SmelDMPs_with_18s_rRNA",
  #"SmelDMPs_with_18s_rRNA_PSBMB",
  "SmelGIFs",
  "SmelGRFs",
  #"Selected_GRF_GIF_Genes",
  "Selected_GRF_GIF_Genes_vAll_GIFs",
  #"Selected_GRF_GIF_Genes_vTwo_GIFs",
  
  # Combined Gene Groups with Control Genes
  #"SmelGIF_with_Cell_Cycle_Control_genes",
  #"SmelGRF_with_Cell_Cycle_Control_genes",
  #"SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Analysis configuration
COUNT_TYPES <- c(
  #"coverage", 
  "fpkm"
  #"tpm"
)

GENE_TYPES <- c(
  #"geneID", 
  "geneName"
)

LABEL_TYPES <- c("SRR", "Organ")

# Normalization schemes configuration
NORMALIZATION_SCHEMES <- c(
  #"raw", 
  #"raw_normalized", 
  #"count_type_normalized", 
  #"zscore", 
  "zscore_scaled_to_ten"
  #"cpm"
)

# Sample label mapping: SRR IDs to developmental stage/organ names
SAMPLE_LABELS <- c(
    # Roots
    "SRR3884675" = "Roots",      # PRJNA328564
    "SRR20722229" = "Roots_2",     # SAMN28540077
    "SRR31755282" = "Roots_3",     # SAMN28540068

    # Stems
    "SRR3884690" = "Stems",      # PRJNA328564
    "SRR20722227" = "Stems_2",     # SAMN28540077
    "SRR20722384" = "Stems_3",     # SAMN28540068

    # Leaves
    "SRR3884689" = "Leaves",     # PRJNA328564
    "SRR20722230" = "Leaves_2",    # SAMN28540077
    "SRR20722386" = "Leaves_3",    # SAMN28540068
    "SRR3884684" = "Senescent_leaves", # PRJNA328564

    # Buds
    "SRR3884686" = "Buds",       # PRJNA328564
    "SRR21010466" = "Buds_2",      # SAMN28540077
    "SRR20722297" = "Buds_3",      # SAMN28540068

    # Opened Buds
    "SRR3884687" = "Opened_Buds", # PRJNA328564

    # Flowers
    "SRR3884597" = "Flowers",    # PRJNA328564
    "SRR20722234" = "Flowers_2",   # SAMN28540077
    "SRR23909863" = "Flowers_3",   # SAMN28540068

    # Fruits
    "SRR3884631" = "Fruits",     # PRJNA328564
    "SRR2072232" = "Fruits_2",     # SAMN28540077
    "SRR20722387" = "Fruits_3",    # SAMN28540068
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


# ===============================================
# OUTPUT DIRECTORY SETUP
# ===============================================

# Create base output directory if it doesn't exist
dir.create(BARGRAPH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# FUNCTIONS
# ===============================================

# Load utility functions
source("9_utility_functions.R")

# Function to calculate coefficient of variation (CV)
calculate_cv <- function(data_matrix) {
  cv_values <- apply(data_matrix, 1, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    if (mean_val == 0 || is.na(mean_val) || is.na(sd_val)) {
      return(0)
    }
    return(sd_val / mean_val)
  })
  return(cv_values)
}

# ===============================================
# MAIN PROCESSING
# ===============================================

# Close any open graphics devices
if (length(dev.list()) > 0) { 
    sapply(dev.list(), dev.off) 
}

# Initialize counters
total_bargraphs <- 0
successful_bargraphs <- 0

# Validate FASTA_GROUPS
if (length(FASTA_GROUPS) == 0) {
    stop("No FASTA groups defined for processing")
}

# Process each FASTA group first
for (group in FASTA_GROUPS) {
  # cat("Processing FASTA group:", group, "\n")
  
  for (normalization_scheme in NORMALIZATION_SCHEMES) {
    # Create normalization-specific output directory within FASTA group
    norm_output_dir <- file.path(BARGRAPH_OUT_DIR, group, normalization_scheme)
    dir.create(norm_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create subdirectories for different sorting methods
    sorted_by_organ_dir <- file.path(norm_output_dir, "sorted_by_organ")
    sorted_by_expression_dir <- file.path(norm_output_dir, "sorted_by_expression")
    dir.create(sorted_by_organ_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sorted_by_expression_dir, recursive = TRUE, showWarnings = FALSE)
    
    # cat("  ✓ Created directories:", norm_output_dir, "\n")
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        # Always use SRR label_type for input file (as Organ is not a file label)
        label_type <- "SRR"
        input_file <- file.path(MATRICES_DIR, group, 
                               paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, MASTER_REFERENCE_SUFFIX, ".tsv"))
        # If that doesn't exist, try without the suffix
        if (!file.exists(input_file)) {
          input_file <- file.path(MATRICES_DIR, group, 
                                 paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
        }
        if (!file.exists(input_file)) next
        input_basename <- tools::file_path_sans_ext(basename(input_file))
        title <- gsub("_", " ", input_basename)
        # Read raw count matrix
        raw_data <- read_count_matrix(input_file)
        if (is.null(raw_data)) {
          next
        }
        # Apply the specified normalization scheme
        processed_data <- apply_normalization(raw_data, normalization_scheme, count_type)
        if (!is.null(processed_data)) {
          # Generate bar graphs sorted by organ (developmental stage order)
          genes_plotted_organ_srr <- generate_gene_bargraphs(processed_data, sorted_by_organ_dir, title, count_type, "SRR", normalization_scheme, sort_by_expression = FALSE)
          genes_plotted_organ_label <- generate_gene_bargraphs(processed_data, sorted_by_organ_dir, title, count_type, "Organ", normalization_scheme, sort_by_expression = FALSE)
          
          # Generate bar graphs sorted by expression level
          genes_plotted_expression_srr <- generate_gene_bargraphs(processed_data, sorted_by_expression_dir, title, count_type, "SRR", normalization_scheme, sort_by_expression = TRUE)
          genes_plotted_expression_label <- generate_gene_bargraphs(processed_data, sorted_by_expression_dir, title, count_type, "Organ", normalization_scheme, sort_by_expression = TRUE)
          
          total_bargraphs <- total_bargraphs + (nrow(processed_data) * 4)  # 2 label types x 2 sorting methods
          successful_bargraphs <- successful_bargraphs + genes_plotted_organ_srr + genes_plotted_expression_srr + genes_plotted_organ_label + genes_plotted_expression_label
        }
      }
    }
  }
}

# ===============================================
# SUMMARY
# ===============================================

while (length(dev.list()) > 0) { dev.off() }

cat("===============================================\n")
cat("BAR GRAPH GENERATION SUMMARY\n")
cat("===============================================\n")
cat("Total bar graphs generated:", successful_bargraphs, "\n")

if (successful_bargraphs > 0) {
  cat("✓ Bar graph generation completed successfully!\n")
  cat("\nOutput folder structure:\n")
  cat(paste0(BARGRAPH_OUT_DIR, "/\n"))
  for (group in FASTA_GROUPS) {
    cat(paste0("├── ", group, "/\n"))
    for (norm_scheme in NORMALIZATION_SCHEMES) {
      group_dir <- file.path(BARGRAPH_OUT_DIR, group, norm_scheme)
      if (dir.exists(group_dir)) {
        organ_dir <- file.path(group_dir, "sorted_by_organ")
        expr_dir <- file.path(group_dir, "sorted_by_expression")
        organ_count <- ifelse(dir.exists(organ_dir), length(list.files(organ_dir, pattern = "*.png")), 0)
        expr_count <- ifelse(dir.exists(expr_dir), length(list.files(expr_dir, pattern = "*.png")), 0)
        cat(paste0("│   ├── ", norm_scheme, "/\n"))
        cat(paste0("│   │   ├── sorted_by_organ/  (", organ_count, " graphs)\n"))
        cat(paste0("│   │   └── sorted_by_expression/  (", expr_count, " graphs)\n"))
      }
    }
  }
} else {
  cat("✗ No bar graphs were generated. Please check input files and paths.\n")
}

# UPDATED STRUCTURE WITH FASTA GROUPS AS PRIMARY ORGANIZATION:
# =============================================================
# - Each combination of (group, normalization_scheme, count_type_normalized, gene_type, label_type) 
#   corresponds to exactly one TSV file input and generates individual bar graphs for each gene
# - Six normalization schemes are applied (UPDATED - Nov 2025):
#   * "raw": Original count data without transformation
#   * "raw_normalized": CPM normalization WITHOUT log transformation (fixed from being identical to raw)
#   * "count_type_normalized": Count-type specific normalization (Coverage: CPM+log2, FPKM/TPM: direct log2)
#   * "zscore": Count-type normalization followed by Z-score standardization
#   * "zscore_scaled_to_ten": Count-type norm + Z-score + PER-GENE min-max scaling to [0,10] (fixed from global scaling)
#   * "cpm": Counts Per Million normalization + log2 (with safeguard against FPKM/TPM misuse)
# - Each gene generates bar graphs showing:
#   * Expression values across developmental stages (x-axis)
#   * Two versions: sorted_by_organ (developmental order) and sorted_by_expression (intensity order)
#   * Color-coded bars by expression level using custom color scale
#   * Professional layout matching reference design
#   * Y-axis labeled according to normalization scheme
# - Output structure: 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_organ/{gene_id}_{normalization_scheme}_bargraph.png
#                   : 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_expression/{gene_id}_{normalization_scheme}_bargraph.png
# - Bar colors represent expression intensity: white (lowest) to violet (highest expression)
# - Individual files for each gene allow easy comparison across conditions
# - Separate folders for different sorting methods facilitate analysis workflows
#
# NORMALIZATION IMPROVEMENTS (Nov 2025):
# ======================================
# 1. raw_normalized: Now performs actual CPM normalization (not log-transformed) instead of being identical to raw
# 2. zscore_scaled_to_ten: Uses GLOBAL min-max scaling across all genes and samples
#    - Entire matrix scaled to [0,10] range for consistent color mapping
#    - Ensures comparability across all genes in the dataset
# 3. cpm: Added validation to prevent inappropriate application to FPKM/TPM data
#    - Warns user and redirects to count_type_normalized when applied to already-normalized data
# 4. Updated y-axis labels to accurately reflect the normalization applied
# 5. Removed duplicate function definitions that were causing the legacy versions to override active ones
# 6. Added edge case handling and warnings for single-sample/single-gene scenarios in Z-score normalizations
# 7. Standardized pseudocount to 1 across all count types for consistency
#    - Previous: coverage used 1, FPKM/TPM used 0.1 (inconsistent)
#    - Updated: all count types now use 1 (standard in RNA-seq analysis, consistent across pipelines)
#    - Rationale: Provides consistency, aligns with common practice, simplifies interpretation
# 8. Enhanced documentation with clear "Use when:" guidelines for each normalization scheme
# 9. Verified mathematical soundness: all schemes follow bioinformatics best practices
#     - Z-score applied to log-transformed data (standard approach for RNA-seq)
#     - Global scaling ensures consistent visualization scale across all genes
#     - Row/column names automatically preserved through all matrix operations
