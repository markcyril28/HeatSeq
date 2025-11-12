# ===============================================
# CV HEATMAP GENERATION WITH COEFFICIENT OF VARIATION
# ===============================================
#
# This script generates ONLY heatmaps with coefficient of variation (CV) annotation.
# Output is saved to "7_Heatmap_with_Covariability_Visualizations" directory.
#
# IMPORTANT: Normalization handling for StringTie metrics
# - Coverage: RAW read depth (NOT library-size normalized) - REQUIRES CPM normalization
# - FPKM/TPM: Already library-size and gene-length normalized - only needs log2
# - CV must ALWAYS be calculated on RAW, untransformed data
#
# Features:
# - Z-score normalized expression heatmaps
# - CV values displayed on the right side of each gene
# - Hierarchical clustering of genes
# - Purple-to-orange color gradient for Z-scores
# - Professional layout matching reference design
# - Saves both PNG heatmaps and TSV CV data files
# - Multiple sorting options:
#   * sorted_by_organ: Samples in original developmental stage order (columns)
#   * sorted_by_expression: Samples sorted by mean expression (lowest to highest from left to right)
#
# CV Calculation: CV = standard deviation / mean for each gene across samples
# High CV = high variability across developmental stages
# Low CV = consistent expression patterns
#
# Output structure: 7_Heatmap_with_Covariability_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_organ/{file}_cv_heatmap.png
#                 : 7_Heatmap_with_Covariability_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_expression/{file}_cv_heatmap.png
#
# Required packages: ComplexHeatmap, circlize, RColorBrewer, dplyr, tibble, grid
# Run test_heatmap_libraries.R first to install missing packages
# ===============================================

# Load required libraries
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(grid)
})

# ===============================================
# CONFIGURATION
# ===============================================

# Legend position configuration
# Options: "right" (vertical legend) or "bottom" (horizontal legend)
LEGEND_POSITION <- "bottom"

# Input and output directories
BASE_DIR <- getwd()
MATRICES_DIR <- file.path("5_stringtie_WD", "b_Method_2_COUNT_MATRICES")
HEATMAP_OUT_DIR <- "7_Heatmap_with_Covariability_Visualizations"

# File naming configuration
MASTER_REFERENCE <- "All_Smel_Genes"
MASTER_REFERENCE_SUFFIX <- paste0("_from_", MASTER_REFERENCE)

# Gene groups
FASTA_GROUPS <- c(
  # Control Gene Groups
  "Best_Cell_Cycle_Associated_Control_Genes",
  "Best_Control_Genes",
  
  # Individual Gene Groups
  "SmelDMPs",
  "SmelDMPs_with_18s_rRNA",
  "SmelDMPs_with_18s_rRNA_PSBMB",
  "SmelGIFs",
  "SmelGRFs",
  "Selected_GRF_GIF_Genes",
  "Selected_GRF_GIF_Genes_vAll_GIFs",
  "Selected_GRF_GIF_Genes_vTwo_GIFs",

  # Combined Gene Groups with Control Genes
  "SmelGIF_with_Cell_Cycle_Control_genes",
  "SmelGRF_with_Cell_Cycle_Control_genes",
  "SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Analysis configuration
COUNT_TYPES <- c(
  #"coverage", 
  "fpkm"
  #"tpm"
)

GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c(
  "SRR", 
  "Organ"
)

# Normalization schemes for gene expression analysis
# "raw" - No transformation (use for initial exploration)
# "count_type_normalized" - Log2-normalized (coverage→log2(CPM+1), FPKM/TPM→log2(value+1))
# "zscore" - Z-score standardization (count_type_normalized → z-score per gene)
# "zscore_scaled_to_ten" - Z-score scaled to [0, 10] range for visualization
# "cpm" - Counts Per Million (coverage only, returns log2(CPM+1))
NORMALIZATION_SCHEMES <- c(
  #"raw", 
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


# Create base output directory if it doesn't exist
dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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
      return(NA_real_)
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
total_heatmaps <- 0
successful_heatmaps <- 0

# Validate FASTA_GROUPS
if (length(FASTA_GROUPS) == 0) {
    stop("No FASTA groups defined for processing")
}

# Process each FASTA group first
for (group in FASTA_GROUPS) {
  # cat("\nProcessing FASTA group:", group, "\n")
  
  # Clear and recreate output directory for this specific group
  group_output_dir <- file.path(HEATMAP_OUT_DIR, group)
  if (dir.exists(group_output_dir)) {
    unlink(group_output_dir, recursive = TRUE)
  }
  dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (normalization_scheme in NORMALIZATION_SCHEMES) {
      # cat("  Processing normalization scheme:", normalization_scheme, "\n")
    
      # Create FASTA group-specific output directory with normalization subdirectory
      norm_output_dir <- file.path(HEATMAP_OUT_DIR, group, normalization_scheme)
      dir.create(norm_output_dir, recursive = TRUE, showWarnings = FALSE)
    
      # Create subdirectories for different orientation and sorting methods
      original_organ_dir <- file.path(norm_output_dir, "Original_RowGene_ColumnOrgan_Sorted_by_Organ")
      original_expr_dir <- file.path(norm_output_dir, "Original_RowGene_ColumnOrgan_Sorted_by_Expression")
      transposed_organ_dir <- file.path(norm_output_dir, "Transposed_RowOrgan_ColumnGene_Sorted_by_Organ")
      transposed_expr_dir <- file.path(norm_output_dir, "Transposed_RowOrgan_ColumnGene_Sorted_by_Expression")
      
      dir.create(original_organ_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(original_expr_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(transposed_organ_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(transposed_expr_dir, recursive = TRUE, showWarnings = FALSE)
    
      for (count_type in COUNT_TYPES) {
        # Skip 'cpm' normalization for FPKM and TPM count types
        if (normalization_scheme == "cpm" && count_type %in% c("fpkm", "tpm")) {
          next
        }
        for (gene_type in GENE_TYPES) {
          # Process both SRR and Organ label types for each input file
          for (label_type in LABEL_TYPES) {
            # Try with master reference suffix first (when QUERY_AGAINST_MASTER_REFERENCE=TRUE in bash script)
            input_file <- file.path(MATRICES_DIR, group, 
                                    paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, MASTER_REFERENCE_SUFFIX, ".tsv"))
          
            # If that doesn't exist, try without the suffix
            if (!file.exists(input_file)) {
              input_file <- file.path(MATRICES_DIR, group, 
                                        paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
            }
          
            if (!file.exists(input_file)) next
          
            # cat("  Processing:", basename(input_file), "\n")
          
            input_basename <- tools::file_path_sans_ext(basename(input_file))
            # Create cleaner title from group name only, avoiding duplication
            title <- gsub("_", " ", group)
          
            # Read raw count matrix
            raw_data <- read_count_matrix(input_file)
            if (is.null(raw_data)) {
              # cat("    Failed to read data from:", input_file, "\n")
              next
            }
          
            # Apply the specified normalization scheme
            processed_data <- apply_normalization(raw_data, normalization_scheme, count_type)
          
            if (!is.null(processed_data)) {
              # Define version configurations for all combinations
              version_configs <- list(
                # Original orientation
                list(dir = original_organ_dir, transpose = FALSE, sort_by_expression = FALSE),
                list(dir = original_expr_dir, transpose = FALSE, sort_by_expression = TRUE),
                # Transposed orientation
                list(dir = transposed_organ_dir, transpose = TRUE, sort_by_expression = FALSE),
                list(dir = transposed_expr_dir, transpose = TRUE, sort_by_expression = TRUE)
              )
              
              # Generate heatmaps with matching label_type
              for (version in version_configs) {
                sort_suffix <- ifelse(version$sort_by_expression, "sorted_by_expression", "sorted_by_organ")
                orientation_suffix <- ifelse(version$transpose, "transposed", "original")
                cv_output <- file.path(version$dir, paste0(input_basename, "_", normalization_scheme, "_", orientation_suffix, "_", sort_suffix, "_cv_heatmap.png"))
                total_heatmaps <- total_heatmaps + 1
                
                if (generate_heatmap_with_cv(processed_data, cv_output, title, count_type, label_type, 
                                            normalization_scheme, transpose = version$transpose, 
                                            sort_by_expression = version$sort_by_expression, 
                                            raw_data_matrix = raw_data)) {
                  successful_heatmaps <- successful_heatmaps + 1
                }
              }
            } else {
              # cat("    ✗ Failed to process data with normalization:", normalization_scheme, "\n")
            }
          }
        }
      }
  }
}

# ===============================================
# SUMMARY
# ===============================================

while (length(dev.list()) > 0) { dev.off() }

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("HEATMAP WITH CV GENERATION SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Total heatmaps attempted:", total_heatmaps, "\n")
cat("Successful heatmaps:", successful_heatmaps, "\n")
cat("Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
cat("Output directory:", HEATMAP_OUT_DIR, "\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

if (successful_heatmaps > 0) {
  cat("✓ Heatmap generation completed successfully!\n")
} else {
  cat("✗ No heatmaps were generated. Please check input files and paths.\n")
}