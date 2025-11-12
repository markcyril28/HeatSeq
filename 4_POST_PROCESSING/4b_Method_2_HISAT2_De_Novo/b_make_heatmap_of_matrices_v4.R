# ===============================================
# HEATMAP GENERATION FOR METHOD2 RESULTS
# ===============================================
#
# Generates heatmaps with ComplexHeatmap package using multiple normalization methods:
# - Raw, Count-Type (log2), Z-score, and Z-score Scaled normalization
# - Professional layout with genes on left, samples on top
# - Saves both PNG heatmaps and TSV data files
#
# Important: Coverage, FPKM, and TPM from StringTie have different normalization states:
# - Coverage: RAW read depth (NOT library-size normalized - REQUIRES CPM normalization)
# - FPKM/TPM: Library-size and gene-length normalized (preferred for cross-gene comparisons)
# - CPM normalization is INAPPROPRIATE for FPKM/TPM (already normalized)
#
# Normalization improvements (v4):
# - Consistent pseudocount (1) for all metrics (appropriate for CPM/FPKM/TPM scale)
# - Per-gene scaling for Z-score_Scaled_to_Ten (preserves biological variance structure)
# - Configurable high-variance gene filtering with user warnings
# - Enhanced documentation and best-practice recommendations
#
# Required packages: ComplexHeatmap, circlize, RColorBrewer, dplyr, tibble, grid
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
HEATMAP_OUT_DIR <- "6_Basic_Heatmap_Visualizations"

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
LABEL_TYPES <- c("SRR", "Organ")

# Gene filtering configuration (for better visualization of large gene sets)
ENABLE_VARIANCE_FILTERING <- TRUE  # Set to FALSE to disable gene filtering
MAX_GENES_FOR_HEATMAP <- 1000      # Maximum genes to display (keeps top by variance)
MIN_GENES_FOR_FILTERING <- 100     # Only filter if gene count exceeds this threshold

# Sample label mapping: SRR IDs to developmental stage/organ names
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

# Create base output directory if it doesn't exist
dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# FUNCTIONS
# ===============================================

# Load utility functions
source("9_utility_functions.R")

# Function to preprocess data based on normalization type
preprocess_data <- function(data_matrix, count_type, norm_type = "default") {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)

  # Handle NAs and infinite values
  data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0

  # Apply normalization based on type
  data_processed <- switch(norm_type,
    "raw" = data_matrix,
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type),
    "default" = preprocess_for_count_type_normalized(data_matrix, count_type)
  )

  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0

  # Filter high-variance genes if too many genes (for better visualization)
  if (ENABLE_VARIANCE_FILTERING && nrow(data_processed) > MIN_GENES_FOR_FILTERING && norm_type != "raw") {
    gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
    gene_vars[is.na(gene_vars)] <- 0
    if (length(gene_vars) > MIN_GENES_FOR_FILTERING) {
      original_gene_count <- nrow(data_processed)
      top_genes <- min(MAX_GENES_FOR_HEATMAP, length(gene_vars))
      keep_genes <- order(gene_vars, decreasing = TRUE)[1:top_genes]
      data_processed <- data_processed[keep_genes, , drop = FALSE]
      
      if (original_gene_count > top_genes) {
        cat("  ⚠️  High-variance filtering: Kept", top_genes, "of", original_gene_count, 
            "genes (", round((top_genes/original_gene_count)*100, 1), "%) for visualization.\n")
        cat("      Filtering criterion: Top genes by variance to improve heatmap clarity.\n")
      }
    }
  }

  return(data_processed)
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

for (group in FASTA_GROUPS) {
  # Clear and recreate output directory for this specific group
  group_output_dir <- file.path(HEATMAP_OUT_DIR, group)
  if (dir.exists(group_output_dir)) {
    unlink(group_output_dir, recursive = TRUE)
  }
  dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (count_type in COUNT_TYPES) {
    for (gene_type in GENE_TYPES) {
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
        
          output_dir <- file.path(HEATMAP_OUT_DIR, group)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          
          raw_data <- read_count_matrix(input_file)
          if (is.null(raw_data)) next
          
          input_basename <- tools::file_path_sans_ext(basename(input_file))
          title <- gsub("_", " ", input_basename)
          
          # Generate all normalization types (CPM removed - inappropriate for pre-normalized data)
          normalization_configs <- list(
            list(type = "raw", data = raw_data, subdir = "raw"),
            list(type = "Count-Type_Normalized", data = preprocess_data(raw_data, count_type, "default"), subdir = "Count-Type_Normalized"),
            list(type = "Z-score_Normalized", data = preprocess_data(raw_data, count_type, "zscore"), subdir = "Z-score_Normalized"),
            list(type = "Z-score_Scaled_to_Ten", data = preprocess_data(raw_data, count_type, "zscore_scaled_to_ten"), subdir = "Z-score_Scaled_to_Ten")
          )
          
          for (config in normalization_configs) {
            if (!is.null(config$data)) {
              # Define all 4 version configurations with subdirectories
              version_configs <- list(
                list(transpose = FALSE, sort_by_expression = FALSE, subdir = "Original_RowGene_ColumnOrgan_Sorted_by_Organ"),
                list(transpose = FALSE, sort_by_expression = TRUE, subdir = "Original_RowGene_ColumnOrgan_Sorted_by_Expression"),
                list(transpose = TRUE, sort_by_expression = FALSE, subdir = "Transposed_RowOrgan_ColumnGene_Sorted_by_Organ"),
                list(transpose = TRUE, sort_by_expression = TRUE, subdir = "Transposed_RowOrgan_ColumnGene_Sorted_by_Expression")
              )
              
              for (version in version_configs) {
                # Create version-specific directory under normalization scheme
                version_dir <- file.path(output_dir, config$subdir, version$subdir)
                dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)
                
                output_path <- file.path(version_dir, paste0(input_basename, "_", gsub("-", "_", tolower(config$type)), "_heatmap.png"))
                total_heatmaps <- total_heatmaps + 1
                if (generate_heatmap_violet(config$data, output_path, title, count_type, label_type, 
                                           config$type, transpose = version$transpose, 
                                           sort_by_expression = version$sort_by_expression)) {
                  successful_heatmaps <- successful_heatmaps + 1
                }
              }
            }
          }
      }
    }
  }
}

# ===============================================
# SUMMARY REPORT
# ===============================================

# Ensure all graphics devices are closed
while (length(dev.list()) > 0) { dev.off() }

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("HEATMAP GENERATION SUMMARY REPORT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Processing Summary:\n")
cat("  • Total heatmaps attempted:", total_heatmaps, "\n")
cat("  • Successful heatmaps generated:", successful_heatmaps, "\n")
cat("  • Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
cat("  • Success rate:", round((successful_heatmaps/total_heatmaps)*100, 1), "%\n")
cat("\nOutput Details:\n")
cat("  • Output directory:", HEATMAP_OUT_DIR, "\n")
cat("  • Normalization types: Raw, Count-Type (Log2), Z-score, Z-score_Scaled_to_Ten\n")
cat("  • File formats: PNG heatmaps + TSV data matrices\n")
cat("  • Coverage normalization: CPM + log2(CPM + 1) - pseudocount = 1\n")
cat("  • FPKM/TPM normalization: log2(value + 1) - pseudocount = 1\n")
cat("  • Z-score_Scaled_to_Ten: Per-gene scaling (preserves variance structure)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

if (successful_heatmaps > 0) {
  cat("✅ Heatmap generation completed successfully!\n")
  cat("📁 Check the output directory for visualization files.\n")
} else {
  cat("❌ No heatmaps were generated.\n")
  cat("🔍 Please check input files and paths.\n")
}

cat("\nAnalysis includes gene groups:", length(FASTA_GROUPS), "groups\n")
cat("Count types processed:", paste(COUNT_TYPES, collapse = ", "), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ===============================================
# NORMALIZATION METHODS SUMMARY
# ===============================================
# 
# This script processes StringTie output metrics and generates heatmaps with 4 normalization schemes:
#
# 1. RAW: Unprocessed values from StringTie
#    - Coverage: Raw read depth (NOT library-size normalized, NOT gene-length normalized)
#    - FPKM/TPM: Fully normalized values
#    - Use case: Quality control, initial data inspection
#
# 2. COUNT-TYPE NORMALIZED (Log2): Count-type specific normalization
#    - Coverage: log2(CPM + 1) where CPM = (counts / library_size) × 1,000,000
#    - FPKM/TPM: log2(value + 1) - already normalized, only variance stabilization
#    - Use case: General expression visualization, reduces dynamic range
#    - Note: Coverage requires CPM normalization as StringTie coverage is NOT library-size normalized
#    - Pseudocount: Consistent value of 1 used for all metrics (appropriate for CPM/FPKM/TPM scale)
#
# 3. Z-SCORE NORMALIZED: Per-gene standardization across samples
#    - First: Apply count-type normalization (see #2 above)
#    - Then: Z-score (mean=0, sd=1 per gene across samples)
#    - Use case: Identifying expression patterns, removes baseline differences between genes
#    - Interpretation: Positive = above-average expression, Negative = below-average
#
# 4. Z-SCORE SCALED TO [0-10]: Z-score with per-gene scaling to fixed range
#    - First: Apply count-type normalization and z-score (see #2 and #3 above)
#    - Then: Per-gene min-max scaling to 0-10 range (each gene scaled independently)
#    - Use case: Consistent color mapping while preserving biological variance structure
#    - Interpretation: For each gene independently, 0 = lowest expression, 10 = highest expression
#    - Note: Per-gene scaling preserves expression patterns better than global scaling
#    - Important: Scale is per-gene, so values are comparable within rows, not between rows
#
# IMPORTANT NOTES:
# - Coverage from StringTie is RAW read depth and REQUIRES CPM normalization for cross-sample comparison
# - FPKM and TPM are already library-size and gene-length normalized - preferred for cross-gene comparisons
# - Pseudocount of 1 is used consistently for all metrics (appropriate for CPM/FPKM/TPM scale: 0.01-10,000+)
# - All z-score normalizations are per-gene (row-wise), preserving sample-to-sample relationships
# - Z-score_Scaled_to_Ten uses PER-GENE scaling (not global) to preserve biological variance structure
# - Optional high-variance filtering can be controlled via ENABLE_VARIANCE_FILTERING parameter
# ===============================================
