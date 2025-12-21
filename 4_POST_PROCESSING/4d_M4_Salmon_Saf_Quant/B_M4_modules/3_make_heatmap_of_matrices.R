#!/usr/bin/env Rscript

# ===============================================
# HEATMAP GENERATION FOR METHOD 4 - SALMON SAF QUANTIFICATION
# ===============================================
# Generates heatmaps from DESeq2/tximport-normalized Salmon data
# Multiple normalization schemes and output formats

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

source("b_modules_for_Method_4/2_utility_functions.R")
source("b_modules_for_Method_4/1_processing_engine.R")

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "bottom"
HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, "I_Basic_Heatmap")

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(HEATMAP_OUT_DIR)

# ===============================================
# PROCESSING CALLBACK
# ===============================================
# Generates basic heatmaps with various orientations and sorting options

process_basic_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
                                  count_type, gene_type, label_type, norm_scheme,
                                  raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  # Initialize counters for tracking success/failure
  local_total <- 0
  local_successful <- 0
  local_skipped <- 0
  
  norm_display <- get_norm_display_name(norm_scheme)
  
  # Generate all versions (transposed/sorted combinations)
  for (orient in get_orientation_options()) {
    for (sorting in get_sorting_options()) {
      
      # Create hierarchical output directory structure
      version_dir <- file.path(
        gene_group_output_dir, processing_level, count_type, gene_type,
        norm_scheme, orient$orient_name, sorting$sort_name
      )
      ensure_output_dir(version_dir)
      
      # Build title and output filename
      title_base <- build_title_base(gene_group, count_type, gene_type, 
                                     label_type, processing_level, norm_scheme)
      output_path <- file.path(version_dir, paste0(title_base, ".png"))
      
      local_total <- local_total + 1
      
      # Check if should skip existing file
      if (should_skip_existing(output_path, overwrite)) {
        cat("      Skipping (already exists):", basename(output_path), "\n")
        local_skipped <- local_skipped + 1
        next
      }
      
      # Generate heatmap with specified parameters
      success <- generate_heatmap_violet(
        data_matrix = normalized_data,
        output_path = output_path,
        title = title_base,
        count_type = count_type,
        label_type = label_type,
        normalization_type = norm_display,
        transpose = orient$transpose,
        sort_by_expression = sorting$sort
      )
      
      if (success) {
        local_successful <- local_successful + 1
        cat("      ✓ Generated:", basename(output_path), "\n")
      } else {
        cat("      ✗ Failed:", basename(output_path), "\n")
      }
    }
  }
  
  list(total = local_total, successful = local_successful, skipped = local_skipped)
}

# ===============================================
# MAIN PROCESSING
# ===============================================

print_config_summary("BASIC HEATMAP GENERATION - METHOD 4 (SALMON SAF)", config)

# Process all combinations
results <- process_all_combinations(
  config = config,
  output_base_dir = HEATMAP_OUT_DIR,
  processing_callback = process_basic_heatmap,
  min_rows = 2
)

# Print summary
print_summary(results$successful, results$total, results$skipped)

if (results$successful > 0) {
  cat("Output directory:", HEATMAP_OUT_DIR, "\n")
  cat("Heatmaps organized by:\n")
  cat("  - Gene group\n")
  cat("  - Processing level (gene/isoform)\n")
  cat("  - Count type (expected_count)\n")
  cat("  - Gene identifier type (geneID/geneName)\n")
  cat("  - Normalization scheme\n")
  cat("  - Orientation (original/transposed)\n")
  cat("  - Sorting (by organ/by expression)\n")
}
