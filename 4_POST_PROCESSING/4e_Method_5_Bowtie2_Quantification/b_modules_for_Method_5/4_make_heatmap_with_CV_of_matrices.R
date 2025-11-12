#!/usr/bin/env Rscript

# ===============================================
# CV HEATMAP GENERATION FOR METHOD 5 (BOWTIE2/RSEM)
# ===============================================
# Generates heatmaps with Coefficient of Variation annotations
# from DESeq2/tximport-normalized RSEM data

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(ggplot2)
})

source("b_modules_for_Method_5/2_utility_functions.R")
source("b_modules_for_Method_5/1_processing_engine.R")

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "top"
CV_HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, "CV_Heatmaps")

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(CV_HEATMAP_OUT_DIR)

# ===============================================
# PROCESSING CALLBACK
# ===============================================

process_cv_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
                               count_type, gene_type, label_type, norm_scheme,
                               raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  local_total <- 0
  local_successful <- 0
  local_skipped <- 0
  
  # Generate versions with different sorting
  # CV heatmaps are only in original orientation (rows=genes, cols=organs)
  for (sorting in get_sorting_options()) {
    
    # Create output directory structure
    version_dir <- file.path(
      gene_group_output_dir, processing_level, count_type, gene_type,
      norm_scheme, "Original_RowGene_ColumnOrgan", sorting$sort_name
    )
    ensure_output_dir(version_dir)
    
    # Build title and output filename
    title_base <- build_title_base(gene_group, count_type, gene_type,
                                   label_type, processing_level, norm_scheme)
    output_path <- file.path(version_dir, paste0(title_base, "_CV.png"))
    
    local_total <- local_total + 1
    
    # Check if should skip
    if (should_skip_existing(output_path, overwrite)) {
      cat("      Skipping (already exists):", basename(output_path), "\n")
      local_skipped <- local_skipped + 1
      next
    }
    
    # Generate CV heatmap
    # CV is ALWAYS calculated on raw, untransformed data
    success <- generate_heatmap_with_cv(
      data_matrix = normalized_data,
      output_path = output_path,
      title = title_base,
      count_type = count_type,
      label_type = label_type,
      normalization_scheme = norm_scheme,
      sort_by_expression = sorting$sort,
      raw_data_matrix = raw_data_matrix  # For CV calculation
    )
    
    if (success) {
      local_successful <- local_successful + 1
    } else {
      cat("      Failed:", basename(output_path), "\n")
    }
  }
  
  list(total = local_total, successful = local_successful, skipped = local_skipped)
}

# ===============================================
# MAIN PROCESSING
# ===============================================

print_config_summary("CV HEATMAP GENERATION - METHOD 5 (BOWTIE2/RSEM)", config)

# Process all combinations
results <- process_all_combinations(
  config = config,
  output_base_dir = CV_HEATMAP_OUT_DIR,
  processing_callback = process_cv_heatmap,
  min_rows = 2
)

# Print summary
print_summary(results$successful, results$total, results$skipped)

if (results$successful > 0) {
  cat("Output directory:", CV_HEATMAP_OUT_DIR, "\n")
  cat("CV heatmaps organized by:\n")
  cat("  - Gene group\n")
  cat("  - Processing level (gene/isoform)\n")
  cat("  - Count type (expected_count)\n")
  cat("  - Gene identifier type (geneID/geneName)\n")
  cat("  - Normalization scheme\n")
  cat("  - Sorting (by organ/by expression)\n")
  cat("\nNote: CV values are calculated on raw, untransformed data\n")
  cat("      while the heatmap displays normalized values for visualization.\n")
}
