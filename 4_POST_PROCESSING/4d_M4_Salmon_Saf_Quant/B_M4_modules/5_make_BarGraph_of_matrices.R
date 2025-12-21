#!/usr/bin/env Rscript

# ===============================================
# BAR GRAPH GENERATION FOR METHOD 4 - SALMON SAF QUANTIFICATION
# ===============================================
# Generates individual bar graphs for each gene
# from DESeq2/tximport-normalized Salmon data

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("b_modules_for_Method_4/2_utility_functions.R")
source("b_modules_for_Method_4/1_processing_engine.R")

# ===============================================
# CONFIGURATION
# ===============================================

BAR_GRAPH_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, "III_Bar_Graphs")

# Load runtime configuration
config <- load_runtime_config()
ensure_output_dir(BAR_GRAPH_OUT_DIR)

# ===============================================
# PROCESSING CALLBACK
# ===============================================
# Generates individual bar graphs for each gene

process_bar_graphs <- function(gene_group, gene_group_output_dir, processing_level,
                               count_type, gene_type, label_type, norm_scheme,
                               raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  # Initialize counters
  local_total <- 0
  local_successful <- 0
  
  # Generate bar graphs with different sorting options
  for (sorting in get_sorting_options()) {
    
    # Create output directory structure
    version_dir <- file.path(
      gene_group_output_dir, processing_level, count_type, gene_type,
      norm_scheme, sorting$sort_name
    )
    ensure_output_dir(version_dir)
    
    # Generate bar graphs (one per gene)
    n_graphs <- generate_gene_bargraphs(
      data_matrix = normalized_data,
      output_dir = version_dir,
      title_prefix = gene_group,
      count_type = count_type,
      label_type = label_type,
      normalization_scheme = norm_scheme,
      sort_by_expression = sorting$sort,
      overwrite = overwrite
    )
    
    local_total <- local_total + nrow(normalized_data)
    local_successful <- local_successful + n_graphs
    
    if (n_graphs > 0) {
      cat("    ✓ Generated", n_graphs, "bar graphs for", processing_level, "|", 
          norm_scheme, "|", sorting$sort_name, "\n")
    }
  }
  
  list(total = local_total, successful = local_successful, skipped = 0)
}

# ===============================================
# MAIN PROCESSING
# ===============================================

print_config_summary("BAR GRAPH GENERATION - METHOD 4 (SALMON SAF)", config)

# Process all combinations (min_rows = 1 for bar graphs)
results <- process_all_combinations(
  config = config,
  output_base_dir = BAR_GRAPH_OUT_DIR,
  processing_callback = process_bar_graphs,
  min_rows = 1
)

# Print summary
print_summary(results$successful, results$total)

if (results$successful > 0) {
  cat("Output directory:", BAR_GRAPH_OUT_DIR, "\n")
  cat("Bar graphs organized by:\n")
  cat("  - Gene group\n")
  cat("  - Processing level (gene/isoform)\n")
  cat("  - Count type (expected_count)\n")
  cat("  - Gene identifier type (geneID/geneName)\n")
  cat("  - Normalization scheme\n")
  cat("  - Sorting (by organ/by expression)\n")
  cat("\nEach gene has its own individual bar graph showing expression across samples.\n")
}
