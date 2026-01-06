#!/usr/bin/env Rscript

# ===============================================
# BASIC HEATMAP GENERATION MODULE
# ===============================================
# Generates standard heatmaps from expression matrices

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", ".")
source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "2_processing_engine.R"))

# ===============================================
# CONFIGURATION
# ===============================================

LEGEND_POSITION <- "bottom"
HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$BASIC_HEATMAP)

# ===============================================
# HEATMAP GENERATION FUNCTION
# ===============================================

generate_heatmap_violet <- function(data_matrix, output_path, title, 
                                     count_type, label_type, normalization_type,
                                     transpose = FALSE, sort_by_expression = FALSE) {
  tryCatch({
    # Prepare data
    if (transpose) {
      data_matrix <- t(data_matrix)
    }
    
    if (sort_by_expression) {
      row_means <- rowMeans(data_matrix, na.rm = TRUE)
      data_matrix <- data_matrix[order(row_means, decreasing = TRUE), , drop = FALSE]
    }
    
    # Color scale
    data_range <- range(data_matrix, na.rm = TRUE)
    color_fun <- colorRamp2(
      seq(data_range[1], data_range[2], length.out = 100),
      get_violet_color_scale(100)
    )
    
    # Legend
    legend_layout <- get_legend_layout(LEGEND_POSITION)
    legend_title <- get_legend_title(normalization_type, count_type)
    
    # Create heatmap
    ht <- Heatmap(
      data_matrix,
      name = legend_title,
      col = color_fun,
      cluster_rows = !sort_by_expression,
      cluster_columns = TRUE,
      show_row_names = nrow(data_matrix) <= 50,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      column_title = title,
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      heatmap_legend_param = list(
        direction = legend_layout$direction,
        legend_height = legend_layout$height,
        legend_width = legend_layout$width
      )
    )
    
    # Save
    png(output_path, width = 1200, height = 800, res = 100)
    draw(ht, heatmap_legend_side = LEGEND_POSITION)
    dev.off()
    
    cat("      Generated:", basename(output_path), "\n")
    return(TRUE)
  }, error = function(e) {
    cat("      Error:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# PROCESSING CALLBACK
# ===============================================

process_basic_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
                                  count_type, gene_type, label_type, norm_scheme,
                                  raw_data_matrix, normalized_data, overwrite, extra_options) {
  
  local_total <- 0
  local_successful <- 0
  local_skipped <- 0
  
  norm_display <- get_norm_display_name(norm_scheme)
  
  for (orient in get_orientation_options()) {
    for (sorting in get_sorting_options()) {
      
      version_dir <- file.path(
        gene_group_output_dir, processing_level, count_type, gene_type,
        norm_scheme, orient$orient_name, sorting$sort_name
      )
      ensure_output_dir(version_dir)
      
      title_base <- build_title_base(gene_group, count_type, gene_type, 
                                     label_type, processing_level, norm_scheme)
      output_path <- file.path(version_dir, paste0(title_base, ".png"))
      
      local_total <- local_total + 1
      
      if (should_skip_existing(output_path, overwrite)) {
        cat("      Skipping (exists):", basename(output_path), "\n")
        local_skipped <- local_skipped + 1
        next
      }
      
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
      
      if (success) local_successful <- local_successful + 1
    }
  }
  
  list(total = local_total, successful = local_successful, skipped = local_skipped)
}

# ===============================================
# MAIN
# ===============================================

run_basic_heatmap <- function(config = NULL, matrices_dir = MATRICES_DIR) {
  if (is.null(config)) config <- load_runtime_config()
  ensure_output_dir(HEATMAP_OUT_DIR)
  
  print_config_summary("BASIC HEATMAP GENERATION", config)
  
  results <- process_all_combinations(
    config = config,
    output_base_dir = HEATMAP_OUT_DIR,
    processing_callback = process_basic_heatmap,
    matrices_dir = matrices_dir
  )
  
  print_summary(results$successful, results$total, results$skipped)
  cat("Output directory:", HEATMAP_OUT_DIR, "\n")
}

# Run if executed directly
if (!interactive() && identical(environment(), globalenv())) {
  run_basic_heatmap()
}
