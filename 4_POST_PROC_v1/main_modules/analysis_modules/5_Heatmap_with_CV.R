#!/usr/bin/env Rscript

# ===============================================
# HEATMAP WITH CV (COEFFICIENT OF VARIATION) MODULE
# ===============================================
# Generates heatmaps with CV annotations per gene

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

CV_HEATMAP_OUT_DIR <- file.path(CONSOLIDATED_BASE_DIR, OUTPUT_SUBDIRS$CV_HEATMAP)

# ===============================================
# CV CALCULATION
# ===============================================

# Calculate Coefficient of Variation: CV = (SD / mean) * 100
# NOTE: CV is most meaningful on raw or linear-scale data, not log-transformed.
#       For log data, consider using SD directly as a variability measure.
calculate_cv <- function(data_matrix, is_log_scale = FALSE) {
  apply(data_matrix, 1, function(row) {
    # For log-transformed data, return SD as variability measure
    if (is_log_scale) {
      row_sd <- sd(row, na.rm = TRUE)
      return(if (is.finite(row_sd)) row_sd else NA)
    }
    # For linear data, calculate true CV
    row_mean <- mean(row, na.rm = TRUE)
    row_sd <- sd(row, na.rm = TRUE)
    if (row_mean > 0 && is.finite(row_mean) && is.finite(row_sd)) {
      return(row_sd / row_mean * 100)
    }
    return(NA)
  })
}

# ===============================================
# CV HEATMAP GENERATION
# ===============================================

generate_heatmap_with_cv <- function(data_matrix, output_path, title,
                                      count_type, normalization_type,
                                      transpose = FALSE, sort_by_expression = FALSE) {
  tryCatch({
    # Determine if data is log-scale based on normalization type
    is_log_scale <- normalization_type %in% c("count_type_normalized", "cpm", 
                                               "deseq2_normalized", "Count-Type_Normalized")
    
    # Calculate CV before any transformation (use raw scale CV if possible)
    cv_values <- calculate_cv(data_matrix, is_log_scale = is_log_scale)
    
    if (transpose) {
      data_matrix <- t(data_matrix)
    }
    
    if (sort_by_expression) {
      row_means <- rowMeans(data_matrix, na.rm = TRUE)
      order_idx <- order(row_means, decreasing = TRUE)
      data_matrix <- data_matrix[order_idx, , drop = FALSE]
      if (!transpose) cv_values <- cv_values[order_idx]
    }
    
    # Color scales
    data_range <- range(data_matrix, na.rm = TRUE)
    color_fun <- colorRamp2(
      seq(data_range[1], data_range[2], length.out = 100),
      get_violet_color_scale(100)
    )
    
    cv_range <- range(cv_values, na.rm = TRUE)
    cv_color_fun <- colorRamp2(
      c(cv_range[1], mean(cv_range), cv_range[2]),
      c("#2166AC", "#F7F7F7", "#B2182B")
    )
    
    # Main heatmap
    ht <- Heatmap(
      data_matrix,
      name = get_legend_title(normalization_type, count_type),
      col = color_fun,
      cluster_rows = !sort_by_expression,
      cluster_columns = TRUE,
      show_row_names = nrow(data_matrix) <= 50,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 10),
      column_title = title,
      column_title_gp = gpar(fontsize = 12, fontface = "bold")
    )
    
    # CV annotation (only if not transposed)
    if (!transpose && length(cv_values) == nrow(data_matrix)) {
      cv_anno <- rowAnnotation(
        CV = cv_values,
        col = list(CV = cv_color_fun),
        annotation_legend_param = list(
          CV = list(title = "CV (%)", at = seq(0, 100, 25))
        )
      )
      ht <- ht + cv_anno
    }
    
    # Save
    png(output_path, width = 1400, height = 800, res = 100)
    draw(ht, heatmap_legend_side = "bottom")
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

process_cv_heatmap <- function(gene_group, gene_group_output_dir, processing_level,
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
      output_path <- file.path(version_dir, paste0(title_base, "_with_CV.png"))
      
      local_total <- local_total + 1
      
      if (should_skip_existing(output_path, overwrite)) {
        local_skipped <- local_skipped + 1
        next
      }
      
      success <- generate_heatmap_with_cv(
        data_matrix = normalized_data,
        output_path = output_path,
        title = paste0(title_base, " (with CV)"),
        count_type = count_type,
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

run_cv_heatmap <- function(config = NULL, matrices_dir = MATRICES_DIR) {
  if (is.null(config)) config <- load_runtime_config()
  ensure_output_dir(CV_HEATMAP_OUT_DIR)
  
  print_config_summary("HEATMAP WITH CV GENERATION", config)
  
  results <- process_all_combinations(
    config = config,
    output_base_dir = CV_HEATMAP_OUT_DIR,
    processing_callback = process_cv_heatmap,
    matrices_dir = matrices_dir
  )
  
  print_summary(results$successful, results$total, results$skipped)
}

if (!interactive() && identical(environment(), globalenv())) {
  run_cv_heatmap()
}
