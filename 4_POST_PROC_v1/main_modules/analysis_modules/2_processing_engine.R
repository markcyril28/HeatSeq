#!/usr/bin/env Rscript

# ===============================================
# PROCESSING ENGINE FOR ALL ANALYSIS MODULES
# ===============================================
# Generic processing loop to avoid repetition across scripts

# ===============================================
# GENERIC PROCESSING FUNCTION
# ===============================================

process_all_combinations <- function(
  config,
  output_base_dir,
  processing_callback,
  min_rows = 2,
  extra_options = NULL,
  matrices_dir = MATRICES_DIR
) {
  counters <- list(total = 0, successful = 0, skipped = 0)
  
  for (gene_group in config$gene_groups) {
    cat("Processing gene group:", gene_group, "\n")
    
    gene_group_output_dir <- file.path(output_base_dir, gene_group)
    ensure_output_dir(gene_group_output_dir, clean = FALSE)
    
    for (processing_level in PROCESSING_LEVELS) {
      for (count_type in COUNT_TYPES) {
        for (gene_type in GENE_TYPES) {
          for (label_type in LABEL_TYPES) {
            
            input_file <- build_input_path(gene_group, processing_level, count_type, 
                                           gene_type, matrices_dir, config$master_reference)
            
            validation <- validate_and_read_matrix(input_file, min_rows)
            
            if (!validation$success) {
              cat("  Skipping (", validation$reason, "):", processing_level, "|", 
                  basename(input_file), "\n", sep = "")
              next
            }
            
            cat("  Processing:", processing_level, "|", count_type, "|", gene_type, "|", 
                label_type, "- Found", validation$n_genes, "genes\n")
            
            for (norm_scheme in NORM_SCHEMES) {
              data_normalized <- apply_normalization(validation$data, norm_scheme, count_type)
              
              if (is.null(data_normalized)) {
                cat("    Failed normalization:", norm_scheme, "\n")
                next
              }
              
              # Apply gene/sample label transformations based on gene_type and label_type
              raw_labeled <- apply_labels(validation$data, gene_group, gene_type, label_type)
              normalized_labeled <- apply_labels(data_normalized, gene_group, gene_type, label_type)
              
              callback_result <- processing_callback(
                gene_group = gene_group,
                gene_group_output_dir = gene_group_output_dir,
                processing_level = processing_level,
                count_type = count_type,
                gene_type = gene_type,
                label_type = label_type,
                norm_scheme = norm_scheme,
                raw_data_matrix = raw_labeled,
                normalized_data = normalized_labeled,
                overwrite = config$overwrite_existing,
                extra_options = extra_options
              )
              
              counters$total <- counters$total + callback_result$total
              counters$successful <- counters$successful + callback_result$successful
              counters$skipped <- counters$skipped + callback_result$skipped
            }
          }
        }
      }
    }
  }
  
  return(counters)
}

# ===============================================
# ORIENTATION AND SORTING OPTIONS
# ===============================================

get_sorting_options <- function() {
  list(
    list(sort = FALSE, sort_name = "Sorted_by_Organ"),
    list(sort = TRUE, sort_name = "Sorted_by_Expression")
  )
}

get_orientation_options <- function() {
  list(
    list(transpose = FALSE, orient_name = "Genes_as_Rows"),
    list(transpose = TRUE, orient_name = "Samples_as_Rows")
  )
}
