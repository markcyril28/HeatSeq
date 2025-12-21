# ===============================================
# PROCESSING ENGINE FOR METHOD 4 HEATMAP GENERATION
# ===============================================
# Purpose: Generic processing loop implementing DRY principle
# Features:
#   - Iterates through all gene groups and processing parameters
#   - Applies normalization schemes
#   - Invokes callback functions for specific visualizations
#   - Tracks success/failure statistics
# Design: Callback pattern for extensibility across different viz types
# ===============================================

source("b_modules_for_Method_4/0_shared_config.R")

#===============================================================================
# MAIN PROCESSING FUNCTION
#===============================================================================
# Generic iterator through all possible parameter combinations
# Args:
#   config: Configuration from load_runtime_config()
#   output_base_dir: Base directory for outputs
#   processing_callback: Function to call for each valid matrix
#   min_rows: Minimum rows required (default 2 for heatmaps, 1 for bar graphs)
#   extra_options: Optional parameters passed to callback
# Returns: List with total, successful, and skipped counts

process_all_combinations <- function(
  config,
  output_base_dir,
  processing_callback,
  min_rows = 2,
  extra_options = NULL
) {
  # Initialize performance counters
  counters <- list(total = 0, successful = 0, skipped = 0)
  
  # Main processing loop: gene groups
  for (gene_group in config$gene_groups) {
    cat("Processing gene group:", gene_group, "\n")
    
    # Create dedicated output directory for this gene group
    gene_group_output_dir <- file.path(output_base_dir, gene_group)
    ensure_output_dir(gene_group_output_dir, clean = TRUE)
    
    # Nested loops through all parameter combinations
    for (processing_level in PROCESSING_LEVELS) {        # gene_level or isoform_level
      for (count_type in COUNT_TYPES) {                  # expected_count
        for (gene_type in GENE_TYPES) {                  # geneID or geneName
          for (label_type in LABEL_TYPES) {              # SRR or Organ
            
            # Construct path to input matrix file
            input_file <- build_input_path(gene_group, processing_level, count_type, gene_type)
            
            # Validate matrix exists and meets minimum requirements
            validation <- validate_and_read_matrix(input_file, min_rows)
            
            if (!validation$success) {
              cat("  Skipping (", validation$reason, "):", processing_level, "|", 
                  basename(input_file), "\n", sep = "")
              next  # Skip to next combination
            }
            
            cat("  Processing:", processing_level, "|", count_type, "|", gene_type, "|", 
                label_type, "- Found", validation$n_genes, "genes\n")
            
            # Apply each normalization scheme
            for (norm_scheme in NORM_SCHEMES) {
              
              # Transform data according to normalization scheme
              data_normalized <- apply_normalization(validation$data, norm_scheme, count_type)
              
              if (is.null(data_normalized)) {
                cat("    Failed normalization:", norm_scheme, "\n")
                next
              }
              
              # Invoke the specific visualization callback function
              # Callback handles actual heatmap/graph generation
              callback_result <- processing_callback(
                gene_group = gene_group,
                gene_group_output_dir = gene_group_output_dir,
                processing_level = processing_level,
                count_type = count_type,
                gene_type = gene_type,
                label_type = label_type,
                norm_scheme = norm_scheme,
                raw_data_matrix = validation$data,      # Original data for CV calculation
                normalized_data = data_normalized,       # Transformed data for visualization
                overwrite = config$overwrite_existing,
                extra_options = extra_options
              )
              
              # Accumulate statistics
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

#===============================================================================
# HEATMAP ORIENTATION AND SORTING OPTIONS
#===============================================================================

# Get sorting configuration options
# Returns: List of sorting configurations
get_sorting_options <- function() {
  list(
    list(sort = FALSE, sort_name = "Sorted_by_Organ"),         # Original sample order
    list(sort = TRUE, sort_name = "Sorted_by_Expression")      # Order by expression level
  )
}

# Get orientation configuration options
# Returns: List of orientation configurations
get_orientation_options <- function() {
  list(
    list(transpose = FALSE, orient_name = "Original_RowGene_ColumnOrgan"),      # Genes as rows
    list(transpose = TRUE, orient_name = "Transposed_RowOrgan_ColumnGene")      # Organs as rows
  )
}

#===============================================================================
# NORMALIZATION DISPLAY NAMES
#===============================================================================

# Convert internal normalization scheme names to display-friendly labels
# Args: norm_scheme - Internal normalization identifier
# Returns: Human-readable label for visualization
get_norm_display_name <- function(norm_scheme) {
  switch(norm_scheme,
    "count_type_normalized" = "Count-Type_Normalized",
    "zscore" = "Z-score_Normalized",
    "zscore_scaled_to_ten" = "Z-score_Scaled_to_Ten",
    norm_scheme  # Fallback to original if no mapping
  )
}
