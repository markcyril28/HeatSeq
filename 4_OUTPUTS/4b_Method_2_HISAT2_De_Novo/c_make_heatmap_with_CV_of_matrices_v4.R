# ===============================================
# CV HEATMAP GENERATION WITH COEFFICIENT OF VARIATION
# ===============================================
#
# This script generates ONLY heatmaps with coefficient of variation (CV) annotation.
# Output is saved to "7_Heatmap_with_Covariability_Visualizations" directory.
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
  #"Best_Cell_Cycle_Associated_Control_Genes",
  #"Best_Control_Genes"
  
  # Individual Gene Groups
  #"SmelDMPs",
  #"SmelGIFs",
  #"SmelGRFs",
  "Selected_GRF_GIF_Genes"
  
  # Combined Gene Groups with Control Genes
  #"SmelGIF_with_Cell_Cycle_Control_genes",
  #"SmelGRF_with_Cell_Cycle_Control_genes",
  #"SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"

)

# Analysis configuration
COUNT_TYPES <- c("coverage", "fpkm", "tpm")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Normalization schemes
NORMALIZATION_SCHEMES <- c("raw", "count_type_normalized", "zscore", "cpm", "zscore_scaled_to_ten") # <-- MODIFIED

# Mapping from SRR IDs to organ names for matrix headers
SAMPLE_LABELS <- c(
    # Roots
    "SRR3884675" = "Roots_1",      # PRJNA328564
    #"SRR20722229" = "Roots_2",      # SAMN28540077
    #"SRR31755282" = "Roots_3",      # SAMN28540068

    # Stems
    "SRR3884690" = "Stems_1",      # PRJNA328564
    #"SRR20722227" = "Stems_2",      # SAMN28540077
    #"SRR20722384" = "Stems_3",      # SAMN28540068

    # Leaves
    "SRR3884689" = "Leaves_1",      # PRJNA328564
    #"SRR20722230" = "Leaves_2",     # SAMN28540077
    #"SRR20722386" = "Leaves_3",     # SAMN28540068
    "SRR3884684" = "Senescent_leaves", # PRJNA328564

    # Buds
    "SRR3884686" = "Buds_1",       # PRJNA328564
    #"SRR21010466" = "Buds_2",       # SAMN28540077
    #"SRR20722297" = "Buds_3",       # SAMN28540068

    # Opened Buds
    "SRR3884687" = "Opened_Buds_1", # PRJNA328564

    # Flowers
    "SRR3884597" = "Flowers_1",     # PRJNA328564
    #"SRR20722234" = "Flowers_2",    # SAMN28540077
    #"SRR23909863" = "Flowers_3",    # SAMN28540068

    # Fruits
    "SRR3884631" = "Fruits_1",      # PRJNA328564
    #"SRR2072232" = "Fruits_2",      # SAMN28540077
    #"SRR20722387" = "Fruits_3",     # SAMN28540068
    "SRR3884608" = "Fruits_1cm",     # PRJNA328564
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

# Clear the output directory at the start of each run
#if (dir.exists(HEATMAP_OUT_DIR)) {
#  #cat("Clearing existing output directory:", HEATMAP_OUT_DIR, "\n")
#  unlink(HEATMAP_OUT_DIR, recursive = TRUE)
#  dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
#}

# ===============================================
# FUNCTIONS
# ===============================================

# ===============================================
# NORMALIZATION FUNCTIONS
# ===============================================

# Function to apply raw normalization (no transformation)
preprocess_for_raw <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # Return raw data without any transformation
  data_processed <- data_matrix
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  return(data_processed)
}

# Function to apply default (count-type specific) normalization
preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  if (count_type == "coverage") {
    # Coverage: Library size normalization (CPM-like) + log2 transformation
    lib_sizes <- colSums(data_matrix, na.rm = TRUE)
    lib_sizes[lib_sizes == 0] <- 1
    data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
    data_processed <- log2(data_normalized + 1)
  } else if (count_type %in% c("fpkm", "tpm")) {
    # FPKM/TPM: Direct log2 transformation
    data_processed <- log2(data_matrix + 0.1)
  } else {
    # Default fallback: log2 transformation
    data_processed <- log2(data_matrix + 1)
  }
  
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  return(data_processed)
}

# Function to apply zscore normalization
preprocess_for_zscore <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # First apply count-type normalization
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (is.null(data_processed) || nrow(data_processed) == 0) return(NULL)
  
  # Then apply Z-score standardization across samples for each gene
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    data_processed[is.na(data_processed)] <- 0
  }
  
  return(data_processed)
}

# We must assume the function preprocess_for_count_type_normalized() 
# exists elsewhere in your environment for this to run.

preprocess_for_zscore_scaled_to_ten <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # 1. Apply count-type normalization (as in your original script)
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (is.null(data_processed) || nrow(data_processed) == 0) return(NULL)
  
  # 2. Apply Z-score standardization across samples (rows)
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    # Replace NaNs (from rows with 0 variance) with 0
    data_processed[is.na(data_processed)] <- 0
  }
  
  # 3. Scale the entire Z-scored matrix to [0, 10] 
  
  # Get the global minimum and maximum of the Z-scored matrix
  min_val <- min(data_processed)
  max_val <- max(data_processed)
  val_range <- max_val - min_val
  
  if (val_range > 0) {
    # Apply standard min-max scaling formula:
    # new_val = ((old_val - min) / (max - min)) * new_range + new_min
    # Here, new_range is 10 and new_min is 0.
    data_scaled <- ((data_processed - min_val) / val_range) * 10
  } else {
    # Handle edge case where all values are identical (e.g., all 0)
    # We can just set all values to 0 (the bottom of our [0, 10] scale)
    data_scaled <- data_processed # Preserve matrix structure
    data_scaled[,] <- 0          # Set all elements to 0
  }
  
  return(data_scaled)
}

# Function to apply CPM normalization
preprocess_for_cpm <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # Calculate library sizes (total counts per sample)
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  
  # Avoid division by zero - set minimum library size to 1
  lib_sizes[lib_sizes == 0] <- 1
  
  # Calculate CPM: (counts / library_size) * 1,000,000
  data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  
  # Handle any remaining NAs or infinite values
  data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
  
  # Apply log2 transformation with pseudocount for visualization
  data_processed <- log2(data_cpm + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  return(data_processed)
}

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

# General preprocessing function that applies the specified normalization scheme
apply_normalization <- function(data_matrix, normalization_scheme, count_type) {
  switch(normalization_scheme,
    "raw" = preprocess_for_raw(data_matrix),
    "count_type_normalized" = preprocess_for_count_type_normalized(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "cpm" = preprocess_for_cpm(data_matrix),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type), # <-- MODIFIED
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}

# Function to generate heatmap with CV annotation
generate_heatmap_with_cv <- function(data_matrix, output_path, title, count_type, label_type, normalization_scheme = "count_type_normalized", sort_by_expression = FALSE) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save matrix data as TSV
  save_matrix_data(data_matrix, output_path)
  
  # Prepare column labels (developmental stages/organs)
  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }
  
  # Set column names for the matrix
  colnames(data_matrix) <- col_labels
  
  # Sort columns by mean expression if requested (for sorted_by_expression version)
  if (sort_by_expression) {
    # Calculate mean expression for each sample (column)
    col_means <- colMeans(data_matrix, na.rm = TRUE)
    # Sort columns by mean expression (lowest to highest)
    col_order <- order(col_means)
    data_matrix <- data_matrix[, col_order, drop = FALSE]
    col_labels <- col_labels[col_order]
    colnames(data_matrix) <- col_labels
  }
  
  # Calculate coefficient of variation (CV) for each gene (using original raw data for CV calculation)
  cv_values <- calculate_cv(data_matrix)
  
  # For CV heatmap visualization, always apply Z-score normalization regardless of the input normalization
  data_scaled <- t(scale(t(data_matrix)))
  data_scaled[is.na(data_scaled)] <- 0
  
  # Create CV annotation data frame
  cv_df <- data.frame(
    Gene = rownames(data_matrix),
    CV = round(cv_values, 3),
    stringsAsFactors = FALSE
  )
  
  # Truncate gene names for display
  row_labels <- truncate_labels(rownames(data_matrix), max_length = 15)
  rownames(data_scaled) <- row_labels
  
  # Define color palette to match the pasted bar exactly (REVERSED):
  # white -> yellow -> lime -> green -> teal -> cyan -> blue -> indigo -> purple (low to high expression)
  heatmap_colors <- colorRamp2(
    seq(-2, 2, length = 10),
    c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
  )
  
  tryCatch({
    # Create CV annotation
    cv_annotation <- rowAnnotation(
      CV = anno_text(
        sprintf("%.3f", cv_values),
        gp = gpar(fontsize = 8, col = "black"),
        just = "left",
        width = unit(1.8, "cm")
      ),
      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
      gap = unit(0.2, "cm"),
      simple_anno_size = unit(1.8, "cm")
    )
    
    # Create the main heatmap
    ht <- Heatmap(
      data_scaled,
      name = "z score",
      col = heatmap_colors,
      
      # Row (gene) settings
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 9),
      row_names_max_width = unit(6, "cm"),
      cluster_rows = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "complete",
      
      # Column (developmental stage) settings
      show_column_names = TRUE,
      column_names_side = "bottom",
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      cluster_columns = FALSE,
      
      # Heatmap body settings
      rect_gp = gpar(col = "white", lwd = 0.3),
      
      # Legend settings
      heatmap_legend_param = list(
        title = "z score",
        legend_direction = "horizontal",
        legend_height = unit(1.5, "cm"),
        legend_width = unit(6, "cm"),
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = unit(0.6, "cm"),
        grid_width = unit(1.2, "cm")
      ),
      
      # Title
      column_title = paste(title, "- CV Heatmap (", normalization_scheme, ifelse(sort_by_expression, " - Sorted by Expression", " - Sorted by Organ"), ")"),
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      
      # Add right annotation for CV
      right_annotation = cv_annotation
    )
    
    # Save the heatmap
    png(output_path, width = 16, height = 4.5, units = "in", res = 300)
    draw(ht, heatmap_legend_side = "bottom", gap = unit(0.5, "cm"))
    dev.off()
    
    # Also save CV data as separate file
    cv_output_path <- gsub("\\.png$", "_CV_data.tsv", output_path)
    write.table(cv_df, file = cv_output_path, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    # cat("Heatmap with CV saved to:", output_path, "\n")
    # cat("CV data saved to:", cv_output_path, "\n")
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap with CV for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# Function to truncate long labels for better readability
truncate_labels <- function(labels, max_length = 25) {
  sapply(labels, function(x) {
    if (nchar(x) > max_length) {
      paste0(substr(x, 1, max_length - 3), "...")
    } else {
      x
    }
  })
}

# Function to save matrix data as TSV
save_matrix_data <- function(data_matrix, output_path) {
  tryCatch({
    # Create matrix path by changing .png to .tsv
    matrix_path <- gsub("\\.png$", ".tsv", output_path)
    
    # Convert matrix to data frame with row names as first column
    matrix_df <- data.frame(
      Gene_ID = rownames(data_matrix),
      data_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Write to TSV file
    write.table(matrix_df, file = matrix_path, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    # cat("Matrix saved to:", matrix_path, "\n")
    return(TRUE)
  }, error = function(e) {
    cat("Error saving matrix:", e$message, "\n")
    return(FALSE)
  })
}

# Function to read and preprocess count matrix
read_count_matrix <- function(file_path) {
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, 
                       na.strings = c("", " ", "NA", "null"))
    
    if (any(duplicated(data[, 1]))) {
      # cat("   ⚠️   Found duplicate row names - making them unique...\n")
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }
    
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]
    
    for (i in 1:ncol(data)) {
      if (!is.numeric(data[, i])) {
        data[, i] <- as.numeric(as.character(data[, i]))
      }
    }
    
    data_matrix <- as.matrix(data)
    data_matrix[is.na(data_matrix)] <- 0
    
    if (!is.numeric(data_matrix)) {
      data_matrix <- apply(data_matrix, 2, as.numeric)
      data_matrix[is.na(data_matrix)] <- 0
    }
    
    # Keep all genes, including those with all zeros
    # row_sums <- rowSums(data_matrix, na.rm = TRUE)
    # data_matrix <- data_matrix[row_sums > 0, , drop = FALSE]  # REMOVED: Now keeping zero-sum rows
    
    # 🔴 Removed the transpose step – keep genes in rows, samples in columns
    
    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
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
  
  for (normalization_scheme in NORMALIZATION_SCHEMES) {
    # cat("  Processing normalization scheme:", normalization_scheme, "\n")
    
    # Create FASTA group-specific output directory with normalization subdirectory
    norm_output_dir <- file.path(HEATMAP_OUT_DIR, group, normalization_scheme)
    dir.create(norm_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create subdirectories for different sorting methods
    sorted_by_organ_dir <- file.path(norm_output_dir, "sorted_by_organ")
    sorted_by_expression_dir <- file.path(norm_output_dir, "sorted_by_expression")
    dir.create(sorted_by_organ_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sorted_by_expression_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (count_type in COUNT_TYPES) {
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
            # Generate CV heatmap sorted by organ (developmental stage order) with SRR labels
            cv_output_organ_srr <- file.path(sorted_by_organ_dir, paste0(input_basename, "_", normalization_scheme, "_SRR_cv_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            
            if (generate_heatmap_with_cv(processed_data, cv_output_organ_srr, title, count_type, "SRR", normalization_scheme, sort_by_expression = FALSE)) {
              successful_heatmaps <- successful_heatmaps + 1
              # cat("    ✓ Successfully generated SRR version (sorted by organ):", basename(cv_output_organ_srr), "\n")
            }
            
            # Generate CV heatmap sorted by organ (developmental stage order) with Organ labels
            cv_output_organ_organ <- file.path(sorted_by_organ_dir, paste0(input_basename, "_", normalization_scheme, "_Organ_cv_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            
            if (generate_heatmap_with_cv(processed_data, cv_output_organ_organ, title, count_type, "Organ", normalization_scheme, sort_by_expression = FALSE)) {
              successful_heatmaps <- successful_heatmaps + 1
              # cat("    ✓ Successfully generated Organ version (sorted by organ):", basename(cv_output_organ_organ), "\n")
            }
            
            # Generate CV heatmap sorted by expression with SRR labels
            cv_output_expr_srr <- file.path(sorted_by_expression_dir, paste0(input_basename, "_", normalization_scheme, "_SRR_cv_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            
            if (generate_heatmap_with_cv(processed_data, cv_output_expr_srr, title, count_type, "SRR", normalization_scheme, sort_by_expression = TRUE)) {
              successful_heatmaps <- successful_heatmaps + 1
              # cat("    ✓ Successfully generated SRR version (sorted by expression):", basename(cv_output_expr_srr), "\n")
            }
            
            # Generate CV heatmap sorted by expression with Organ labels
            cv_output_expr_organ <- file.path(sorted_by_expression_dir, paste0(input_basename, "_", normalization_scheme, "_Organ_cv_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            
            if (generate_heatmap_with_cv(processed_data, cv_output_expr_organ, title, count_type, "Organ", normalization_scheme, sort_by_expression = TRUE)) {
              successful_heatmaps <- successful_heatmaps + 1
              # cat("    ✓ Successfully generated Organ version (sorted by expression):", basename(cv_output_expr_organ), "\n")
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

# cat("\n", paste(rep("=", 50), collapse = ""), "\n")
# cat("HEATMAP GENERATION SUMMARY\n")
# cat(paste(rep("=", 50), collapse = ""), "\n")
# cat("Total heatmaps attempted:", total_heatmaps, "\n")
# cat("Successful heatmaps:", successful_heatmaps, "\n")
# cat("Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
# cat("Output directory:", HEATMAP_OUT_DIR, "\n")
# cat(paste(rep("=", 50), collapse = ""), "\n")

 if (successful_heatmaps > 0) {
    # cat("Heatmap generation completed successfully!\n")
 } else {
   cat("No heatmaps were generated. Please check input files and paths.\n")
 }