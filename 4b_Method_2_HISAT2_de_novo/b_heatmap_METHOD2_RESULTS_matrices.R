# ===============================================
# HEATMAP GENERATION FOR METHOD2 RESULTS
# ===============================================

# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)

# ===============================================
# CONFIGURATION
# ===============================================

# Input and output directories
BASE_DIR <- getwd()
MATRICES_DIR <- "5_stringtie/a_Method2_RESULTS_matrices"
HEATMAP_OUT_DIR <- "6_Visualization/a_Method2_RESULTS_matrices_heatmaps"

# Clear and create output directory
if (dir.exists(HEATMAP_OUT_DIR)) {
  unlink(HEATMAP_OUT_DIR, recursive = TRUE)
}
dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Gene groups and versions
FASTA_GROUPS <- c(
  "TEST",
  "SmelDMP_CDS_Control_Best", 
  "SmelGIF_with_Best_Control_Cyclo",
  "SmelGRF_with_Best_Control_Cyclo",
  "SmelGRF-GIF_with_Best_Control_Cyclo"
)

VERSIONS <- c("v1", "v2")
COUNT_TYPES <- c("coverage", "fpkm", "tpm")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Sample metadata for better labeling (ordered to match bash script)
SAMPLE_LABELS <- c(
  "SRR3884631" = "Fruits_6cm",
  "SRR3884677" = "Cotyledons", 
  "SRR3884679" = "Pistils",
  "SRR3884597" = "Flowers",
  "SRR3884687" = "Buds_Opened",
  "SRR3884686" = "Buds_0.7cm",
  "SRR3884689" = "Leaves",
  "SRR3884690" = "Stems", 
  "SRR3884685" = "Radicles",
  "SRR3884675" = "Roots"
)

# ===============================================
# FUNCTIONS
# ===============================================

# Function to read and preprocess count matrix
read_count_matrix <- function(file_path) {
  tryCatch({
    # Read with na.strings to handle empty cells properly
    data <- read.table(file_path, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, 
                      na.strings = c("", " ", "NA", "null"))
    
    # Handle duplicate row names by making them unique
    if (any(duplicated(data[, 1]))) {
      cat("  ⚠️  Found duplicate row names - making them unique...\n")
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }
    
    # Set gene IDs as row names
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]
    
    # Check if data contains non-numeric columns and convert
    for (i in 1:ncol(data)) {
      if (!is.numeric(data[, i])) {
        # Convert to numeric, NA values will become NA
        data[, i] <- as.numeric(as.character(data[, i]))
      }
    }
    
    # Convert to numeric matrix
    data_matrix <- as.matrix(data)
    
    # Replace all NA values with 0
    data_matrix[is.na(data_matrix)] <- 0

    
    # Ensure all values are numeric
    if (!is.numeric(data_matrix)) {
      data_matrix <- apply(data_matrix, 2, as.numeric)
      data_matrix[is.na(data_matrix)] <- 0
    }
    
    # Remove rows that are all zeros (genes with no expression anywhere)
    row_sums <- rowSums(data_matrix, na.rm = TRUE)
    data_matrix <- data_matrix[row_sums > 0, , drop = FALSE]
    
    # Transpose the matrix to swap rows and columns
    data_matrix <- t(data_matrix)
    
    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
}

# Function to preprocess data for heatmap
preprocess_for_heatmap <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    return(NULL)
  }
  
  # Handle different count types
  if (count_type %in% c("fpkm", "tpm")) {
    # For FPKM/TPM, add small pseudocount and log transform
    data_processed <- log2(data_matrix + 0.1)
  } else {
    # For coverage, add pseudocount and log transform
    data_processed <- log2(data_matrix + 0.1)
  }
  
  # Replace any NA, Inf, or -Inf values with 0
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  # Check for rows with zero variance (all same values) which cause scaling issues
  row_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  row_vars[is.na(row_vars)] <- 0
  
  # Remove rows with zero variance to prevent scaling issues
  if (any(row_vars == 0)) {
    cat("Removing", sum(row_vars == 0), "genes with zero variance\n")
    data_processed <- data_processed[row_vars > 0, , drop = FALSE]
  }
  
  # Filter genes with low variance only if we have enough genes
  if (nrow(data_processed) > 10) {
    gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
    gene_vars[is.na(gene_vars)] <- 0
    
    if (length(gene_vars) > 100) {
      # Keep top 1000 most variable genes or all if less than 1000
      top_genes <- min(1000, length(gene_vars))
      keep_genes <- order(gene_vars, decreasing = TRUE)[1:top_genes]
      data_processed <- data_processed[keep_genes, ]
    }
  }
  
  return(data_processed)
}

# Preprocessing for raw data (no normalization, just basic cleaning)
preprocess_for_raw_data <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    return(NULL)
  }
  
  # Just clean the data without normalization
  data_processed <- data_matrix
  
  # Replace any NA, Inf, or -Inf values with 0
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  # Remove genes with zero variance across all samples
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  if (any(gene_vars == 0)) {
    data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  }
  
  return(data_processed)
}

# Preprocessing for count-type normalized data
preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    return(NULL)
  }
  
  # Handle different count types with appropriate normalization
  if (count_type == "coverage") {
    # For raw coverage counts: library size normalization + log transform
    lib_sizes <- colSums(data_matrix, na.rm = TRUE)
    lib_sizes[lib_sizes == 0] <- 1  # Avoid division by zero
    
    # Normalize by library size (CPM-like normalization)
    data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
    
    # Log2 transform with pseudocount
    data_processed <- log2(data_normalized + 1)
    
  } else if (count_type == "fpkm") {
    # FPKM is already normalized for library size and gene length
    data_processed <- log2(data_matrix + 0.1)
    
  } else if (count_type == "tpm") {
    # TPM is already normalized - just log transform
    data_processed <- log2(data_matrix + 0.1)
    
  } else {
    # Default: simple log2 transform
    data_processed <- log2(data_matrix + 1)
  }
  
  # Replace any NA, Inf, or -Inf values with 0
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  # Remove genes with zero variance across all samples
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  if (any(gene_vars == 0)) {
    data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  }
  
  return(data_processed)
}

# Preprocessing for Z-score normalized data
preprocess_for_zscore_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    return(NULL)
  }
  
  # First apply count-type normalization
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  
  if (is.null(data_processed) || nrow(data_processed) == 0) {
    return(NULL)
  }
  
  # Apply Z-score normalization across samples for each gene
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    data_processed[is.na(data_processed)] <- 0
  }
  
  return(data_processed)
}

# Function to generate heatmap
generate_heatmap <- function(data_matrix, output_path, title, count_type, label_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    cat("Skipping heatmap - no data for:", title, "\n")
    return(FALSE)
  }
  
  # Check if we have enough data for clustering
  if (nrow(data_matrix) < 2) {
    cat("Skipping heatmap - not enough genes (", nrow(data_matrix), ") for:", title, "\n")
    return(FALSE)
  }
  
  # Additional check for data quality
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    cat("Warning: Data contains NA or infinite values, cleaning...\n")
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  
  # Check if all values are the same (would cause scaling issues)
  if (length(unique(as.vector(data_matrix))) == 1) {
    cat("Skipping heatmap - all values are identical for:", title, "\n")
    return(FALSE)
  }
  
  # Define custom eggplant violet color palette
  # Gradient: white (low expression) -> deep eggplant violet (high expression)
  eggplant_colors <- colorRampPalette(c(
    "#FFFFFF",  # White (low/no expression)
    "#F3E5F5",  # Very light violet
    "#CE93D8",  # Pale violet
    "#AB47BC",  # Soft violet
    "#8E24AA",  # Light violet
    "#6A1B9A",  # Medium violet
    "#4A148C",  # Rich purple
    "#2F1B69"   # Deep eggplant (high expression)
  ))(100)
  
  # Create heatmap with original labels
  tryCatch({
    pheatmap(
      data_matrix,
      color = eggplant_colors,
      scale = "row",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = ifelse(ncol(data_matrix) > 50, FALSE, TRUE),
      legend = TRUE,
      main = title,
      fontsize = 12,
      fontsize_row = 10,
      fontsize_col = 10,
      filename = output_path,
      width = 16,
      height = 12,
      dpi = 300,
      units = "in",
      res = 300
    )
    
    # Close any open graphics devices
    if (length(dev.list()) > 0) {
      dev.off()
    }
    
    # cat("Heatmap saved:", output_path, "\n")
    return(TRUE)
  }, error = function(e) {
    # Try without scaling if row scaling fails
    cat("Row scaling failed, trying without scaling...\n")
    tryCatch({
      pheatmap(
        data_matrix,
        color = eggplant_colors,
        scale = "none",
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = ifelse(ncol(data_matrix) > 50, FALSE, TRUE),
        legend = TRUE,
        main = paste(title, "(No Scaling)"),
        fontsize = 12,
        fontsize_row = 10,
        fontsize_col = 10,
        filename = output_path,
        width = 16,
        height = 12,
        dpi = 300,
        units = "in",
        res = 300
      )
      
      # Close any open graphics devices
      if (length(dev.list()) > 0) {
        dev.off()
      }
      
      # cat("Heatmap saved (no scaling):", output_path, "\n")
      return(TRUE)
    }, error = function(e2) {
      cat("Error generating heatmap for:", title, "\n")
      cat("Error message:", e2$message, "\n")
      return(FALSE)
    })
  })
}

# Function to generate normalized heatmap with different normalization types
generate_normalized_heatmap <- function(data_matrix, output_path, title, count_type, label_type, normalization_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    cat("Skipping heatmap - no data for:", title, "\n")
    return(FALSE)
  }
  
  # Check if we have enough data for clustering
  if (nrow(data_matrix) < 2) {
    cat("Skipping heatmap - not enough genes (", nrow(data_matrix), ") for:", title, "\n")
    return(FALSE)
  }
  
  # Additional check for data quality
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    cat("Warning: Data contains NA or infinite values, cleaning...\n")
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  
  # Check if all values are the same (would cause scaling issues)
  if (length(unique(as.vector(data_matrix))) == 1) {
    cat("Skipping heatmap - all values are identical for:", title, "\n")
    return(FALSE)
  }
  
  # Define custom eggplant violet color palette
  eggplant_colors <- colorRampPalette(c(
    "#FFFFFF",  # White (low/no expression)
    "#F3E5F5",  # Very light violet
    "#CE93D8",  # Pale violet
    "#AB47BC",  # Soft violet
    "#8E24AA",  # Light violet
    "#6A1B9A",  # Medium violet
    "#4A148C",  # Rich purple
    "#2F1B69"   # Deep eggplant (high expression)
  ))(100)
  
  # Create better column labels if using SRR codes
  col_labels <- colnames(data_matrix)
  if (label_type == "SRR" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                        SAMPLE_LABELS[colnames(data_matrix)],
                        colnames(data_matrix))
  }
  
  # Determine scaling method based on normalization type
  scale_method <- switch(normalization_type,
                        "raw_normalized" = "row",
                        "Count-Type_Normalized" = "none", 
                        "Z-score_Normalized" = "none",
                        "row")  # default
  
  # Create normalized heatmap (no clustering or sorting)
  tryCatch({
    pheatmap(
      data_matrix,
      color = eggplant_colors,
      scale = scale_method,
      cluster_rows = FALSE,     # Disable clustering for genes
      cluster_cols = FALSE,     # Disable clustering for samples
      show_rownames = ifelse(nrow(data_matrix) > 50, FALSE, TRUE),
      show_colnames = TRUE,
      labels_col = col_labels,
      legend = TRUE,
      main = paste(title, "-", normalization_type),
      fontsize = 12,
      fontsize_row = 8,
      fontsize_col = 10,
      filename = output_path,
      width = 18,
      height = 14,
      dpi = 300,
      units = "in",
      res = 300,
      border_color = NA
    )
    
    # Close any open graphics devices
    if (length(dev.list()) > 0) {
      dev.off()
    }
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating normalized heatmap for:", title, "\n")
    cat("Error message:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# MAIN PROCESSING
# ===============================================

# Close any existing graphics devices before starting
while (length(dev.list()) > 0) {
  dev.off()
}

total_heatmaps <- 0
successful_heatmaps <- 0

for (group in FASTA_GROUPS) {
  for (version in VERSIONS) {
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        for (label_type in LABEL_TYPES) {
          # Define input file path based on bash script naming convention
          input_file <- file.path(MATRICES_DIR, group, version, 
                                 paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
          # Check if input file exists
          if (!file.exists(input_file)) {
            next
          }
          
          # Mirror the exact input folder structure in output
          output_dir <- file.path(HEATMAP_OUT_DIR, group, version)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Create additional empty subfolders
          dir.create(file.path(output_dir, "raw_copy"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "raw_normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Count-Type_Normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Z-score_Normalized"), showWarnings = FALSE)
          
          # Use input filename as base for output filename
          input_basename <- tools::file_path_sans_ext(basename(input_file))
          output_file <- file.path(output_dir, paste0(input_basename, "_heatmap.png"))
          
          # Create title based on input filename
          title <- gsub("_", " ", input_basename)
          
          # Read and preprocess data
          raw_data <- read_count_matrix(input_file)
          processed_data <- preprocess_for_heatmap(raw_data, count_type)
          
          # Generate heatmap
          total_heatmaps <- total_heatmaps + 1
          if (generate_heatmap(processed_data, output_file, title, count_type, label_type)) {
            successful_heatmaps <- successful_heatmaps + 1
          }
          
          # Create a copy of the original heatmap in raw_copy folder
          raw_copy_output <- file.path(output_dir, "raw_copy", paste0(input_basename, "_raw_copy_heatmap.png"))
          total_heatmaps <- total_heatmaps + 1
          if (generate_heatmap(processed_data, raw_copy_output, title, count_type, label_type)) {
            successful_heatmaps <- successful_heatmaps + 1
          }
          
          # Generate additional normalized and sorted heatmaps
          
          # 1. Raw normalized heatmap (no normalization, original gene order)
          raw_data_processed <- preprocess_for_raw_data(raw_data)
          if (!is.null(raw_data_processed)) {
            raw_output <- file.path(output_dir, "raw_normalized", paste0(input_basename, "_raw_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            if (generate_normalized_heatmap(raw_data_processed, raw_output, title, count_type, label_type, "raw_normalized")) {
              successful_heatmaps <- successful_heatmaps + 1
            }
          }
          
          # 2. Count-type normalized heatmap
          count_normalized_data <- preprocess_for_count_type_normalized(raw_data, count_type)
          if (!is.null(count_normalized_data)) {
            count_normalized_output <- file.path(output_dir, "Count-Type_Normalized", paste0(input_basename, "_count_normalized_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            if (generate_normalized_heatmap(count_normalized_data, count_normalized_output, title, count_type, label_type, "Count-Type_Normalized")) {
              successful_heatmaps <- successful_heatmaps + 1
            }
          }
          
          # 3. Z-score normalized heatmap
          zscore_normalized_data <- preprocess_for_zscore_normalized(raw_data, count_type)
          if (!is.null(zscore_normalized_data)) {
            zscore_normalized_output <- file.path(output_dir, "Z-score_Normalized", paste0(input_basename, "_zscore_normalized_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            if (generate_normalized_heatmap(zscore_normalized_data, zscore_normalized_output, title, count_type, label_type, "Z-score_Normalized")) {
              successful_heatmaps <- successful_heatmaps + 1
            }
          }
        }
      }
    }
  }
}

# ===============================================
# SUMMARY
# ===============================================

# Close any remaining graphics devices
while (length(dev.list()) > 0) {
  dev.off()
}

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("HEATMAP GENERATION SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Total heatmaps attempted:", total_heatmaps, "\n")
cat("Successful heatmaps:", successful_heatmaps, "\n")
cat("Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
cat("Output directory:", HEATMAP_OUT_DIR, "\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

if (successful_heatmaps > 0) {
  cat("Heatmap generation completed successfully!\n")
} else {
  cat("No heatmaps were generated. Please check input files and paths.\n")
}

# The nested loops already ensure one heatmap per TSV file:
# - Each combination of (group, version, count_type, gene_type, label_type) 
#   corresponds to exactly one TSV file
# - Each TSV file generates exactly one heatmap
# - The transposed matrix ensures samples are in rows, genes in columns
# - Output files use the input filename as base name

# Current structure processes:
# SmelDMP_CDS_Control_Best: 24 TSV files → 24 heatmaps
# SmelGIF_with_Best_Control_Cyclo: 24 TSV files → 24 heatmaps  
# SmelGRF-GIF_with_Best_Control_Cyclo: 24 TSV files → 24 heatmaps
# SmelGRF_with_Best_Control_Cyclo: 24 TSV files → 24 heatmaps
# TEST: 24 TSV files → 24 heatmaps (if present)
# Total: Up to 120 heatmaps for 120 TSV files
