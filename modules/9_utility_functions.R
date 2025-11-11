# ===============================================
# UTILITY FUNCTIONS FOR HEATSEQ R SCRIPTS
# ===============================================
#
# This script contains shared functions used by:
# - b_make_heatmap_of_matrices_v4.R
# - c_make_heatmap_with_CV_of_matrices_v4.R
# - d_make_BarGraph_of_matrices_v4.R
#
# ===============================================

# ===============================================
# FILE I/O AND DATA PREPARATION
# ===============================================

# Function to read and preprocess count matrix
# Inputs:
#  - file_path: path to a tab-delimited file where the first column is gene IDs
# Outputs:
#  - numeric matrix with gene IDs as rownames, NA/Inf coerced to 0
# Notes:
#  - duplicates in the first column are made unique using make.unique()
#  - non-numeric columns (after the first) are coerced to numeric
read_count_matrix <- function(file_path) {
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE,
                      na.strings = c("", " ", "NA", "null"))

    # Make duplicate gene IDs unique (Gene, Gene.1, ...) to avoid rowname collisions
    if (any(duplicated(data[, 1]))) {
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }

    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]

    # Coerce each column to numeric where possible (preserves NAs which are handled below)
    for (i in seq_len(ncol(data))) {
      if (!is.numeric(data[, i])) {
        data[, i] <- as.numeric(as.character(data[, i]))
      }
    }

  data_matrix <- as.matrix(data)
  # Replace NA values with 0 (common for count matrices)
  data_matrix[is.na(data_matrix)] <- 0

    # Ensure storage mode is numeric
    if (!is.numeric(data_matrix)) {
      data_matrix <- apply(data_matrix, 2, as.numeric)
      data_matrix[is.na(data_matrix)] <- 0
    }

    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
}

# Function to save matrix data as TSV
# Inputs:
#  - data_matrix: numeric matrix with rownames representing gene IDs
#  - output_path: target filename (extension will be replaced with .tsv)
# Output: TRUE on success, FALSE on failure
save_matrix_data <- function(data_matrix, output_path) {
  tryCatch({
    # Create matrix path by changing extension
    matrix_path <- paste0(tools::file_path_sans_ext(output_path), ".tsv")

    # Convert matrix to data frame with row names as first column
    # Convert matrix to data frame with gene IDs in first column
    matrix_df <- data.frame(
      Gene_ID = rownames(data_matrix),
      data_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Write to TSV file
    write.table(matrix_df, file = matrix_path, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    return(TRUE)
  }, error = function(e) {
    cat("Error saving matrix:", e$message, "\n")
    return(FALSE)
  })
}

# Function to truncate long labels for better readability
# Inputs:
#  - labels: vector of strings (can be NULL)
#  - max_length: maximum allowed length; longer labels will be truncated with '...'
# Output: character vector of same length with truncated labels
truncate_labels <- function(labels, max_length = 25) {
  # Ensure labels is a vector
  if (!is.vector(labels) && !is.null(labels)) {
    labels <- as.character(labels)
  }
  
  # Handle NULL or empty input
  if (is.null(labels) || length(labels) == 0) {
    return(character(0))
  }
  
  # Convert to character vector
  labels <- as.character(labels)
  
  # Apply truncation to each label; preserve order and return character vector
  sapply(labels, function(x) {
    if (is.na(x) || is.null(x)) {
      return("")
    }
    x <- as.character(x)
    if (nchar(x) > max_length) {
      paste0(substr(x, 1, max_length - 3), "...")
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

# ===============================================
# NORMALIZATION FUNCTIONS
# ===============================================

# Function to apply raw normalization (no transformation)
# Inputs: numeric matrix
# Output: same matrix with NA/Inf replaced by 0
preprocess_for_raw <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- data_matrix
  # Replace NA and infinite values with zero to avoid downstream errors
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

# Function to apply raw_normalized (CPM without log)
# Converts counts to CPM (counts per million) but does NOT log-transform
# Inputs: numeric count matrix
# Output: CPM matrix (same shape) with NA/Inf -> 0
preprocess_for_raw_normalized <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  # Use column sums as library sizes; guard against zero-sum columns
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- 1
  data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_normalized[is.na(data_normalized) | is.infinite(data_normalized)] <- 0
  return(data_normalized)
}

# Function to apply CPM normalization (for coverage)
# Inputs:
#  - data_matrix: numeric matrix
#  - count_type: expected to be "coverage" for CPM
# Output: log2(CPM + 1) matrix with NA/Inf -> 0
preprocess_for_cpm <- function(data_matrix, count_type = "coverage") {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  if (count_type != "coverage") {
    stop("ERROR: CPM normalization can only be applied to 'coverage' count type.")
  }
  # Calculate library sizes and convert to CPM
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- 1
  data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
  # Use log2(CPM + 1) to stabilize variance
  data_processed <- log2(data_cpm + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

# Function for count-type specific normalization
# - For 'coverage' uses CPM (with log); otherwise applies log2(x + 1)
# Returns a numeric matrix with NA/Inf replaced by 0
preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  if (count_type == "coverage") {
    data_processed <- preprocess_for_cpm(data_matrix, count_type)
  } else {
    # For FPKM/TPM or other count types, use log2(x+1)
    data_processed <- log2(data_matrix + 1)
  }
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

# Function for z-score normalization (DATASET-WIDE)
# Computes z-scores using global mean and SD across the entire matrix (not per-gene).
# Inputs:
#  - data_matrix: numeric matrix
#  - count_type: forwarded to count-type normalization to ensure consistent preprocessing
# Output: matrix of z-scores (global) with NA/Inf -> 0
preprocess_for_zscore <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_norm <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
    # DATASET-WIDE z-score: subtract global mean and divide by global SD
    global_mean <- mean(data_norm, na.rm = TRUE)
    global_sd <- sd(as.vector(data_norm), na.rm = TRUE)
    
    if (global_sd > 0 && is.finite(global_sd)) {
      data_zscore <- (data_norm - global_mean) / global_sd
    } else {
      # If SD is zero/invalid, just center by mean
      data_zscore <- data_norm - global_mean
    }
    
    data_zscore[is.na(data_zscore) | is.infinite(data_zscore)] <- 0
    return(data_zscore)
  } else {
    # If matrix too small to z-score, return the normalized matrix
    return(data_norm)
  }
}

# Function for z-score scaled to [0, 10] with GLOBAL scaling
# Steps:
# 1) Apply count-type normalization (log or CPM)
# 2) Compute dataset-wide z-score (global mean and SD)
# 3) Min-max scale to the [0, 10] range globally
# Useful when a bounded visualization range is required
preprocess_for_zscore_scaled_to_ten <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    # DATASET-WIDE z-score: center and scale using global mean/SD
    global_mean <- mean(data_processed, na.rm = TRUE)
    global_sd <- sd(as.vector(data_processed), na.rm = TRUE)
    
    if (global_sd > 0 && is.finite(global_sd)) {
      data_processed <- (data_processed - global_mean) / global_sd
    } else {
      data_processed <- data_processed - global_mean
    }
    data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  }

  # Apply GLOBAL min-max scaling to the range [0, 10]
  min_val <- min(data_processed, na.rm = TRUE)
  max_val <- max(data_processed, na.rm = TRUE)
  val_range <- max_val - min_val
  
  if (val_range > 0 && is.finite(val_range)) {
    data_scaled <- ((data_processed - min_val) / val_range) * 10
  } else {
    # If constant matrix, return zeros of same shape
    data_scaled <- matrix(0, nrow = nrow(data_processed), ncol = ncol(data_processed))
  }

  data_scaled[is.na(data_scaled) | is.infinite(data_scaled)] <- 0
  return(data_scaled)
}


# General preprocessing function
# Convenience wrapper that selects the requested normalization scheme.
# Inputs:
#  - data_matrix: numeric matrix
#  - normalization_scheme: one of the supported keys (see switch cases)
#  - count_type: forwarded where needed (e.g., 'coverage')
apply_normalization <- function(data_matrix, normalization_scheme, count_type) {
  switch(normalization_scheme,
    "raw" = preprocess_for_raw(data_matrix),
    "raw_normalized" = preprocess_for_raw_normalized(data_matrix),
    "count_type_normalized" = preprocess_for_count_type_normalized(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type),
    "cpm" = preprocess_for_cpm(data_matrix, count_type),
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}


# ===============================================
# BASIC HEATMAP GENERATION
# ===============================================

# Function to generate heatmap (unified for all normalization types)
# Inputs:
#  - data_matrix: numeric matrix (genes x samples)
#  - output_path: path to save PNG
#  - title: plot title
#  - count_type: used for legend labelling and preprocessing choices
#  - label_type: 'Organ' or 'SRR' affects column labeling
#  - normalization_type: textual label used to choose palettes/legends
# Returns: TRUE on success, FALSE on failure
generate_heatmap_violet <- function(data_matrix, output_path, title, count_type, label_type, normalization_type = "raw") {
  # Basic sanity checks
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || nrow(data_matrix) < 2) return(FALSE)
  # Replace problematic numeric values
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  # If matrix is constant, heatmap not informative
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save matrix data as TSV for downstream inspection
  save_matrix_data(data_matrix, output_path)
  
  # Prepare column labels (replace with SAMPLE_LABELS if provided and label_type == 'Organ')
  col_labels <- colnames(data_matrix)
  if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
    colnames(data_matrix) <- col_labels
  }
  
  # Choose palette: a special case for Z-score scaled to ten (0..10) vs others
  if (normalization_type == "Z-score_Scaled_to_Ten") {
    # For Z-score scaled to ten: use fixed breaks in [0,10]
    eggplant_colors <- colorRamp2(
      seq(0, 10, length = 8),
      c("#dab3ddff","#d9afe0ff","#c57fd1ff","#ac44beff",
        "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
    )
  } else {
    # For other schemes, map colors across observed data range
    eggplant_colors <- colorRamp2(
      seq(min(data_matrix), max(data_matrix), length = 8),
      c("#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
        "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
    )
  }
  
  # Legend title mapping (note: callers control normalization_type string)
  legend_title <- switch(normalization_type,
    "raw" = paste0("Raw ", toupper(count_type)),
    "Count-Type_Normalized" = "Log2 Expression",
    "Z-score_Normalized" = "Z-score",
    "Z-score_Scaled_to_Ten" = "Z-Score Scaled to [0-10]",
    "Expression"  # Fallback
  )
  
  # Configure legend placement parameters based on global LEGEND_POSITION
  if (exists("LEGEND_POSITION") && LEGEND_POSITION == "bottom") {
    legend_direction <- "horizontal"
    legend_height <- unit(2, "cm")
    legend_width <- unit(10, "cm")
    grid_height <- unit(0.8, "cm")
    grid_width <- unit(2.5, "cm")
  } else {  # "right" (default)
    legend_direction <- "vertical"
    legend_height <- unit(10, "cm")
    legend_width <- unit(2, "cm")
    grid_height <- unit(2.5, "cm")
    grid_width <- unit(0.8, "cm")
  }
  
  tryCatch({
    # Truncate row names for display if many rows
    max_length <- if (nrow(data_matrix) > 50) 20 else 30
    
    # Ensure rownames exist
    if (is.null(rownames(data_matrix)) || length(rownames(data_matrix)) == 0) {
      rownames(data_matrix) <- paste0("Gene_", seq_len(nrow(data_matrix)))
    }
    
    row_labels <- truncate_labels(rownames(data_matrix), max_length = max_length)
    rownames(data_matrix) <- row_labels
    
    # Build the Heatmap object (ComplexHeatmap)
    ht <- Heatmap(
      data_matrix,
      name = legend_title,
      col = eggplant_colors,
      
      # Row (gene) settings
      show_row_names = ifelse(nrow(data_matrix) > 50, FALSE, TRUE),
      row_names_side = "left",
      row_names_gp = gpar(fontsize = if (nrow(data_matrix) > 50) 7 else 8),
      row_names_max_width = unit(10, "cm"),
      cluster_rows = FALSE,
      
      # Column (sample) settings
      show_column_names = TRUE,
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      cluster_columns = FALSE,
      
      # Heatmap body settings
      rect_gp = gpar(col = "white", lwd = if (nrow(data_matrix) > 50) 0.3 else 0.5),
      
      # Legend settings (dynamic based on LEGEND_POSITION)
      heatmap_legend_param = list(
        title = legend_title,
        legend_direction = legend_direction,
        legend_height = legend_height,
        legend_width = legend_width,
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        grid_height = grid_height,
        grid_width = grid_width,
        at = if (normalization_type == "Z-score_Scaled_to_Ten") seq(0, 10, by = 2) else NULL,
        just = c("left", "top")
      ),
      
      # Title
      column_title = if (normalization_type == "raw") title else paste(title, "-", normalization_type),
      column_title_gp = gpar(fontsize = 10, fontface = "bold")
    )
    
    # Calculate dynamic image dimensions based on matrix size
    n_rows <- nrow(data_matrix)
    n_cols <- ncol(data_matrix)
    base_cell_size <- 0.3
    base_width <- n_cols * base_cell_size
    base_height <- n_rows * base_cell_size
    padding <- 1
    final_width <- max(8, min(20, base_width + 2 * padding))
    final_height <- max(6, min(16, base_height + 2 * padding))
    
    # Write PNG using fixed size (commented alternative uses computed sizes)
    png(output_path, width = 14, height = 10.5, units = "in", res = 500)
    #png(output_path, width = final_width, height = final_height, units = "in", res = 300)
    draw(ht, 
         heatmap_legend_side = if (exists("LEGEND_POSITION")) LEGEND_POSITION else "right", 
         annotation_legend_side = if (exists("LEGEND_POSITION")) LEGEND_POSITION else "right",
         merge_legends = TRUE,
         padding = unit(c(padding, padding, padding, padding), "inches"))
    dev.off()
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# HEATMAP WITH CV GENERATION
# ===============================================

## Function to generate heatmap with coefficient of variation annotation
## Inputs:
##  - data_matrix: matrix used for the heatmap (may be transformed for viz)
##  - output_path: PNG path to write
##  - title: plot title
##  - count_type: used for legend/labeling
##  - label_type: how to label columns ('SRR' or 'Organ')
##  - normalization_scheme: string describing how data_matrix was normalized
##  - sort_by_expression: logical; if TRUE, columns are ordered by mean expression
##  - raw_data_matrix: raw/untransformed counts REQUIRED for CV calculation
## Returns: TRUE on success, FALSE on failure
generate_heatmap_with_cv <- function(data_matrix, output_path, title, count_type, label_type, normalization_scheme = "count_type_normalized", sort_by_expression = FALSE, raw_data_matrix = NULL) {
  # Basic checks
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save the matrix for traceability
  save_matrix_data(data_matrix, output_path)
  
  # Prepare column labels
  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }
  colnames(data_matrix) <- col_labels
  if (!is.null(raw_data_matrix)) {
    colnames(raw_data_matrix) <- col_labels
  }
  
  # Optional sorting of columns by mean expression
  if (sort_by_expression) {
    col_means <- colMeans(data_matrix, na.rm = TRUE)
    col_order <- order(col_means)
    data_matrix <- data_matrix[, col_order, drop = FALSE]
    col_labels <- col_labels[col_order]
    colnames(data_matrix) <- col_labels
    if (!is.null(raw_data_matrix)) {
      raw_data_matrix <- raw_data_matrix[, col_order, drop = FALSE]
      colnames(raw_data_matrix) <- col_labels
    }
  }
  
  # CV MUST be computed on raw, untransformed counts
  if (is.null(raw_data_matrix)) {
    stop("ERROR: Cannot calculate CV without raw data.\n",
         "       raw_data_matrix parameter is required for all normalizations.\n",
         "       CV must be calculated on untransformed counts for valid interpretation.\n",
         "       Current normalization scheme: ", normalization_scheme)
  }
  # Compute CV per gene: sd/mean (NA if mean 0)
  cv_values <- apply(raw_data_matrix, 1, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    if (mean_val == 0 || is.na(mean_val) || is.na(sd_val)) {
      return(NA_real_)
    }
    return(sd_val / mean_val)
  })
  
  # Prepare data for visualization: avoid double-normalization
  if (normalization_scheme %in% c("zscore", "zscore_scaled_to_ten")) {
    # Already z-scored or scaled; use as-is
    data_scaled <- data_matrix
  } else {
    # For other schemes, z-score (per-gene) to make heatmap comparable across genes
    data_scaled <- t(scale(t(data_matrix)))
    data_scaled[is.na(data_scaled)] <- 0
  }
  
  # Make CV dataframe for output
  cv_df <- data.frame(
    Gene = rownames(data_matrix),
    CV = round(cv_values, 3),
    stringsAsFactors = FALSE
  )
  
  # Ensure row names exist and truncate for display
  if (is.null(rownames(data_matrix)) || length(rownames(data_matrix)) == 0) {
    rownames(data_matrix) <- paste0("Gene_", seq_len(nrow(data_matrix)))
  }
  row_labels <- truncate_labels(rownames(data_matrix), max_length = 15)
  rownames(data_scaled) <- row_labels
  
  # Choose color mapping depending on normalization scheme
  if (normalization_scheme == "zscore") {
    max_abs_zscore <- max(abs(data_scaled), na.rm = TRUE)
    if (is.finite(max_abs_zscore) && max_abs_zscore > 0) {
      zscore_limit <- max(2, ceiling(max_abs_zscore))
      heatmap_colors <- colorRamp2(
        seq(-zscore_limit, zscore_limit, length = 10),
        c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
      )
    } else {
      heatmap_colors <- colorRamp2(c(0, 1), c("#FFFFFF", "#4A148C"))
    }
  } else if (normalization_scheme == "zscore_scaled_to_ten") {
    heatmap_colors <- colorRamp2(
      seq(0, 10, length = 10),
      c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
    )
  } else {
    data_range <- range(data_scaled, na.rm = TRUE, finite = TRUE)
    if (length(data_range) == 2 && data_range[1] != data_range[2]) {
      heatmap_colors <- colorRamp2(
        seq(data_range[1], data_range[2], length = 10),
        c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
      )
    } else {
      heatmap_colors <- colorRamp2(c(0, 1), c("#FFFFFF", "#4A148C"))
    }
  }
  
  tryCatch({
    # Create a textual CV annotation on the right side of the heatmap
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
    
    # Legend title mapping
    legend_title <- switch(normalization_scheme,
      "raw" = "Raw Counts (Z-scored for viz)",
      "count_type_normalized" = "Log2 Expression (Z-scored for viz)",
      "zscore" = "Z-score",
      "zscore_scaled_to_ten" = "ZScore Scaled to [0-10]",
      "cpm" = "Log2(CPM+1) (Z-scored for viz)",
      "Expression"
    )
    
    # Legend ticks and layout
    legend_at <- if (normalization_scheme == "zscore_scaled_to_ten") {
      seq(0, 10, by = 2)
    } else if (normalization_scheme == "zscore") {
      max_abs_zscore <- max(abs(data_scaled), na.rm = TRUE)
      if (is.finite(max_abs_zscore) && max_abs_zscore > 0) {
        zscore_limit <- max(2, ceiling(max_abs_zscore))
        seq(-zscore_limit, zscore_limit, length.out = 10)
      } else {
        NULL
      }
    } else {
      NULL
    }
    
    if (exists("LEGEND_POSITION") && LEGEND_POSITION == "bottom") {
      legend_direction <- "horizontal"
      legend_height <- unit(1.5, "cm")
      legend_width <- unit(6, "cm")
      grid_height <- unit(0.6, "cm")
      grid_width <- unit(1.2, "cm")
    } else {
      legend_direction <- "vertical"
      legend_height <- unit(8, "cm")
      legend_width <- unit(2, "cm")
      grid_height <- unit(1.8, "cm")
      grid_width <- unit(0.8, "cm")
    }
    
    # Construct heatmap object
    ht <- Heatmap(
      data_scaled,
      name = legend_title,
      col = heatmap_colors,
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 9),
      row_names_max_width = unit(6, "cm"),
      cluster_rows = FALSE,
      show_column_names = TRUE,
      column_names_side = "bottom",
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      cluster_columns = FALSE,
      rect_gp = gpar(col = "white", lwd = 0.3),
      heatmap_legend_param = list(
        title = legend_title,
        legend_direction = legend_direction,
        legend_height = legend_height,
        legend_width = legend_width,
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = grid_height,
        grid_width = grid_width,
        at = legend_at
      ),
      column_title = paste(title, "- CV Heatmap (", normalization_scheme, ifelse(sort_by_expression, " - Sorted by Expression", " - Sorted by Organ"), ")"),
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      right_annotation = cv_annotation
    )
    
    # Compute dynamic image size and render
    n_rows <- nrow(data_scaled)
    n_cols <- ncol(data_scaled)
    base_cell_size <- 0.5
    min_width <- 8
    min_height <- 6
    max_width <- 24
    max_height <- 20
    label_width <- 4
    label_height <- 3
    legend_space <- if (exists("LEGEND_POSITION") && LEGEND_POSITION == "bottom") 2 else 3
    final_width <- max(min_width, min(max_width, n_cols * base_cell_size + label_width + legend_space))
    final_height <- max(min_height, min(max_height, n_rows * base_cell_size + label_height + legend_space))
    
    png(output_path, width = final_width, height = final_height, units = "in", res = 300)
    draw(ht, heatmap_legend_side = ifelse(exists("LEGEND_POSITION"), LEGEND_POSITION, "right"), gap = unit(0.5, "cm"))
    dev.off()
    
    # Save CV table alongside the image for downstream use
    cv_output_path <- gsub("\\.png$", "_CV_data.tsv", output_path)
    write.table(cv_df, file = cv_output_path, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap with CV for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# BAR GRAPH GENERATION
# ===============================================

# Function to generate individual bar graphs for each gene
# Inputs:
#  - data_matrix: genes x samples numeric matrix
#  - output_dir: directory to save PNGs
#  - title_prefix: prefix applied to plot titles/files
#  - count_type, label_type, normalization_scheme: used for labels
#  - sort_by_expression: if TRUE, bars are ordered by expression
# Returns: number of successfully saved plots
generate_gene_bargraphs <- function(data_matrix, output_dir, title_prefix, count_type, label_type, normalization_scheme = "count_type_normalized", sort_by_expression = FALSE) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(0)
  if (ncol(data_matrix) == 0) return(0)
  
  # Prepare column labels
  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }
    
  successful_plots <- 0
  
  # Loop through genes and create a bar graph per gene
  for (i in seq_len(nrow(data_matrix))) {
    gene_name <- rownames(data_matrix)[i]
    gene_values <- as.numeric(data_matrix[i, ])
    
    # Skip genes with all zeros
    if (all(gene_values == 0)) next
    
    # Create data frame for ggplot
    plot_data <- data.frame(
      Stage = col_labels,
      Expression = gene_values,
      stringsAsFactors = FALSE
    )
    
    # Control ordering of bars
    if (sort_by_expression) {
      plot_data <- plot_data[order(plot_data$Expression), ]
      plot_data$Stage <- factor(plot_data$Stage, levels = plot_data$Stage)
    } else {
      plot_data$Stage <- factor(plot_data$Stage, levels = col_labels)
    }
    
    # Y-axis label mapping by normalization scheme
    y_label <- switch(normalization_scheme,
      "raw" = "Raw counts",
      "raw_normalized" = "CPM (not log-transformed)", 
      "count_type_normalized" = ifelse(count_type == "coverage", "Log2(CPM + 1)", 
                        ifelse(count_type %in% c("fpkm", "tpm"), "Log2(Expression + 1)", "Log2(Expression + 1)")),
      "zscore" = "Z-score",
      "zscore_scaled_to_ten" = "Z-Score Scaled to [0-10]",
      "cpm" = "Log2(CPM + 1)",
      "Normalized Expression"
    )
    
    # Color ramp for fill aesthetic
    custom_colors <- c("#eeecd6ff", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
    
    # Build the ggplot bar graph
    p <- ggplot(plot_data, aes(x = .data$Stage, y = .data$Expression, fill = .data$Expression)) +
      geom_bar(stat = "identity", color = "black", size = 0.3) +
      labs(title = gene_name, x = "", y = y_label) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      scale_fill_gradientn(colors = custom_colors, name = "Expression") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    
    # Save the plot using a safe filename
    clean_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene_name)
    output_file <- file.path(output_dir, paste0(clean_gene_name, "_", normalization_scheme, "_bargraph.png"))
    
    tryCatch({
      ggsave(output_file, plot = p, width = 8, height = 6, units = "in", dpi = 300, bg = "white")
      successful_plots <- successful_plots + 1
    }, error = function(e) {
      cat("Error saving bar graph for gene:", gene_name, "\n")
      cat("Error details:", e$message, "\n")
    })
  }
  
  return(successful_plots)
}
