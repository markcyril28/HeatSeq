#!/usr/bin/env Rscript

# ===============================================
# UTILITY FUNCTIONS FOR POST-PROCESSING ANALYSES
# ===============================================
# Comprehensive utility functions for all analysis modules

source(file.path(dirname(sys.frame(1)$ofile), "0_shared_config.R"))

# ===============================================
# GPU-ACCELERATED COMPUTATION FUNCTIONS
# ===============================================
# These functions use GPU when available, with automatic CPU fallback

# GPU-accelerated correlation matrix computation
# WGCNA cor() is the bottleneck for large datasets
gpu_cor <- function(x, method = "pearson") {
  if (!GPU_AVAILABLE || nrow(x) < 100) {
    # Use CPU for small matrices or when GPU unavailable
    return(cor(x, method = method, use = "pairwise.complete.obs"))
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      # Use torch for GPU correlation
      x_tensor <- torch::torch_tensor(as.matrix(x), device = "cuda")
      # Center the data
      x_centered <- x_tensor - x_tensor$mean(dim = 1, keepdim = TRUE)
      # Compute covariance
      cov_matrix <- torch::torch_mm(x_centered, x_centered$t()) / (x_tensor$size(2) - 1)
      # Compute standard deviations
      std_dev <- torch::torch_sqrt(torch::torch_diag(cov_matrix))
      # Compute correlation
      cor_matrix <- cov_matrix / torch::torch_outer(std_dev, std_dev)
      result <- as.matrix(cor_matrix$cpu())
      rownames(result) <- rownames(x)
      colnames(result) <- rownames(x)
      return(result)
    } else if (GPU_BACKEND == "gpuR") {
      # Use gpuR for GPU correlation
      x_gpu <- gpuR::vclMatrix(as.matrix(x), type = "float")
      result <- as.matrix(gpuR::cov(x_gpu))
      # Convert covariance to correlation
      std_dev <- sqrt(diag(result))
      result <- result / outer(std_dev, std_dev)
      rownames(result) <- rownames(x)
      colnames(result) <- rownames(x)
      return(result)
    }
  }, error = function(e) {
    message("[GPU] Correlation failed, falling back to CPU: ", e$message)
  })
  
  # Fallback to CPU
  return(cor(x, method = method, use = "pairwise.complete.obs"))
}

# GPU-accelerated PCA
gpu_prcomp <- function(x, center = TRUE, scale. = FALSE, rank. = NULL) {
  if (!GPU_AVAILABLE || nrow(x) < 50) {
    return(prcomp(x, center = center, scale. = scale., rank. = rank.))
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      x_mat <- as.matrix(x)
      
      # Center and scale on CPU first (small operation)
      if (center) {
        col_means <- colMeans(x_mat, na.rm = TRUE)
        x_mat <- sweep(x_mat, 2, col_means, "-")
      }
      if (scale.) {
        col_sds <- apply(x_mat, 2, sd, na.rm = TRUE)
        col_sds[col_sds == 0] <- 1
        x_mat <- sweep(x_mat, 2, col_sds, "/")
      }
      
      # SVD on GPU
      x_tensor <- torch::torch_tensor(x_mat, device = "cuda", dtype = torch::torch_float32())
      svd_result <- torch::linalg_svd(x_tensor, full_matrices = FALSE)
      
      # Extract results
      n <- nrow(x_mat)
      sdev <- as.numeric(svd_result[[2]]$cpu()) / sqrt(max(1, n - 1))
      rotation <- as.matrix(svd_result[[3]]$t()$cpu())
      x_scores <- as.matrix(torch::torch_mm(x_tensor, svd_result[[3]]$t()$t())$cpu())
      
      # Build prcomp-compatible result
      result <- list(
        sdev = sdev,
        rotation = rotation,
        x = x_scores,
        center = if (center) col_means else FALSE,
        scale = if (scale.) col_sds else FALSE
      )
      class(result) <- "prcomp"
      rownames(result$x) <- rownames(x)
      colnames(result$x) <- paste0("PC", seq_len(ncol(result$x)))
      colnames(result$rotation) <- paste0("PC", seq_len(ncol(result$rotation)))
      rownames(result$rotation) <- colnames(x)
      
      return(result)
    }
  }, error = function(e) {
    message("[GPU] PCA failed, falling back to CPU: ", e$message)
  })
  
  # Fallback to CPU
  return(prcomp(x, center = center, scale. = scale., rank. = rank.))
}

# GPU-accelerated matrix multiplication (for large TOM calculations in WGCNA)
gpu_matmult <- function(A, B) {
  if (!GPU_AVAILABLE || nrow(A) < 100) {
    return(A %*% B)
  }
  
  tryCatch({
    if (GPU_BACKEND == "torch") {
      A_tensor <- torch::torch_tensor(as.matrix(A), device = "cuda", dtype = torch::torch_float32())
      B_tensor <- torch::torch_tensor(as.matrix(B), device = "cuda", dtype = torch::torch_float32())
      result <- as.matrix(torch::torch_mm(A_tensor, B_tensor)$cpu())
      rownames(result) <- rownames(A)
      colnames(result) <- colnames(B)
      return(result)
    } else if (GPU_BACKEND == "gpuR") {
      A_gpu <- gpuR::vclMatrix(as.matrix(A), type = "float")
      B_gpu <- gpuR::vclMatrix(as.matrix(B), type = "float")
      result <- as.matrix(A_gpu %*% B_gpu)
      rownames(result) <- rownames(A)
      colnames(result) <- colnames(B)
      return(result)
    }
  }, error = function(e) {
    message("[GPU] Matrix multiplication failed, falling back to CPU: ", e$message)
  })
  
  return(A %*% B)
}

# ===============================================
# FILE I/O FUNCTIONS
# ===============================================

read_count_matrix <- function(file_path) {
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE,
                      na.strings = c("", " ", "NA", "null"))
    if (any(duplicated(data[, 1]))) {
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]
    for (i in seq_len(ncol(data))) {
      if (!is.numeric(data[, i])) {
        data[, i] <- suppressWarnings(as.numeric(as.character(data[, i])))
      }
    }
    data_matrix <- as.matrix(data)
    data_matrix[is.na(data_matrix)] <- 0
    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "-", e$message, "\n")
    return(NULL)
  })
}

save_matrix_data <- function(data_matrix, output_path, metadata = NULL) {
  tryCatch({
    base_path <- tools::file_path_sans_ext(output_path)
    matrix_df <- data.frame(
      Gene_ID = rownames(data_matrix),
      data_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    write.table(matrix_df, file = paste0(base_path, ".tsv"), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    return(TRUE)
  }, error = function(e) {
    cat("Error saving matrix:", e$message, "\n")
    return(FALSE)
  })
}

# Read gene list from file (supports CSV and plain text formats)
# CSV format expected: Gene_ID,Shortened_Name,... (uses Gene_ID column)
# Plain text format: one gene ID per line
read_gene_list_from_file <- function(gene_list_file) {
  if (!file.exists(gene_list_file)) return(NULL)
  
  file_ext <- tolower(tools::file_ext(gene_list_file))
  
  tryCatch({
    if (file_ext == "csv") {
      gene_df <- read.csv(gene_list_file, stringsAsFactors = FALSE, header = TRUE)
      if ("Gene_ID" %in% colnames(gene_df)) {
        gene_list <- gene_df$Gene_ID
      } else {
        gene_list <- gene_df[[1]]  # Fallback to first column
      }
    } else {
      gene_list <- suppressWarnings(readLines(gene_list_file))
      gene_list <- gene_list[!grepl("^#|^Gene_ID", gene_list, ignore.case = TRUE) & nzchar(gene_list)]
    }
    return(trimws(gene_list))
  }, error = function(e) {
    cat("Error reading gene list:", e$message, "\n")
    return(NULL)
  })
}

truncate_labels <- function(labels, max_length = 25) {
  sapply(as.character(labels), function(x) {
    if (is.na(x) || is.null(x)) return("")
    if (nchar(x) > max_length) paste0(substr(x, 1, max_length - 3), "...")
    else x
  }, USE.NAMES = FALSE)
}

# ===============================================
# NORMALIZATION FUNCTIONS
# ===============================================

preprocess_for_raw <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  return(data_matrix)
}

preprocess_for_cpm <- function(data_matrix, count_type = "expected_count") {
  # Counts Per Million (CPM): library-size normalization
  # Formula: CPM = (count / total_library_count) * 1e6
  # Note: CPM accounts for sequencing depth differences between samples
  # Log2(CPM + 1) applied for variance stabilization (pseudocount of 1)
  # Suitable for visualization; for DEA, use DESeq2's internal normalization
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  lib_sizes <- colSums(data_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- 1  # Avoid division by zero
  data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
  data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
  return(log2(data_cpm + 1))  # +1 pseudocount before log
}

preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- log2(data_matrix + 1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  return(data_processed)
}

preprocess_for_zscore <- function(data_matrix, count_type) {
  # Z-score normalization: row-wise (per gene) standardization
  # Each gene scaled to mean=0, sd=1 across samples
  # Appropriate for comparing relative expression patterns
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_norm <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
    # Row-wise z-score (per gene normalization)
    row_means <- rowMeans(data_norm, na.rm = TRUE)
    row_sds <- apply(data_norm, 1, sd, na.rm = TRUE)
    row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1  # Avoid division by zero
    data_zscore <- sweep(data_norm, 1, row_means, "-")
    data_zscore <- sweep(data_zscore, 1, row_sds, "/")
    data_zscore[is.na(data_zscore) | is.infinite(data_zscore)] <- 0
    return(data_zscore)
  }
  return(data_norm)
}

preprocess_for_zscore_scaled_to_ten <- function(data_matrix, count_type) {
  data_processed <- preprocess_for_zscore(data_matrix, count_type)
  if (is.null(data_processed)) return(NULL)
  min_val <- min(data_processed, na.rm = TRUE)
  max_val <- max(data_processed, na.rm = TRUE)
  val_range <- max_val - min_val
  if (val_range > 0 && is.finite(val_range)) {
    data_scaled <- ((data_processed - min_val) / val_range) * 10
  } else {
    data_scaled <- matrix(0, nrow = nrow(data_processed), ncol = ncol(data_processed))
  }
  data_scaled[is.na(data_scaled) | is.infinite(data_scaled)] <- 0
  return(data_scaled)
}

preprocess_for_deseq2_normalized <- function(data_matrix, count_type = "expected_count") {
  # DESeq2-style median-of-ratios normalization
  # Computes size factors based on geometric mean of each gene across samples
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # Replace zeros/NAs with small value for geometric mean calculation
  data_clean <- data_matrix
  data_clean[data_clean == 0 | is.na(data_clean)] <- 0.5
  
  # Calculate geometric mean per gene (row)
  log_data <- log(data_clean)
  geo_means <- exp(rowMeans(log_data, na.rm = TRUE))
  
  # Remove genes with zero geometric mean
  valid_genes <- geo_means > 0 & is.finite(geo_means)
  if (sum(valid_genes) == 0) {
    # Fallback to simple log2 if no valid genes
    return(preprocess_for_count_type_normalized(data_matrix, count_type))
  }
  
  # Calculate size factors (median of ratios for each sample)
  ratios <- sweep(data_clean[valid_genes, , drop = FALSE], 1, geo_means[valid_genes], FUN = "/")
  size_factors <- apply(ratios, 2, median, na.rm = TRUE)
  size_factors[size_factors == 0 | !is.finite(size_factors)] <- 1
  
  # Normalize counts by size factors
  data_normalized <- sweep(data_matrix, 2, size_factors, FUN = "/")
  
  # Apply log2 transformation
  data_normalized <- log2(data_normalized + 1)
  data_normalized[is.na(data_normalized) | is.infinite(data_normalized)] <- 0
  
  return(data_normalized)
}

apply_normalization <- function(data_matrix, normalization_scheme, count_type) {
  switch(normalization_scheme,
    "raw" = preprocess_for_raw(data_matrix),
    "count_type_normalized" = preprocess_for_count_type_normalized(data_matrix, count_type),
    "deseq2_normalized" = preprocess_for_deseq2_normalized(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "zscore_scaled_to_ten" = preprocess_for_zscore_scaled_to_ten(data_matrix, count_type),
    "cpm" = preprocess_for_cpm(data_matrix, count_type),
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}

# ===============================================
# HEATMAP LEGEND HELPERS
# ===============================================

get_legend_layout <- function(position = "bottom") {
  if (position == "bottom") {
    list(direction = "horizontal", height = unit(2, "cm"), width = unit(10, "cm"),
         grid_height = unit(0.8, "cm"), grid_width = unit(2.5, "cm"))
  } else {
    list(direction = "vertical", height = unit(10, "cm"), width = unit(2, "cm"),
         grid_height = unit(2.5, "cm"), grid_width = unit(0.8, "cm"))
  }
}

get_legend_title <- function(normalization_type, count_type = NULL) {
  switch(normalization_type,
    "raw" = if (!is.null(count_type)) paste0("Raw ", toupper(count_type)) else "Raw Counts",
    "count_type_normalized" = "Log2 Expression",
    "deseq2_normalized" = "DESeq2 Norm. Expression",
    "zscore" = "Z-score",
    "zscore_scaled_to_ten" = "Z-Score [0-10]",
    "cpm" = "Log2(CPM+1)",
    "Expression"
  )
}

get_norm_display_name <- function(norm_scheme) {
  switch(norm_scheme,
    "raw" = "Raw",
    "count_type_normalized" = "Count-Type_Normalized",
    "deseq2_normalized" = "DESeq2_Normalized",
    "zscore" = "Z-score_Normalized",
    "zscore_scaled_to_ten" = "Z-score_Scaled_to_Ten",
    "cpm" = "CPM_Normalized",
    norm_scheme
  )
}

# ===============================================
# SAMPLE LABEL CONVERSION
# ===============================================

convert_to_organ_labels <- function(counts_matrix) {
  colnames_organ <- SAMPLE_LABELS[colnames(counts_matrix)]
  colnames_organ[is.na(colnames_organ)] <- colnames(counts_matrix)[is.na(colnames_organ)]
  result <- counts_matrix
  colnames(result) <- colnames_organ
  result
}

# ===============================================
# GENE NAME CONVERSION
# ===============================================

# Load gene name mapping from gene_groups CSV files
# CSV format: Gene_ID,Shortened_Name,...
load_gene_name_mapping <- function(gene_group, gene_groups_dir = GENE_GROUPS_DIR) {
  csv_file <- file.path(gene_groups_dir, paste0(gene_group, ".csv"))
  if (!file.exists(csv_file)) return(NULL)
  
  tryCatch({
    df <- read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE)
    if ("Gene_ID" %in% colnames(df) && "Shortened_Name" %in% colnames(df)) {
      return(setNames(df$Shortened_Name, df$Gene_ID))
    }
    return(NULL)
  }, error = function(e) NULL)
}

# Convert Gene_ID to Shortened_Name in row names
convert_to_shortened_names <- function(counts_matrix, gene_group) {
  mapping <- load_gene_name_mapping(gene_group)
  if (is.null(mapping)) return(counts_matrix)
  
  new_rownames <- mapping[rownames(counts_matrix)]
  # Keep original name if no mapping found
  new_rownames[is.na(new_rownames)] <- rownames(counts_matrix)[is.na(new_rownames)]
  result <- counts_matrix
  rownames(result) <- new_rownames
  result
}

# Apply label transformations based on gene_type and label_type
apply_labels <- function(counts_matrix, gene_group, gene_type, label_type) {
  result <- counts_matrix
  
  # Apply row label transformation (gene names)
  if (gene_type == "Shortened_Name") {
    result <- convert_to_shortened_names(result, gene_group)
  }
  # gene_type == "Gene_ID" means keep original row names
  
  # Apply column label transformation (sample names)
  if (label_type == "Organ") {
    result <- convert_to_organ_labels(result)
  }
  # label_type == "SRR_ID" means keep original column names
  
  return(result)
}

# ===============================================
# TISSUE GROUP MAPPING (DRY - used by PCA, Correlation, DEA)
# ===============================================
# Maps organ names from CSV to biological tissue groups
# Update this function when CSV organ names change

map_tissue_to_group <- function(tissue_name) {
  if (grepl("Root|Stem|Leaf|Leaves|Senescent", tissue_name, ignore.case = TRUE)) {
    return("Vegetative")
  }
  if (grepl("Flower|Bud|Pistil", tissue_name, ignore.case = TRUE)) {
    return("Reproductive")
  }
  if (grepl("Fruit|peduncle", tissue_name, ignore.case = TRUE)) {
    return("Fruit")
  }
  if (grepl("Radicle|Cotyledon", tissue_name, ignore.case = TRUE)) {
    return("Seedling")
  }
  return("Other")
}

# Vectorized version for sapply
get_tissue_groups <- function(tissue_names) {
  sapply(tissue_names, map_tissue_to_group, USE.NAMES = FALSE)
}

# ===============================================
# COLOR FUNCTIONS
# ===============================================

get_violet_color_scale <- function(n_breaks = 100) {
  colorRampPalette(c("#F5F5F5", "#E6D9F2", "#C9A6E5", "#9B59B6", "#6C3483", "#4A235A"))(n_breaks)
}

get_blue_red_color_scale <- function(n_breaks = 100) {
  colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(n_breaks)
}
