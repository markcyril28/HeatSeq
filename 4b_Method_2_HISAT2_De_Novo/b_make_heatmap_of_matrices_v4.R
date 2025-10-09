# ===============================================
# HEATMAP GENERATION FOR METHOD2 RESULTS
# ===============================================
#
# Generates heatmaps with ComplexHeatmap package using multiple normalization methods:
# - Raw, Count-Type, Z-score, and CPM normalization
# - Professional layout with genes on left, samples on top
# - Saves both PNG heatmaps and TSV data files
#
# Required packages: ComplexHeatmap, circlize, RColorBrewer, dplyr, tibble, grid
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
HEATMAP_OUT_DIR <- "6_Heatmap_Visualizations"

# File naming configuration
MASTER_REFERENCE <- "All_Smel_Genes"
MASTER_REFERENCE_SUFFIX <- paste0("_from_", MASTER_REFERENCE)

# Gene groups
FASTA_GROUPS <- c(
  # Control Gene Groups
  "Best_Cell_Cycle_Associated_Control_Genes",
  "Best_Control_Genes",
  
  # Individual Gene Groups
  "SmelDMPs",
  "SmelGIFs",
  "SmelGRFs",
  
  # Combined Gene Groups with Control Genes
  "SmelGIF_with_Cell_Cycle_Control_genes",
  "SmelGRF_with_Cell_Cycle_Control_genes",
  "SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"

)

# Analysis configuration
COUNT_TYPES <- c("coverage", "fpkm", "tpm")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Mapping from SRR IDs to organ names for matrix headers
SAMPLE_LABELS <- c(
  # Roots
  "SRR3884675" = "Roots_1",       # PRJNA328564
  "SRR20722229" = "Roots_2",      # SAMN28540077
  "SRR31755282" = "Roots_3",      # SAMN28540068
  
  # Stems
  "SRR3884690" = "Stems_1",       # PRJNA328564
  "SRR20722227" = "Stems_2",      # SAMN28540077
  "SRR20722384" = "Stems_3",      # SAMN28540068
  
  # Leaves
  "SRR3884689" = "Leaves_1",      # PRJNA328564
  "SRR20722230" = "Leaves_2",     # SAMN28540077
  "SRR20722386" = "Leaves_3",     # SAMN28540068
  
  # Opened Buds
  "SRR3884687" = "Opened_Buds_1", # PRJNA328564

  # Buds
  "SRR3884686" = "Buds_1",        # PRJNA328564
  "SRR21010466" = "Buds_2",       # SAMN28540077
  "SRR20722297" = "Buds_3",       # SAMN28540068
  
  # Flowers
  "SRR3884597" = "Flowers_1",     # PRJNA328564
  "SRR20722234" = "Flowers_2",    # SAMN28540077
  "SRR23909863" = "Flowers_3",    # SAMN28540068
  
  # Fruits
  "SRR3884631" = "Fruits_1",      # PRJNA328564
  "SRR2072232" = "Fruits_2",      # SAMN28540077
  "SRR20722387" = "Fruits_3"      # SAMN28540068
)


# Create base output directory if it doesn't exist
dir.create(HEATMAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===============================================
# FUNCTIONS
# ===============================================

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
      cat("  ⚠️  Found duplicate row names - making them unique...\n")
      data[, 1] <- make.unique(as.character(data[, 1]), sep = "_")
    }
    
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]
    
    for (i in seq_len(ncol(data))) {
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

# Function to preprocess data based on normalization type
preprocess_data <- function(data_matrix, count_type, norm_type = "default") {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  # Handle NAs and infinite values
  data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  
  # Apply normalization based on type
  data_processed <- switch(norm_type,
    "raw" = data_matrix,
    "cpm" = {
      lib_sizes <- colSums(data_matrix, na.rm = TRUE)
      lib_sizes[lib_sizes == 0] <- 1
      data_cpm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
      data_cpm[is.na(data_cpm) | is.infinite(data_cpm)] <- 0
      log2(data_cpm + 1)
    },
    "zscore" = {
      # First apply count-type normalization
      if (count_type == "coverage") {
        lib_sizes <- colSums(data_matrix, na.rm = TRUE)
        lib_sizes[lib_sizes == 0] <- 1
        data_norm <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
        data_norm <- log2(data_norm + 1)
      } else if (count_type %in% c("fpkm", "tpm")) {
        data_norm <- log2(data_matrix + 0.1)
      } else {
        data_norm <- log2(data_matrix + 1)
      }
      # Then apply z-score normalization
      data_norm[is.na(data_norm) | is.infinite(data_norm)] <- 0
      if (nrow(data_norm) > 1 && ncol(data_norm) > 1) {
        data_zscore <- t(scale(t(data_norm), center = TRUE, scale = TRUE))
        data_zscore[is.na(data_zscore)] <- 0
        data_zscore
      } else {
        data_norm
      }
    },
    {  # Default case (count-type normalization)
      if (count_type == "coverage") {
        lib_sizes <- colSums(data_matrix, na.rm = TRUE)
        lib_sizes[lib_sizes == 0] <- 1
        data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
        log2(data_normalized + 1)
      } else if (count_type %in% c("fpkm", "tpm")) {
        log2(data_matrix + 0.1)
      } else {
        log2(data_matrix + 1)
      }
    }
  )
  
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  # Filter high-variance genes if too many genes
  if (nrow(data_processed) > 100 && norm_type != "raw") {
    gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
    gene_vars[is.na(gene_vars)] <- 0
    if (length(gene_vars) > 100) {
      top_genes <- min(1000, length(gene_vars))
      keep_genes <- order(gene_vars, decreasing = TRUE)[1:top_genes]
      data_processed <- data_processed[keep_genes, ]
    }
  }
  
  return(data_processed)
}

# Function to generate heatmap (unified for all normalization types)
generate_heatmap <- function(data_matrix, output_path, title, count_type, label_type, normalization_type = "raw") {
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save matrix data as TSV
  save_matrix_data(data_matrix, output_path)
  
  # Prepare column labels
  col_labels <- colnames(data_matrix)
  if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
    colnames(data_matrix) <- col_labels
  }
  
  # Define color palette and legend title
  eggplant_colors <- colorRamp2(
    seq(min(data_matrix), max(data_matrix), length = 8),
    c("#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
      "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
  )
  
  legend_title <- switch(normalization_type,
    "raw" = "Raw Counts",
    "raw_normalized" = "Raw Counts", 
    "Count-Type_Normalized" = "Log2 Expression",
    "Z-score_Normalized" = "Z-score",
    "CPM_Normalized" = "Log2(CPM + 1)",
    "Expression"
  )
  
  tryCatch({
    # Truncate row names for better display
    max_length <- if (nrow(data_matrix) > 50) 20 else 30
    row_labels <- truncate_labels(rownames(data_matrix), max_length = max_length)
    rownames(data_matrix) <- row_labels
    
    # Create the heatmap
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
      
      # Legend settings
      heatmap_legend_param = list(
        title = legend_title,
        legend_direction = "vertical",
        legend_height = unit(8, "cm"),
        legend_width = unit(2, "cm"),
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        grid_height = unit(1.8, "cm"),
        grid_width = unit(0.8, "cm"),
        just = c("left", "top")
      ),
      
      # Title
      column_title = if (normalization_type == "raw") title else paste(title, "-", normalization_type),
      column_title_gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Save the heatmap
    png(output_path, width = if (normalization_type == "raw") 16 else 18, 
        height = if (normalization_type == "raw") 12 else 14, units = "in", res = 300)
    draw(ht, 
         heatmap_legend_side = "right", 
         annotation_legend_side = "right",
         merge_legends = TRUE,
         gap = unit(7, "mm"))
    dev.off()
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
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

for (group in FASTA_GROUPS) {
  # Clear and recreate output directory for this specific group
  group_output_dir <- file.path(HEATMAP_OUT_DIR, group)
  if (dir.exists(group_output_dir)) {
    unlink(group_output_dir, recursive = TRUE)
  }
  
  for (count_type in COUNT_TYPES) {
    for (gene_type in GENE_TYPES) {
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
        
        output_dir <- file.path(HEATMAP_OUT_DIR, group)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          dir.create(file.path(output_dir, "raw"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "raw_normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Count-Type_Normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Z-score_Normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "CPM_Normalized"), showWarnings = FALSE)
          
          raw_data <- read_count_matrix(input_file)
          if (is.null(raw_data)) next
          
          input_basename <- tools::file_path_sans_ext(basename(input_file))
          title <- gsub("_", " ", input_basename)
          
          # Generate all normalization types
          normalization_configs <- list(
            list(type = "raw", data = raw_data, subdir = "raw"),
            list(type = "raw_normalized", data = preprocess_data(raw_data, count_type, "raw"), subdir = "raw_normalized"),
            list(type = "Count-Type_Normalized", data = preprocess_data(raw_data, count_type, "default"), subdir = "Count-Type_Normalized"),
            list(type = "Z-score_Normalized", data = preprocess_data(raw_data, count_type, "zscore"), subdir = "Z-score_Normalized"),
            list(type = "CPM_Normalized", data = preprocess_data(raw_data, count_type, "cpm"), subdir = "CPM_Normalized")
          )
          
          for (config in normalization_configs) {
            if (!is.null(config$data)) {
              output_path <- file.path(output_dir, config$subdir, 
                                     paste0(input_basename, "_", gsub("-", "_", tolower(config$type)), "_heatmap.png"))
              total_heatmaps <- total_heatmaps + 1
              if (generate_heatmap(config$data, output_path, title, count_type, label_type, config$type)) {
                successful_heatmaps <- successful_heatmaps + 1
              }
            }
          }
      }
    }
  }
}

# ===============================================
# SUMMARY REPORT
# ===============================================

# Ensure all graphics devices are closed
while (length(dev.list()) > 0) { dev.off() }

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("HEATMAP GENERATION SUMMARY REPORT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Processing Summary:\n")
cat("  • Total heatmaps attempted:", total_heatmaps, "\n")
cat("  • Successful heatmaps generated:", successful_heatmaps, "\n")
cat("  • Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
cat("  • Success rate:", round((successful_heatmaps/total_heatmaps)*100, 1), "%\n")
cat("\nOutput Details:\n")
cat("  • Output directory:", HEATMAP_OUT_DIR, "\n")
cat("  • Normalization types: Raw, Raw_Normalized, Count-Type, Z-score, CPM\n")
cat("  • File formats: PNG heatmaps + TSV data matrices\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

if (successful_heatmaps > 0) {
  cat("✅ Heatmap generation completed successfully!\n")
  cat("📁 Check the output directory for visualization files.\n")
} else {
  cat("❌ No heatmaps were generated.\n")
  cat("🔍 Please check input files and paths.\n")
}

cat("\nAnalysis includes gene groups:", length(FASTA_GROUPS), "groups\n")
cat("Count types processed:", paste(COUNT_TYPES, collapse = ", "), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Current structure processes:
# - Each combination of (group, count_type, gene_type, label_type) 
#   corresponds to exactly one TSV file
# - Each TSV file generates multiple heatmaps under different normalization types:
#   * raw: Raw count data without any transformation or scaling
#   * raw_normalized: Raw count data without any normalization or scaling
#   * Count-Type_Normalized: Type-specific normalization (coverage->CPM-like, fpkm/tpm->log2)
#   * Z-score_Normalized: Count-type normalized then Z-score standardized
#   * CPM_Normalized: Counts Per Million normalization with log2 transformation
# - Orientation: genes in rows, samples in columns (no transposition)
# - Output files use the input filename as base name
# - TSV matrices are saved alongside PNG heatmaps for each normalization type
