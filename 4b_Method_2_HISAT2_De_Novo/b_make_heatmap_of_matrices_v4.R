# ===============================================
# HEATMAP GENERATION FOR METHOD2 RESULTS
# ===============================================
#
# UPDATES in this version:
# - Replaced pheatmap with ComplexHeatmap for better label control
# - Row labels (genes) now positioned on the LEFT side
# - Column labels (samples) now positioned at the TOP
# - Added 45-degree rotation for column labels for better readability
# - Implemented label truncation for long gene names
# - Enhanced color palette with proper color mapping
# - Improved legend positioning and formatting
# - Better error handling and debugging output
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
MATRICES_DIR <- "5_stringtie_WD/b_Method_2_COUNT_MATRICES"
HEATMAP_OUT_DIR <- "6_Heatmap_Visualizations"

# Gene groups
FASTA_GROUPS <- c(
  #"TEST",
  #"All_Smel_Genes",
  "SmelDMP_CDS_Control_Best"
  #"SmelGIF_with_Best_Control_Cyclo"
  #"SmelGRF_with_Best_Control_Cyclo",
  #"SmelGRF-GIF_with_Best_Control_Cyclo",
  #"SmelGIF_with_Cell_Cycle_Control_genes",
  #"SmelGRF_with_Cell_Cycle_Control_genes"
)

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
    
    cat("Matrix saved to:", matrix_path, "\n")
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

# Function to preprocess data for heatmap
preprocess_for_heatmap <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  data_processed <- log2(data_matrix + 0.1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  row_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  row_vars[is.na(row_vars)] <- 0
  
  # Keep genes with zero variance by adding a small pseudocount for heatmap visualization
  # if (any(row_vars == 0)) {
  #   cat("Removing", sum(row_vars == 0), "genes with zero variance\n")
  #   data_processed <- data_processed[row_vars > 0, , drop = FALSE]
  # }
  # MODIFIED: Now keeping all genes, including those with zero variance
  cat("Keeping all", nrow(data_processed), "genes (including", sum(row_vars == 0), "with zero variance)\n")
  
  if (nrow(data_processed) > 10) {
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

# Preprocessing for raw data
preprocess_for_raw_data <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- data_matrix
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  # Keep all genes, including those with zero variance
  # if (any(gene_vars == 0)) {
  #   data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  # }
  # MODIFIED: Now keeping all genes for complete representation
  
  return(data_processed)
}

# Preprocessing for count-type normalized data
preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  
  if (count_type == "coverage") {
    lib_sizes <- colSums(data_matrix, na.rm = TRUE)
    lib_sizes[lib_sizes == 0] <- 1
    data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
    data_processed <- log2(data_normalized + 1)
  } else if (count_type %in% c("fpkm", "tpm")) {
    data_processed <- log2(data_matrix + 0.1)
  } else {
    data_processed <- log2(data_matrix + 1)
  }
  
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  gene_vars[is.na(gene_vars)] <- 0
  
  # Keep all genes, including those with zero variance
  # if (any(gene_vars == 0)) {
  #   data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  # }
  # MODIFIED: Now keeping all genes for complete representation
  
  return(data_processed)
}

# Preprocessing for Z-score normalized data
preprocess_for_zscore_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (is.null(data_processed) || nrow(data_processed) == 0) return(NULL)
  
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    data_processed[is.na(data_processed)] <- 0
  }
  return(data_processed)
}

# Function to generate heatmap
generate_heatmap <- function(data_matrix, output_path, title, count_type, label_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save matrix data as TSV
  save_matrix_data(data_matrix, output_path)
  
  # Define eggplant color palette
  eggplant_colors <- colorRamp2(
    seq(min(data_matrix), max(data_matrix), length = 8),
    c("#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
      "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
  )
  
  tryCatch({
    # Scale by rows (genes)
    data_scaled <- t(scale(t(data_matrix)))
    data_scaled[is.na(data_scaled)] <- 0
    
    # Truncate row names for better display
    row_labels <- truncate_labels(rownames(data_matrix), max_length = 30)
    rownames(data_scaled) <- row_labels
    
    # Create the heatmap
    ht <- Heatmap(
      data_scaled,
      name = "Expression",
      col = eggplant_colors,
      
      # Row (gene) settings - labels on the LEFT
      show_row_names = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      row_names_max_width = unit(10, "cm"),
      cluster_rows = FALSE,
      
      # Column (sample) settings - labels at the TOP
      show_column_names = TRUE,
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      cluster_columns = FALSE,
      
      # Heatmap body settings
      rect_gp = gpar(col = "white", lwd = 0.5),
      
      # Legend settings - larger and positioned at top right
      heatmap_legend_param = list(
        title = "Z-score",
        legend_direction = "vertical",
        legend_height = unit(8, "cm"),
        legend_width = unit(2.0, "cm"),
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        grid_height = unit(1.8, "cm"),
        grid_width = unit(0.8, "cm"),
        just = c("left", "top")
      ),
      
      # Title
      column_title = title,
      column_title_gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Save the heatmap
    png(output_path, width = 16, height = 12, units = "in", res = 300)
    draw(ht, 
         heatmap_legend_side = "right", 
         annotation_legend_side = "right",
         merge_legends = TRUE,
         gap = unit(7, "mm"))
    dev.off()
    
    return(TRUE)
  }, error = function(e) {
    cat("Row scaling failed, trying without scaling...\n")
    tryCatch({
      # Truncate row names for better display
      row_labels <- truncate_labels(rownames(data_matrix), max_length = 30)
      rownames(data_matrix) <- row_labels
      
      # Create heatmap without scaling
      ht <- Heatmap(
        data_matrix,
        name = "Expression",
        col = eggplant_colors,
        
        # Row (gene) settings - labels on the LEFT
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        row_names_max_width = unit(10, "cm"),
        cluster_rows = FALSE,
        
        # Column (sample) settings - labels at the TOP
        show_column_names = TRUE,
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        cluster_columns = FALSE,
        
        # Heatmap body settings
        rect_gp = gpar(col = "white", lwd = 0.5),
        
        # Legend settings - larger and positioned at top right
        heatmap_legend_param = list(
          title = "Raw Values",
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
        column_title = paste(title, "(No Scaling)"),
        column_title_gp = gpar(fontsize = 14, fontface = "bold")
      )
      
      # Save the heatmap
      png(output_path, width = 16, height = 12, units = "in", res = 300)
      draw(ht, 
           heatmap_legend_side = "right", 
           annotation_legend_side = "right",
           merge_legends = TRUE,
           gap = unit(7, "mm"))
      dev.off()
      
      return(TRUE)
    }, error = function(e2) {
      cat("Error generating heatmap for:", title, "\n")
      cat("Error details:", e2$message, "\n")
      return(FALSE)
    })
  })
}

# Function to generate normalized heatmap
generate_normalized_heatmap <- function(data_matrix, output_path, title, count_type, label_type, normalization_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (nrow(data_matrix) < 2) return(FALSE)
  if (any(is.na(data_matrix)) || any(is.infinite(data_matrix))) {
    data_matrix[is.na(data_matrix) | is.infinite(data_matrix)] <- 0
  }
  if (length(unique(as.vector(data_matrix))) == 1) return(FALSE)
  
  # Save matrix data as TSV
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
  
  # Set column names for the matrix
  colnames(data_matrix) <- col_labels
  
  # Define eggplant color palette
  eggplant_colors <- colorRamp2(
    seq(min(data_matrix), max(data_matrix), length = 8),
    c("#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
      "#8E24AA","#6A1B9A","#4A148C","#2F1B69")
  )
  
  # Determine scaling method
  scale_method <- switch(normalization_type,
                         "raw_normalized" = "row",
                         "Count-Type_Normalized" = "none", 
                         "Z-score_Normalized" = "none",
                         "row")
  
  tryCatch({
    # Apply scaling if needed
    if (scale_method == "row") {
      data_scaled <- t(scale(t(data_matrix)))
      data_scaled[is.na(data_scaled)] <- 0
      legend_title <- "Z-score"
    } else {
      data_scaled <- data_matrix
      legend_title <- "Expression"
    }
    
    # Truncate row names for better display (only if many genes)
    if (nrow(data_matrix) > 50) {
      row_labels <- truncate_labels(rownames(data_scaled), max_length = 20)
    } else {
      row_labels <- truncate_labels(rownames(data_scaled), max_length = 30)
    }
    rownames(data_scaled) <- row_labels
    
    # Create the heatmap
    ht <- Heatmap(
      data_scaled,
      name = legend_title,
      col = eggplant_colors,
      
      # Row (gene) settings - labels on the LEFT
      show_row_names = ifelse(nrow(data_matrix) > 50, FALSE, TRUE),
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 7),
      row_names_max_width = unit(10, "cm"),
      cluster_rows = FALSE,
      
      # Column (sample) settings - labels at the TOP
      show_column_names = TRUE,
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      cluster_columns = FALSE,
      
      # Heatmap body settings
      rect_gp = gpar(col = "white", lwd = 0.3),
      
      # Legend settings - larger and positioned at top right
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
      column_title = paste(title, "-", normalization_type),
      column_title_gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Save the heatmap
    png(output_path, width = 18, height = 14, units = "in", res = 300)
    draw(ht, 
         heatmap_legend_side = "right", 
         annotation_legend_side = "right",
         merge_legends = TRUE,
         gap = unit(7, "mm"))
    dev.off()
    
    return(TRUE)
  }, error = function(e) {
    cat("Error generating normalized heatmap for:", title, "\n")
    cat("Error details:", e$message, "\n")
    return(FALSE)
  })
}

# ===============================================
# MAIN PROCESSING
# ===============================================

while (length(dev.list()) > 0) { dev.off() }

total_heatmaps <- 0
successful_heatmaps <- 0

for (group in FASTA_GROUPS) {
  # Clear and recreate output directory for this specific group
  group_output_dir <- file.path(HEATMAP_OUT_DIR, group)
  if (dir.exists(group_output_dir)) {
    unlink(group_output_dir, recursive = TRUE)
  }
  
  for (count_type in COUNT_TYPES) {
    for (gene_type in GENE_TYPES) {
      for (label_type in LABEL_TYPES) {
        input_file <- file.path(MATRICES_DIR, group, 
                               paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
        if (!file.exists(input_file)) next
        
        output_dir <- file.path(HEATMAP_OUT_DIR, group)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          dir.create(file.path(output_dir, "raw_copy"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "raw_normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Count-Type_Normalized"), showWarnings = FALSE)
          dir.create(file.path(output_dir, "Z-score_Normalized"), showWarnings = FALSE)
          
          input_basename <- tools::file_path_sans_ext(basename(input_file))
          output_file <- file.path(output_dir, paste0(input_basename, "_heatmap.png"))
          title <- gsub("_", " ", input_basename)
          
          raw_data <- read_count_matrix(input_file)
          processed_data <- preprocess_for_heatmap(raw_data, count_type)
          
          total_heatmaps <- total_heatmaps + 1
          if (generate_heatmap(processed_data, output_file, title, count_type, label_type)) {
            successful_heatmaps <- successful_heatmaps + 1
          }
          
          raw_copy_output <- file.path(output_dir, "raw_copy", paste0(input_basename, "_raw_copy_heatmap.png"))
          total_heatmaps <- total_heatmaps + 1
          if (generate_heatmap(processed_data, raw_copy_output, title, count_type, label_type)) {
            successful_heatmaps <- successful_heatmaps + 1
          }
          
          raw_data_processed <- preprocess_for_raw_data(raw_data)
          if (!is.null(raw_data_processed)) {
            raw_output <- file.path(output_dir, "raw_normalized", paste0(input_basename, "_raw_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            if (generate_normalized_heatmap(raw_data_processed, raw_output, title, count_type, label_type, "raw_normalized")) {
              successful_heatmaps <- successful_heatmaps + 1
            }
          }
          
          count_normalized_data <- preprocess_for_count_type_normalized(raw_data, count_type)
          if (!is.null(count_normalized_data)) {
            count_normalized_output <- file.path(output_dir, "Count-Type_Normalized", paste0(input_basename, "_count_normalized_heatmap.png"))
            total_heatmaps <- total_heatmaps + 1
            if (generate_normalized_heatmap(count_normalized_data, count_normalized_output, title, count_type, label_type, "Count-Type_Normalized")) {
              successful_heatmaps <- successful_heatmaps + 1
            }
          }
          
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

# ===============================================
# SUMMARY
# ===============================================

while (length(dev.list()) > 0) { dev.off() }

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

# Current structure processes:
# - Each combination of (group, count_type, gene_type, label_type) 
#   corresponds to exactly one TSV file
# - Each TSV file generates multiple heatmaps under different normalization types
# - Orientation: genes in rows, samples in columns (no transposition)
# - Output files use the input filename as base name
