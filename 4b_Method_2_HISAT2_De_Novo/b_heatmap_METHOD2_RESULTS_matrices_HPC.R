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
BASE_DIR <- "/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/4b_Method_2_HISAT2_De_Novo"
MATRICES_DIR <- file.path(BASE_DIR, "5_stringtie/a_Method_2_Results_matrices_post-processed")
HEATMAP_OUT_DIR <- file.path(BASE_DIR, "6_Heatmap_Visualization/a_Method_2_Results_matriced_heatmapped")

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
    
    row_sums <- rowSums(data_matrix, na.rm = TRUE)
    data_matrix <- data_matrix[row_sums > 0, , drop = FALSE]
    
    return(data_matrix)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
}

# Preprocessing functions
preprocess_for_heatmap <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- log2(data_matrix + 0.1)
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  row_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  row_vars[is.na(row_vars)] <- 0
  data_processed <- data_processed[row_vars > 0, , drop = FALSE]
  if (nrow(data_processed) > 1000) {
    gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
    keep_genes <- order(gene_vars, decreasing = TRUE)[1:1000]
    data_processed <- data_processed[keep_genes, ]
  }
  return(data_processed)
}

preprocess_for_raw_data <- function(data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  data_processed <- data_matrix
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  return(data_processed)
}

preprocess_for_count_type_normalized <- function(data_matrix, count_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(NULL)
  if (count_type == "coverage") {
    lib_sizes <- colSums(data_matrix, na.rm = TRUE)
    lib_sizes[lib_sizes == 0] <- 1
    data_normalized <- sweep(data_matrix, 2, lib_sizes/1e6, FUN = "/")
    data_processed <- log2(data_normalized + 1)
  } else {
    data_processed <- log2(data_matrix + 0.1)
  }
  data_processed[is.na(data_processed) | is.infinite(data_processed)] <- 0
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  data_processed <- data_processed[gene_vars > 0, , drop = FALSE]
  return(data_processed)
}

preprocess_for_zscore_normalized <- function(data_matrix, count_type) {
  data_processed <- preprocess_for_count_type_normalized(data_matrix, count_type)
  if (is.null(data_processed) || nrow(data_processed) == 0) return(NULL)
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    data_processed[is.na(data_processed)] <- 0
  }
  return(data_processed)
}

# Function to generate heatmap
generate_heatmap <- function(data_matrix, output_path, title, col_labels) {
  if (is.null(data_matrix) || nrow(data_matrix) < 2) return(FALSE)
  
  eggplant_colors <- colorRampPalette(c(
    "#FFFFFF","#F3E5F5","#CE93D8","#AB47BC",
    "#8E24AA","#6A1B9A","#4A148C","#2F1B69"
  ))(100)
  
  tryCatch({
    png(output_path, width = 16, height = 12, units = "in", res = 300)
    pheatmap(
      data_matrix,
      color = eggplant_colors,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      labels_col = col_labels,
      legend = TRUE,
      main = title,
      fontsize = 12,
      fontsize_row = 10,
      fontsize_col = 10,
      angle_col = 135,   # 🔴 Rotate column labels 270°
      border_color = NA
    )
    dev.off()
    return(TRUE)
  }, error = function(e) {
    cat("Error generating heatmap for:", title, "\n")
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
  for (version in VERSIONS) {
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        for (label_type in LABEL_TYPES) {
          
          input_file <- file.path(MATRICES_DIR, group, version, 
                                 paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
          if (!file.exists(input_file)) next
          
          output_dir <- file.path(HEATMAP_OUT_DIR, group, version)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          
          input_basename <- tools::file_path_sans_ext(basename(input_file))
          output_file <- file.path(output_dir, paste0(input_basename, "_heatmap.png"))
          title <- gsub("_", " ", input_basename)
          
          raw_data <- read_count_matrix(input_file)
          processed_data <- preprocess_for_heatmap(raw_data, count_type)
          
          # Apply sample labels if SRR
          col_labels <- colnames(processed_data)
          if (label_type == "SRR" && exists("SAMPLE_LABELS")) {
            col_labels <- ifelse(colnames(processed_data) %in% names(SAMPLE_LABELS),
                                 SAMPLE_LABELS[colnames(processed_data)],
                                 colnames(processed_data))
          }
          
          total_heatmaps <- total_heatmaps + 1
          if (generate_heatmap(processed_data, output_file, title, col_labels)) {
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
