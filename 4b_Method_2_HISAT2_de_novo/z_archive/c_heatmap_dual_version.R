# ===============================================
# DUAL VERSION HEATMAP GENERATION
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
HEATMAP_OUT_DIR <- "6_Visualization"

# Create output directory
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
MATRIX_TYPES <- c("geneID", "geneName")

# Sample metadata for better labeling
SAMPLE_LABELS <- c(
  "SRR3884631" = "Fruits_6cm",
  "SRR3884677" = "Cotyledons", 
  "SRR3884679" = "Pistils",
  "SRR3884597" = "Flowers",
  "SRR3884686" = "Buds_0.7cm",
  "SRR3884687" = "Buds_Opened",
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
                      stringsAsFactors = FALSE, row.names = 1)
    
    # Convert to numeric matrix
    data_matrix <- as.matrix(data)
    
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
    data_processed <- log2(data_matrix + 1)
  } else {
    data_processed <- log2(data_matrix + 1)
  }
  
  # Filter genes with low variance (optional)
  gene_vars <- apply(data_processed, 1, var, na.rm = TRUE)
  if (length(gene_vars) > 100) {
    top_genes <- min(1000, length(gene_vars))
    keep_genes <- order(gene_vars, decreasing = TRUE)[1:top_genes]
    data_processed <- data_processed[keep_genes, ]
  }
  
  return(data_processed)
}

# Function to create sample annotation
create_sample_annotation <- function(sample_names, use_organ_names = FALSE) {
  if (use_organ_names) {
    # Sample names are already organ names
    tissue_types <- sapply(sample_names, function(x) {
      if (grepl("Fruits|Buds", x)) return("Reproductive")
      if (grepl("Flowers|Pistils", x)) return("Floral")  
      if (grepl("Leaves|Stems", x)) return("Vegetative")
      if (grepl("Roots|Radicles", x)) return("Root")
      if (grepl("Cotyledons", x)) return("Seed")
      return("Other")
    })
  } else {
    # Sample names are SRR IDs, convert to organ names first
    tissue_types <- sapply(sample_names, function(x) {
      if (x %in% names(SAMPLE_LABELS)) {
        label <- SAMPLE_LABELS[x]
      } else if (x %in% SAMPLE_LABELS) {
        label <- x
      } else {
        label <- x
      }
      
      if (grepl("Fruits|Buds", label)) return("Reproductive")
      if (grepl("Flowers|Pistils", label)) return("Floral")  
      if (grepl("Leaves|Stems", label)) return("Vegetative")
      if (grepl("Roots|Radicles", label)) return("Root")
      if (grepl("Cotyledons", label)) return("Seed")
      return("Other")
    })
  }
  
  annotation_df <- data.frame(
    Tissue_Type = tissue_types,
    row.names = sample_names
  )
  
  return(annotation_df)
}

# Function to generate heatmap
generate_dual_heatmaps <- function(data_matrix, output_path_srr, output_path_organ, title_base, count_type, matrix_type) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) {
    cat("Skipping heatmap - no data for:", title_base, "\n")
    return(c(FALSE, FALSE))
  }
  
  # Define colors
  if (count_type %in% c("fpkm", "tpm")) {
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
  } else {
    colors <- colorRampPalette(c("white", "red"))(100)
  }
  
  # Annotation colors
  ann_colors <- list(
    Tissue_Type = c(
      "Reproductive" = "#FF6B6B",
      "Floral" = "#4ECDC4", 
      "Vegetative" = "#45B7D1",
      "Root" = "#96CEB4",
      "Seed" = "#FFEAA7",
      "Other" = "#DDA0DD"
    )
  )
  
  success_results <- c(FALSE, FALSE)
  
  # Generate SRR vs Gene heatmap
  data_srr <- data_matrix
  sample_annotation_srr <- create_sample_annotation(colnames(data_srr), use_organ_names = FALSE)
  title_srr <- paste(title_base, "- SRR vs", ifelse(matrix_type == "geneID", "Gene ID", "Gene Name"))
  
  tryCatch({
    pheatmap(
      data_srr,
      color = colors,
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean", 
      clustering_method = "complete",
      annotation_col = sample_annotation_srr["Tissue_Type"],
      annotation_colors = ann_colors,
      show_rownames = ifelse(nrow(data_srr) > 50, FALSE, TRUE),
      show_colnames = TRUE,
      main = title_srr,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      filename = output_path_srr,
      width = 12,
      height = 10
    )
    
    cat("SRR Heatmap saved:", output_path_srr, "\n")
    success_results[1] <- TRUE
  }, error = function(e) {
    cat("Error generating SRR heatmap for:", title_base, "\n")
    cat("Error message:", e$message, "\n")
  })
  
  # Generate Organ vs Gene heatmap
  data_organ <- data_matrix
  # Rename columns to organ names
  colnames(data_organ) <- sapply(colnames(data_organ), function(x) {
    label <- SAMPLE_LABELS[x]
    if (!is.na(label)) return(label) else return(x)
  })
  
  sample_annotation_organ <- create_sample_annotation(colnames(data_organ), use_organ_names = TRUE)
  title_organ <- paste(title_base, "- Organ vs", ifelse(matrix_type == "geneID", "Gene ID", "Gene Name"))
  
  tryCatch({
    pheatmap(
      data_organ,
      color = colors,
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean", 
      clustering_method = "complete",
      annotation_col = sample_annotation_organ["Tissue_Type"],
      annotation_colors = ann_colors,
      show_rownames = ifelse(nrow(data_organ) > 50, FALSE, TRUE),
      show_colnames = TRUE,
      main = title_organ,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 10,
      filename = output_path_organ,
      width = 12,
      height = 10
    )
    
    cat("Organ Heatmap saved:", output_path_organ, "\n")
    success_results[2] <- TRUE
  }, error = function(e) {
    cat("Error generating Organ heatmap for:", title_base, "\n")
    cat("Error message:", e$message, "\n")
  })
  
  return(success_results)
}

# ===============================================
# MAIN PROCESSING
# ===============================================

cat("Starting dual heatmap generation...\n")

total_heatmaps <- 0
successful_heatmaps <- 0

for (group in FASTA_GROUPS) {
  for (version in VERSIONS) {
    for (count_type in COUNT_TYPES) {
      for (matrix_type in MATRIX_TYPES) {
        
        # Define input file path
        input_file <- file.path(MATRICES_DIR, group, version, 
                               paste0(group, "_", count_type, "_counts_", matrix_type, ".tsv"))
        
        # Check if input file exists
        if (!file.exists(input_file)) {
          cat("File not found:", input_file, "\n")
          next
        }
        
        # Create output directory structure
        output_dir <- file.path(HEATMAP_OUT_DIR, group, version, matrix_type)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Define output file paths
        output_file_srr <- file.path(output_dir, 
                                    paste0(group, "_", version, "_", count_type, "_", matrix_type, "_SRR_heatmap.png"))
        output_file_organ <- file.path(output_dir, 
                                      paste0(group, "_", version, "_", count_type, "_", matrix_type, "_Organ_heatmap.png"))
        
        # Create title base
        title_base <- paste("Gene Expression:", group, version, toupper(count_type))
        
        cat("\nProcessing:", group, "-", version, "-", count_type, "-", matrix_type, "\n")
        
        # Read and preprocess data
        raw_data <- read_count_matrix(input_file)
        processed_data <- preprocess_for_heatmap(raw_data, count_type)
        
        # Generate dual heatmaps
        total_heatmaps <- total_heatmaps + 2  # Two heatmaps per iteration
        success_results <- generate_dual_heatmaps(processed_data, output_file_srr, output_file_organ, 
                                                  title_base, count_type, matrix_type)
        successful_heatmaps <- successful_heatmaps + sum(success_results)
      }
    }
  }
}

# ===============================================
# SUMMARY
# ===============================================

cat("\n", "="*50, "\n")
cat("DUAL HEATMAP GENERATION SUMMARY\n")
cat("="*50, "\n")
cat("Total heatmaps attempted:", total_heatmaps, "\n")
cat("Successful heatmaps:", successful_heatmaps, "\n")
cat("Failed heatmaps:", total_heatmaps - successful_heatmaps, "\n")
cat("Output directory:", HEATMAP_OUT_DIR, "\n")
cat("="*50, "\n")

if (successful_heatmaps > 0) {
  cat("Dual heatmap generation completed successfully!\n")
  cat("Generated both SRR vs Gene and Organ vs Gene versions\n")
} else {
  cat("No heatmaps were generated. Please check input files and paths.\n")
}
