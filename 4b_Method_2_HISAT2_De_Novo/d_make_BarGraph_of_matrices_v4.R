# ===============================================
# BAR GRAPH GENERATION FOR GENE EXPRESSION
# ===============================================
#
# This script generates bar graphs showing individual gene expression across developmental stages.
# Output is saved to "8_Bar_Graph_Visualizations" directory.
#
# Features:
# - Individual bar plots for each gene
# - Expression values across developmental stages (organs)
# - Color-coded by expression level using custom color scale
# - Professional layout matching reference design
# - Saves both PNG bar graphs and TSV data files
# - Multiple normali# - Each gene generates TWO bar graphs:
#   * sorted_by_organ: Expression values in original developmental stage order (x-axis)
#   * sorted_by_expression: Expression values sorted by intensity (lowest to highest from left to right)
#   * Color-coded bars by expression level using custom color scale
#   * Professional layout matching reference design
#   * Y-axis labeled according to normalization scheme
# - Output structure: 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_organ/{gene_id}_{normalization_scheme}_bargraph.png
#                   : 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_expression/{gene_id}_{normalization_scheme}_bargraph.png
#
# Required packages: ggplot2, dplyr, RColorBrewer, grid, scales
# ===============================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(grid)
  library(scales)
  library(gridExtra)
})

# ===============================================
# CONFIGURATION
# ===============================================

# Input and output directories
BASE_DIR <- getwd()
MATRICES_DIR <- file.path("5_stringtie_WD", "b_Method_2_COUNT_MATRICES")
BARGRAPH_OUT_DIR <- "8_Bar_Graph_Visualizations"

# File naming configuration
MASTER_REFERENCE <- "All_Smel_Genes"
MASTER_REFERENCE_SUFFIX <- paste0("_from_", MASTER_REFERENCE)

# Gene groups
FASTA_GROUPS <- c(
  # Control Gene Groups
  "Best_Cell_Cycle_Associated_Control_Genes",
  "Best_Control_Genes"
  
  # Individual Gene Groups
  #"SmelDMPs",
  #"SmelGIFs",
  #"SmelGRFs"
  
  # Combined Gene Groups with Control Genes
  #"SmelGIF_with_Cell_Cycle_Control_genes",
  #"SmelGRF_with_Cell_Cycle_Control_genes",
  #"SmelGRF-GIF_with_Best_Cell_Cycle_Control_Genes"
)

# Analysis configuration
COUNT_TYPES <- c("coverage", "fpkm", "tpm")
GENE_TYPES <- c("geneID", "geneName")
LABEL_TYPES <- c("SRR", "Organ")

# Normalization schemes configuration
NORMALIZATION_SCHEMES <- c("raw", "raw_normalized", "default", "zscore", "cpm")

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


# ===============================================
# OUTPUT DIRECTORY SETUP
# ===============================================

# Create base output directory if it doesn't exist
dir.create(BARGRAPH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Clean up the output directory at the start of each run
# cat("===============================================\n")
# cat("PREPARING OUTPUT DIRECTORY\n")
# cat("===============================================\n")

if (dir.exists(BARGRAPH_OUT_DIR)) {
  # cat("Cleaning up existing output directory:", BARGRAPH_OUT_DIR, "\n")
  unlink(BARGRAPH_OUT_DIR, recursive = TRUE, force = TRUE)
  # cat("✓ Previous output directory removed\n")
}

dir.create(BARGRAPH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# cat("✓ Fresh output directory created:", BARGRAPH_OUT_DIR, "\n")
# cat("\n")



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

# Function to apply raw_normalized (identical to raw)
preprocess_for_raw_normalized <- function(data_matrix) {
  # Identical to raw normalization
  return(preprocess_for_raw(data_matrix))
}

# Function to apply default (count-type specific) normalization
preprocess_for_default <- function(data_matrix, count_type) {
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
  data_processed <- preprocess_for_default(data_matrix, count_type)
  if (is.null(data_processed) || nrow(data_processed) == 0) return(NULL)
  
  # Then apply Z-score standardization across samples for each gene
  if (nrow(data_processed) > 1 && ncol(data_processed) > 1) {
    data_processed <- t(scale(t(data_processed), center = TRUE, scale = TRUE))
    data_processed[is.na(data_processed)] <- 0
  }
  
  return(data_processed)
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
    "raw_normalized" = preprocess_for_raw_normalized(data_matrix),
    "default" = preprocess_for_default(data_matrix, count_type),
    "zscore" = preprocess_for_zscore(data_matrix, count_type),
    "cpm" = preprocess_for_cpm(data_matrix),
    stop("Unknown normalization scheme: ", normalization_scheme)
  )
}

# ===============================================
# BAR GRAPH GENERATION FUNCTIONS
# ===============================================

# Function to create individual gene bar plots
generate_gene_bargraphs <- function(data_matrix, output_dir, title_prefix, count_type, label_type, normalization_scheme = "default", sort_by_expression = FALSE) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0) return(FALSE)
  if (ncol(data_matrix) == 0) return(FALSE)
  
  # Prepare column labels (developmental stages/organs)
  col_labels <- colnames(data_matrix)
  if (label_type == "SRR") {
    col_labels <- colnames(data_matrix)
  } else if (label_type == "Organ" && exists("SAMPLE_LABELS")) {
    col_labels <- ifelse(colnames(data_matrix) %in% names(SAMPLE_LABELS),
                         SAMPLE_LABELS[colnames(data_matrix)],
                         colnames(data_matrix))
  }
  

  
  successful_plots <- 0
  
  # Create individual bar plots for each gene
  for (i in seq_len(nrow(data_matrix))) {
    gene_name <- rownames(data_matrix)[i]
    gene_values <- as.numeric(data_matrix[i, ])
    
    # Skip genes with all zeros
    if (all(gene_values == 0)) next
    
    # Create data frame for plotting
    plot_data <- data.frame(
      Stage = col_labels,
      Expression = gene_values,
      stringsAsFactors = FALSE
    )
    
    # Set factor levels to control ordering
    if (sort_by_expression) {
      # Sort by expression level (highest on the right)
      plot_data <- plot_data[order(plot_data$Expression), ]
      plot_data$Stage <- factor(plot_data$Stage, levels = plot_data$Stage)
    } else {
      # Keep original developmental stage ordering
      plot_data$Stage <- factor(plot_data$Stage, levels = col_labels)
    }
    
    # Determine y-axis label based on normalization scheme
    y_label <- switch(normalization_scheme,
      "raw" = "Raw counts",
      "raw_normalized" = "Raw counts", 
      "default" = ifelse(count_type == "coverage", "Log2(CPM + 1)", 
                        ifelse(count_type %in% c("fpkm", "tpm"), "Log2(Expression + 0.1)", "Log2(Expression + 1)")),
      "zscore" = "Z-score",
      "cpm" = "Log2(CPM + 1)",
      "normalized counts (CPM)"
    )
    
    # Create the bar plot
    # Define custom color scheme with violet as highest expression
    custom_colors <- c("#FFFFFF", "#FFEB3B", "#CDDC39", "#8BC34A", "#4CAF50", "#009688", "#00BCD4", "#3F51B5", "#673AB7", "#4A148C")
    
    # Use a continuous color scale for fill based on expression
    p <- ggplot(plot_data, aes(x = .data$Stage, y = .data$Expression, fill = .data$Expression)) +
      geom_bar(stat = "identity", color = "black", size = 0.3) +
      labs(
        title = gene_name,
        x = "",
        y = y_label
      ) +
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
    
    # Save the plot
    clean_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene_name)
    output_file <- file.path(output_dir, paste0(clean_gene_name, "_", normalization_scheme, "_bargraph.png"))
    
    tryCatch({
      ggsave(output_file, plot = p, width = 8, height = 6, units = "in", dpi = 300, 
             bg = "white")
      successful_plots <- successful_plots + 1
    }, error = function(e) {
      cat("Error saving bar graph for gene:", gene_name, "\n")
      cat("Error details:", e$message, "\n")
    })
  }
  
  return(successful_plots)
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
      # cat("  ⚠️  Found duplicate row names - making them unique...\n")
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

# Preprocessing for CPM (Counts Per Million) normalized data
preprocess_for_cpm_normalized <- function(data_matrix) {
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

# ===============================================
# MAIN PROCESSING
# ===============================================

# Close any open graphics devices
if (length(dev.list()) > 0) { 
    sapply(dev.list(), dev.off) 
}

# Initialize counters
total_bargraphs <- 0
successful_bargraphs <- 0

# Validate FASTA_GROUPS
if (length(FASTA_GROUPS) == 0) {
    stop("No FASTA groups defined for processing")
}

# Process each FASTA group first
for (group in FASTA_GROUPS) {
  # cat("Processing FASTA group:", group, "\n")
  
  for (normalization_scheme in NORMALIZATION_SCHEMES) {
    # Create normalization-specific output directory within FASTA group
    norm_output_dir <- file.path(BARGRAPH_OUT_DIR, group, normalization_scheme)
    dir.create(norm_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create subdirectories for different sorting methods
    sorted_by_organ_dir <- file.path(norm_output_dir, "sorted_by_organ")
    sorted_by_expression_dir <- file.path(norm_output_dir, "sorted_by_expression")
    dir.create(sorted_by_organ_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(sorted_by_expression_dir, recursive = TRUE, showWarnings = FALSE)
    
    # cat("  ✓ Created directories:", norm_output_dir, "\n")
    for (count_type in COUNT_TYPES) {
      for (gene_type in GENE_TYPES) {
        # Always use SRR label_type for input file (as Organ is not a file label)
        label_type <- "SRR"
        input_file <- file.path(MATRICES_DIR, group, 
                               paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, MASTER_REFERENCE_SUFFIX, ".tsv"))
        # If that doesn't exist, try without the suffix
        if (!file.exists(input_file)) {
          input_file <- file.path(MATRICES_DIR, group, 
                                 paste0(group, "_", count_type, "_counts_", gene_type, "_", label_type, ".tsv"))
        }
        if (!file.exists(input_file)) next
        input_basename <- tools::file_path_sans_ext(basename(input_file))
        title <- gsub("_", " ", input_basename)
        # Read raw count matrix
        raw_data <- read_count_matrix(input_file)
        if (is.null(raw_data)) {
          next
        }
        # Apply the specified normalization scheme
        processed_data <- apply_normalization(raw_data, normalization_scheme, count_type)
        if (!is.null(processed_data)) {
          # Generate bar graphs sorted by organ (developmental stage order)
          # SRR labels, sorted by developmental stage
          genes_plotted_dev_srr <- generate_gene_bargraphs(processed_data, sorted_by_organ_dir, title, count_type, "SRR", normalization_scheme, sort_by_expression = FALSE)
          # Organ labels, sorted by developmental stage  
          genes_plotted_dev_org <- generate_gene_bargraphs(processed_data, sorted_by_organ_dir, title, count_type, "Organ", normalization_scheme, sort_by_expression = FALSE)
          
          # Generate bar graphs sorted by expression level
          # SRR labels, sorted by expression
          genes_plotted_expr_srr <- generate_gene_bargraphs(processed_data, sorted_by_expression_dir, title, count_type, "SRR", normalization_scheme, sort_by_expression = TRUE)
          # Organ labels, sorted by expression
          genes_plotted_expr_org <- generate_gene_bargraphs(processed_data, sorted_by_expression_dir, title, count_type, "Organ", normalization_scheme, sort_by_expression = TRUE)
          
          total_bargraphs <- total_bargraphs + (nrow(processed_data) * 4)  # 2 label types x 2 sortings
          successful_bargraphs <- successful_bargraphs + genes_plotted_dev_srr + genes_plotted_expr_srr + genes_plotted_dev_org + genes_plotted_expr_org
        }
      }
    }
  }
}

# ===============================================
# SUMMARY
# ===============================================

while (length(dev.list()) > 0) { dev.off() }

# cat("===============================================\n")
# cat("BAR GRAPH GENERATION SUMMARY\n")
# cat("===============================================\n")
# cat("Total bar graphs generated:", successful_bargraphs, "\n")

# if (successful_bargraphs > 0) {
#   cat("✓ Bar graph generation completed successfully!\n")
#   cat("\nOutput folder structure:\n")
#   cat(paste0(BARGRAPH_OUT_DIR, "/\n"))
#   for (group in FASTA_GROUPS) {
#     cat(paste0("├── ", group, "/\n"))
#     for (norm_scheme in NORMALIZATION_SCHEMES) {
#       group_dir <- file.path(BARGRAPH_OUT_DIR, group, norm_scheme)
#       if (dir.exists(group_dir)) {
#         organ_dir <- file.path(group_dir, "sorted_by_organ")
#         expr_dir <- file.path(group_dir, "sorted_by_expression")
#         organ_count <- ifelse(dir.exists(organ_dir), length(list.files(organ_dir, pattern = "*.png")), 0)
#         expr_count <- ifelse(dir.exists(expr_dir), length(list.files(expr_dir, pattern = "*.png")), 0)
#         cat(paste0("│   ├── ", norm_scheme, "/\n"))
#         cat(paste0("│   │   ├── sorted_by_organ/  (", organ_count, " graphs)\n"))
#         cat(paste0("│   │   └── sorted_by_expression/  (", expr_count, " graphs)\n"))
#       }
#     }
#   }
# } else {
#   cat("✗ No bar graphs were generated. Please check input files and paths.\n")
# }

# UPDATED STRUCTURE WITH FASTA GROUPS AS PRIMARY ORGANIZATION:
# =============================================================
# - Each combination of (group, normalization_scheme, count_type, gene_type, label_type) 
#   corresponds to exactly one TSV file input and generates individual bar graphs for each gene
# - Five normalization schemes are applied:
#   * "raw": Original count data without transformation
#   * "raw_normalized": Identical to raw (for comparison)
#   * "default": Count-type specific normalization (Coverage: CPM+log2, FPKM/TPM: direct log2)
#   * "zscore": Count-type normalization followed by Z-score standardization
#   * "cpm": Counts Per Million normalization + log2 transformation
# - Each gene generates bar graphs showing:
#   * Expression values across developmental stages (x-axis)
#   * Two versions: sorted_by_organ (developmental order) and sorted_by_expression (intensity order)
#   * Color-coded bars by expression level using custom color scale
#   * Professional layout matching reference design
#   * Y-axis labeled according to normalization scheme
# - Output structure: 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_organ/{gene_id}_{normalization_scheme}_bargraph.png
#                   : 8_Bar_Graph_Visualizations/{gene_group}/{normalization_scheme}/sorted_by_expression/{gene_id}_{normalization_scheme}_bargraph.png
# - Bar colors represent expression intensity: white (lowest) to violet (highest expression)
# - Individual files for each gene allow easy comparison across conditions
# - Separate folders for different sorting methods facilitate analysis workflows
