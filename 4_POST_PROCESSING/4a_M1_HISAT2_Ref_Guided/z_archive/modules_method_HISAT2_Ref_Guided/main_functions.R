# ---
# Title: "Main Analysis Functions"
# ---

# Load libraries
library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)

# Source utility functions and configurations
source("modules_method_HISAT2_Ref_Guided/utility_functions.R")
source("modules_method_HISAT2_Ref_Guided/inputs_configs_outputs.R")

# --- DESeq2 Analysis ---
run_deseq2 <- function(count_matrix, coldata) {
  log_message("Running DESeq2...")
  
  # Ensure row names of coldata match column names of count_matrix
  if (!all(rownames(coldata) == colnames(count_matrix))) {
    log_message("Row names of coldata do not match column names of count matrix. Reordering coldata.", level = "WARNING")
    coldata <- coldata[colnames(count_matrix), , drop = FALSE]
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq(dds)
  log_message("DESeq2 analysis complete.")
  return(dds)
}

# --- Heatmap Generation ---
generate_basic_heatmap <- function(dds, gene_group_name, count_type, norm_scheme, organ_label_version, output_dir) {
  log_message(paste("Generating basic heatmap for", gene_group_name, "with", norm_scheme, "normalization..."))
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Get gene group
  gene_group_file <- file.path(gene_groups_dir, paste0(gene_group_name, ".txt"))
  if (!file.exists(gene_group_file)) {
    log_message(paste("Gene group file not found:", gene_group_file), level = "ERROR")
    return()
  }
  gene_group <- read.table(gene_group_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  # Filter counts for gene group
  plot_data <- norm_counts[rownames(norm_counts) %in% gene_group, ]
  
  # Apply normalization scheme
  if (norm_scheme == "zscore") {
    plot_data <- t(scale(t(plot_data)))
  }
  
  # Create heatmap
  output_file <- file.path(output_dir, paste0(gene_group_name, "_", count_type, "_geneName_", organ_label_version, "_", norm_scheme, ".png"))
  
  # Define color palette
  colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ht <- Heatmap(
    plot_data,
    name = "Expression",
    col = colors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = paste(gene_group_name, norm_scheme, "normalized"),
    row_names_gp = gpar(fontsize = 8)
  )
  
  png(output_file, width = 8, height = 10, units = "in", res = 300)
  draw(ht)
  dev.off()
  
  log_message(paste("Heatmap saved to", output_file))
}

# --- CV Heatmap Generation ---
generate_cv_heatmap <- function(dds, gene_group_name, count_type, norm_scheme, organ_label_version, output_dir) {
  log_message(paste("Generating CV heatmap for", gene_group_name, "..."))
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Get gene group
  gene_group_file <- file.path(gene_groups_dir, paste0(gene_group_name, ".txt"))
  if (!file.exists(gene_group_file)) {
    log_message(paste("Gene group file not found:", gene_group_file), level = "ERROR")
    return()
  }
  gene_group <- read.table(gene_group_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  # Filter counts for gene group
  plot_data <- norm_counts[rownames(norm_counts) %in% gene_group, ]
  
  # Calculate CV
  cv_values <- apply(plot_data, 1, function(x) sd(x) / mean(x))
  
  # Create heatmap with CV annotation
  output_file <- file.path(output_dir, paste0(gene_group_name, "_", count_type, "_geneName_", organ_label_version, "_cv_heatmap.png"))
  
  # Define color palette
  colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  # Row annotation
  ha <- rowAnnotation(CV = anno_barplot(cv_values))
  
  ht <- Heatmap(
    t(scale(t(plot_data))),
    name = "Z-score",
    col = colors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = paste(gene_group_name, "CV Heatmap"),
    right_annotation = ha,
    row_names_gp = gpar(fontsize = 8)
  )
  
  png(output_file, width = 10, height = 10, units = "in", res = 300)
  draw(ht)
  dev.off
  
  log_message(paste("CV Heatmap saved to", output_file))
}

# --- Bar Graph Generation ---
generate_bar_graphs <- function(dds, gene_group_name, output_dir) {
  log_message(paste("Generating bar graphs for", gene_group_name, "..."))
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Get gene group
  gene_group_file <- file.path(gene_groups_dir, paste0(gene_group_name, ".txt"))
  if (!file.exists(gene_group_file)) {
    log_message(paste("Gene group file not found:", gene_group_file), level = "ERROR")
    return()
  }
  gene_group <- read.table(gene_group_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  # Filter counts for gene group
  plot_data <- norm_counts[rownames(norm_counts) %in% gene_group, ]
  
  # Create bar graphs for each gene
  for (i in 1:nrow(plot_data)) {
    gene_name <- rownames(plot_data)[i]
    df <- data.frame(
      sample = colnames(plot_data),
      counts = plot_data[i, ]
    )
    
    p <- ggplot(df, aes(x = sample, y = counts, fill = sample)) +
      geom_bar(stat = "identity") +
      labs(title = gene_name, x = "Sample", y = "Normalized Counts") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    output_file <- file.path(output_dir, paste0(gene_group_name, "_", gene_name, "_bargraph.png"))
    ggsave(output_file, plot = p, width = 8, height = 6, units = "in")
  }
  
  log_message(paste("Bar graphs saved to", output_dir))
}