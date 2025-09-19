#!/usr/bin/env Rscript

# ===============================================
# Package Installation and Loading
# ===============================================
# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      
      # Check if it's a Bioconductor package
      if (pkg %in% c("DESeq2", "edgeR", "limma", "GenomicFeatures", "GenomicAlignments")) {
        # Install BiocManager if not available
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cran.r-project.org/")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        # Regular CRAN package
        install.packages(pkg, repos = "https://cran.r-project.org/")
      }
    }
    
    # Load the package
    library(pkg, character.only = TRUE)
    cat("✓ Loaded:", pkg, "\n")
  }
}

# Define required packages
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "DESeq2",       # Differential expression analysis
  "pheatmap",     # Heatmap generation
  "RColorBrewer", # Color palettes
  "viridis"       # Color palettes
)

cat("Checking and installing required packages...\n")
cat("Required packages:", paste(required_packages, collapse = ", "), "\n\n")

# Install and load packages
install_and_load(required_packages)

cat("\n✓ All packages loaded successfully!\n\n")

# ===============================================
# Libraries (now redundant but kept for compatibility)
# ===============================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
})

# ===============================================
# Configuration
# ===============================================
count_dir <- "5_stringtie/method_2.2_count_matrices"   # <-- output of bash script
out_dir   <- "6_Visualization/heatmaps_method_2.2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================================
# Helper function
# ===============================================
process_group <- function(file) {
  # Extract group name from filename (removes _counts.tsv suffix)
  group_name <- gsub("_counts\.tsv$", "", basename(file))
  
  # Create a clean name for output files (replace special characters)
  clean_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)
  
  message("\n=== Processing: ", group_name, " ===")
  message("Input file: ", file)
  message("Clean output name: ", clean_name)

  # Read count matrix
  counts <- read.delim(file, row.names = 1, check.names = FALSE)

  # Remove rows with all zeros
  counts <- counts[rowSums(counts) > 0, ]

  # Build metadata - attempt to infer conditions from sample names
  # You should customize this logic based on your experimental design
  samples <- colnames(counts)
  
  # Improved tissue/condition classification based on SRR IDs
  condition <- case_when(
    grepl("SRR3884631", samples) ~ "Fruits_6cm",
    grepl("SRR3884677", samples) ~ "Cotyledons", 
    grepl("SRR3884679", samples) ~ "Pistils",
    grepl("SRR3884597", samples) ~ "Flowers",
    grepl("SRR3884686", samples) ~ "Buds_0.7cm",
    grepl("SRR3884687", samples) ~ "Buds_Opened",
    grepl("SRR3884689", samples) ~ "Leaves",
    grepl("SRR3884690", samples) ~ "Stems",
    grepl("SRR3884685", samples) ~ "Radicles",
    grepl("SRR3884675", samples) ~ "Roots",
    TRUE ~ "Other"
  )
  
  # If all samples have same condition, create artificial groups for visualization
  if(length(unique(condition)) == 1) {
    condition <- rep(c("GroupA", "GroupB"), length.out = length(samples))
    message("Warning: All samples have same condition. Creating artificial groups for visualization.")
  }
  
  coldata <- data.frame(row.names = samples, condition = factor(condition))

  # DESeq2 analysis
  mat <- NULL
  normalized_counts <- NULL
  
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData   = coldata,
                                  design    = ~ condition)
    
    # Filter low count genes (at least 10 reads total)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    
    message("Filtered to ", nrow(dds), " genes with sufficient counts")
    
    # Run DESeq2 analysis if we have multiple conditions
    if (length(unique(coldata$condition)) > 1) {
      dds <- DESeq(dds)
      
      # Get normalized counts
      normalized_counts <- counts(dds, normalized = TRUE)
      
      # Variance-stabilizing transform for visualization
      vsd <- vst(dds, blind = FALSE)
      mat <- assay(vsd)
      
      message("DESeq2 analysis completed successfully")
    } else {
      # If only one condition, use rlog transformation
      vsd <- vst(dds, blind = TRUE)
      mat <- assay(vsd)
      normalized_counts <- counts(dds, normalized = FALSE)
    }
    
  }, error = function(e) {
    message("Error in DESeq2 analysis for ", group_name, ": ", e$message)
    message("Using fallback log2 transformation")
    
    # Fallback: simple log2 transformation with pseudocount
    filtered_counts <- counts[rowSums(counts) >= 10, ]
    mat <- log2(filtered_counts + 1)
    normalized_counts <- filtered_counts
  })

  # Calculate gene variability
  gene_vars <- rowVars(mat)
  
  # Select top variable genes (different numbers for different views)
  top50_genes <- head(order(gene_vars, decreasing = TRUE), 50)
  top100_genes <- head(order(gene_vars, decreasing = TRUE), 100)
  
  # Create color annotation for samples
  annotation_colors <- list(
    condition = rainbow(length(unique(coldata$condition)))
  )
  names(annotation_colors$condition) <- unique(coldata$condition)
  
  # Generate multiple heatmaps
  
  # Generate timestamp for unique naming
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create output file names with clean naming
  base_name <- paste0(clean_name, "_", timestamp)
  
  # 1. Top 50 most variable genes - detailed view
  top50_file <- file.path(out_dir, paste0(base_name, "_top50_genes_heatmap.jpg"))
  jpeg(top50_file, width = 1400, height = 1000, res = 150, quality = 95)
  pheatmap(t(mat[top50_genes, ]),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_row = coldata,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste("Top 50 Variable Genes -", group_name),
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
  
  # 2. Top 100 most variable genes - overview
  top100_file <- file.path(out_dir, paste0(base_name, "_top100_genes_heatmap.jpg"))
  jpeg(top100_file, width = 1600, height = 1000, res = 150, quality = 95)
  pheatmap(t(mat[top100_genes, ]),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_row = coldata,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = FALSE,
           fontsize = 10,
           fontsize_row = 8,
           main = paste("Top 100 Variable Genes -", group_name),
           color = viridis(100))
  dev.off()
  
  # 3. Sample correlation heatmap (samples vs samples)
  sample_cor <- cor(mat, use = "complete.obs")
  correlation_file <- file.path(out_dir, paste0(base_name, "_sample_correlation.jpg"))
  jpeg(correlation_file, width = 1000, height = 1000, res = 150, quality = 95)
  pheatmap(sample_cor,
           annotation_col = coldata,
           annotation_row = coldata,
           annotation_colors = annotation_colors,
           main = paste("Sample Correlation -", group_name),
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
  
  # Save normalized count matrix
  if (!is.null(normalized_counts)) {
    norm_counts_file <- file.path(out_dir, paste0(base_name, "_normalized_counts.csv"))
    write.csv(normalized_counts, norm_counts_file, quote = FALSE)
    message("Saved normalized counts: ", basename(norm_counts_file))
  }
  
  # Save top variable genes list
  top_genes_df <- data.frame(
    GeneID = rownames(mat)[top100_genes],
    Variance = gene_vars[top100_genes],
    Rank = 1:100
  )
  top_genes_file <- file.path(out_dir, paste0(base_name, "_top_variable_genes.csv"))
  write.csv(top_genes_df, top_genes_file, quote = FALSE, row.names = FALSE)
  
  message("\n✓ Successfully generated outputs for: ", group_name)
  message("  Files created:")
  message("    - ", basename(top50_file))
  message("    - ", basename(top100_file)) 
  message("    - ", basename(correlation_file))
  if (!is.null(normalized_counts)) {
    message("    - ", basename(norm_counts_file))
  }
  message("    - ", basename(top_genes_file))
  message("  Total genes analyzed: ", nrow(mat))
  message("  Samples in analysis: ", ncol(mat))
  message("=== Completed: ", group_name, " ===")
}

# ===============================================
# Main
# ===============================================
# Check if count directory exists
if (!dir.exists(count_dir)) {
  stop("Count directory not found: ", count_dir, 
       "\nMake sure the bash script has been run successfully.")
}

files <- list.files(count_dir, pattern = "_counts.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No *_counts.tsv files found in ", count_dir)
}

# Display found TSV files
message("Found ", length(files), " TSV count matrices to process:")
for (i in seq_along(files)) {
  message("  ", i, ". ", basename(files[i]))
}
message("\nStarting heatmap generation...\n")

# Process each TSV file
lapply(files, process_group)

# Generate summary report
summary_file <- file.path(out_dir, "analysis_summary.txt")
cat("DESeq2 Heatmap Analysis Summary\n", file = summary_file)
cat("==============================\n", file = summary_file, append = TRUE)
cat("Analysis Date:", Sys.time(), "\n", file = summary_file, append = TRUE)
cat("Input Directory:", count_dir, "\n", file = summary_file, append = TRUE)
cat("Output Directory:", out_dir, "\n", file = summary_file, append = TRUE)
cat("Number of groups processed:", length(files), "\n", file = summary_file, append = TRUE)
cat("\nProcessed Groups:\n", file = summary_file, append = TRUE)

for (file in files) {
  group_name <- gsub("_counts.tsv$", "", basename(file))
  cat("- ", group_name, "\n", file = summary_file, append = TRUE)
}

cat("\nOutput Files Generated (per group):\n", file = summary_file, append = TRUE)
cat("- [group]_top50_heatmap.jpg\n", file = summary_file, append = TRUE)
cat("- [group]_top100_heatmap.jpg\n", file = summary_file, append = TRUE)
cat("- [group]_sample_correlation.jpg\n", file = summary_file, append = TRUE)
cat("- [group]_normalized_counts.csv\n", file = summary_file, append = TRUE)
cat("- [group]_top_variable_genes.csv\n", file = summary_file, append = TRUE)

message("Analysis complete! Summary saved to: ", summary_file)
message("All heatmaps and data saved in: ", out_dir)
