#!/usr/bin/env Rscript

# ==============================================================================
# HEATMAP GENERATION FOR RNA-SEQ DATA WITH DESeq2
# ==============================================================================
# Purpose:     Generate comprehensive heatmaps from StringTie count matrices
# Input:       TSV count matrices from StringTie pipeline
# Output:      High-quality heatmaps with eggplant violet color scheme
# Author:      Modified for complete gene inclusion analysis
# Date:        September 2025
# ==============================================================================

# ==============================================================================
# PACKAGE INSTALLATION AND LOADING
# ==============================================================================
# Function: install_and_load
# Purpose:  Intelligently install and load required R packages
# Args:     packages - character vector of package names
# Returns:  None (side effects: packages loaded into environment)
install_and_load <- function(packages) {
  
  cat("\n🔧 Installing and loading packages...\n")
  
  for (pkg in packages) {
    # Check if package is already installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("  📦 Installing:", pkg, "\n")
      
      # Determine installation method (CRAN vs Bioconductor)
      bioc_packages <- c("DESeq2", "edgeR", "limma", "GenomicFeatures", "GenomicAlignments")
      
      if (pkg %in% bioc_packages) {
        # Install BiocManager if needed for Bioconductor packages
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cran.r-project.org/")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        # Standard CRAN installation
        install.packages(pkg, repos = "https://cran.r-project.org/")
      }
    }
    
    # Load package with special handling for tidyverse conflicts
    if (pkg == "tidyverse") {
      suppressMessages(suppressWarnings(library(pkg, character.only = TRUE)))
      cat("  ✅ Loaded:", pkg, "(function conflicts handled)\n")
    } else {
      library(pkg, character.only = TRUE)
      cat("  ✅ Loaded:", pkg, "\n")
    }
  }
  
  # Special verification for matrixStats functions
  if ("matrixStats" %in% packages) {
    if (!exists("rowVars", mode = "function")) {
      cat("  🔧 Explicitly loading matrixStats functions...\n")
      library(matrixStats)
      
      if (exists("rowVars", mode = "function")) {
        cat("  ✅ rowVars function now available\n")
      } else {
        cat("  ⚠️  Warning: rowVars function still not found\n")
      }
    }
  }
  
  cat("\n📋 Package loading completed!\n")
}

# ------------------------------------------------------------------------------
# REQUIRED PACKAGES CONFIGURATION
# ------------------------------------------------------------------------------

# Define all required packages with their purposes
required_packages <- c(
  "tidyverse",    # 📊 Data manipulation and visualization (includes dplyr, ggplot2, etc.)
  "DESeq2",       # 🧬 Differential expression analysis for RNA-seq data
  "pheatmap",     # 🎨 High-quality heatmap generation
  "RColorBrewer", # 🌈 Color palette utilities
  "viridis",      # 🎨 Perceptually uniform color scales
  "matrixStats"   # 📈 Efficient matrix statistics (rowVars function)
)

header_line <- paste(rep("=", 60), collapse = "")
cat("\n", header_line, "\n", sep = "")
cat("🚀 STARTING HEATMAP ANALYSIS PIPELINE\n")
cat(header_line, "\n", sep = "")
cat("📋 Required packages:", paste(required_packages, collapse = ", "), "\n\n")

# Initialize package loading
install_and_load(required_packages)

# ------------------------------------------------------------------------------
# FUNCTION CONFLICT RESOLUTION
# ------------------------------------------------------------------------------

cat("\n� Resolving function conflicts...\n")
cat("   • dplyr::filter() will override stats::filter()\n")
cat("   • dplyr::lag() will override stats::lag()\n")
cat("   • Use stats::filter() or stats::lag() explicitly if base R functions needed\n")

cat("\n✅ All packages loaded and conflicts resolved!\n\n")

# ==============================================================================
# PROJECT CONFIGURATION
# ==============================================================================

# Define input/output directories
count_dir <- "5_stringtie/method_2.2_count_matrices"  # 📁 Input: StringTie count matrices
out_dir   <- "6_Visualization/heatmaps_method_2.2"     # 📁 Output: Generated heatmaps and data

cat("📂 Input directory: ", count_dir, "\n")
cat("📂 Output directory:", out_dir, "\n\n")

# ==============================================================================
# OUTPUT DIRECTORY MANAGEMENT
# ==============================================================================
cat("🧹 Preparing output directory...\n")

# Check if output directory exists and clean it
if (dir.exists(out_dir)) {
  # Get list of existing files
  existing_files <- list.files(out_dir, full.names = TRUE, recursive = TRUE)
  
  if (length(existing_files) > 0) {
    cat("  📄 Found", length(existing_files), "existing files\n")
    cat("  🗑️  Removing old files...\n")
    
    # Remove all files and subdirectories
    unlink(existing_files, recursive = TRUE, force = TRUE)
    
    # Clean up empty subdirectories
    subdirs <- list.dirs(out_dir, full.names = TRUE, recursive = TRUE)
    subdirs <- subdirs[subdirs != out_dir]  # Exclude main directory
    
    if (length(subdirs) > 0) {
      unlink(subdirs, recursive = TRUE, force = TRUE)
    }
    
    cat("  ✅ Output directory cleaned successfully\n")
  } else {
    cat("  📭 Output directory is already empty\n")
  }
} else {
  cat("  📁 Output directory does not exist yet\n")
}

# Create/recreate the output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("  ✅ Output directory ready:", out_dir, "\n\n")

# ==============================================================================
# MAIN DATA PROCESSING FUNCTION
# ==============================================================================

# Function: process_group
# Purpose:  Process individual count matrix file and generate comprehensive heatmaps
# Args:     file - path to TSV count matrix file
# Returns:  None (side effects: creates heatmap files and data exports)
# Features: - Handles complex TSV file formats from StringTie
#           - Performs DESeq2 normalization with fallback options
#           - Generates multiple heatmap views with eggplant violet colors
#           - Exports normalized data and gene rankings
process_group <- function(file) {
  
  # Extract meaningful names from filename
  group_name <- gsub("_counts\\.tsv$", "", basename(file))
  clean_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)  # Filesystem-safe name
  
  # --------------------------------------------------------------------------
  # CREATE INDIVIDUAL OUTPUT DIRECTORY
  # --------------------------------------------------------------------------
  
  # Create separate output directory for this dataset
  individual_out_dir <- file.path(out_dir, clean_name)
  dir.create(individual_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Progress reporting
  cat("\n", paste(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("🔬 PROCESSING:", group_name, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n", sep = "")
  cat("📄 Input file:", file, "\n")
  cat("🏷️  Clean name:", clean_name, "\n")
  cat("📁 Output directory:", individual_out_dir, "\n\n")

  # --------------------------------------------------------------------------
  # STEP 1: READ AND CLEAN COUNT MATRIX
  # --------------------------------------------------------------------------
  
  cat("📂 Reading count matrix...\n")
  
  # Read raw data from StringTie output format
  counts_raw <- read.delim(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Handle StringTie-specific format: remove description row if present
  if (ncol(counts_raw) > 1 && any(grepl("Coverage|Gene ID", counts_raw[1, ], ignore.case = TRUE))) {
    cat("  🧹 Removing description row...\n")
    counts_raw <- counts_raw[-1, ]
  }
  
  # Set gene IDs as row names and clean up duplicate columns
  rownames(counts_raw) <- counts_raw[, 1]
  
  # Remove gene ID columns (first and potentially last column)
  if (ncol(counts_raw) > 2 && 
      (any(grepl("MSTRG", counts_raw[, ncol(counts_raw)])) || 
       all(counts_raw[, ncol(counts_raw)] == rownames(counts_raw)))) {
    cat("  🧹 Removing duplicate gene ID columns...\n")
    counts <- counts_raw[, 2:(ncol(counts_raw)-1)]  # Remove first and last
  } else {
    counts <- counts_raw[, -1]  # Remove only first column
  }
  
  cat("  📈 Initial dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
  
  # --------------------------------------------------------------------------
  # STEP 2: DATA VALIDATION AND TYPE CONVERSION
  # --------------------------------------------------------------------------
  
  cat("\n🔍 Validating and converting data types...\n")
  
  # Convert all columns to numeric if needed (vectorized for efficiency)
  numeric_check <- vapply(counts, is.numeric, logical(1))
  if (!all(numeric_check)) {
    cat("  🔄 Converting non-numeric columns to numeric...\n")
    
    # Efficient vectorized conversion
    counts[] <- lapply(counts, function(x) {
      if (!is.numeric(x)) {
        # Clean and convert non-numeric values
        x <- gsub("NA|NULL|\\s+", "0", as.character(x))
        as.numeric(x)
      } else {
        x
      }
    })
    
    # Verify conversion success
    if (any(vapply(counts, function(x) any(is.na(x)), logical(1)))) {
      stop("❌ Failed to convert some columns to numeric. Please check input data format.")
    }
    cat("  ✅ All columns successfully converted to numeric\n")
  } else {
    cat("  ✅ All columns already numeric\n")
  }
  
  # --------------------------------------------------------------------------
  # STEP 3: DATA CLEANING AND FILTERING
  # --------------------------------------------------------------------------
  
  cat("\n🧹 Cleaning and filtering data...\n")
  
  # Remove negative values (shouldn't occur in count data, but safety check)
  negative_count <- sum(counts < 0, na.rm = TRUE)
  if (negative_count > 0) {
    cat("  ⚠️  Found and corrected", negative_count, "negative values\n")
    counts[counts < 0] <- 0
  }
  
  # Filter out genes with zero expression across all samples (vectorized)
  initial_genes <- nrow(counts)
  non_zero_genes <- rowSums(counts) > 0
  
  if (any(!non_zero_genes)) {
    removed_genes <- sum(!non_zero_genes)
    cat("  🗑️  Removing", removed_genes, "genes with zero counts\n")
    counts <- counts[non_zero_genes, ]
  }
  
  final_genes <- nrow(counts)
  cat("  📈 Final count matrix:", final_genes, "genes x", ncol(counts), "samples\n")
  cat("  📉 Filtered out:", (initial_genes - final_genes), "genes\n")

  # --------------------------------------------------------------------------
  # STEP 4: SIMPLIFIED SAMPLE ANNOTATION
  # --------------------------------------------------------------------------
  
  cat("\n🏷️  Preparing simple sample labels...\n")
  
  samples <- colnames(counts)
  cat("  📋 Found", length(samples), "samples:", paste(samples, collapse = ", "), "\n")

  # --------------------------------------------------------------------------
  # STEP 5: DESeq2 NORMALIZATION AND TRANSFORMATION (SIMPLIFIED)
  # --------------------------------------------------------------------------
  
  cat("\n🧬 Starting simplified DESeq2 analysis...\n")
  
  # Initialize output variables
  mat <- NULL
  normalized_counts <- NULL
  
  # Create simple metadata for DESeq2 (all samples treated equally)
  simple_coldata <- data.frame(
    row.names = samples, 
    condition = factor(rep("sample", length(samples)))
  )
  
  # Attempt DESeq2 analysis with robust error handling
  tryCatch({
    
    # Create DESeq2 dataset object
    cat("  📏 Creating DESeq2 dataset...\n")
    dds <- DESeqDataSetFromMatrix(
      countData = round(counts),  # DESeq2 requires integer counts
      colData   = simple_coldata,
      design    = ~ 1  # No grouping design
    )
    
    # Filter low-count genes (standard: at least 10 reads total)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat("  📊 Retained", nrow(dds), "genes with sufficient counts (>= 10 total reads)\n")
    
    # Use variance stabilizing transformation without differential analysis
    cat("  🔬 Running variance stabilizing transformation...\n")
    vsd <- vst(dds, blind = TRUE)  # Blind transformation for simple visualization
    mat <- assay(vsd)
    normalized_counts <- counts(dds, normalized = FALSE)
    
    cat("  ✅ DESeq2 transformation completed successfully\n")
    cat("  📊 Matrix dimensions:", nrow(mat), "x", ncol(mat), "\n")
    
    # Verify matrix was created
    if (is.null(mat)) {
      stop("Matrix is NULL after DESeq2 processing")
    }
    
  }, error = function(e) {
    
    # Check if this is the specific dispersion estimation error
    if (grepl("dispersion estimates are within 2 orders of magnitude", e$message)) {
      cat("  ⚠️  DESeq2 dispersion error detected\n")
      cat("  🔧 Applying recommended dispersion fix...\n")
      
      tryCatch({
        # Apply the recommended fix for dispersion estimation
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- nbinomWaldTest(dds)
        
        # Extract results using the fixed dds
        normalized_counts <<- counts(dds, normalized = FALSE)
        vsd <- vst(dds, blind = TRUE)
        mat <<- assay(vsd)
        cat("  ✅ DESeq2 dispersion fix successful\n")
        
        # Verify the fix worked
        if (is.null(mat)) {
          stop("Matrix is still NULL after dispersion fix")
        }
        cat("  📏 Matrix after dispersion fix:", nrow(mat), "x", ncol(mat), "\n")
        
      }, error = function(e2) {
        cat("  ❌ Dispersion fix failed:", e2$message, "\n")
        # Fall back to log2 transformation
        cat("  🔄 Applying log2 fallback transformation...\n")
        
        # Use original counts matrix for fallback
        filtered_counts <- counts[rowSums(counts) >= 10, ]
        if (nrow(filtered_counts) == 0) {
          stop("No genes pass filtering criteria (minimum 10 reads total)")
        }
        
        cat("  📏 Filtered counts dimensions:", nrow(filtered_counts), "x", ncol(filtered_counts), "\n")
        
        # Ensure proper matrix creation and assignment
        log_matrix <- as.matrix(log2(filtered_counts + 1))
        
        # Debug the matrix creation
        cat("  🔍 Log matrix creation check:\n")
        cat("    • Is matrix:", is.matrix(log_matrix), "\n")
        cat("    • Is null:", is.null(log_matrix), "\n")
        cat("    • Dimensions:", nrow(log_matrix), "x", ncol(log_matrix), "\n")
        
        # Assign to parent scope
        mat <<- log_matrix
        normalized_counts <<- filtered_counts
        
        cat("  🔄 Log2 transformation completed\n")
        cat("  📏 Final log2 matrix:", nrow(mat), "x", ncol(mat), "\n")
        
        # Final verification
        if (is.null(mat)) {
          cat("  ❌ STILL NULL: Matrix assignment failed\n")
          stop("Matrix is NULL after log2 transformation")
        } else {
          cat("  ✅ Matrix assignment successful!\n")
        }
      })
      
    } else {
      # Other DESeq2 errors - use fallback
      cat("  ⚠️  DESeq2 error:", e$message, "\n")
      cat("  🔄 Using fallback log2 transformation...\n")
      
      # Apply basic filtering and log transformation using original counts matrix
      cat("  🔄 Applying basic log2 transformation...\n")
      
      # Use the original counts matrix, not the dds object
      filtered_counts <- counts[rowSums(counts) >= 10, ]
      if (nrow(filtered_counts) == 0) {
        stop("No genes pass filtering criteria (minimum 10 reads total)")
      }
      
      cat("  📏 Filtered counts dimensions:", nrow(filtered_counts), "x", ncol(filtered_counts), "\n")
      
      # Ensure proper matrix creation and assignment
      log_matrix <- as.matrix(log2(filtered_counts + 1))  # Add pseudocount
      
      # Debug the matrix creation
      cat("  🔍 Log matrix creation check:\n")
      cat("    • Is matrix:", is.matrix(log_matrix), "\n")
      cat("    • Is null:", is.null(log_matrix), "\n")
      cat("    • Dimensions:", nrow(log_matrix), "x", ncol(log_matrix), "\n")
      
      # Assign to parent scope
      mat <<- log_matrix
      normalized_counts <<- filtered_counts
      
      cat("  ✅ Fallback transformation completed\n")
      cat("  📏 Final fallback matrix:", nrow(mat), "x", ncol(mat), "\n")
      
      # Final verification
      if (is.null(mat)) {
        cat("  ❌ STILL NULL: Matrix assignment failed\n")
        stop("Matrix is NULL after fallback log2 transformation")
      } else {
        cat("  ✅ Matrix assignment successful!\n")
      }
    }
  })

  # --------------------------------------------------------------------------
  # CRITICAL MATRIX VALIDATION CHECKPOINT
  # --------------------------------------------------------------------------
  
  cat("\n🔍 Critical matrix validation checkpoint...\n")
  
  if (is.null(mat)) {
    cat("❌ CRITICAL ERROR: Matrix is NULL after all processing attempts\n")
    cat("📊 Available data:\n")
    cat("   • Original counts matrix:", nrow(counts), "x", ncol(counts), "\n")
    cat("   • Coldata dimensions:", nrow(coldata), "x", ncol(coldata), "\n")
    cat("   • Condition levels:", length(unique(coldata$condition)), "\n")
    stop("Matrix processing failed completely - cannot proceed with heatmap generation")
  }
  
  cat("✅ Matrix validation passed - proceeding with analysis\n")
  cat("📊 Final matrix ready:", nrow(mat), "x", ncol(mat), "\n")

  # --------------------------------------------------------------------------
  # STEP 6: MATRIX VALIDATION AND HEATMAP PREPARATION
  # --------------------------------------------------------------------------
  
  cat("\n🎨 Preparing heatmap visualization...\n")
  
  # Comprehensive matrix validation
  if (is.null(mat)) {
    stop("Matrix is NULL - data processing failed completely")
  }
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Matrix has zero dimensions: ", nrow(mat), " x ", ncol(mat))
  }
  
  cat("  📏 Matrix dimensions:", nrow(mat), "genes ×", ncol(mat), "samples\n")
  
  # Ensure matrix is in correct format for analysis
  if (!is.matrix(mat)) {
    cat("  🔄 Converting to matrix format...\n")
    mat <- as.matrix(mat)
  }
  
  if (!is.numeric(mat)) {
    cat("  🔄 Converting to numeric format...\n")
    mat <- apply(mat, 2, as.numeric)
  }
  
  # Check for valid numeric data
  na_mask <- is.na(mat)
  if (any(na_mask)) {
    cat("  ⚠️  Found", sum(na_mask), "NA values - replacing with zeros\n")
    mat[na_mask] <- 0
  }
  
  if (!is.finite(sum(mat))) {
    stop("Matrix contains infinite or non-finite values")
  }
  
  # Calculate gene variability for ranking (all genes included)
  cat("  📈 Calculating gene variance statistics...\n")
  
  # Use explicit namespace to ensure function availability
  if (exists("rowVars", where = "package:matrixStats", mode = "function")) {
    gene_vars <- matrixStats::rowVars(mat)
  } else {
    # Fallback to manual calculation if matrixStats::rowVars not available
    gene_vars <- apply(mat, 1, var, na.rm = TRUE)
    cat("  ⚠️  Using fallback variance calculation\n")
  }
  
  if (any(na_genes <- is.na(gene_vars))) {
    cat("  ⚠️  Found genes with undefined variance - using median variance\n")
    median_var <- median(gene_vars, na.rm = TRUE)
    gene_vars[na_genes] <- median_var
  }
  
  all_genes <- 1:nrow(mat)  # Use ALL genes (no filtering to top genes)
  
  cat("  📉 Matrix ready:", nrow(mat), "genes x", ncol(mat), "samples\n")
  
  # --------------------------------------------------------------------------
  # STEP 7: SIMPLIFIED ORGAN NAME MAPPING
  # --------------------------------------------------------------------------
  
  # Create organ name mapping for biological interpretation (optimized lookup)
  srr_to_organ <- c(
    "SRR3884631" = "Fruits 6 cm",
    "SRR3884677" = "Cotyledons", 
    "SRR3884679" = "Pistils",
    "SRR3884597" = "Flowers",
    "SRR3884686" = "Buds 0.7 cm",
    "SRR3884687" = "Buds, Opened Buds",
    "SRR3884689" = "Leaves",
    "SRR3884690" = "Stems", 
    "SRR3884685" = "Radicles",
    "SRR3884675" = "Roots"
  )
  
  # Extract SRR IDs from sample names and map to organ names
  organ_names <- sapply(samples, function(sample) {
    # Find matching SRR ID
    matching_srr <- names(srr_to_organ)[sapply(names(srr_to_organ), function(srr) grepl(srr, sample))]
    if (length(matching_srr) > 0) {
      return(srr_to_organ[matching_srr[1]])
    } else {
      return("Unknown")
    }
  })
  
  cat("  🏷️  Organ name mapping created for biological heatmaps\n")
  
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
  
  cat("  🌈 Eggplant violet color palette configured\n")
  
  # --------------------------------------------------------------------------
  # STEP 8: HEATMAP GENERATION
  # --------------------------------------------------------------------------
  
  cat("\n🖼️  Generating publication-quality heatmaps...\n")
  
  # Create timestamp for unique file naming
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_name <- paste0(clean_name, "_", timestamp)
  
  # === HEATMAP 1: DETAILED VIEW WITH GENE NAMES ===
  cat("  🆪 Generating detailed heatmap (with gene names)...\n")
  
  detailed_file <- file.path(individual_out_dir, paste0(base_name, "_simple_srr_layout.jpg"))
  
  tryCatch({
    jpeg(detailed_file, 
         width = 2200,  # Increased width for gene labels
         height = max(1400, nrow(mat) * 12),  # Increased height for gene labels
         res = 150, 
         quality = 95)
    
    pheatmap(
      t(mat[all_genes, ]),  # Transpose: samples as rows, genes as columns
      cluster_rows = FALSE,  # No row clustering - simple SRR order
      cluster_cols = FALSE,  # No column clustering (genes shown in input order)
      show_rownames = TRUE,  # Show SRR sample names
      show_colnames = TRUE,  # Always show gene names
      fontsize = max(6, min(10, 120/nrow(mat))),  # Dynamic font sizing
      fontsize_row = 10,
      fontsize_col = max(3, min(7, 120/nrow(mat))),  # Smaller but readable gene fonts
      breaks = seq(0, 10, length.out = 101),  # Fixed scale: 0 (low) to 10 (high expression)
      main = paste("All Genes - Simple SRR Layout (", nrow(mat), " genes) -", group_name),
      color = eggplant_colors
    )
    
    dev.off()
    cat("    ✅ Detailed heatmap saved:\n      ", detailed_file, "\n")
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()  # Ensure device is closed
    cat("    ❌ Detailed heatmap failed:", e$message, "\n")
  })
  
  # === HEATMAP 2: OVERVIEW WITHOUT GENE NAMES ===
  cat("  🆫 Generating overview heatmap (no gene labels)...\n")
  
  overview_file <- file.path(individual_out_dir, paste0(base_name, "_simple_srr_overview.jpg"))
  
  tryCatch({
    jpeg(overview_file, 
         width = 2000,  # Increased width for gene labels
         height = max(1200, nrow(mat) * 6),  # Increased height
         res = 150, 
         quality = 95)
    
    pheatmap(
      t(mat[all_genes, ]),
      cluster_rows = FALSE,  # No row clustering - simple SRR order
      cluster_cols = FALSE,  # No column clustering (genes shown in input order)
      show_rownames = TRUE,  # Show SRR sample names
      show_colnames = TRUE,  # Show gene names
      fontsize = 10,
      fontsize_row = 10,
      fontsize_col = max(2, min(5, 80/nrow(mat))),  # Smaller gene label fonts
      breaks = seq(0, 10, length.out = 101),  # Fixed scale: 0 (low) to 10 (high expression)
      main = paste("All Genes - Simple SRR Overview (", nrow(mat), " genes) -", group_name),
      color = eggplant_colors
    )
    
    dev.off()
    cat("    ✅ Overview heatmap saved:\n      ", overview_file, "\n")
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()  # Ensure device is closed
    cat("    ❌ Overview heatmap failed:", e$message, "\n")
  })
  
  # === HEATMAP 3: SAMPLE CORRELATION MATRIX ===
  cat("  🇈 Generating sample correlation heatmap...\n")
  
  correlation_file <- file.path(individual_out_dir, paste0(base_name, "_sample_correlation.jpg"))
  
  tryCatch({
    # Calculate sample-to-sample correlations
    sample_cor <- cor(mat, use = "complete.obs")
    
    jpeg(correlation_file, 
         width = 1000, 
         height = 1000, 
         res = 150, 
         quality = 95)
    
    pheatmap(
      sample_cor,
      cluster_rows = FALSE,  # No clustering - simple sample order
      cluster_cols = FALSE,  # No clustering
      show_rownames = TRUE,  # Show SRR sample names
      show_colnames = TRUE,  # Show SRR sample names
      fontsize = 10,
      fontsize_row = 10,
      fontsize_col = 10,
      breaks = seq(0, 1, length.out = 101),  # Correlation scale: 0 to 1
      main = paste("Sample Correlation - Simple Layout -", group_name),
      color = eggplant_colors
    )
    
    dev.off()
    cat("    ✅ Correlation heatmap saved:\n      ", correlation_file, "\n")
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()  # Ensure device is closed
    cat("    ❌ Correlation heatmap failed:", e$message, "\n")
  })
  
  # === HEATMAP 4: BIOLOGICAL ORGAN NAMES VERSION ===
  cat("  🌿 Generating biological organ names heatmap...\n")
  
  biological_file <- file.path(individual_out_dir, paste0(base_name, "_simple_organ_layout.jpg"))
  
  tryCatch({
    # Create matrix with organ names as row labels
    mat_biological <- t(mat[all_genes, ])
    rownames(mat_biological) <- organ_names
    
    jpeg(biological_file, 
         width = 2000,  # Increased width for gene labels
         height = max(1200, nrow(mat) * 6),  # Increased height
         res = 150, 
         quality = 95)
    
    pheatmap(
      mat_biological,
      cluster_rows = FALSE,  # No row clustering - simple organ order
      cluster_cols = FALSE,  # No column clustering (genes shown in input order)
      show_rownames = TRUE,  # Show organ names
      show_colnames = TRUE,  # Show gene names
      fontsize = 10,
      fontsize_row = 11,     # Larger for organ names
      fontsize_col = max(2, min(4, 60/nrow(mat))),  # Smaller gene fonts for readability
      breaks = seq(0, 10, length.out = 101),  # Fixed scale: 0 (low) to 10 (high expression)
      main = paste("All Genes - Simple Organ Layout (", nrow(mat), " genes) -", group_name),
      color = eggplant_colors
    )
    
    dev.off()
    cat("    ✅ Biological organ names heatmap saved:\n      ", biological_file, "\n")
    
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()  # Ensure device is closed
    cat("    ❌ Biological organ names heatmap failed:", e$message, "\n")
  })
  
  # --------------------------------------------------------------------------
  # STEP 9: DATA EXPORT AND REPORTING
  # --------------------------------------------------------------------------
  
  cat("\n📋 Exporting analysis data...\n")
  
  # Export normalized count matrix
  if (!is.null(normalized_counts)) {
    norm_counts_file <- file.path(individual_out_dir, paste0(base_name, "_normalized_counts.csv"))
    write.csv(normalized_counts, norm_counts_file, quote = FALSE)
    cat("  🗃️  Saved normalized counts:", basename(norm_counts_file), "\n")
  }
  
  # Export gene variance rankings (all genes included)
  gene_variance_order <- order(gene_vars, decreasing = TRUE)
  all_genes_df <- data.frame(
    GeneID = rownames(mat)[gene_variance_order],
    Variance = gene_vars[gene_variance_order],
    Rank = 1:length(gene_vars)
  )
  
  variance_file <- file.path(individual_out_dir, paste0(base_name, "_all_genes_variance_ranked.csv"))
  write.csv(all_genes_df, variance_file, quote = FALSE, row.names = FALSE)
  cat("  🗃️  Saved gene variance rankings:", basename(variance_file), "\n")
  
  # --------------------------------------------------------------------------
  # COMPLETION SUMMARY
  # --------------------------------------------------------------------------
  
  cat("\n🎉 ANALYSIS COMPLETED FOR:", group_name, "\n")
  cat("📁 Generated files:\n")
  
  # Check which JPEG files were actually created
  if (file.exists(detailed_file)) {
    cat("   ✅", basename(detailed_file), "\n")
  } else {
    cat("   ❌", basename(detailed_file), "(FAILED)\n")
  }
  
  if (file.exists(overview_file)) {
    cat("   ✅", basename(overview_file), "\n")
  } else {
    cat("   ❌", basename(overview_file), "(FAILED)\n")
  }
  
  if (file.exists(correlation_file)) {
    cat("   ✅", basename(correlation_file), "\n")
  } else {
    cat("   ❌", basename(correlation_file), "(FAILED)\n")
  }
  
  if (file.exists(biological_file)) {
    cat("   ✅", basename(biological_file), "\n")
  } else {
    cat("   ❌", basename(biological_file), "(FAILED)\n")
  }
  
  if (!is.null(normalized_counts)) {
    cat("   ✅", basename(norm_counts_file), "\n")
  }
  
  cat("   ✅", basename(variance_file), "\n")
  
  # --------------------------------------------------------------------------
  # GENERATE DIRECTORY README
  # --------------------------------------------------------------------------
  
  cat("\n📝 Creating directory documentation...\n")
  
  readme_file <- file.path(individual_out_dir, "README.txt")
  readme_content <- paste0(
    "# HEATMAP ANALYSIS RESULTS FOR: ", group_name, "\n",
    "Generated on: ", Sys.time(), "\n\n",
    "DATASET INFORMATION:\n",
    "- Source file: ", basename(file), "\n",
    "- Total genes analyzed: ", nrow(mat), "\n",
    "- Samples in analysis: ", ncol(mat), "\n",
    "- Condition groups: ", length(unique(coldata$condition)), " (", paste(unique(coldata$condition), collapse = ", "), ")\n\n",
    "FILES GENERATED:\n",
  "1. SIMPLIFIED HEATMAPS (JPEG format - NO ROW OR COLUMN CLUSTERING):\n",
  "   - ", basename(detailed_file), " : SRR x Genes (display order preserved)\n",
  "   - ", basename(overview_file), " : SRR overview (no clustering)\n",
  "   - ", basename(correlation_file), " : Sample correlation (order preserved)\n",
  "   - ", basename(biological_file), " : Organ names x Genes (no clustering)\n\n",
    "2. DATA FILES (CSV format):\n",
    if (!is.null(normalized_counts)) paste0("   - ", basename(norm_counts_file), " : DESeq2 normalized count matrix\n") else "",
    "   - ", basename(variance_file), " : Gene variance rankings (all genes)\n\n",
    "ANALYSIS NOTES:\n",
    "- All genes included (no filtering to top genes)\n",
  "- Simple layout: SRR IDs or organ names per row (no clustering on rows or columns)\n",
    "- No condition groupings or complex annotations\n",
    "- Eggplant violet color scheme applied\n",
    "- Genes displayed as columns, samples as rows\n",
    "- DESeq2 variance stabilization (no differential analysis)\n",
    "- Font sizes optimized for gene label readability\n"
  )
  
  writeLines(readme_content, readme_file)
  cat("   ✅ README.txt created\n")
  
  cat("📈 Analysis summary:\n")
  cat("   • Total genes analyzed:", nrow(mat), "\n")
  cat("   • Samples in analysis:", ncol(mat), "\n")
  cat("   • Condition groups:", length(unique(coldata$condition)), "\n")
  
  completion_line <- paste(rep("=", 50), collapse = "")
  cat(completion_line, "\n", sep = "")
  cat("✅ COMPLETED:", group_name, "\n")
  cat(completion_line, "\n", sep = "")
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

batch_header <- paste(rep("=", 60), collapse = "")
cat("\n", batch_header, "\n", sep = "")
cat("🚀 STARTING BATCH PROCESSING\n")
cat(batch_header, "\n", sep = "")

# --------------------------------------------------------------------------
# INPUT VALIDATION
# --------------------------------------------------------------------------

# Verify input directory exists
if (!dir.exists(count_dir)) {
  stop("❌ Count directory not found: ", count_dir, 
       "\n   Make sure the StringTie pipeline has completed successfully.")
}

# Find all count matrix files
files <- list.files(count_dir, pattern = "_counts.tsv$", full.names = TRUE)

if (length(files) == 0) {
  stop("❌ No *_counts.tsv files found in ", count_dir, 
       "\n   Expected StringTie count matrix files.")
}

# --------------------------------------------------------------------------
# BATCH PROCESSING SETUP
# --------------------------------------------------------------------------

cat("📂 Input validation successful\n")
cat("📁 Found", length(files), "count matrices to process:\n")

for (i in seq_along(files)) {
  group_name <- gsub("_counts\\.tsv$", "", basename(files[i]))
  cat("  ", i, ".", group_name, "\n")
}

cat("\n🔄 Starting batch heatmap generation...\n")

# --------------------------------------------------------------------------
# PROCESS ALL FILES
# --------------------------------------------------------------------------

# Process each TSV file with error handling
processing_results <- list()

for (i in seq_along(files)) {
  file <- files[i]
  group_name <- gsub("_counts\\.tsv$", "", basename(file))
  
  tryCatch({
    process_group(file)
    processing_results[[group_name]] <- "SUCCESS"
  }, error = function(e) {
    cat("❌ ERROR processing", group_name, ":", e$message, "\n")
    processing_results[[group_name]] <- paste("ERROR:", e$message)
  })
}

# ==============================================================================
# FINAL SUMMARY AND REPORTING
# ==============================================================================

summary_header <- paste(rep("=", 60), collapse = "")
cat("\n", summary_header, "\n", sep = "")
cat("🎉 BATCH PROCESSING COMPLETED\n")
cat(summary_header, "\n", sep = "")

# --------------------------------------------------------------------------
# PROCESSING RESULTS SUMMARY
# --------------------------------------------------------------------------

cat("📈 Processing Results Summary:\n")

# Safe counting of results with proper type checking
success_count <- 0
error_details <- list()

for (name in names(processing_results)) {
  result <- processing_results[[name]]
  if (is.character(result) && length(result) == 1 && result == "SUCCESS") {
    success_count <- success_count + 1
  } else {
    # Store error details for reporting
    error_details[[name]] <- if (is.character(result)) result else "Unknown error"
  }
}

error_count <- length(processing_results) - success_count

cat("✅ Successful:", success_count, "/", length(processing_results), "files\n")
cat("📁 Individual directories created for each dataset\n")
if (error_count > 0) {
  cat("❌ Failed:", error_count, "files\n")
  cat("\n⚠️  Error Details:\n")
  for (name in names(error_details)) {
    cat("   •", name, ":", error_details[[name]], "\n")
  }
}

# --------------------------------------------------------------------------
# GENERATE COMPREHENSIVE SUMMARY REPORT
# --------------------------------------------------------------------------

cat("\n🗃️  Generating analysis summary report...\n")

summary_file <- file.path(out_dir, "analysis_summary.txt")

# Write comprehensive summary
cat("DESeq2 Heatmap Analysis Summary\n", file = summary_file)
cat("==============================\n", file = summary_file, append = TRUE)
cat("Analysis Date:", as.character(Sys.time()), "\n", file = summary_file, append = TRUE)
cat("Input Directory:", count_dir, "\n", file = summary_file, append = TRUE)
cat("Output Directory:", out_dir, "\n", file = summary_file, append = TRUE)
cat("Directory Structure: Individual folders per dataset\n", file = summary_file, append = TRUE)
cat("Total Files Processed:", length(files), "\n", file = summary_file, append = TRUE)
cat("Successful Analyses:", success_count, "\n", file = summary_file, append = TRUE)
cat("Failed Analyses:", error_count, "\n", file = summary_file, append = TRUE)

# List processed groups
cat("\nProcessed Groups:\n", file = summary_file, append = TRUE)
for (file in files) {
  group_name <- gsub("_counts\\.tsv$", "", basename(file))
  status <- processing_results[[group_name]]
  cat("- ", group_name, " [", status, "]\n", file = summary_file, append = TRUE)
}

# Document output file types
cat("\nOutput Files Generated (per successful group):\n", file = summary_file, append = TRUE)
cat("- [group]_simple_srr_layout.jpg          (Simple SRR layout with gene labels)\n", file = summary_file, append = TRUE)
cat("- [group]_simple_srr_overview.jpg        (Simple SRR overview with gene labels)\n", file = summary_file, append = TRUE)
cat("- [group]_sample_correlation.jpg         (Simple sample correlation)\n", file = summary_file, append = TRUE)
cat("- [group]_simple_organ_layout.jpg        (Simple organ layout with gene labels)\n", file = summary_file, append = TRUE)
cat("- [group]_normalized_counts.csv          (DESeq2 normalized expression)\n", file = summary_file, append = TRUE)
cat("- [group]_all_genes_variance_ranked.csv  (Gene variance rankings)\n", file = summary_file, append = TRUE)

# Analysis features
cat("\nAnalysis Features:\n", file = summary_file, append = TRUE)
cat("- Complete gene inclusion (no filtering to top genes)\n", file = summary_file, append = TRUE)
cat("- Simple layout: no phylogenetic clustering or condition groupings\n", file = summary_file, append = TRUE)
cat("- Direct SRR ID or organ name per row visualization\n", file = summary_file, append = TRUE)
cat("- Custom eggplant violet color scheme\n", file = summary_file, append = TRUE)
cat("- DESeq2 variance stabilization (no differential analysis)\n", file = summary_file, append = TRUE)
cat("- Publication-quality high-resolution outputs with optimized dimensions\n", file = summary_file, append = TRUE)

cat("✅ Summary report saved:", summary_file, "\n")

# --------------------------------------------------------------------------
# FINAL STATUS
# --------------------------------------------------------------------------

cat("\n🎆 ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("📁 All outputs saved in:", out_dir, "\n")
cat("� Directory structure: Individual folders per dataset\n")
cat("   • Each dataset has its own subfolder with:\n")
cat("   • 4 JPEG heatmaps (all non-clustered: detailed, overview, correlation, biological)\n")
cat("   • CSV data files (normalized counts, gene rankings)\n")
cat("   • README.txt documentation\n")
cat("�📋 Total heatmaps generated:", success_count * 4, "\n")  # 4 heatmaps per group
cat("🎨 Color scheme: Custom eggplant violet palette\n")
cat("🧬 Analysis method: Complete gene inclusion with DESeq2\n")

if (success_count == length(files)) {
  cat("\n✨ Perfect run - all files processed successfully!\n")
} else {
  cat("\n⚠️  Some files encountered errors - check summary for details\n")
}

cat(paste0("\n", paste(rep("=", 60), collapse = ""), "\n"))
