# ===============================================
# Complete R Libraries Installation for HeatSeq Project
# ===============================================
# 
# This script installs all required R packages for:
# - b_make_heatmap_of_matrices_v4.R (ComplexHeatmap visualizations)
# - c_make_heatmap_with_CV_of_matrices_v4.R (CV heatmaps)
# - d_make_BarGraph_of_matrices_v4.R (ggplot2 bar graphs)
# ===============================================

# Function to install packages if not available
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package, "...\n")
    
    # Handle Bioconductor packages
    if (package %in% c("ComplexHeatmap")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(package, update = FALSE, ask = FALSE)
    } else {
      # CRAN packages
      install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    }
    
    suppressPackageStartupMessages(library(package, character.only = TRUE))
    cat("✓ Successfully installed and loaded", package, "\n")
  } else {
    cat("✓ Package", package, "is already available\n")
  }
}

# Required packages for all R scripts in the project
required_packages <- c(
  # Bioconductor packages
  "ComplexHeatmap",  # Main heatmap package
  
  # CRAN packages for heatmaps
  "circlize",        # Color mapping
  "RColorBrewer",    # Color palettes
  "dplyr",           # Data manipulation
  "tibble",          # Data frames
  "grid",            # Graphics
  
  # Additional packages for bar graphs
  "ggplot2",         # Grammar of graphics plotting
  "scales",          # Scale functions for ggplot2
  "gridExtra",       # Arranging multiple plots
  
  # Additional utility packages
  "readr",           # Fast reading of delimited files
  "tidyr",           # Data tidying functions
  "stringr"          # String manipulation
)

cat("===============================================\n")
cat("Installing Required Libraries for All R Scripts\n")
cat("===============================================\n")

# Install all required packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("\n===============================================\n")
cat("Testing Installation\n")
cat("===============================================\n")

# Test basic functionality
tryCatch({
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(dplyr)
    library(tibble)
    library(grid)
    library(ggplot2)
    library(scales)
    library(gridExtra)
    library(readr)
    library(tidyr)
    library(stringr)
  })
  
  # Quick test for heatmap functionality
  test_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(test_matrix) <- paste0("Gene_", 1:4)
  colnames(test_matrix) <- paste0("Sample_", 1:5)
  
  # Test color mapping
  colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  # Quick test for ggplot2 functionality
  test_df <- data.frame(x = 1:5, y = rnorm(5))
  test_plot <- ggplot(test_df, aes(x = x, y = y)) + geom_point()
  
  cat("✓ All packages loaded and tested successfully!\n")
  cat("✓ Ready to run heatmap, bar graph, and data processing scripts.\n")
  
}, error = function(e) {
  cat("✗ Error during testing:", e$message, "\n")
  cat("Please check the installation and try again.\n")
})

# From the Bash script. 
# Install packages that need special handling
packages_to_install <- c()

# Check and install ComplexHeatmap if not available
if (!require('ComplexHeatmap', quietly = TRUE)) {
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
        install.packages('BiocManager', repos = 'https://cloud.r-project.org/')
    }
    BiocManager::install('ComplexHeatmap', ask = FALSE, update = FALSE)
    packages_to_install <- c(packages_to_install, 'ComplexHeatmap')
}

# Check and install circlize if not available
if (!require('circlize', quietly = TRUE)) {
    install.packages('circlize', repos = 'https://cloud.r-project.org/')
    packages_to_install <- c(packages_to_install, 'circlize')
}

# Check and install grid (usually comes with base R)
if (!require('grid', quietly = TRUE)) {
    install.packages('grid', repos = 'https://cloud.r-project.org/')
    packages_to_install <- c(packages_to_install, 'grid')
}

if (length(packages_to_install) > 0) {
    cat('Successfully installed:', paste(packages_to_install, collapse = ', '), '\n')
} else {
    cat('All required packages are already available\n')
}

# Test that everything works
cat('Testing package loading...\n')
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(grid)
cat('All packages loaded successfully!\n')
"