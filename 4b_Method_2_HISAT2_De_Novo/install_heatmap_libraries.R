# ===============================================
# Simplified Heatmap Libraries Installation
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
    
    library(package, character.only = TRUE)
    cat("✓ Successfully installed and loaded", package, "\n")
  } else {
    cat("✓ Package", package, "is already available\n")
  }
}

# Required packages for heatmap generation
required_packages <- c(
  "ComplexHeatmap",  # Main heatmap package
  "circlize",        # Color mapping
  "RColorBrewer",    # Color palettes
  "dplyr",           # Data manipulation
  "tibble",          # Data frames
  "grid"             # Graphics
)

cat("===============================================\n")
cat("Installing Required Heatmap Libraries\n")
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
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(grid)
  
  # Quick test
  test_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(test_matrix) <- paste0("Gene_", 1:4)
  colnames(test_matrix) <- paste0("Sample_", 1:5)
  
  # Test color mapping
  colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  cat("✓ All packages loaded and tested successfully!\n")
  cat("✓ Ready to run heatmap generation scripts.\n")
  
}, error = function(e) {
  cat("✗ Error during testing:", e$message, "\n")
  cat("Please check the installation and try again.\n")
})