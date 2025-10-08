# Test script to check if required libraries are available
# and install them if needed

# Function to install packages if not available
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package, "...\n")
    if (package %in% c("ComplexHeatmap", "circlize")) {
      # These are Bioconductor packages
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package, update = FALSE, ask = FALSE)
    } else {
      # CRAN packages
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
    cat("Successfully installed and loaded", package, "\n")
  } else {
    cat("Package", package, "is already available\n")
  }
}

# List of required packages
required_packages <- c("ComplexHeatmap", "circlize", "RColorBrewer", "dplyr", "tibble", "grid")

cat("Checking and installing required packages...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("All required packages are now available!\n")

# Test basic functionality
cat("\nTesting basic ComplexHeatmap functionality...\n")
library(ComplexHeatmap)
library(circlize)

# Create a simple test matrix
test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(test_matrix) <- paste0("Gene_", 1:10)
colnames(test_matrix) <- paste0("Sample_", 1:10)

# Test color mapping
test_colors <- colorRamp2(
  seq(min(test_matrix), max(test_matrix), length = 5),
  c("#FFFFFF", "#CE93D8", "#8E24AA", "#4A148C", "#2F1B69")
)

cat("ComplexHeatmap is working correctly!\n")
cat("You can now run the main heatmap script.\n")