#!/usr/bin/env Rscript
# Install WGCNA and dependencies for coexpression analysis
# Usage: Rscript install_wgcna.R

required_packages <- c(
  "WGCNA",
  "dynamicTreeCut",
  "fastcluster",
  "ggplot2",
  "RColorBrewer"
)

cat("Installing WGCNA and dependencies...\n\n")

# Install CRAN packages
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
    } else {
      cat(pkg, "already installed\n")
    }
  }
}

install_if_missing(required_packages)

cat("\nVerifying installations...\n")
success <- TRUE
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "\n")
  } else {
    cat("✗", pkg, "FAILED\n")
    success <- FALSE
  }
}

if (success) {
  cat("\n✓ All packages installed successfully!\n")
  cat("You can now run the WGCNA coexpression analysis.\n")
} else {
  cat("\n✗ Some packages failed to install.\n")
  cat("Try installing manually with: install.packages('WGCNA')\n")
}
