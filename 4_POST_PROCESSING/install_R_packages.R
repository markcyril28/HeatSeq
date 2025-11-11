#!/usr/bin/env Rscript
# Install necessary CRAN and Bioconductor packages used by the pipeline
# Usage: Rscript install_R_packages.R

required_cran <- c(
  "tidyverse",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "RColorBrewer",
  "circlize",
  "getopt",
  "grid"
)

required_bioc <- c(
  "DESeq2",
  "ComplexHeatmap",
  "ballgown",
  "tximport",
  "tximeta",
  "AnnotationDbi",
  "org.Mm.eg.db"
)

install_if_missing_cran <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"]) ]
  if (length(to_install) > 0) {
    message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  } else {
    message("All CRAN packages already installed.")
  }
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"]) ]
  if (length(to_install) > 0) {
    message("Installing Bioconductor packages: ", paste(to_install, collapse = ", "))
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages already installed.")
  }
}

message("Checking / installing CRAN packages...")
install_if_missing_cran(required_cran)

message("Checking / installing Bioconductor packages...")
install_if_missing_bioc(required_bioc)

message("Done. You can now load required packages in R (e.g. library(ballgown)).")
