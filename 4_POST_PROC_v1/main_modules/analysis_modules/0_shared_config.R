#!/usr/bin/env Rscript

# ===============================================
# SHARED CONFIGURATION FOR POST-PROCESSING ANALYSES
# ===============================================
# Centralized configuration to avoid duplication across scripts
# This file is sourced by all analysis modules

# ===============================================
# GLOBAL CONFIGURATION (Override via environment variables)
# ===============================================

# Method being processed (set by bash wrapper)
CURRENT_METHOD <- Sys.getenv("CURRENT_METHOD", unset = "M5_RSEM_Bowtie2")

# Master reference genome/transcriptome
MASTER_REFERENCE <- Sys.getenv("MASTER_REFERENCE", unset = "Eggplant_V4.1_transcripts.function")

# GPU acceleration flag
ENABLE_GPU <- as.logical(Sys.getenv("ENABLE_GPU", unset = "FALSE"))

# Thread count
THREADS <- as.integer(Sys.getenv("THREADS", unset = "8"))

# Available RAM (GB) - used to determine if memory-intensive operations are safe
AVAILABLE_RAM_GB <- as.integer(Sys.getenv("AVAILABLE_RAM_GB", unset = "24"))

# Memory-aware settings
# With 24GB+ RAM: prioritize accuracy over memory conservation
HIGH_MEMORY_MODE <- AVAILABLE_RAM_GB >= 16

# ===============================================
# GPU DETECTION AND CONFIGURATION
# ===============================================
# Check for GPU availability and compatible packages
# Falls back to CPU if GPU unavailable or packages not installed
# Note: R torch 0.16.x only supports CUDA 11.6-11.8
# For CUDA 12.x systems, GPU ops fall back to CPU

GPU_AVAILABLE <- FALSE
GPU_BACKEND <- "cpu"  # "cpu", "cuda", or "torch"

# Get CUDA version - prefer conda CUDA over system
get_cuda_version <- function() {
  # First check for conda CUDA environment variable
  conda_cuda <- Sys.getenv("TORCH_CUDA_VERSION", unset = "")
  if (nzchar(conda_cuda)) {
    return(list(version = conda_cuda, source = "conda"))
  }
  
  # Check nvcc first (more accurate for installed toolkit)
  nvcc_version <- tryCatch({
    output <- system("nvcc --version 2>/dev/null", intern = TRUE)
    ver_line <- grep("release", output, value = TRUE)
    if (length(ver_line) > 0) {
      version <- gsub(".*release ([0-9]+\\.[0-9]+).*", "\\1", ver_line[1])
      return(list(version = version, source = "nvcc"))
    }
    NULL
  }, error = function(e) NULL)
  
  if (!is.null(nvcc_version)) return(nvcc_version)
  
  # Fallback: detect from nvidia-smi (driver's max supported version)
  tryCatch({
    cuda_output <- system("nvidia-smi 2>/dev/null", intern = TRUE)
    cuda_line <- grep("CUDA Version", cuda_output, value = TRUE)
    if (length(cuda_line) > 0) {
      version <- gsub(".*CUDA Version: ([0-9]+\\.[0-9]+).*", "\\1", cuda_line[1])
      return(list(version = version, source = "nvidia-smi"))
    }
    return(list(version = NULL, source = "none"))
  }, error = function(e) list(version = NULL, source = "none"))
}

# Detect GPU and compatible packages
detect_gpu <- function() {
  if (!ENABLE_GPU) {
    return(list(available = FALSE, backend = "cpu", message = "GPU disabled by configuration"))
  }
  
  # Check for GPU hardware via nvidia-smi
  gpu_available <- tryCatch({
    result <- system("nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null", intern = TRUE)
    length(result) > 0 && !grepl("error|fail", result[1], ignore.case = TRUE)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  
  if (!gpu_available) {
    return(list(available = FALSE, backend = "cpu", message = "No NVIDIA GPU detected"))
  }
  
  # Check CUDA version
  cuda_info <- get_cuda_version()
  cuda_version <- if (!is.null(cuda_info$version)) as.numeric(cuda_info$version) else 0
  
  # R torch 0.16.x only supports CUDA 11.6-11.8
  torch_supported <- cuda_version >= 11.6 && cuda_version < 12.0
  
  if (!is.null(cuda_info$version)) {
    message("[GPU] CUDA ", cuda_info$version, " detected (", cuda_info$source, ")")
    if (!torch_supported) {
      message("[GPU] Note: R torch 0.16.x requires CUDA 11.6-11.8, found ", cuda_info$version)
      message("[GPU] GPU matrix operations will use CPU. This does not affect pipeline results.")
    }
  }
  
  # Check for torch (preferred for matrix operations)
  if (requireNamespace("torch", quietly = TRUE)) {
    # Try to check CUDA availability without triggering the version check error
    cuda_works <- tryCatch({
      # This will error if lantern isn't installed or CUDA version unsupported
      torch::cuda_is_available()
    }, error = function(e) {
      # Lantern not installed or CUDA version check failed
      FALSE
    })
    
    if (cuda_works) {
      return(list(available = TRUE, backend = "torch", 
                  message = paste0("GPU available via torch (", torch::cuda_device_count(), " device(s))")))
    } else if (torch_supported) {
      message("[GPU] torch package found but CUDA backend not installed")
      message("[GPU] Run: Rscript -e 'torch::install_torch(type=\"cuda\")'")
    }
  }
  
  # Check for gpuR (alternative for matrix operations)
  if (requireNamespace("gpuR", quietly = TRUE)) {
    tryCatch({
      gpuR::detectGPUs()
      return(list(available = TRUE, backend = "gpuR", message = "GPU available via gpuR"))
    }, error = function(e) NULL)
  }
  
  return(list(available = FALSE, backend = "cpu", 
              message = "GPU detected but R GPU acceleration not available (using CPU - this is fine)"))
}

# Initialize GPU detection at load time
.gpu_info <- detect_gpu()
GPU_AVAILABLE <- .gpu_info$available
GPU_BACKEND <- .gpu_info$backend

# Log GPU status (quieter - only log if GPU is used or explicitly requested)
if (GPU_AVAILABLE) {
  message("[GPU] ", .gpu_info$message)
} else if (ENABLE_GPU && interactive()) {
  message("[GPU] ", .gpu_info$message)
}

# ===============================================
# DIRECTORY CONSTANTS
# ===============================================

# Base directories (relative to method folder)
# Different methods have different quantification output structures:
#   M1/M2 (HISAT2+StringTie): 5_stringtie_WD/
#   M4 (Salmon): 5_Salmon_Quant_WD/
#   M5 (RSEM): 5_RSEM_Quant_WD/
MATRICES_DIR <- "6_matrices"
CONSOLIDATED_BASE_DIR <- "7_Figures_Outputs"

# Gene groups directory - use environment variable if set, otherwise compute from script location
# The bash wrapper exports GENE_GROUPS_DIR as absolute path
GENE_GROUPS_DIR <- Sys.getenv("GENE_GROUPS_DIR", unset = "")
if (GENE_GROUPS_DIR == "") {
  # Fallback: compute relative to ANALYSIS_MODULES_DIR
  ANALYSIS_MODULES_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", unset = ".")
  GENE_GROUPS_DIR <- file.path(dirname(ANALYSIS_MODULES_DIR), "gene_groups_csv")
}

# Output subdirectories
OUTPUT_SUBDIRS <- list(
  MATRIX_CREATION = "0_Matrix_Creation",
  BASIC_HEATMAP = "I_Basic_Heatmap",
  CV_HEATMAP = "II_Heatmap_with_CV",
  BAR_GRAPH = "III_Bar_Graphs",
  WGCNA = "IV_Coexpression_WGCNA",
  DEA = "V_Differential_Expression",
  GSEA = "VI_Gene_Set_Enrichment",
  DIM_REDUCTION = "VII_Dimensionality_Reduction",
  CORRELATION = "VIII_Sample_Correlation",
  TISSUE_SPEC = "IX_Tissue_Specificity"
)

# ===============================================
# ANALYSIS PARAMETERS
# ===============================================

# Count and gene type options
# NOTE: Available count types depend on the quantification method:
#   - Salmon (M4): NumReads, TPM (length-corrected, bias-corrected)
#   - RSEM (M5): expected_count, TPM, FPKM
#   - StringTie (M1/M2): coverage, FPKM, TPM
# For DESeq2: use expected_count/NumReads (raw counts)
# For visualization: TPM is preferred for cross-sample comparison
COUNT_TYPES <- c(
  #"coverage",
  "expected_count",
  #"fpkm",
  "tpm"
)

# Gene name display options
# These correspond to column names in the gene_groups_csv/*.csv files
# The processing engine will apply label transformations based on selection
GENE_TYPES <- c(
  #"Gene_ID",       # Column 1: Original gene identifier (e.g., SMEL4.1_01g005840)
  "Shortened_Name" # Column 2: Human-readable short name (e.g., SmelDMP01.840)
)

# Sample label display options  
# These correspond to column names in the SRR_csv/*.csv files
# The processing engine will apply label transformations based on selection
LABEL_TYPES <- c(
  #"SRR_ID", # Column 1: SRA accession number (e.g., SRR3884686)
  "Organ"   # Column 2: Tissue/organ description (e.g., Flower_Buds)
)

PROCESSING_LEVELS <- c(
  "gene_level"
  #,"isoform_level"
)

# Minimum sample/gene thresholds for various analyses
MIN_SAMPLES <- 3      # Minimum samples for correlation/clustering
MIN_GENES_DEA <- 10   # Minimum genes for differential expression
MIN_GENES_WGCNA <- 20 # Minimum genes for WGCNA network

# ===============================================
# METHOD TYPE DETECTION
# ===============================================

# Detect quantification method type from CURRENT_METHOD
get_method_type <- function(method = CURRENT_METHOD) {
  if (grepl("HISAT2|StringTie|M1|M2", method, ignore.case = TRUE)) {
    return("stringtie")
  } else if (grepl("Salmon|M4", method, ignore.case = TRUE)) {
    return("salmon")
  } else if (grepl("RSEM|M5", method, ignore.case = TRUE)) {
    return("rsem")
  } else if (grepl("STAR|M3", method, ignore.case = TRUE)) {
    return("star")
  }
  return("unknown")
}

# Get method-specific quant directory
get_quant_dir <- function(method = CURRENT_METHOD) {
  method_type <- get_method_type(method)
  switch(method_type,
    "stringtie" = "5_stringtie_WD",
    "salmon" = "5_Salmon_Quant_WD",
    "rsem" = "5_RSEM_Quant_WD",
    "star" = "5_STAR_WD",
    "5_quant_WD"  # default
  )
}

# Normalization schemes
# Available options:
#   - "raw": No transformation (raw counts)
#   - "count_type_normalized": Log2(counts + 1) transformation
#   - "deseq2_normalized": DESeq2-style median-of-ratios + log2
#   - "zscore": Z-score normalization (after log2 transformation)
#   - "zscore_scaled_to_ten": Z-score scaled to 0-10 range
#   - "cpm": Counts Per Million (library-size normalized) + log2
# NOTE: For DEA, use raw counts with DESeq2 internal normalization.
#       For visualization (heatmaps, PCA), use log2 or VST-transformed data.
NORM_SCHEMES <- c(
  "count_type_normalized",
  #"raw",
  #"deseq2_normalized",
  #"zscore",
  "zscore_scaled_to_ten"
  #,"cpm"
)

# ===============================================
# SAMPLE LABELS LOADING FROM CSV
# ===============================================
# Load sample labels dynamically from SRR CSV files
# CSV format: SRR_ID,Organ,Notes

SRR_CSV_DIR <- Sys.getenv("SRR_CSV_DIR", unset = "")
if (SRR_CSV_DIR == "") {
  ANALYSIS_MODULES_DIR_TMP <- Sys.getenv("ANALYSIS_MODULES_DIR", unset = ".")
  SRR_CSV_DIR <- file.path(dirname(ANALYSIS_MODULES_DIR_TMP), "SRR_csv")
}

# Load sample labels from all CSV files in SRR_csv directory
load_sample_labels_from_csv <- function(srr_csv_dir = SRR_CSV_DIR) {
  labels <- c()
  if (!dir.exists(srr_csv_dir)) return(labels)
  
  csv_files <- list.files(srr_csv_dir, pattern = "\\.csv$", full.names = TRUE)
  for (csv_file in csv_files) {
    tryCatch({
      df <- read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE, comment.char = "#")
      if ("SRR_ID" %in% colnames(df) && "Organ" %in% colnames(df)) {
        # Filter out empty or invalid SRR_IDs
        df <- df[!is.na(df$SRR_ID) & nzchar(trimws(df$SRR_ID)), ]
        new_labels <- setNames(df$Organ, df$SRR_ID)
        labels <- c(labels, new_labels[!names(new_labels) %in% names(labels)])
      }
    }, error = function(e) NULL)
  }
  return(labels)
}

# Initialize SAMPLE_LABELS from CSV files
# Uses first column (SRR_ID) and second column (Organ) from SRR_csv/*.csv files
SAMPLE_LABELS <- load_sample_labels_from_csv()

if (length(SAMPLE_LABELS) == 0) {
  warning("No sample labels loaded from CSV files in: ", SRR_CSV_DIR,
          "\nEnsure CSV files have SRR_ID and Organ columns.")
}

SAMPLE_IDS <- names(SAMPLE_LABELS)

# ===============================================
# RUNTIME CONFIGURATION LOADING
# ===============================================

read_config_file <- function(file_path, default_value, is_boolean = FALSE) {
  if (!file.exists(file_path)) return(default_value)
  value <- trimws(readLines(file_path, warn = FALSE))
  if (is_boolean) {
    return(tolower(value[1]) == "true")
  } else {
    return(value[nzchar(value)])
  }
}

load_runtime_config <- function(method_modules_dir = ".") {
  gene_groups <- read_config_file(
    file.path(method_modules_dir, ".gene_groups_temp.txt"),
    default_value = c("SmelDMPs", "SmelGRF-GIFs")
  )
  
  master_reference <- read_config_file(
    file.path(method_modules_dir, ".master_reference_temp.txt"),
    default_value = MASTER_REFERENCE
  )
  if (length(master_reference) > 1) master_reference <- master_reference[1]
  
  overwrite <- read_config_file(
    file.path(method_modules_dir, ".overwrite_temp.txt"),
    default_value = TRUE,
    is_boolean = TRUE
  )
  
  list(
    gene_groups = gene_groups,
    master_reference = master_reference,
    overwrite_existing = overwrite
  )
}

# ===============================================
# UTILITY FUNCTIONS
# ===============================================

print_separator <- function(char = "=", width = 60) {
  cat("\n", paste(rep(char, width), collapse = ""), "\n")
}

print_config_summary <- function(title, config) {
  print_separator()
  cat(title, "\n")
  print_separator()
  cat("\nConfiguration:\n")
  cat("  • Method:", CURRENT_METHOD, "\n")
  cat("  • Master Reference:", config$master_reference, "\n")
  cat("  • Overwrite existing:", config$overwrite_existing, "\n")
  cat("  • Gene groups:", paste(config$gene_groups, collapse = ", "), "\n")
  cat("  • Threads:", THREADS, "\n")
  # GPU status
  if (GPU_AVAILABLE) {
    cat("  • GPU:", "ENABLED (", GPU_BACKEND, ")\n", sep = "")
  } else if (ENABLE_GPU) {
    cat("  • GPU: REQUESTED but unavailable (using CPU)\n")
  } else {
    cat("  • GPU: disabled\n")
  }
  cat("\n")
}

print_summary <- function(successful, total, skipped = NULL) {
  print_separator()
  cat("SUMMARY:", successful, "/", total, "items generated")
  if (!is.null(skipped) && skipped > 0) {
    cat(" (", skipped, " skipped)\n", sep = "")
  } else {
    cat("\n")
  }
  print_separator()
}

ensure_output_dir <- function(dir_path, clean = FALSE) {
  if (clean && dir.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)
  }
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

build_input_path <- function(gene_group, processing_level, count_type, gene_type, 
                             matrices_dir = MATRICES_DIR, master_ref = MASTER_REFERENCE) {
  if (gene_group == master_ref) {
    file.path(matrices_dir, master_ref, processing_level,
              paste0(master_ref, "_", count_type, "_", gene_type,
                     "_from_", master_ref, "_", processing_level, ".tsv"))
  } else {
    file.path(matrices_dir, master_ref, processing_level, gene_group,
              paste0(gene_group, "_", count_type, "_", gene_type,
                     "_from_", master_ref, "_", processing_level, ".tsv"))
  }
}

build_title_base <- function(gene_group, count_type, gene_type, label_type, 
                             processing_level, norm_scheme, master_ref = MASTER_REFERENCE) {
  paste0(gene_group, "_", count_type, "_", gene_type, "_",
         label_type, "_from_", master_ref, "_", processing_level, "_", norm_scheme)
}

validate_and_read_matrix <- function(input_file, min_rows = 2) {
  if (!file.exists(input_file)) {
    return(list(success = FALSE, reason = "file not found"))
  }
  matrix_data <- read_count_matrix(input_file)
  if (is.null(matrix_data)) {
    return(list(success = FALSE, reason = "failed to read"))
  }
  if (nrow(matrix_data) < min_rows) {
    return(list(success = FALSE, reason = paste0("need ≥", min_rows, " rows")))
  }
  list(success = TRUE, data = matrix_data, n_genes = nrow(matrix_data))
}

should_skip_existing <- function(output_path, overwrite) {
  !overwrite && file.exists(output_path)
}
