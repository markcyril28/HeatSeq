# ---
# Title: "Utility Functions"
# ---

# --- Logging ---
log_message <- function(message, level = "INFO") {
  cat(paste0("[ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ] [ ", level, " ] ", message, "\n"))
  # Append to log file
  log_dir <- "logs_HISAT2_Ref_Guided"
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  log_file <- file.path(log_dir, "run_all.log")
  cat(paste0("[ ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ] [ ", level, " ] ", message, "\n"), file = log_file, append = TRUE)
}

# --- File Handling ---
read_sample_metadata <- function(metadata_file) {
  if (!file.exists(metadata_file)) {
    log_message(paste("Sample metadata file not found:", metadata_file), level = "ERROR")
    return(NULL)
  }
  
  metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Validate required columns
  if (!all(c("sample", "condition") %in% colnames(metadata))) {
    log_message("Sample metadata must contain 'sample' and 'condition' columns", level = "ERROR")
    return(NULL)
  }
  
  log_message(paste("Loaded metadata for", nrow(metadata), "samples"))
  return(metadata)
}

# --- Data Manipulation ---
# Add data manipulation functions here if needed
