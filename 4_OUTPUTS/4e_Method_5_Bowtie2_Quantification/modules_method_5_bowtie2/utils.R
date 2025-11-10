# --- 
# Title: "Utility Functions for Method 5"
# --- 

# --- Logging ---
log_message <- function(message, level = "INFO") {
  cat(paste0("[", Sys.time(), "][", level, "] ", message, "\n"))
}
