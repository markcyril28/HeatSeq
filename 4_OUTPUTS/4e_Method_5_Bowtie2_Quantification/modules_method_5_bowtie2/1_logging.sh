#!/bin/bash

# ===============================================
# LOGGING SYSTEM FOR METHOD 5
# Based on modules/utils.sh
# ===============================================

# Logging configuration
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-9_logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_full_log.log}"

# Timestamp function
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# Core logging functions
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; }
log_error() { log ERROR "$@"; }
log_step() { log INFO "=============== $* ==============="; }

# Setup logging with tee
setup_logging() {
    local clear_logs="${1:-false}"  # Default to keeping logs if not specified
    
    mkdir -p "$LOG_DIR"
    
    # Clear previous logs if requested
    if [ "$clear_logs" = true ]; then
        if [ -d "$LOG_DIR" ] && [ "$(ls -A "$LOG_DIR" 2>/dev/null)" ]; then
            log_info "Clearing previous log files from $LOG_DIR"
            rm -f "$LOG_DIR"/*.log
            log_info "Previous logs cleared"
        fi
    else
        log_info "Keeping previous log files (CLEAR_LOGS_ON_RUN=false)"
    fi
    
    exec > >(tee -a "$LOG_FILE") 2>&1
    log_info "Logging to: $LOG_FILE"
}

# Error handling
trap 'log_error "Command failed at line $LINENO: ${BASH_COMMAND}"; exit 1' ERR
trap 'log_info "Script finished. Log: $LOG_FILE"' EXIT
