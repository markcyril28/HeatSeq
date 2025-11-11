#!/bin/bash

# ==============================================================================
# TIME METRICS ANALYZER
# ==============================================================================
# Description: Analyze pipeline performance metrics from CSV time logs
# Author: Mark Cyril R. Mercado
# Date: November 2025
# ==============================================================================

set -euo pipefail

TIME_DIR="logs/time_logs"
SPACE_DIR="logs/space_logs"
COMBINED_DIR="logs/space_time_logs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${CYAN}================================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}================================================${NC}"
}

print_section() {
    echo -e "\n${BLUE}>>> $1${NC}"
}

format_time() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local mins=$(( (seconds % 3600) / 60 ))
    local secs=$((seconds % 60))
    printf "%02d:%02d:%02d" $hours $mins $secs
}

format_memory() {
    local kb=$1
    local mb=$((kb / 1024))
    local gb=$((mb / 1024))
    if [ $gb -gt 0 ]; then
        echo "${gb} GB"
    elif [ $mb -gt 0 ]; then
        echo "${mb} MB"
    else
        echo "${kb} KB"
    fi
}

# Find the most recent time metrics file
find_latest_metrics() {
    local latest=$(ls -t "$TIME_DIR"/pipeline_*_time_metrics.csv 2>/dev/null | head -1)
    if [ -z "$latest" ]; then
        echo -e "${RED}Error: No time metrics files found in $TIME_DIR${NC}" >&2
        exit 1
    fi
    echo "$latest"
}

# Main analysis function
analyze_metrics() {
    local csv_file="$1"
    
    if [ ! -f "$csv_file" ]; then
        echo -e "${RED}Error: File not found: $csv_file${NC}" >&2
        exit 1
    fi
    
    print_header "PIPELINE PERFORMANCE ANALYSIS"
    echo -e "File: ${GREEN}$csv_file${NC}\n"
    
    # Count total commands
    local total_cmds=$(tail -n +2 "$csv_file" | wc -l)
    print_section "Summary Statistics"
    echo "Total commands executed: $total_cmds"
    
    # Calculate total elapsed time
    local total_time=$(tail -n +2 "$csv_file" | awk -F',' '{sum+=$3} END {print sum}')
    local total_time_int=${total_time%.*}
    echo "Total elapsed time: $(format_time $total_time_int) (${total_time} seconds)"
    
    # Calculate total CPU time
    local total_cpu=$(tail -n +2 "$csv_file" | awk -F',' '{sum+=$6+$7} END {print sum}')
    echo "Total CPU time: ${total_cpu} seconds"
    
    # Average CPU utilization
    local avg_cpu=$(tail -n +2 "$csv_file" | awk -F',' '{sum+=$4; count++} END {if(count>0) print sum/count; else print 0}')
    printf "Average CPU utilization: %.2f%%\n" "$avg_cpu"
    
    # Peak memory usage
    local max_mem=$(tail -n +2 "$csv_file" | awk -F',' '{if($5>max) max=$5} END {print max}')
    echo "Peak memory usage: $(format_memory ${max_mem%.*})"
    
    # Count failures
    local failures=$(tail -n +2 "$csv_file" | awk -F',' '$8!=0 {count++} END {print count+0}')
    if [ "$failures" -gt 0 ]; then
        echo -e "${RED}Failed commands: $failures${NC}"
    else
        echo -e "${GREEN}Failed commands: 0${NC}"
    fi
    
    # Top 5 slowest commands
    print_section "Top 5 Slowest Commands"
    echo -e "${YELLOW}Rank | Time (sec) | Command${NC}"
    echo "-----------------------------------------------------------"
    tail -n +2 "$csv_file" | sort -t',' -k3 -n -r | head -5 | \
        awk -F',' '{printf "%4d | %10.2f | %.60s...\n", NR, $3, $2}'
    
    # Top 5 memory-intensive commands
    print_section "Top 5 Memory-Intensive Commands"
    echo -e "${YELLOW}Rank | Memory     | Command${NC}"
    echo "-----------------------------------------------------------"
    tail -n +2 "$csv_file" | sort -t',' -k5 -n -r | head -5 | \
        awk -F',' '{
            mem_kb = $5
            if (mem_kb > 1048576) printf "%4d | %7.2f GB | %.60s...\n", NR, mem_kb/1048576, $2
            else if (mem_kb > 1024) printf "%4d | %7.2f MB | %.60s...\n", NR, mem_kb/1024, $2
            else printf "%4d | %7d KB | %.60s...\n", NR, mem_kb, $2
        }'
    
    # Commands with low CPU efficiency (< 50%)
    print_section "Low CPU Efficiency Commands (< 50%)"
    local low_cpu_count=$(tail -n +2 "$csv_file" | awk -F',' '$4<50 {count++} END {print count+0}')
    if [ "$low_cpu_count" -gt 0 ]; then
        echo -e "${YELLOW}CPU% | Time (sec) | Command${NC}"
        echo "-----------------------------------------------------------"
        tail -n +2 "$csv_file" | awk -F',' '$4<50 {printf "%4s%% | %10.2f | %.60s...\n", $4, $3, $2}' | head -10
    else
        echo -e "${GREEN}All commands have good CPU efficiency!${NC}"
    fi
    
    # Failed commands (if any)
    if [ "$failures" -gt 0 ]; then
        print_section "Failed Commands"
        echo -e "${RED}Exit Code | Command${NC}"
        echo "-----------------------------------------------------------"
        tail -n +2 "$csv_file" | awk -F',' '$8!=0 {printf "%9s | %.70s...\n", $8, $2}'
    fi
    
    # Time distribution by command prefix
    print_section "Time Distribution by Tool"
    echo -e "${YELLOW}Tool          | Total Time (sec) | Avg Time (sec) | Count${NC}"
    echo "-----------------------------------------------------------"
    tail -n +2 "$csv_file" | awk -F',' '
        {
            # Extract first word (tool name) from command
            match($2, /[a-zA-Z0-9_-]+/)
            tool = substr($2, RSTART, RLENGTH)
            time[tool] += $3
            count[tool]++
        }
        END {
            for (t in time) {
                printf "%-13s | %16.2f | %14.2f | %5d\n", t, time[t], time[t]/count[t], count[t]
            }
        }
    ' | sort -t'|' -k2 -n -r | head -10
    
    echo ""
}

# Export metrics to simple text summary
export_summary() {
    local csv_file="$1"
    local output="${csv_file%.csv}_summary.txt"
    
    {
        echo "PIPELINE PERFORMANCE SUMMARY"
        echo "Generated: $(date)"
        echo "Source: $csv_file"
        echo "================================"
        echo ""
        
        local total_cmds=$(tail -n +2 "$csv_file" | wc -l)
        echo "Total commands: $total_cmds"
        
        local total_time=$(tail -n +2 "$csv_file" | awk -F',' '{sum+=$3} END {print sum}')
        echo "Total time: ${total_time} seconds"
        
        local max_mem=$(tail -n +2 "$csv_file" | awk -F',' '{if($5>max) max=$5} END {print max}')
        echo "Peak memory: ${max_mem} KB"
        
        echo ""
        echo "Top 10 slowest commands:"
        tail -n +2 "$csv_file" | sort -t',' -k3 -n -r | head -10 | \
            awk -F',' '{printf "  %10.2f sec | %s\n", $3, $2}'
    } > "$output"
    
    echo -e "${GREEN}Summary exported to: $output${NC}"
}

# Main script
main() {
    local csv_file=""
    local export_flag=false
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -f|--file)
                csv_file="$2"
                shift 2
                ;;
            -e|--export)
                export_flag=true
                shift
                ;;
            -h|--help)
                echo "Usage: $0 [-f FILE] [-e] [-h]"
                echo ""
                echo "Options:"
                echo "  -f, --file FILE    Specify CSV file to analyze (default: most recent)"
                echo "  -e, --export       Export summary to text file"
                echo "  -h, --help         Show this help message"
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}"
                exit 1
                ;;
        esac
    done
    
    # Use most recent file if not specified
    if [ -z "$csv_file" ]; then
        csv_file=$(find_latest_metrics)
    fi
    
    # Run analysis
    analyze_metrics "$csv_file"
    
    # Export if requested
    if [ "$export_flag" = true ]; then
        export_summary "$csv_file"
    fi
}

main "$@"
