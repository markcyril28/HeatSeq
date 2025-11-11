#!/bin/bash

# ==============================================================================
# QUICK DEMO: Four-Version Logging System
# ==============================================================================
# Demonstrates the four logging outputs: full logs, time, space, and combined
# ==============================================================================

echo "==================================================="
echo "  FOUR-VERSION LOGGING SYSTEM - QUICK DEMO"
echo "==================================================="
echo ""

# Check if directories exist
if [ ! -d "logs/log_files" ] || [ ! -d "logs/time_logs" ] || [ ! -d "logs/space_logs" ] || [ ! -d "logs/space_time_logs" ]; then
    echo "Creating logging directories..."
    mkdir -p logs/log_files logs/time_logs logs/space_logs logs/space_time_logs
fi

# Check for existing logs
echo "1. Checking for existing logs..."
LOG_COUNT=$(ls logs/log_files/pipeline_*.log 2>/dev/null | wc -l)
TIME_COUNT=$(ls logs/time_logs/pipeline_*.csv 2>/dev/null | wc -l)
SPACE_COUNT=$(ls logs/space_logs/pipeline_*.csv 2>/dev/null | wc -l)
COMBINED_COUNT=$(ls logs/space_time_logs/pipeline_*.csv 2>/dev/null | wc -l)

echo "   - Full logs found: $LOG_COUNT"
echo "   - Time metrics found: $TIME_COUNT"
echo "   - Space metrics found: $SPACE_COUNT"
echo "   - Combined metrics found: $COMBINED_COUNT"
echo ""

if [ $TIME_COUNT -gt 0 ]; then
    echo "2. Most recent logs:"
    LATEST_TIME=$(ls -t logs/time_logs/pipeline_*_time_metrics.csv 2>/dev/null | head -1)
    LATEST_SPACE=$(ls -t logs/space_logs/pipeline_*_space_metrics.csv 2>/dev/null | head -1)
    LATEST_COMBINED=$(ls -t logs/space_time_logs/pipeline_*_combined_metrics.csv 2>/dev/null | head -1)
    echo "   Time: $LATEST_TIME"
    echo "   Space: $LATEST_SPACE"
    echo "   Combined: $LATEST_COMBINED"
    echo ""
    
    if [ -f "$LATEST_TIME" ]; then
        echo "3. Time metrics preview (first 5 entries):"
        echo "   -------------------------------------------"
        head -6 "$LATEST_TIME" | column -t -s','
        echo ""
    fi
    
    if [ -f "$LATEST_SPACE" ]; then
        echo "4. Space metrics preview (first 5 entries):"
        echo "   -------------------------------------------"
        head -6 "$LATEST_SPACE" | column -t -s','
        echo ""
    fi
    
    if [ -f "$LATEST_COMBINED" ]; then
        echo "5. Combined metrics preview (first 5 entries):"
        echo "   -------------------------------------------"
        head -6 "$LATEST_COMBINED" | column -t -s','
        echo ""
    fi
else
    echo "2. No logs found yet."
    echo "   Run the pipeline to generate all four log versions"
    echo ""
fi

