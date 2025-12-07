#!/usr/bin/env python3
"""
Complexity Analysis Module
==========================
Measures and reports runtime and space complexity for alignment algorithms.
Generates graphs and comprehensive reports.

Only basic Python: variables, loops, conditionals, functions.
"""

import time
import os
import sys

# =============================================================================
# STEP 1: RUNTIME MEASUREMENT FUNCTIONS
# =============================================================================

def measure_time(func, *args, **kwargs):
    """
    Measure execution time of a function.
    
    Returns:
        Tuple of (result, elapsed_time_seconds)
    """
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    
    elapsed = end_time - start_time
    return result, elapsed


def measure_memory_usage():
    """
    Get current memory usage in MB.
    Uses /proc/self/status on Linux.
    """
    try:
        with open('/proc/self/status', 'r') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    # VmRSS is resident set size in kB
                    parts = line.split()
                    return int(parts[1]) / 1024  # Convert to MB
    except:
        pass
    return 0


# =============================================================================
# STEP 2: COMPLEXITY TRACKING
# =============================================================================

def create_complexity_tracker():
    """Create a new complexity tracker dictionary."""
    return {
        'measurements': [],
        'algorithm': '',
        'input_sizes': [],
        'runtimes': [],
        'memory_usages': [],
        'operations_count': []
    }


def add_measurement(tracker, input_size, runtime, memory_mb, operations=0, label=""):
    """Add a measurement to the tracker."""
    measurement = {
        'input_size': input_size,
        'runtime': runtime,
        'memory_mb': memory_mb,
        'operations': operations,
        'label': label
    }
    tracker['measurements'].append(measurement)
    tracker['input_sizes'].append(input_size)
    tracker['runtimes'].append(runtime)
    tracker['memory_usages'].append(memory_mb)
    tracker['operations_count'].append(operations)


# =============================================================================
# STEP 3: REPORT GENERATION (TEXT)
# =============================================================================

def generate_text_report(tracker, output_file):
    """Generate a text-based complexity report."""
    lines = []
    lines.append("=" * 70)
    lines.append("COMPLEXITY ANALYSIS REPORT")
    lines.append("Algorithm: " + tracker['algorithm'])
    lines.append("=" * 70)
    lines.append("")
    
    # Summary statistics
    if len(tracker['runtimes']) > 0:
        avg_runtime = sum(tracker['runtimes']) / len(tracker['runtimes'])
        max_runtime = max(tracker['runtimes'])
        min_runtime = min(tracker['runtimes'])
        
        lines.append("RUNTIME STATISTICS:")
        lines.append("  Average: " + str(round(avg_runtime, 4)) + " seconds")
        lines.append("  Maximum: " + str(round(max_runtime, 4)) + " seconds")
        lines.append("  Minimum: " + str(round(min_runtime, 4)) + " seconds")
        lines.append("")
    
    if len(tracker['memory_usages']) > 0:
        avg_memory = sum(tracker['memory_usages']) / len(tracker['memory_usages'])
        max_memory = max(tracker['memory_usages'])
        
        lines.append("MEMORY STATISTICS:")
        lines.append("  Average: " + str(round(avg_memory, 2)) + " MB")
        lines.append("  Maximum: " + str(round(max_memory, 2)) + " MB")
        lines.append("")
    
    # Detailed measurements
    lines.append("DETAILED MEASUREMENTS:")
    lines.append("-" * 70)
    lines.append("{:<15} {:>15} {:>15} {:>15}".format(
        "Input Size", "Runtime (s)", "Memory (MB)", "Operations"))
    lines.append("-" * 70)
    
    for m in tracker['measurements']:
        lines.append("{:<15} {:>15.4f} {:>15.2f} {:>15}".format(
            m['input_size'], m['runtime'], m['memory_mb'], m['operations']))
    
    lines.append("-" * 70)
    lines.append("")
    
    # Complexity estimation
    lines.append("COMPLEXITY ESTIMATION:")
    complexity = estimate_complexity(tracker)
    lines.append("  Time Complexity: " + complexity['time'])
    lines.append("  Space Complexity: " + complexity['space'])
    lines.append("")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    return '\n'.join(lines)


def estimate_complexity(tracker):
    """
    Estimate time and space complexity based on measurements.
    Simple heuristic approach.
    """
    if len(tracker['input_sizes']) < 2:
        return {'time': 'Unknown (insufficient data)', 'space': 'Unknown'}
    
    # Calculate growth ratios
    n_values = tracker['input_sizes']
    t_values = tracker['runtimes']
    m_values = tracker['memory_usages']
    
    # Time complexity estimation
    time_ratios = []
    for i in range(1, len(n_values)):
        if n_values[i-1] > 0 and t_values[i-1] > 0:
            n_ratio = n_values[i] / n_values[i-1]
            t_ratio = t_values[i] / t_values[i-1] if t_values[i-1] > 0 else 1
            if n_ratio > 0:
                time_ratios.append(t_ratio / n_ratio)
    
    avg_time_ratio = sum(time_ratios) / len(time_ratios) if len(time_ratios) > 0 else 1
    
    if avg_time_ratio < 0.5:
        time_complexity = "O(1) or O(log n)"
    elif avg_time_ratio < 1.5:
        time_complexity = "O(n)"
    elif avg_time_ratio < 3:
        time_complexity = "O(n log n)"
    elif avg_time_ratio < 10:
        time_complexity = "O(n^2)"
    else:
        time_complexity = "O(n^2) or higher"
    
    # Space complexity estimation
    space_ratios = []
    for i in range(1, len(n_values)):
        if n_values[i-1] > 0 and m_values[i-1] > 0:
            n_ratio = n_values[i] / n_values[i-1]
            m_ratio = m_values[i] / m_values[i-1] if m_values[i-1] > 0 else 1
            if n_ratio > 0:
                space_ratios.append(m_ratio / n_ratio)
    
    avg_space_ratio = sum(space_ratios) / len(space_ratios) if len(space_ratios) > 0 else 1
    
    if avg_space_ratio < 0.5:
        space_complexity = "O(1)"
    elif avg_space_ratio < 1.5:
        space_complexity = "O(n)"
    else:
        space_complexity = "O(n) or higher"
    
    return {'time': time_complexity, 'space': space_complexity}


# =============================================================================
# STEP 4: ASCII GRAPH GENERATION
# =============================================================================

def generate_ascii_graph(values, labels, title, width=60, height=15):
    """
    Generate an ASCII bar graph.
    
    Args:
        values: List of numeric values
        labels: List of labels for each value
        title: Graph title
        width: Graph width in characters
        height: Graph height in lines
    
    Returns:
        String containing the ASCII graph
    """
    if len(values) == 0:
        return "No data to display"
    
    lines = []
    lines.append("")
    lines.append(title)
    lines.append("=" * len(title))
    lines.append("")
    
    max_val = max(values) if max(values) > 0 else 1
    max_label_len = max(len(str(l)) for l in labels)
    
    # Generate bars
    for i in range(len(values)):
        bar_len = int((values[i] / max_val) * (width - max_label_len - 15))
        bar_len = max(1, bar_len)
        
        label = str(labels[i]).ljust(max_label_len)
        bar = "#" * bar_len
        value_str = str(round(values[i], 4))
        
        lines.append(label + " |" + bar + " " + value_str)
    
    lines.append("")
    return '\n'.join(lines)


def generate_runtime_graph(tracker):
    """Generate ASCII graph for runtime measurements."""
    if len(tracker['measurements']) == 0:
        return "No measurements available"
    
    labels = []
    values = []
    
    for m in tracker['measurements']:
        labels.append(str(m['input_size']))
        values.append(m['runtime'])
    
    return generate_ascii_graph(
        values, labels, 
        "RUNTIME vs INPUT SIZE (seconds)"
    )


def generate_memory_graph(tracker):
    """Generate ASCII graph for memory measurements."""
    if len(tracker['measurements']) == 0:
        return "No measurements available"
    
    labels = []
    values = []
    
    for m in tracker['measurements']:
        labels.append(str(m['input_size']))
        values.append(m['memory_mb'])
    
    return generate_ascii_graph(
        values, labels,
        "MEMORY USAGE vs INPUT SIZE (MB)"
    )


# =============================================================================
# STEP 5: CSV EXPORT
# =============================================================================

def export_to_csv(tracker, output_file):
    """Export measurements to CSV file for external graphing."""
    lines = []
    lines.append("input_size,runtime_sec,memory_mb,operations,label")
    
    for m in tracker['measurements']:
        line = ",".join([
            str(m['input_size']),
            str(m['runtime']),
            str(m['memory_mb']),
            str(m['operations']),
            '"' + m['label'] + '"'
        ])
        lines.append(line)
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))


# =============================================================================
# STEP 6: MATPLOTLIB GRAPH GENERATION (IF AVAILABLE)
# =============================================================================

def generate_matplotlib_graphs(tracker, output_dir):
    """
    Generate PNG graphs using matplotlib if available.
    Includes: runtime, memory, throughput, efficiency, scalability analysis.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import math
        
        has_matplotlib = True
    except ImportError:
        print("matplotlib not available, using ASCII graphs only")
        return False
    
    if len(tracker['measurements']) == 0:
        return False
    
    # Prepare data
    input_sizes = [m['input_size'] for m in tracker['measurements']]
    runtimes = [m['runtime'] for m in tracker['measurements']]
    memories = [m['memory_mb'] for m in tracker['measurements']]
    operations = [m['operations'] for m in tracker['measurements']]
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Runtime complexity graph
    plt.figure(figsize=(10, 6))
    plt.plot(input_sizes, runtimes, '-o', color='#E63946', linewidth=2, markersize=8)
    plt.xlabel('Input Size (reads)')
    plt.ylabel('Runtime (seconds)')
    plt.title('Runtime Complexity: ' + tracker['algorithm'])
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'runtime_complexity.png'), dpi=150)
    plt.close()
    
    # 2. Memory complexity graph
    plt.figure(figsize=(10, 6))
    plt.plot(input_sizes, memories, '-o', color='#2A9D8F', linewidth=2, markersize=8)
    plt.xlabel('Input Size (reads)')
    plt.ylabel('Memory Usage (MB)')
    plt.title('Space Complexity: ' + tracker['algorithm'])
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'memory_complexity.png'), dpi=150)
    plt.close()
    
    # 3. Throughput analysis (reads per second)
    throughput = []
    for i in range(len(input_sizes)):
        if runtimes[i] > 0:
            throughput.append(input_sizes[i] / runtimes[i])
        else:
            throughput.append(0)
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(input_sizes)), throughput, color='#9B5DE5', alpha=0.8)
    plt.xticks(range(len(input_sizes)), [str(s) for s in input_sizes])
    plt.xlabel('Input Size (reads)')
    plt.ylabel('Throughput (reads/second)')
    plt.title('Throughput Analysis: ' + tracker['algorithm'])
    plt.grid(True, alpha=0.3, axis='y')
    for i, v in enumerate(throughput):
        plt.text(i, v + max(throughput)*0.02, str(int(v)), ha='center', fontsize=9)
    plt.savefig(os.path.join(output_dir, 'throughput_analysis.png'), dpi=150)
    plt.close()
    
    # 4. Memory efficiency (reads processed per MB)
    mem_efficiency = []
    for i in range(len(input_sizes)):
        if memories[i] > 0:
            mem_efficiency.append(input_sizes[i] / memories[i])
        else:
            mem_efficiency.append(0)
    
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(input_sizes)), mem_efficiency, color='#F4D03F', alpha=0.8)
    plt.xticks(range(len(input_sizes)), [str(s) for s in input_sizes])
    plt.xlabel('Input Size (reads)')
    plt.ylabel('Memory Efficiency (reads/MB)')
    plt.title('Memory Efficiency: ' + tracker['algorithm'])
    plt.grid(True, alpha=0.3, axis='y')
    plt.savefig(os.path.join(output_dir, 'memory_efficiency.png'), dpi=150)
    plt.close()
    
    # 5. Scalability analysis (log-log plot for complexity estimation)
    if len(input_sizes) >= 2 and min(input_sizes) > 0 and min(runtimes) > 0:
        plt.figure(figsize=(10, 6))
        log_sizes = [math.log10(s) if s > 0 else 0 for s in input_sizes]
        log_times = [math.log10(t) if t > 0 else 0 for t in runtimes]
        plt.plot(log_sizes, log_times, 'b-o', linewidth=2, markersize=8, label='Measured')
        
        # Add reference lines for O(n), O(n log n), O(n^2)
        if len(log_sizes) >= 2:
            x_range = [min(log_sizes), max(log_sizes)]
            # O(n) reference
            y_on = [log_times[0] + (x - log_sizes[0]) for x in x_range]
            plt.plot(x_range, y_on, '--', color='green', alpha=0.5, label='O(n)')
            # O(n^2) reference
            y_on2 = [log_times[0] + 2*(x - log_sizes[0]) for x in x_range]
            plt.plot(x_range, y_on2, '--', color='red', alpha=0.5, label='O(n²)')
        
        plt.xlabel('log₁₀(Input Size)')
        plt.ylabel('log₁₀(Runtime)')
        plt.title('Scalability Analysis (Log-Log): ' + tracker['algorithm'])
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(output_dir, 'scalability_analysis.png'), dpi=150)
        plt.close()
    
    # 6. Runtime per read (amortized cost)
    runtime_per_read = []
    for i in range(len(input_sizes)):
        if input_sizes[i] > 0:
            runtime_per_read.append(runtimes[i] / input_sizes[i] * 1000)  # ms per read
        else:
            runtime_per_read.append(0)
    
    plt.figure(figsize=(10, 6))
    plt.plot(input_sizes, runtime_per_read, 'o-', color='#d62728', linewidth=2, markersize=8)
    plt.xlabel('Input Size (reads)')
    plt.ylabel('Runtime per Read (ms)')
    plt.title('Amortized Cost Analysis: ' + tracker['algorithm'])
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(output_dir, 'amortized_cost.png'), dpi=150)
    plt.close()
    
    # 7. Combined 4-panel dashboard
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Runtime
    axes[0, 0].plot(input_sizes, runtimes, 'b-o', linewidth=2, markersize=6)
    axes[0, 0].set_xlabel('Input Size')
    axes[0, 0].set_ylabel('Runtime (s)')
    axes[0, 0].set_title('Runtime Complexity')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Memory
    axes[0, 1].plot(input_sizes, memories, 'r-o', linewidth=2, markersize=6)
    axes[0, 1].set_xlabel('Input Size')
    axes[0, 1].set_ylabel('Memory (MB)')
    axes[0, 1].set_title('Space Complexity')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Throughput
    axes[1, 0].bar(range(len(input_sizes)), throughput, color='#2ca02c', alpha=0.8)
    axes[1, 0].set_xticks(range(len(input_sizes)))
    axes[1, 0].set_xticklabels([str(s) for s in input_sizes])
    axes[1, 0].set_xlabel('Input Size')
    axes[1, 0].set_ylabel('Reads/sec')
    axes[1, 0].set_title('Throughput')
    axes[1, 0].grid(True, alpha=0.3, axis='y')
    
    # Amortized cost
    axes[1, 1].plot(input_sizes, runtime_per_read, 'o-', color='#d62728', linewidth=2, markersize=6)
    axes[1, 1].set_xlabel('Input Size')
    axes[1, 1].set_ylabel('ms/read')
    axes[1, 1].set_title('Amortized Cost')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.suptitle(tracker['algorithm'] + ' - Performance Dashboard', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'performance_dashboard.png'), dpi=150)
    plt.close()
    
    return True


# =============================================================================
# STEP 7: FULL REPORT GENERATION
# =============================================================================

def generate_full_report(tracker, output_dir):
    """
    Generate complete complexity report with text, CSV, and graphs.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    report_lines = []
    
    # Text report
    text_file = os.path.join(output_dir, 'complexity_report.txt')
    text_report = generate_text_report(tracker, text_file)
    report_lines.append(text_report)
    
    # ASCII graphs
    runtime_graph = generate_runtime_graph(tracker)
    memory_graph = generate_memory_graph(tracker)
    
    # Append graphs to text report
    with open(text_file, 'a') as f:
        f.write("\n" + runtime_graph + "\n")
        f.write("\n" + memory_graph + "\n")
    
    report_lines.append(runtime_graph)
    report_lines.append(memory_graph)
    
    # CSV export
    csv_file = os.path.join(output_dir, 'complexity_data.csv')
    export_to_csv(tracker, csv_file)
    
    # Try to generate matplotlib graphs
    generate_matplotlib_graphs(tracker, output_dir)
    
    print("Reports generated in:", output_dir)
    print("  - complexity_report.txt")
    print("  - complexity_data.csv")
    print("  - *.png (if matplotlib available)")
    
    return '\n'.join(report_lines)


# =============================================================================
# STEP 8: COMBINED COMPARISON GRAPHS FOR ALL ALGORITHMS
# =============================================================================

def generate_combined_comparison(trackers, output_dir):
    """
    Generate combined comparison graphs for multiple algorithms.
    
    Args:
        trackers: List of tracker dictionaries, each with 'algorithm' key
        output_dir: Directory to save combined graphs
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate ASCII combined report
    lines = []
    lines.append("=" * 70)
    lines.append("COMBINED ALGORITHM COMPARISON REPORT")
    lines.append("=" * 70)
    lines.append("")
    
    for tracker in trackers:
        lines.append("-" * 50)
        lines.append(tracker['algorithm'])
        lines.append("-" * 50)
        if len(tracker['runtimes']) > 0:
            lines.append("  Avg Runtime: " + str(round(sum(tracker['runtimes'])/len(tracker['runtimes']), 4)) + "s")
            lines.append("  Max Runtime: " + str(round(max(tracker['runtimes']), 4)) + "s")
        if len(tracker['memory_usages']) > 0:
            lines.append("  Avg Memory: " + str(round(sum(tracker['memory_usages'])/len(tracker['memory_usages']), 2)) + " MB")
        lines.append("")
    
    # Write text report
    text_file = os.path.join(output_dir, 'combined_comparison.txt')
    with open(text_file, 'w') as f:
        f.write('\n'.join(lines))
    
    # Generate matplotlib combined graphs
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Rainbow color palette - distinct hues for algorithms
        colors = [
            '#E63946',  # Red - HISAT2
            '#2A9D8F',  # Teal - Bowtie2
            '#9B5DE5',  # Purple - Salmon
            '#F4A261',  # Orange
            '#00CED1',  # Cyan
            '#FFD700',  # Gold/Yellow
        ]
        
        # Runtime comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        for i, tracker in enumerate(trackers):
            if len(tracker['measurements']) > 0:
                sizes = [m['input_size'] for m in tracker['measurements']]
                times = [m['runtime'] for m in tracker['measurements']]
                ax1.plot(sizes, times, '-o', label=tracker['algorithm'], 
                        color=colors[i % len(colors)], linewidth=2, markersize=6)
        
        ax1.set_xlabel('Input Size (reads)')
        ax1.set_ylabel('Runtime (seconds)')
        ax1.set_title('Runtime Comparison')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Memory comparison
        for i, tracker in enumerate(trackers):
            if len(tracker['measurements']) > 0:
                sizes = [m['input_size'] for m in tracker['measurements']]
                mems = [m['memory_mb'] for m in tracker['measurements']]
                ax2.plot(sizes, mems, '-o', label=tracker['algorithm'],
                        color=colors[i % len(colors)], linewidth=2, markersize=6)
        
        ax2.set_xlabel('Input Size (reads)')
        ax2.set_ylabel('Memory Usage (MB)')
        ax2.set_title('Memory Comparison')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle('Algorithm Comparison: HISAT2 vs Bowtie2 vs Salmon')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'combined_algorithms_comparison.png'), dpi=150)
        plt.close()
        
        # Bar chart comparison for final test size
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        names = [t['algorithm'] for t in trackers if len(t['runtimes']) > 0]
        final_runtimes = [t['runtimes'][-1] for t in trackers if len(t['runtimes']) > 0]
        final_memories = [t['memory_usages'][-1] for t in trackers if len(t['memory_usages']) > 0]
        
        bars1 = ax1.bar(names, final_runtimes, color=colors[:len(names)])
        ax1.set_ylabel('Runtime (seconds)')
        ax1.set_title('Final Test Runtime')
        for bar, val in zip(bars1, final_runtimes):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                    str(round(val, 3))+'s', ha='center', va='bottom', fontsize=9)
        
        bars2 = ax2.bar(names, final_memories, color=colors[:len(names)])
        ax2.set_ylabel('Memory (MB)')
        ax2.set_title('Final Test Memory Usage')
        for bar, val in zip(bars2, final_memories):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    str(round(val, 1))+'MB', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'algorithm_bar_comparison.png'), dpi=150)
        plt.close()
        
        # Throughput comparison
        fig, ax = plt.subplots(figsize=(12, 6))
        x = range(len(trackers[0]['input_sizes'])) if trackers else []
        width = 0.25
        
        for i, tracker in enumerate(trackers):
            if len(tracker['measurements']) > 0:
                throughputs = []
                for j in range(len(tracker['input_sizes'])):
                    if tracker['runtimes'][j] > 0:
                        throughputs.append(tracker['input_sizes'][j] / tracker['runtimes'][j])
                    else:
                        throughputs.append(0)
                offset = (i - len(trackers)/2 + 0.5) * width
                bars = ax.bar([xi + offset for xi in x], throughputs, width, 
                             label=tracker['algorithm'], color=colors[i % len(colors)], alpha=0.8)
        
        ax.set_xlabel('Test Index')
        ax.set_ylabel('Throughput (reads/second)')
        ax.set_title('Throughput Comparison Across Algorithms')
        ax.set_xticks(x)
        ax.set_xticklabels([str(s) for s in trackers[0]['input_sizes']] if trackers else [])
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'throughput_comparison.png'), dpi=150)
        plt.close()
        
        # Combined 6-panel dashboard
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        
        # Runtime line plot
        for i, tracker in enumerate(trackers):
            if len(tracker['measurements']) > 0:
                axes[0, 0].plot(tracker['input_sizes'], tracker['runtimes'], '-o',
                               label=tracker['algorithm'], color=colors[i % len(colors)], linewidth=2)
        axes[0, 0].set_xlabel('Input Size')
        axes[0, 0].set_ylabel('Runtime (s)')
        axes[0, 0].set_title('Runtime Comparison')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Memory line plot
        for i, tracker in enumerate(trackers):
            if len(tracker['measurements']) > 0:
                axes[0, 1].plot(tracker['input_sizes'], tracker['memory_usages'], '-o',
                               label=tracker['algorithm'], color=colors[i % len(colors)], linewidth=2)
        axes[0, 1].set_xlabel('Input Size')
        axes[0, 1].set_ylabel('Memory (MB)')
        axes[0, 1].set_title('Memory Comparison')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Final runtime bar
        axes[0, 2].bar(names, final_runtimes, color=colors[:len(names)])
        axes[0, 2].set_ylabel('Runtime (s)')
        axes[0, 2].set_title('Final Runtime')
        axes[0, 2].grid(True, alpha=0.3, axis='y')
        
        # Final memory bar
        axes[1, 0].bar(names, final_memories, color=colors[:len(names)])
        axes[1, 0].set_ylabel('Memory (MB)')
        axes[1, 0].set_title('Final Memory')
        axes[1, 0].grid(True, alpha=0.3, axis='y')
        
        # Average throughput bar
        avg_throughputs = []
        for tracker in trackers:
            if len(tracker['runtimes']) > 0:
                total_reads = sum(tracker['input_sizes'])
                total_time = sum(tracker['runtimes'])
                avg_throughputs.append(total_reads / total_time if total_time > 0 else 0)
        axes[1, 1].bar(names, avg_throughputs, color=colors[:len(names)])
        axes[1, 1].set_ylabel('Reads/second')
        axes[1, 1].set_title('Average Throughput')
        axes[1, 1].grid(True, alpha=0.3, axis='y')
        
        # Efficiency score (higher = better)
        efficiency_scores = []
        for tracker in trackers:
            if len(tracker['runtimes']) > 0 and len(tracker['memory_usages']) > 0:
                avg_time = sum(tracker['runtimes']) / len(tracker['runtimes'])
                avg_mem = sum(tracker['memory_usages']) / len(tracker['memory_usages'])
                score = 1 / (avg_time * avg_mem + 0.001) * 100
                efficiency_scores.append(score)
            else:
                efficiency_scores.append(0)
        axes[1, 2].bar(names, efficiency_scores, color=colors[:len(names)])
        axes[1, 2].set_ylabel('Efficiency Score')
        axes[1, 2].set_title('Overall Efficiency')
        axes[1, 2].grid(True, alpha=0.3, axis='y')
        
        plt.suptitle('Algorithm Comparison Dashboard: HISAT2 vs Bowtie2 vs Salmon', fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'combined_dashboard.png'), dpi=150)
        plt.close()
        
        print("Combined graphs saved to:", output_dir)
        
    except ImportError:
        print("matplotlib not available for combined graphs")
    
    return '\n'.join(lines)


# =============================================================================
# STEP 9: EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    print("Complexity Analysis Module Test")
    print("=" * 50)
    
    # Create tracker
    tracker = create_complexity_tracker()
    tracker['algorithm'] = 'Test Algorithm'
    
    # Add sample measurements
    add_measurement(tracker, 100, 0.01, 10, 1000, "100 reads")
    add_measurement(tracker, 500, 0.05, 25, 5000, "500 reads")
    add_measurement(tracker, 1000, 0.12, 50, 10000, "1000 reads")
    add_measurement(tracker, 2000, 0.28, 95, 20000, "2000 reads")
    add_measurement(tracker, 5000, 0.85, 220, 50000, "5000 reads")
    
    # Generate report
    report = generate_full_report(tracker, "test_output")
    print(report)
