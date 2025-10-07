#!/usr/bin/env python
"""
get_max_dmu.py - Extract maximum moment difference from fort.9006

Purpose: Analyzes moment difference data (dmu) from MaxEnt output to find 
         the maximum value and compute the average. The dmu quantity represents
         the difference between input spectral moments and those recovered from 
         the MaxEnt reconstructed function.

Usage: get_max_dmu.py <fort.9006_file>

Output format: <iteration_of_max> <max_value> <average_value>

The fort.9006 file is expected to have space-separated columns:
  Column 0: iteration index
  Column 3: dmu (moment difference) value
"""

import sys
import os


def analyze_dmu(filename):
    """
    Read dmu data and find maximum and average values.
    
    dmu represents the difference between input spectral moments and those
    recovered from the MaxEnt reconstructed spectral function.
    
    Args:
        filename: Path to fort.9006 file
    
    Returns:
        Tuple of (max_iteration, max_value, average_value)
    """
    if not os.path.isfile(filename):
        print(f"Error: File '{filename}' not found", file=sys.stderr)
        sys.exit(1)
    
    max_value = -1.0
    max_iteration = -1
    total_sum = 0.0
    count = 0
    
    try:
        with open(filename, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                # Skip empty lines
                line = line.strip()
                if not line:
                    continue
                
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                try:
                    fields = line.split()
                    
                    # Check if we have enough columns
                    if len(fields) < 4:
                        print(f"Warning: Line {line_num} has insufficient columns, "
                              f"skipping", file=sys.stderr)
                        continue
                    
                    iteration = int(fields[0])
                    dmu = float(fields[3])
                    
                    # Update statistics
                    total_sum += dmu
                    count += 1
                    
                    # Track maximum
                    if dmu > max_value:
                        max_value = dmu
                        max_iteration = iteration
                
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line {line_num}: {e}", 
                          file=sys.stderr)
                    continue
    
    except IOError as e:
        print(f"Error: Could not read file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if count == 0:
        print("Error: No valid data found in file", file=sys.stderr)
        sys.exit(1)
    
    average = total_sum / count
    
    return max_iteration, max_value, average


def main():
    """Main entry point for the script."""
    if len(sys.argv) != 2:
        print("Usage: get_max_dmu.py <fort.9006_file>", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    
    # Analyze the file
    max_iter, max_val, avg_val = analyze_dmu(filename)
    
    # Format output (scientific notation)
    max_str = f"{max_val:.2e}"
    avg_str = f"{avg_val:.2e}"
    
    # Output: iteration max_value average_value
    print(f"{max_iter}  {max_str}  {avg_str}")


if __name__ == "__main__":
    main()
