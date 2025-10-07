#!/usr/bin/env python
"""
maxentAverage.py - Average Gamma(E) over a range of M values

Purpose: Computes averaged spectral function Gamma(E) from multiple MaxEnt 
         runs with different numbers of moments M. Also calculates standard 
         deviations and estimates width and lifetime at E=0.

Usage: maxentAverage.py <func.all> <M_low> <M_high>

Example: maxentAverage.py func.all 03 20
         Averages columns corresponding to M=03 through M=20

Input file format:
  Column 0: Energy E
  Columns 1,5,9,...: Gamma values for M=0,1,2,... (every 4th column)

Output:
  For each energy: E average(Gamma) std_dev std_error
  Summary lines at end: Gamma(0), tau (lifetime) with uncertainties
"""

import sys
import math
import os


def parse_arguments():
    """Parse and validate command line arguments."""
    if len(sys.argv) != 4:
        print("Usage: maxentAverage.py <func.all> <M_low> <M_high>", 
              file=sys.stderr)
        print("\nExample: maxentAverage.py func.all 03 20", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    
    try:
        m_low = int(sys.argv[2])
        m_high = int(sys.argv[3])
    except ValueError:
        print("Error: M_low and M_high must be integers", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isfile(filename):
        print(f"Error: File '{filename}' not found", file=sys.stderr)
        sys.exit(1)
    
    if m_low > m_high:
        print("Error: M_low must be <= M_high", file=sys.stderr)
        sys.exit(1)
    
    if m_low < 0:
        print("Error: M values must be non-negative", file=sys.stderr)
        sys.exit(1)
    
    return filename, m_low, m_high


def extract_gamma_values(fields, m_low, m_high):
    """
    Extract Gamma values for specified M range from line fields.
    
    The file has Gamma values at columns 1, 5, 9, 13, ... for M=0, 1, 2, 3, ...
    
    Args:
        fields: List of string fields from data line
        m_low: Lower M index
        m_high: Upper M index
    
    Returns:
        List of Gamma values (as floats)
    """
    values = []
    
    for m in range(m_low, m_high + 1):
        # Gamma for M=m is at column index 1 + m*4
        col_index = 1 + m * 4
        
        if col_index < len(fields):
            try:
                values.append(float(fields[col_index]))
            except ValueError:
                print(f"Warning: Could not parse value at column {col_index}", 
                      file=sys.stderr)
        #else:
        #    print(f"Warning: M={m} (column {col_index}) not available in data", 
        #          file=sys.stderr)
    
    return values


def compute_statistics(values):
    """
    Compute mean, standard deviation, and standard error.
    
    Args:
        values: List of numerical values
    
    Returns:
        Tuple of (mean, std_dev, std_error)
    """
    if not values:
        return 0.0, 0.0, 0.0
    
    n = len(values)
    mean = sum(values) / n
    
    # Compute standard deviation
    variance = sum((x - mean)**2 for x in values) / n
    std_dev = math.sqrt(variance)
    
    # Standard error of the mean
    std_error = std_dev / math.sqrt(n)
    
    return mean, std_dev, std_error


def interpolate_at_zero(e_below, gamma_below, sigma_below, 
                        e_above, gamma_above, sigma_above):
    """
    Linear interpolation to estimate Gamma(0) and uncertainty.
    
    Args:
        e_below: Energy below zero
        gamma_below: Gamma at e_below
        sigma_below: Uncertainty at e_below
        e_above: Energy above zero
        gamma_above: Gamma at e_above
        sigma_above: Uncertainty at e_above
    
    Returns:
        Tuple of (gamma_at_zero, sigma_at_zero)
    """
    if e_above == e_below:
        return gamma_below, sigma_below
    
    # Linear interpolation weight
    de = e_above - e_below
    weight = -e_below / de
    
    gamma_0 = gamma_below + (gamma_above - gamma_below) * weight
    sigma_0 = sigma_below + (sigma_above - sigma_below) * weight
    
    return gamma_0, sigma_0


def main():
    """Main entry point for the script."""
    filename, m_low, m_high = parse_arguments()
    
    # Track points bracketing E=0 for interpolation
    e_below = -1.0
    gamma_below = -1.0
    sigma_below = -1.0
    
    e_above = -1.0
    gamma_above = -1.0
    sigma_above = -1.0
    
    found_bracketing = False
    
    # Process the data file
    try:
        with open(filename, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split()
                
                if len(fields) < 2:
                    print(f"Warning: Line {line_num} has too few fields", 
                          file=sys.stderr)
                    continue
                
                try:
                    energy = float(fields[0])
                except ValueError:
                    print(f"Warning: Could not parse energy at line {line_num}", 
                          file=sys.stderr)
                    continue
                
                # Extract Gamma values for selected M range
                gamma_values = extract_gamma_values(fields, m_low, m_high)
                
                if not gamma_values:
                    print(f"Warning: No valid Gamma values at E={energy}", 
                          file=sys.stderr)
                    continue
                
                # Compute statistics
                mean, std_dev, std_error = compute_statistics(gamma_values)
                
                # Track bracketing points around E=0
                if energy < 0.0:
                    e_below = energy
                    gamma_below = mean
                    sigma_below = std_dev
                elif not found_bracketing and energy >= 0.0:
                    e_above = energy
                    gamma_above = mean
                    sigma_above = std_dev
                    found_bracketing = True
                
                # Output: Energy, Mean, StdDev, StdError
                print(f"{fields[0]}  {mean}  {std_dev}  {std_error}")
    
    except IOError as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # ========================================================================
    # Compute Gamma(0) and lifetime tau by interpolation
    # ========================================================================
    
    if e_below < 0 and e_above >= 0:
        gamma_0, sigma_0 = interpolate_at_zero(
            e_below, gamma_below, sigma_below,
            e_above, gamma_above, sigma_above
        )
        
        # Convert to lifetime (hbar/Gamma, with conversion factor)
        # tau [fs] = 658.21189916 / Gamma [meV]
        HBAR_EV_FS = 658.21189916  # hbar in meV·fs
        
        if gamma_0 > 0:
            tau = HBAR_EV_FS / gamma_0
            # Error propagation: sigma_tau = tau^2 / hbar * sigma_gamma
            sigma_tau = (HBAR_EV_FS / gamma_0**2) * sigma_0
            
            print(f"# Gamma(0) = {gamma_0:.6f} +/- {sigma_0:.6f} meV")
            print(f"# tau = {tau:.6f} +/- {sigma_tau:.6f} fs")
        else:
            print(f"# Gamma(0) = {gamma_0:.6f} +/- {sigma_0:.6f} meV")
            print("# Warning: Cannot compute lifetime (Gamma <= 0)", 
                  file=sys.stderr)
    else:
        print("# Warning: Could not bracket E=0 for interpolation", 
              file=sys.stderr)


if __name__ == "__main__":
    main()
