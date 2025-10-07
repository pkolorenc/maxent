#!/usr/bin/env python
"""
inp_subs.py - Input file variable substitution utility

Purpose: Substitutes a variable value in an input file.
         Useful for modifying Fortran namelist-style input files.

Usage: inp_subs.py <filename> <variable_name> <new_value>

Example: inp_subs.py entropy.input NT 05

The script reads the input file, finds lines with "variable_name = value" format,
and replaces the value, writing output to <filename>.new
"""

import sys
import os


def substitute_variable(filename, var_name, new_value):
    """
    Read input file and substitute variable value.
    
    Args:
        filename: Path to input file
        var_name: Variable name to search for
        new_value: New value to assign
    
    Returns:
        Number of substitutions made
    """
    if not os.path.isfile(filename):
        print(f"Error: Input file '{filename}' not found", file=sys.stderr)
        sys.exit(1)
    
    output_filename = f"{filename}.new"
    substitutions = 0
    
    try:
        with open(filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line in infile:
                # Check if this line contains the variable name
                if var_name in line and '=' in line:
                    # Split by '=' and clean whitespace
                    parts = [part.strip() for part in line.split('=')]
                    
                    # Check if the left side matches our variable exactly
                    if len(parts) >= 2 and parts[0] == var_name:
                        # Write the substituted line
                        outfile.write(f" {parts[0]} = {new_value}\n")
                        substitutions += 1
                    else:
                        # Variable name appears but doesn't match exactly
                        outfile.write(line)
                else:
                    # Line doesn't contain variable, write as-is
                    outfile.write(line)
    
    except IOError as e:
        print(f"Error: Could not read/write files: {e}", file=sys.stderr)
        sys.exit(1)
    
    return substitutions


def main():
    """Main entry point for the script."""
    if len(sys.argv) != 4:
        print("Usage: inp_subs.py <filename> <variable_name> <new_value>", 
              file=sys.stderr)
        print("\nExample: inp_subs.py entropy.input NT 05", file=sys.stderr)
        sys.exit(1)
    
    filename = sys.argv[1]
    var_name = sys.argv[2]
    new_value = sys.argv[3]
    
    # Perform substitution
    num_subs = substitute_variable(filename, var_name, new_value)
    
    if num_subs == 0:
        print(f"Warning: Variable '{var_name}' not found in {filename}", 
              file=sys.stderr)
    else:
        print(f"Substituted {num_subs} occurrence(s) of '{var_name}' "
              f"with '{new_value}' in {filename}")
        print(f"Output written to {filename}.new")


if __name__ == "__main__":
    main()
