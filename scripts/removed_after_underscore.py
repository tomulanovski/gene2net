#!/usr/bin/env python3
"""
Newick Tree Taxon Name Cleaner

This script processes Newick format tree files and removes underscores and everything after them
in taxon names. It can handle multiple trees in a single file.

Usage:
    python newick_cleaner.py input_file output_file
"""

import sys
import re


def clean_newick_tree(newick_str):
    """
    Clean taxon names in a Newick tree string with special handling for 'sect_nov_X' patterns:
    - Keep 'sect_nov_A' and 'sect_nov_B' intact
    - For longer names like 'sect_nov_B_C_1', keep up to the second underscore ('sect_nov_B')
    - For all other taxa, remove underscore and everything after it
    
    Args:
        newick_str (str): A tree in Newick format
        
    Returns:
        str: Cleaned Newick tree string
    """
    # Function to handle individual taxon name replacements
    def replacement(match):
        full_name = match.group(0)
        
        # Special case for sect_nov patterns
        if full_name.startswith('sect_nov_'):
            parts = full_name.split('_')
            if len(parts) <= 3:  # sect_nov_A or sect_nov_B - keep as is
                return full_name
            else:  # sect_nov_B_C_1 - keep up to the second underscore
                return f"sect_nov_{parts[2]}"
        else:
            # Regular case - remove first underscore and everything after
            return full_name.split('_')[0]
    
    # Pattern to match taxon names (any alphanumeric sequence that might include underscores)
    pattern = r'[A-Za-z0-9\.\-]+(?:_[^,:\(\)]+)*'
    
    # Apply the replacement function to each match
    cleaned_newick = re.sub(pattern, replacement, newick_str)
    
    return cleaned_newick


def process_newick_file(input_file, output_file):
    """
    Process a file containing one or multiple Newick trees.
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
    """
    try:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()
        
        cleaned_lines = []
        for line in lines:
            if line.strip():  # Skip empty lines
                cleaned_line = clean_newick_tree(line)
                cleaned_lines.append(cleaned_line)
        
        with open(output_file, 'w') as outfile:
            outfile.writelines(cleaned_lines)
            
        print(f"Successfully processed {len(cleaned_lines)} tree(s) from '{input_file}' to '{output_file}'")
            
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python newick_cleaner.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_newick_file(input_file, output_file)