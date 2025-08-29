#!/usr/bin/env python3
"""
Script to reformat gene trees for GRAMPA using ETE3 Toolkit.
This script adds a gene identifier to each tip label in Newick formatted gene trees.

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/polyphest_to_grampa_reformat.py" <input_file> <output_file>
"""

import sys
import argparse
import re
from collections import defaultdict

def parse_newick_leaves(newick_string):
    """Extract all leaf names from a Newick string using regex."""
    # This regex finds leaf names that are proper identifiers
    # Must start with a letter and can contain letters, numbers, underscores
    pattern = r'([a-zA-Z][a-zA-Z0-9_]*)'
    
    # Find all potential leaf names
    matches = re.findall(pattern, newick_string)
    
    # Filter out obvious non-leaf elements
    leaves = []
    for match in matches:
        # Skip common internal node indicators and very short names
        if match in ['internal', 'node'] or len(match) < 2:
            continue
        leaves.append(match)
    
    return leaves

def replace_leaf_names_in_order(newick_string, leaves):
    """Replace leaf names in the order they appear, adding copy numbers."""
    species_counts = defaultdict(int)
    result = ""
    i = 0
    
    while i < len(newick_string):
        # Check if we're at the start of a potential leaf name
        found_leaf = False
        
        # Try to match the longest possible leaf name first
        for leaf in sorted(set(leaves), key=len, reverse=True):
            if newick_string[i:].startswith(leaf):
                # Check if this is a complete word boundary
                # Before: should be (, or start of string
                before_ok = i == 0 or newick_string[i-1] in "(,"
                # After: should be :, ,, ), ; or end of string
                after_ok = (i + len(leaf) >= len(newick_string) or 
                           newick_string[i + len(leaf)] in ",:);")
                
                if before_ok and after_ok:
                    species_counts[leaf] += 1
                    copy_number = species_counts[leaf]
                    new_name = f"{copy_number}_{leaf}"
                    result += new_name
                    i += len(leaf)
                    found_leaf = True
                    break
        
        if not found_leaf:
            result += newick_string[i]
            i += 1
    
    return result

def reformat_tree(newick_string):
    """Reformat a single tree by adding copy numbers to duplicate species."""
    # Extract all leaf names
    leaves = parse_newick_leaves(newick_string)
    
    # Replace leaf names in order they appear
    result = replace_leaf_names_in_order(newick_string, leaves)
    
    return result

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Reformat gene trees for GRAMPA by adding copy numbers to duplicate species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python grampa_reformatter.py input_trees.tre output_trees.tre
    python grampa_reformatter.py /path/to/gene_trees.txt /path/to/grampa_formatted.tre
        """
    )
    
    parser.add_argument('input_file', 
                       help='Input file containing gene trees (one tree per line in Newick format)')
    parser.add_argument('output_file', 
                       help='Output file for GRAMPA-formatted trees')
    parser.add_argument('-v', '--verbose', 
                       action='store_true', 
                       help='Print detailed progress information')
    
    # Parse arguments
    args = parser.parse_args()
    
    try:
        # Open output file
        with open(args.output_file, 'w') as outfile:
            # Open and process input file
            with open(args.input_file, 'r') as infile:
                for i, line in enumerate(infile, 1):
                    line = line.strip()
                    if not line:
                        continue  # Skip empty lines
                    
                    try:
                        # Reformat the tree
                        reformatted_tree = reformat_tree(line)
                        
                        # Write the reformatted tree to the output file
                        outfile.write(reformatted_tree + '\n')
                        
                        # Print progress
                        if args.verbose and i % 100 == 0:
                            print(f"Processed {i} trees...")
                        elif i % 1000 == 0:
                            print(f"Processed {i} trees...")
                    
                    except Exception as e:
                        print(f"Warning: Could not process tree {i}, skipping. Error: {e}")
                        continue
            
            print(f"Successfully processed {i} trees. Reformatted trees saved to {args.output_file}")
    
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()