#!/usr/bin/env python3
"""
Gene Tree Taxa Copy Number Adder

This script reads a Newick format gene trees file and creates a new file
where each repeated taxon name is appended with a sequential copy number
(e.g., RS04_1, RS04_2, etc.)

usage : python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/add_copy_number.py" -i input_trees.tre -o output_trees.tre
"""

import argparse
import os
from ete3 import Tree
from collections import defaultdict

def add_copy_numbers_to_tree(tree):
    """
    Modifies a tree by adding sequential copy numbers to repeated taxa.
    
    Args:
        tree (ete3.Tree): The tree to modify
        
    Returns:
        ete3.Tree: The modified tree
    """
    # Count occurrences of each taxon
    taxa_counts = defaultdict(int)
    taxa_dict = {}
    
    # First pass: count occurrences
    for leaf in tree.get_leaves():
        if leaf.name:
            taxa_counts[leaf.name] += 1
    
    # Second pass: rename leaves with copy numbers
    for leaf in tree.get_leaves():
        if leaf.name:
            if taxa_counts[leaf.name] > 1:
                # This taxon appears multiple times
                if leaf.name not in taxa_dict:
                    taxa_dict[leaf.name] = 1
                
                # Add copy number and increment counter
                leaf.name = f"{leaf.name}_{taxa_dict[leaf.name]}"
                taxa_dict[leaf.name.split('_')[0]] += 1
    
    return tree

def process_trees_file(input_file, output_file, tree_mode="multiple"):
    """
    Reads trees from input file, adds copy numbers to taxa, and writes to output file.
    
    Args:
        input_file (str): Path to input Newick trees file
        output_file (str): Path to output Newick trees file
        tree_mode (str): "multiple" for file with multiple trees or "single" for one tree
    
    Returns:
        int: Number of trees processed
    """
    try:
        # Read the entire file content
        with open(input_file, 'r') as file:
            content = file.read()
        
        # Open output file
        with open(output_file, 'w') as out_file:
            tree_count = 0
            
            if tree_mode == "single":
                # Process a single tree
                try:
                    tree = Tree(content)
                    modified_tree = add_copy_numbers_to_tree(tree)
                    out_file.write(modified_tree.write(format=9))
                    tree_count = 1
                except Exception as e:
                    print(f"Error processing tree: {e}")
            else:
                # Process multiple trees
                # Split by ';' to get individual trees
                tree_strings = [ts.strip() for ts in content.split(';') if ts.strip()]
                
                # Process each tree
                for tree_str in tree_strings:
                    try:
                        # Parse tree and add ';' back as needed for Newick format
                        tree = Tree(tree_str + ';')
                        
                        # Add copy numbers
                        modified_tree = add_copy_numbers_to_tree(tree)
                        
                        # Write the modified tree to the output file
                        out_file.write(modified_tree.write(format=9))
                        out_file.write('\n')  # Add newline for readability
                        
                        tree_count += 1
                        
                    except Exception as e:
                        print(f"Error processing tree #{tree_count}: {e}")
                        print(f"Problematic tree string: {tree_str[:50]}...")
        
        return tree_count
        
    except Exception as e:
        print(f"Error opening or reading file {input_file}: {e}")
        return 0

def main():
    """
    Main function to parse arguments and process the trees file.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Add sequential copy numbers to repeated taxa in gene trees."
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Path to the input Newick trees file"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Path to the output Newick trees file (default: input_copynum.tre)"
    )
    
    parser.add_argument(
        "-m", "--mode",
        choices=["multiple", "single"],
        default="multiple",
        help="File mode: 'multiple' for file with multiple trees, 'single' for one tree (default: multiple)"
    )
    
    args = parser.parse_args()
    
    # Set default output file if not specified
    if not args.output:
        base_name = os.path.splitext(args.input)[0]
        args.output = f"{base_name}_copynum.tre"
    
    # Process the file
    print(f"Reading trees from {args.input} in {args.mode} tree mode...")
    tree_count = process_trees_file(args.input, args.output, args.mode)
    
    # Print results
    if tree_count > 0:
        print(f"Successfully processed {tree_count} trees.")
        print(f"Modified trees written to {args.output}")
    else:
        print("No trees were processed. Check the input file or error messages.")

if __name__ == "__main__":
    main()