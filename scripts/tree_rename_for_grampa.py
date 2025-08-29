#!/usr/bin/env python3
"""
Newick Tree Node Renamer

This script processes a file containing multiple Newick trees and renames nodes
by adding 'genei_' prefix to each node name, where 'i' is the tree counter.

Usage:
    python newick_renamer.py input_file.nwk output_file.nwk
"""

import sys
import re

def rename_nodes_in_tree(newick_string, tree_number):
    """
    Add 'genei_' prefix to all node names in a Newick tree string,
    where 'i' is the tree number.
    
    Args:
        newick_string (str): A Newick format tree string
        tree_number (int): The current tree number
        
    Returns:
        str: Modified Newick tree with renamed nodes
    """
    prefix = f"gene{tree_number}_"
    
    # This pattern specifically targets node names but avoids branch lengths
    # It looks for names that are:
    # 1. At the start of the string or after a bracket or comma
    # 2. Followed by a colon (branch length), bracket, comma, or end of string
    # 3. Handles quoted names properly
    pattern = r'((?<=\()|(?<=,)|^)([A-Za-z0-9_\-\.]+|\'.+?\')(?=:|,|\)|$)'
    
    # Function to add prefix to each matched node name
    def add_prefix(match):
        # The full match includes the lookbehind which is zero-width
        # Group 2 contains just the node name
        node_name = match.group(2)
        
        # Check if the node name is in quotes
        if node_name.startswith("'") and node_name.endswith("'"):
            # Keep the quotes but add prefix to the name inside
            return match.group(1) + f"'{prefix}{node_name[1:-1]}'"
        else:
            # Add prefix to unquoted name
            return match.group(1) + f"{prefix}{node_name}"
    
    # Replace all node names with prefixed versions
    return re.sub(pattern, add_prefix, newick_string)

def process_tree_file(input_file, output_file):
    """
    Process a file containing multiple Newick trees and create a new file with renamed nodes.
    
    Args:
        input_file (str): Path to the input file with multiple Newick trees
        output_file (str): Path to the output file to write renamed trees
    """
    try:
        # Read the input file
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Split content by semicolons to get individual trees
        # Need to add semicolons back since they're delimiters
        trees = [tree.strip() for tree in content.split(";") if tree.strip()]
        
        # Process each tree
        renamed_trees = []
        for i, tree in enumerate(trees, start=1):
            renamed_tree = rename_nodes_in_tree(tree, i)
            renamed_trees.append(renamed_tree)
        
        # Write all renamed trees to the output file
        with open(output_file, 'w') as f:
            for tree in renamed_trees:
                f.write(tree + ";\n")
            
        print(f"Processed {len(trees)} trees from {input_file}")
        print(f"Output written to {output_file}")
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")

def main():
    """Main function to process command line arguments and handle tree files"""
    if len(sys.argv) != 3:
        print("Usage: python newick_renamer.py input_file.nwk output_file.nwk")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Process the input file
    process_tree_file(input_file, output_file)

if __name__ == "__main__":
    main()