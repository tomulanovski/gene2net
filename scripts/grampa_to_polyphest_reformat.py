#!/usr/bin/env python3
"""
Script to process multiple Newick-formatted gene trees from a single file 
and rename nodes to use only the portion after the underscore.
Example: Qlo03p072170_Qlo will be renamed to Qlo

python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/grampa_to_polyphest_reformat.py" input.tre output.tre
"""

import sys
import re

def rename_nodes_in_tree(tree_string):
    """
    Rename all nodes in a Newick tree string to use only the text after the underscore.
    
    Args:
        tree_string (str): Newick-formatted tree string
    
    Returns:
        str: Modified tree string with renamed nodes
    """
    # Regular expression to find node names with underscores
    # This pattern matches anything that's not a parenthesis, comma, colon, or semicolon
    # followed by an underscore and then captures the part after the underscore
    pattern = r'([^(),;:]+)_([^(),;:]+)(?=[:;,)]|$)'
    
    # Replace each match with just the part after the underscore
    renamed_tree = re.sub(pattern, r'\2', tree_string)
    
    return renamed_tree

def process_tree_file(input_file, output_file):
    """
    Process all trees in the input file and write renamed trees to the output file.
    Trees are identified by semicolon (;) endings, as per Newick format.
    
    Args:
        input_file (str): Path to the input file containing Newick trees
        output_file (str): Path for the output file
    
    Returns:
        int: Number of trees processed
    """
    try:
        with open(input_file, 'r') as f_in:
            content = f_in.read()
        
        # Split the content by semicolons to separate trees
        # We add the semicolon back to each tree after splitting
        tree_strings = [t.strip() + ';' for t in content.split(';') if t.strip()]
        
        processed_trees = []
        for tree_string in tree_strings:
            renamed_tree = rename_nodes_in_tree(tree_string)
            processed_trees.append(renamed_tree)
        
        with open(output_file, 'w') as f_out:
            for tree in processed_trees:
                f_out.write(tree + '\n')
        
        return len(processed_trees)
    
    except Exception as e:
        print(f"Error processing trees: {str(e)}")
        return 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python rename_tree_nodes.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    num_trees = process_tree_file(input_file, output_file)
    print(f"Successfully processed {num_trees} trees from {input_file} to {output_file}")

if __name__ == "__main__":
    main()