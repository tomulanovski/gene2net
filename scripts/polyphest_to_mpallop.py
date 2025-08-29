#!/usr/bin/env python3
"""
Newick to NEXUS Converter (Semicolon Delimiter)
This script converts a file containing Newick-formatted gene trees separated by semicolons
to a NEXUS file format with a TREES block and a PHYLONET block.
Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/polyphest_to_mpallop.py" gene_trees.txt output.nex
"""
import sys
import re

def create_nexus_file(input_file, output_file):
    """Convert Newick trees to NEXUS format, separating trees by semicolons."""
    try:
        # Read the entire input file
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Remove all whitespace for easier processing
        content = content.replace('\n', ' ').replace('\r', '')
        
        # Find all trees by matching patterns ending with semicolons
        # This regex matches content between parentheses ending with semicolons
        tree_pattern = re.compile(r'(\([^;]*;)')
        trees = tree_pattern.findall(content)
        
        # Clean up any extra whitespace in trees
        trees = [tree.strip() for tree in trees if tree.strip()]
        
        if not trees:
            print("Warning: No valid Newick trees found in the input file.")
            return False
        
        # Write the NEXUS file
        with open(output_file, 'w') as out:
            # Write NEXUS header
            out.write("#NEXUS\n")
            
            # Write TREES block
            out.write("BEGIN TREES;\n")
            for i, tree in enumerate(trees):
                out.write(f"Tree geneTree{i+1} = {tree}\n")
            out.write("END;\n")
            
            # Add PHYLONET block
            out.write("\nBEGIN PHYLONET;\n")
            out.write("InferNetwork_MP_Allopp (all) 1 -x 1 -pl 12 -di;\n")
            out.write("END;\n")
        
        print(f"Successfully converted {len(trees)} trees to NEXUS format in {output_file}")
        
        # Display the first few lines of the output file
        with open(output_file, 'r') as f:
            preview_lines = f.readlines()
            preview = "".join(preview_lines[:min(10, len(preview_lines))])
            print("\nPreview of the NEXUS file:")
            print(preview + "...\n")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return False
    except Exception as e:
        print(f"Error: {str(e)}")
        return False
    
    return True

def main():
    """Main function to handle command-line arguments."""
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = create_nexus_file(input_file, output_file)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()