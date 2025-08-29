#!/usr/bin/env python3
"""
Simple script to resolve polytomies in gene trees
Takes a file with multiple Newick trees and outputs resolved trees

python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/resolve_polytomies.py" input file out file
"""

import sys
from ete3 import Tree

def resolve_polytomies(input_file, output_file):
    """
    Read trees from input file, resolve polytomies, write to output file
    Each tree is on a separate line
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for i, line in enumerate(infile, 1):
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                
                # Try parsing with different formats
                tree_processed = False
                for fmt in [1, 0]:
                    try:
                        # Parse tree
                        tree = Tree(line, format=fmt)
                        # Resolve polytomies
                        tree.resolve_polytomy(recursive=True)
                        # Write resolved tree
                        outfile.write(tree.write(format=fmt) + '\n')
                        tree_processed = True
                        print(f"Successfully processed tree {i} with format {fmt}")
                        break
                        
                    except Exception as e:
                        if fmt == 0:  # Last format attempt failed
                            print(f"Error processing tree {i}: {e}")
                            print(f"Tree string (first 100 chars): {line[:100]}...")
                
                if not tree_processed:
                    print(f"Failed to process tree {i} with any format")
                    
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied accessing files")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python resolve_polytomies.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    resolve_polytomies(input_file, output_file)
    print(f"Resolved trees written to {output_file}")