#!/usr/bin/env python3
"""
simplify_tree.py

Simplifies Newick tree files using ete3:
- Keeps leaf names
- Removes internal node labels
- Removes branch lengths and support values

Usage:
    python simplify_tree.py <input_file> <output_file>
"""

from ete3 import Tree
import sys

def simplify_tree_stream(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            t = Tree(line, format=1)
            for n in t.traverse():
                if not n.is_leaf():
                    n.name = ""  # remove internal labels
                n.dist = 0.0     # set branch lengths to 0 internally
            # format=9 ? clean Newick, no branch lengths or internal labels, keeps leaf names
            outfile.write(t.write(format=9) + "\n")
    print(f"Simplified trees written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python simplify_tree.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    simplify_tree_stream(input_file, output_file)
