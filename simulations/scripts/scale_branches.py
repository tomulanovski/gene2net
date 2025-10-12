#!/usr/bin/env python3
"""
Scale branch lengths in a Nexus tree file and output SimPhy-compatible Nexus.

Usage:
    python scale_branches.py input.nex output.nex scale_factor
"""

import sys
from io import StringIO
from Bio import Phylo

def scale_tree(tree, factor):
    """Multiply all branch lengths by factor."""
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length *= factor
    return tree

def write_simphy_nexus(tree, output_file):
    """Write tree in SimPhy Nexus format with no root branch length."""
    handle = StringIO()
    Phylo.write([tree], handle, "newick")
    newick_str = handle.getvalue().strip().replace("\n", "")

    # Ensure it ends with a semicolon
    if not newick_str.endswith(";"):
        newick_str += ";"

    # Remove root branch length if present (SimPhy requirement)
    newick_str = newick_str.replace("):0.0;", ");")
    newick_str = newick_str.replace("):0.00000;", ");")

    with open(output_file, "w") as out:
        out.write("#NEXUS\n")
        out.write("begin trees;\n")
        out.write(f"tree 1={newick_str}\n")
        out.write("end;\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python scale_branches.py input.nex output.nex scale_factor")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    factor = float(sys.argv[3])

    # Read first tree from input
    trees = list(Phylo.parse(input_file, "nexus"))
    if not trees:
        print("No trees found in input file.")
        sys.exit(1)

    tree = trees[0]
    scaled_tree = scale_tree(tree, factor)
    write_simphy_nexus(scaled_tree, output_file)

if __name__ == "__main__":
    main()
