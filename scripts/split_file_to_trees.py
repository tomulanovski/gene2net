#!/usr/bin/env python3
"""
Split a file containing multiple Newick trees into separate files.

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/split_file_to_trees.py" input.tre output_dir
"""

import sys
import os

def main():
    if len(sys.argv) != 3:
        print("Usage: python split_trees.py input.tre output_dir")
        sys.exit(1)

    infile, outdir = sys.argv[1], sys.argv[2]

    # Make sure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Read trees (ignoring empty lines)
    with open(infile, "r") as f:
        trees = [line.strip() for line in f if line.strip()]

    # Write each tree to a separate file
    for i, tree in enumerate(trees, start=1):
        outpath = os.path.join(outdir, f"tree{i}.tre")
        with open(outpath, "w") as f:
            f.write(tree + "\n")

    print(f"? Wrote {len(trees)} trees to {outdir}")

if __name__ == "__main__":
    main()