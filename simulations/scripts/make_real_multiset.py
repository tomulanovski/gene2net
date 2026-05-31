#!/usr/bin/env python3
"""
Generate a REAL (ground-truth) Polyphest multiset from a MUL-tree.

The network .tre files in simulations/networks/ are MUL-trees: each species
appears as a leaf exactly its true-copy-number times. Counting leaf-name
occurrences therefore gives the true ploidy of each species, with no inference.

Output format matches copies_smoothing_with_multiset.py's multi_set.txt:
one species name per line, repeated once per copy, species sorted.

Usage:
    python make_real_multiset.py -i networks/Marcussen_2012.tre -m multi_set_real.txt
"""

import argparse
import sys
from collections import Counter
from ete3 import Tree


def load_mul_tree(path):
    """Load a single MUL-tree, trying a few Newick flavors for robustness."""
    last_err = None
    for fmt in (1, 0, 5):
        try:
            return Tree(path, format=fmt)
        except Exception as e:  # noqa: BLE001 - report after exhausting formats
            last_err = e
    raise SystemExit(f"ERROR: could not parse tree {path}: {last_err}")


def count_copies(tree):
    """Return a Counter of leaf-name -> occurrence count (the true ploidy)."""
    counts = Counter()
    for leaf in tree.get_leaves():
        if leaf.name:
            counts[leaf.name] += 1
    return counts


def write_multiset(counts, output_file):
    """Write each species name `count` times, sorted for stable output."""
    with open(output_file, "w") as f:
        for species in sorted(counts):
            for _ in range(counts[species]):
                f.write(f"{species}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
                        help="MUL-tree Newick file (simulations/networks/<network>.tre)")
    parser.add_argument("-m", "--multiset", required=True,
                        help="Output multi_set_real.txt path")
    args = parser.parse_args()

    tree = load_mul_tree(args.input)
    counts = count_copies(tree)

    if not counts:
        print(f"ERROR: no leaves found in {args.input}", file=sys.stderr)
        sys.exit(1)

    write_multiset(counts, args.multiset)

    polyploids = {n: c for n, c in counts.items() if c > 1}
    total_copies = sum(counts.values())
    print(f"Real multiset written to {args.multiset}")
    print(f"  Species: {len(counts)}  |  Polyploid: {len(polyploids)}  |  "
          f"Total leaves: {total_copies}  |  Max copies: {max(counts.values())}")
    for name, count in sorted(polyploids.items()):
        print(f"  {name}: {count}")


if __name__ == "__main__":
    main()
