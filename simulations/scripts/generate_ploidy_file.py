#!/usr/bin/env python3
"""
Generate a ploidy prior file for GRANDMA_SPLIT from a Polyphest multiset file.

The multiset lists one species name per line; duplicate lines = multiple copies.
Output: one line per unique name with its occurrence count.

If a taxa_map file is provided (original TAB replacement, one per line), species
names are renamed from original → replacement before writing, so the ploidy file
matches the substring-fixed gene tree labels used by GRANDMA_SPLIT.

Usage:
    python generate_ploidy_file.py MULTISET_FILE PLOIDY_FILE
    python generate_ploidy_file.py MULTISET_FILE PLOIDY_FILE --taxa-map taxa_map.txt
"""

import sys
import argparse
from collections import Counter


def load_taxa_map(taxa_map_file):
    """Load original->replacement mapping from taxa_map.txt."""
    mapping = {}
    with open(taxa_map_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) == 2:
                original, replacement = parts
                mapping[original] = replacement
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('multiset_file', help='Polyphest multi_set.txt')
    parser.add_argument('ploidy_file',   help='Output split_ploidies.txt')
    parser.add_argument('--taxa-map',    default=None,
                        help='taxa_map.txt (original TAB replacement); applies substring fix')
    args = parser.parse_args()

    with open(args.multiset_file) as f:
        names = [line.strip() for line in f if line.strip()]

    if not names:
        print(f"ERROR: multiset file is empty: {args.multiset_file}", file=sys.stderr)
        sys.exit(1)

    counts = Counter(names)

    # Apply taxa map rename if provided
    if args.taxa_map:
        mapping = load_taxa_map(args.taxa_map)
        if mapping:
            renamed = {}
            for name, count in counts.items():
                new_name = mapping.get(name, name)
                renamed[new_name] = renamed.get(new_name, 0) + count
            counts = renamed
            print(f"Applied taxa map: {len(mapping)} renames")

    with open(args.ploidy_file, 'w') as f:
        for name, count in sorted(counts.items()):
            f.write(f"{name} {count}\n")

    polyploids = {n: c for n, c in counts.items() if c > 1}
    print(f"Written {args.ploidy_file}: {len(counts)} species, {len(polyploids)} polyploid")
    for name, count in sorted(polyploids.items()):
        print(f"  {name}: {count}")


if __name__ == '__main__':
    main()
