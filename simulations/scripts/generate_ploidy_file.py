#!/usr/bin/env python3
"""
Generate a ploidy prior file for GRANDMA_SPLIT from a Polyphest multiset file.

The multiset lists one species name per line; duplicate lines = multiple copies.
Output: one line per unique name with its occurrence count.

Usage:
    python generate_ploidy_file.py MULTISET_FILE PLOIDY_FILE
    python generate_ploidy_file.py multi_set.txt split_ploidies.txt
"""

import sys
from collections import Counter


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} MULTISET_FILE PLOIDY_FILE", file=sys.stderr)
        sys.exit(1)

    multiset_file = sys.argv[1]
    ploidy_file   = sys.argv[2]

    with open(multiset_file) as f:
        names = [line.strip() for line in f if line.strip()]

    if not names:
        print(f"ERROR: multiset file is empty: {multiset_file}", file=sys.stderr)
        sys.exit(1)

    counts = Counter(names)

    with open(ploidy_file, 'w') as f:
        for name, count in sorted(counts.items()):
            f.write(f"{name} {count}\n")

    polyploids = {n: c for n, c in counts.items() if c > 1}
    print(f"Written {ploidy_file}: {len(counts)} species, {len(polyploids)} polyploid")
    for name, count in sorted(polyploids.items()):
        print(f"  {name}: {count}")


if __name__ == '__main__':
    main()
