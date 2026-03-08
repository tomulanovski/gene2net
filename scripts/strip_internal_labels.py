#!/usr/bin/env python3
"""
Strip Internal Node Labels from Newick Trees

Removes internal node labels added by rooting tools (e.g. r98, n110).
These labels appear after ')' and interfere with GRAMPA and GRANDMA_SPLIT.

Leaf labels, branch lengths, and numeric bootstrap values are preserved.

Usage:
    python strip_internal_labels.py input.tre output.tre
    python strip_internal_labels.py input.tre          # writes to input_stripped.tre
    python strip_internal_labels.py input.tre --dry-run
"""

import re
import sys
import argparse


# Matches a label after ')' that starts with a letter.
# Stops before ':', ',', ')', '(', ';' (i.e. branch lengths and tree structure).
# Purely numeric bootstrap values (e.g. )98) are NOT matched and thus preserved.
INTERNAL_LABEL_PATTERN = re.compile(r'\)([A-Za-z][A-Za-z0-9_.-]*)')


def strip_internal_labels(newick: str) -> str:
    return INTERNAL_LABEL_PATTERN.sub(')', newick)


def count_labels(newick: str) -> int:
    return len(INTERNAL_LABEL_PATTERN.findall(newick))


def process_file(input_file: str, output_file: str, dry_run: bool = False) -> bool:
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"ERROR: File not found: {input_file}")
        return False

    trees = [l.strip() for l in lines if l.strip()]
    print(f"Read {len(trees)} tree(s) from {input_file}")

    total_removed = 0
    stripped_trees = []
    for i, tree in enumerate(trees, 1):
        n = count_labels(tree)
        total_removed += n
        stripped_trees.append(strip_internal_labels(tree))
        if n > 0:
            print(f"  Tree {i}: removed {n} internal label(s)")

    print(f"Total internal labels removed: {total_removed}")

    if dry_run:
        print("Dry run - no file written.")
        return True

    try:
        with open(output_file, 'w') as f:
            for tree in stripped_trees:
                f.write(tree + '\n')
        print(f"Written to: {output_file}")
        return True
    except Exception as e:
        print(f"ERROR writing {output_file}: {e}")
        return False


def default_output(input_file: str) -> str:
    for ext in ('.tre', '.nwk', '.nex', '.txt'):
        if input_file.endswith(ext):
            return input_file[:-len(ext)] + '_stripped' + ext
    return input_file + '_stripped'


def main():
    parser = argparse.ArgumentParser(
        description="Strip internal node labels from Newick tree files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python strip_internal_labels.py grampa_trees.tre grampa_trees_stripped.tre
    python strip_internal_labels.py grampa_trees.tre
    python strip_internal_labels.py grampa_trees.tre --dry-run
        """
    )
    parser.add_argument('input_file', help='Input Newick tree file')
    parser.add_argument('output_file', nargs='?', help='Output file (optional)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Report labels found without writing output')
    args = parser.parse_args()

    output_file = args.output_file or default_output(args.input_file)
    success = process_file(args.input_file, output_file, dry_run=args.dry_run)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
