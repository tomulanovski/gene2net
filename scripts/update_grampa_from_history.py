#!/usr/bin/env python3
"""
Update GRAMPA network.tre files from grandma_split history.json.

Takes the (0, 0) iteration's best_mt from history.json (which is GRAMPA's
first iteration result), cleans it up, and saves it as the grampa network.tre.

Cleaning:
- Removes * and + suffixes from leaf names
- Removes <label> internal node annotations
- Strips the genus prefix if other methods use species-only names

Usage:
    python update_grampa_from_history.py /path/to/papers [--dry-run] [--dataset DATASET]
"""

import os
import re
import sys
import json
import shutil
import argparse


def clean_grampa_tree(tree_str):
    """Clean a GRAMPA best_mt string to produce a standard Newick tree.

    - Remove * and + suffixes from leaf/node names
    - Remove <label> internal node annotations (e.g., <1>, <P*>, <14>)
    """
    # Remove <...> internal node labels
    cleaned = re.sub(r'<[^>]*>', '', tree_str)

    # Remove * and + suffixes from names (but not from within parentheses/structure)
    # Match word characters followed by * or + before a Newick delimiter
    cleaned = re.sub(r'([A-Za-z0-9_])([*+]+)(?=[,):;])', r'\1', cleaned)

    # Clean up any double commas or empty nodes that might result
    # e.g., (A,,B) -> (A,B)
    while ',,' in cleaned:
        cleaned = cleaned.replace(',,', ',')

    # Remove trailing/leading commas in groups: (,A) -> (A), (A,) -> (A)
    cleaned = re.sub(r'\(,', '(', cleaned)
    cleaned = re.sub(r',\)', ')', cleaned)

    return cleaned.strip()


def main():
    parser = argparse.ArgumentParser(
        description='Update GRAMPA network.tre from grandma_split history.json')
    parser.add_argument('papers_dir', help='Path to papers directory')
    parser.add_argument('--dry-run', action='store_true',
                        help='Only show what would be done')
    parser.add_argument('--dataset', help='Process only this dataset')
    args = parser.parse_args()

    if not os.path.isdir(args.papers_dir):
        print(f"Error: {args.papers_dir} is not a directory")
        sys.exit(1)

    datasets = sorted(d for d in os.listdir(args.papers_dir)
                      if os.path.isdir(os.path.join(args.papers_dir, d))
                      and (not args.dataset or d == args.dataset))

    updated = 0
    skipped = 0

    for dataset in datasets:
        history_path = os.path.join(
            args.papers_dir, dataset, 'networks', 'grandma_split', 'history.json')
        grampa_dir = os.path.join(
            args.papers_dir, dataset, 'networks', 'grampa')
        grampa_tree_path = os.path.join(grampa_dir, 'network.tre')

        if not os.path.isfile(history_path):
            continue

        # Read history.json
        with open(history_path) as f:
            history = json.load(f)

        # Get (0, 0) iteration
        first_iter = history.get('(0, 0)')
        if not first_iter:
            print(f"  {dataset}: No (0, 0) iteration in history.json — skipping")
            skipped += 1
            continue

        best_mt = first_iter.get('best_mt', '')
        if not best_mt:
            print(f"  {dataset}: Empty best_mt in (0, 0) — skipping")
            skipped += 1
            continue

        # Clean the tree
        cleaned = clean_grampa_tree(best_mt)

        # Ensure it ends with semicolon
        if not cleaned.endswith(';'):
            cleaned += ';'

        print(f"\n{'='*60}")
        print(f"Dataset: {dataset}")
        print(f"{'='*60}")
        print(f"  Source: {history_path}")
        print(f"  Target: {grampa_tree_path}")
        print(f"  h_leaves (iteration 0,0): {first_iter.get('h_leaves', [])}")
        print(f"  passed: {first_iter.get('passed', 'N/A')}")

        if args.dry_run:
            print(f"  [DRY RUN] Would write cleaned tree ({len(cleaned)} chars)")
            # Show first 200 chars
            preview = cleaned[:200] + ('...' if len(cleaned) > 200 else '')
            print(f"  Preview: {preview}")
        else:
            # Backup existing grampa network.tre
            os.makedirs(grampa_dir, exist_ok=True)
            if os.path.isfile(grampa_tree_path):
                bak_path = grampa_tree_path + '.bak'
                if not os.path.exists(bak_path):
                    shutil.copy2(grampa_tree_path, bak_path)
                    print(f"  Backed up: {bak_path}")

            with open(grampa_tree_path, 'w') as f:
                f.write(cleaned)
            print(f"  Written: {grampa_tree_path}")

        updated += 1

    print(f"\n{'='*60}")
    action = "Would update" if args.dry_run else "Updated"
    print(f"{action}: {updated} datasets, Skipped: {skipped}")


if __name__ == '__main__':
    main()
