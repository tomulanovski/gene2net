#!/usr/bin/env python3
"""
Update GRAMPA grampa_result.tre files from grandma_split history.json.

Takes the (0, 0) iteration's best_mt from history.json (which is GRAMPA's
first iteration result), cleans it up, and saves it as grampa_result.tre.

This replaces the GRAMPA output that was previously extracted from
grampa-scores.txt by postprocess_results.py.

Cleaning:
- Removes * and + suffixes from leaf names
- Removes <label> internal node annotations

Usage:
    # Process all configs, networks, replicates
    python update_grampa_from_grandma_history.py

    # Process specific config(s)
    python update_grampa_from_grandma_history.py --config conf_ils_low_10M conf_ils_high_10M

    # Process specific network
    python update_grampa_from_grandma_history.py --config conf_ils_low_10M --network Bendiksby_2011

    # Dry run
    python update_grampa_from_grandma_history.py --config conf_ils_low_10M --dry-run
"""

import os
import re
import sys
import json
import shutil
import argparse
from pathlib import Path


BASE_DIR = Path("/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")

CONFIGS = [
    "conf_ils_low_10M",
    "conf_ils_medium_10M",
    "conf_ils_high_10M",
    "conf_dup_loss_low_10M",
    "conf_dup_loss_medium_10M",
    "conf_dup_loss_high_10M",
    "conf_dup_loss_low_10M_ne1M",
    "conf_dup_loss_medium_10M_ne1M",
    "conf_dup_loss_high_10M_ne1M",
]

NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019",
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011",
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014",
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014",
]

NUM_REPLICATES = 5


def clean_grampa_tree(tree_str):
    """Clean a GRAMPA best_mt string to produce a standard Newick tree.

    - Remove <...> internal node labels (e.g., <1>, <P*>, <14>)
    - Remove * and + suffixes from names
    """
    # Remove <...> internal node labels
    cleaned = re.sub(r'<[^>]*>', '', tree_str)

    # Remove * and + suffixes from names before Newick delimiters
    cleaned = re.sub(r'([A-Za-z0-9_])([*+]+)(?=[,):;])', r'\1', cleaned)

    # Clean up any double commas or empty nodes
    while ',,' in cleaned:
        cleaned = cleaned.replace(',,', ',')

    # Remove trailing/leading commas in groups: (,A) -> (A), (A,) -> (A)
    cleaned = re.sub(r'\(,', '(', cleaned)
    cleaned = re.sub(r',\)', ')', cleaned)

    return cleaned.strip()


def main():
    parser = argparse.ArgumentParser(
        description='Update GRAMPA grampa_result.tre from grandma_split history.json')
    parser.add_argument('--config', nargs='+', default=None,
                        help='Configuration(s) to process (default: all)')
    parser.add_argument('--network', help='Process only this network')
    parser.add_argument('--dry-run', action='store_true',
                        help='Only show what would be done')
    args = parser.parse_args()

    configs = args.config if args.config else CONFIGS
    networks = [args.network] if args.network else NETWORKS

    updated = 0
    skipped = 0
    missing = 0

    for config in configs:
        print(f"\n{'='*70}")
        print(f"Config: {config}")
        print(f"{'='*70}")

        for network in networks:
            for replicate in range(1, NUM_REPLICATES + 1):
                label = f"{network}/replicate_{replicate}"

                # Source: grandma_split history.json
                history_path = (BASE_DIR / network / "results" / config /
                                "grandma_split" / f"replicate_{replicate}" / "history.json")

                # Target: grampa result tree
                grampa_dir = (BASE_DIR / network / "results" / config /
                              "grampa" / f"replicate_{replicate}")
                grampa_tree_path = grampa_dir / "grampa_result.tre"

                if not history_path.is_file():
                    missing += 1
                    continue

                # Read history.json
                try:
                    with open(history_path) as f:
                        history = json.load(f)
                except Exception as e:
                    print(f"  {label}: Failed to read history.json — {e}")
                    skipped += 1
                    continue

                # Get (0, 0) iteration
                first_iter = history.get('(0, 0)')
                if not first_iter:
                    print(f"  {label}: No (0, 0) iteration in history.json — skipping")
                    skipped += 1
                    continue

                best_mt = first_iter.get('best_mt', '')
                if not best_mt:
                    print(f"  {label}: Empty best_mt in (0, 0) — skipping")
                    skipped += 1
                    continue

                # Clean the tree
                cleaned = clean_grampa_tree(best_mt)

                # Ensure it ends with semicolon
                if not cleaned.endswith(';'):
                    cleaned += ';'

                if args.dry_run:
                    preview = cleaned[:120] + ('...' if len(cleaned) > 120 else '')
                    print(f"  {label}: [DRY RUN] Would write ({len(cleaned)} chars) — {preview}")
                else:
                    # Backup existing grampa_result.tre
                    grampa_dir.mkdir(parents=True, exist_ok=True)
                    if grampa_tree_path.is_file():
                        bak_path = grampa_tree_path.with_suffix('.tre.bak')
                        if not bak_path.exists():
                            shutil.copy2(grampa_tree_path, bak_path)

                    with open(grampa_tree_path, 'w') as f:
                        f.write(cleaned)
                        if not cleaned.endswith('\n'):
                            f.write('\n')

                    print(f"  {label}: Updated grampa_result.tre")

                updated += 1

    print(f"\n{'='*70}")
    action = "Would update" if args.dry_run else "Updated"
    print(f"{action}: {updated}, Skipped: {skipped}, No history.json: {missing}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
