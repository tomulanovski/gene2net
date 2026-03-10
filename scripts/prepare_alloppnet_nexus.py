#!/usr/bin/env python3
"""
Prepare NEXUS files for AlloppNET by renaming taxa from Species@ID format
to Species_A / Species_B format.

For diploids:    Species@ID  -> Species_A
For tetraploids: Species@ID1 -> Species_A  (sorted alphabetically)
                 Species@ID2 -> Species_B

Also outputs taxa_copies.json with copy numbers (1=diploid, 2=tetraploid)
for use with create_taxa_table.py.

Usage:
    python prepare_alloppnet_nexus.py --nexus-dir <dir> --ploidy-file <json> --output-dir <dir>
"""

import os
import json
import shutil
import argparse
from collections import defaultdict


def extract_species(taxon_name):
    """Extract species name from Species@ID format."""
    if '@' in taxon_name:
        return taxon_name.split('@')[0]
    return taxon_name


def parse_nexus_taxa(filepath):
    """Extract all @-format taxon names from the MATRIX section of a NEXUS file."""
    taxa = []
    in_matrix = False

    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            upper = stripped.upper()

            if 'MATRIX' in upper and not stripped.startswith('['):
                in_matrix = True
                continue
            if stripped.startswith(';') or (upper.startswith('END') and in_matrix):
                in_matrix = False
                continue

            if in_matrix and stripped and not stripped.startswith('['):
                parts = stripped.split()
                if parts and '@' in parts[0]:
                    taxa.append(parts[0])

    return taxa


def build_rename_map(taxa_in_file, copy_num_dict):
    """
    Map original taxa names to Species_A / Species_B names.
    Within each species, sort taxa alphabetically before assigning A/B
    for deterministic results across all gene files.
    """
    species_taxa = defaultdict(list)
    for t in taxa_in_file:
        species = extract_species(t)
        if t not in species_taxa[species]:
            species_taxa[species].append(t)

    for sp in species_taxa:
        species_taxa[sp].sort()

    rename_map = {}
    suffixes = ['_A', '_B', '_C', '_D']

    for species, taxon_list in species_taxa.items():
        for i, t in enumerate(taxon_list):
            suffix = suffixes[i] if i < len(suffixes) else f'_{chr(65 + i)}'
            rename_map[t] = f"{species}{suffix}"

    return rename_map


def apply_rename(filepath, rename_map, output_path):
    """Apply renaming to NEXUS file content (longest names first to avoid partial matches)."""
    with open(filepath) as f:
        content = f.read()

    for old, new in sorted(rename_map.items(), key=lambda x: len(x[0]), reverse=True):
        content = content.replace(old, new)

    with open(output_path, 'w') as f:
        f.write(content)


def main():
    parser = argparse.ArgumentParser(
        description='Rename taxa in NEXUS files for AlloppNET (Species@ID -> Species_A/B)'
    )
    parser.add_argument('--nexus-dir', required=True,
                        help='Directory containing input .nex files')
    parser.add_argument('--ploidy-file', required=True,
                        help='JSON with ploidy levels (2=diploid, 4=tetraploid)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for renamed .nex files')
    args = parser.parse_args()

    # Load ploidy levels and convert to copy numbers (2->1, 4->2)
    with open(args.ploidy_file) as f:
        raw_ploidy = json.load(f)
    copy_num_dict = {
        sp: (1 if level <= 2 else level // 2)
        for sp, level in raw_ploidy.items()
    }

    os.makedirs(args.output_dir, exist_ok=True)

    nex_files = sorted(f for f in os.listdir(args.nexus_dir) if f.endswith('.nex'))
    if not nex_files:
        print(f"No .nex files found in {args.nexus_dir}")
        return

    print(f"Processing {len(nex_files)} NEXUS files...")

    no_taxa = 0
    for nf in nex_files:
        inp = os.path.join(args.nexus_dir, nf)
        out = os.path.join(args.output_dir, nf)

        taxa = parse_nexus_taxa(inp)
        if not taxa:
            no_taxa += 1
            shutil.copy(inp, out)
            continue

        rename_map = build_rename_map(taxa, copy_num_dict)
        apply_rename(inp, rename_map, out)

    if no_taxa:
        print(f"  Warning: {no_taxa} files had no @-format taxa (copied as-is)")

    # Write taxa_copies.json for create_taxa_table.py
    copies_path = os.path.join(args.output_dir, 'taxa_copies.json')
    with open(copies_path, 'w') as f:
        json.dump(copy_num_dict, f, indent=4)

    print(f"Done. Renamed files written to: {args.output_dir}")
    print(f"taxa_copies.json written to: {copies_path}")
    print(f"\nCopy numbers used:")
    for sp, cn in sorted(copy_num_dict.items()):
        label = 'diploid' if cn == 1 else 'tetraploid'
        print(f"  {sp}: {cn} ({label})")


if __name__ == '__main__':
    main()
