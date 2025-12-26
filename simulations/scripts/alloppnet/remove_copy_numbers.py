#!/usr/bin/env python3
"""
remove_copy_numbers.py - Remove copy number suffixes from AlloppNET consensus trees

AlloppNET adds _0 and _1 suffixes to tetraploid species for homeolog tracking.
This script removes these suffixes to produce final MUL-trees for comparison.

Usage:
    python remove_copy_numbers.py <input_tree> <output_tree>

Example:
    python remove_copy_numbers.py alloppnet_consensus.tre alloppnet_final.tre
"""

import argparse
import sys
import re
import json
import os
from Bio import Phylo


def load_ploidy_json(ploidy_json_path, verbose=False):
    """
    Load ploidy JSON file and return mapping of species names to ploidy levels.

    Args:
        ploidy_json_path (str): Path to ploidy_level.json file
        verbose (bool): Print progress messages

    Returns:
        dict: Mapping of species name -> ploidy level (2 or 4)
    """
    if not os.path.exists(ploidy_json_path):
        if verbose:
            print(f"WARNING: Ploidy JSON not found: {ploidy_json_path}")
        return {}

    try:
        with open(ploidy_json_path, 'r') as f:
            ploidy_data = json.load(f)

        if verbose:
            polyploid_count = sum(1 for p in ploidy_data.values() if p == 4)
            diploid_count = sum(1 for p in ploidy_data.values() if p == 2)
            print(f"  Loaded ploidy data: {polyploid_count} polyploids, {diploid_count} diploids")

        return ploidy_data
    except Exception as e:
        if verbose:
            print(f"WARNING: Failed to load ploidy JSON: {e}")
        return {}


def detect_ploidy_json_path(input_tree_path, verbose=False):
    """
    Auto-detect ploidy JSON path from input tree file path.

    Path structure:
    - Input: ${BASE_DIR}/${network}/results/${CONFIG}/alloppnet/replicate_${REPLICATE}/alloppnet_consensus.tre
    - Output: ${BASE_DIR}/${network}/processed/${CONFIG}/alloppnet_input/replicate_${REPLICATE}/ploidy_level.json

    Args:
        input_tree_path (str): Path to input tree file
        verbose (bool): Print progress messages

    Returns:
        str: Path to ploidy JSON file, or None if detection fails
    """
    try:
        # Get absolute path
        abs_path = os.path.abspath(input_tree_path)
        dir_path = os.path.dirname(abs_path)

        # Check if path contains expected structure
        if '/results/' in dir_path and '/alloppnet/replicate_' in dir_path:
            # Replace results/.../alloppnet/replicate_N with processed/.../alloppnet_input/replicate_N
            ploidy_path = dir_path.replace('/results/', '/processed/')
            ploidy_path = ploidy_path.replace('/alloppnet/replicate_', '/alloppnet_input/replicate_')
            ploidy_path = os.path.join(ploidy_path, 'ploidy_level.json')

            if os.path.exists(ploidy_path):
                if verbose:
                    print(f"  Auto-detected ploidy JSON: {ploidy_path}")
                return ploidy_path
            elif verbose:
                print(f"  Auto-detection failed: {ploidy_path} does not exist")

        return None
    except Exception as e:
        if verbose:
            print(f"  Auto-detection failed: {e}")
        return None


def remove_annotations(tree):
    """
    Remove all annotations (comments) and branch lengths from tree nodes.
    
    BEAST/TreeAnnotator adds annotations like [&length_range={...}, tti="X", ...]
    This function removes all such annotations and branch lengths to produce a clean
    tree with only topology and species names.
    
    Args:
        tree: Bio.Phylo tree object
        
    Returns:
        Modified tree with annotations and branch lengths removed
    """
    for node in tree.find_clades():
        # Clear all comments (annotations)
        node.comment = None
        # Clear branch lengths (keep only topology)
        node.branch_length = None
    
    return tree


def remove_copy_suffixes(tree, ploidy_data=None, verbose=False):
    """
    Remove copy number suffixes from terminal nodes using ploidy information.

    AlloppNET adds 0/1 suffixes to polyploid species (ploidy=4) for homeolog tracking.
    This function only removes suffixes from species that are polyploids according to
    the ploidy JSON file.

    Args:
        tree: Bio.Phylo tree object
        ploidy_data (dict): Optional mapping of species name -> ploidy level
                           If None, falls back to regex-based removal
        verbose (bool): Print progress messages

    Returns:
        Modified tree with suffixes removed
    """
    # Build set of polyploid species names
    polyploid_species = set()
    if ploidy_data:
        polyploid_species = {name for name, ploidy in ploidy_data.items() if ploidy == 4}
        if verbose:
            print(f"  Found {len(polyploid_species)} polyploid species in ploidy data")
    else:
        if verbose:
            print("  WARNING: No ploidy data available. No suffixes will be removed.")
            print("           Provide --ploidy-json to enable suffix removal.")

    changes_made = []
    for leaf in tree.get_terminals():
        if not leaf.name:
            continue

        original_name = leaf.name

        # Try to remove suffix 0 or 1
        if leaf.name.endswith('0'):
            candidate_name = leaf.name[:-1]  # Remove last '0'
            # Check if candidate is a polyploid species
            if candidate_name in polyploid_species:
                leaf.name = candidate_name
                if verbose and original_name != leaf.name:
                    changes_made.append(f"    {original_name} → {leaf.name}")
        elif leaf.name.endswith('1'):
            candidate_name = leaf.name[:-1]  # Remove last '1'
            # Check if candidate is a polyploid species
            if candidate_name in polyploid_species:
                leaf.name = candidate_name
                if verbose and original_name != leaf.name:
                    changes_made.append(f"    {original_name} → {leaf.name}")

    if verbose and changes_made:
        print(f"  Made {len(changes_made)} suffix removals")
        if len(changes_made) <= 10:
            print("  Example changes:")
            for change in changes_made[:10]:
                print(change)
        else:
            print("  Example changes (showing first 10):")
            for change in changes_made[:10]:
                print(change)

    return tree


def process_tree_file(input_file, output_file, ploidy_json_path=None, verbose=False):
    """
    Process tree file to remove copy number suffixes.

    Args:
        input_file (str): Input tree file (Newick or NEXUS format)
        output_file (str): Output tree file (Newick format)
        ploidy_json_path (str): Optional path to ploidy_level.json file.
                                If None, will attempt auto-detection.
        verbose (bool): Print progress messages
    """
    if verbose:
        print(f"Reading tree from: {input_file}")

    # Detect file format by checking file content
    try:
        with open(input_file, 'r') as f:
            content_start = f.read(1000)  # Read first 1000 chars
            is_nexus = '#NEXUS' in content_start.upper() or 'BEGIN TREES' in content_start.upper()
    except Exception as e:
        print(f"ERROR: Could not read file: {e}")
        sys.exit(1)

    # Read tree based on format
    try:
        if is_nexus:
            # NEXUS format (from TreeAnnotator)
            if verbose:
                print(f"  Detected NEXUS format")
            trees = list(Phylo.parse(input_file, 'nexus'))
            if not trees:
                print(f"ERROR: No trees found in NEXUS file")
                sys.exit(1)
            tree = trees[0]  # Get first (and only) tree
        else:
            # Newick format
            if verbose:
                print(f"  Detected Newick format")
            tree = Phylo.read(input_file, 'newick')
    except Exception as e:
        print(f"ERROR: Failed to read tree file: {e}")
        if is_nexus:
            print(f"       Attempted to parse as NEXUS format")
        else:
            print(f"       Attempted to parse as Newick format")
        sys.exit(1)

    if verbose:
        original_tips = [leaf.name for leaf in tree.get_terminals()]
        print(f"  Original tree has {len(original_tips)} tips")

    # Load ploidy data if available
    ploidy_data = None
    if ploidy_json_path:
        if verbose:
            print(f"  Loading ploidy data from: {ploidy_json_path}")
        ploidy_data = load_ploidy_json(ploidy_json_path, verbose)
    else:
        # Try auto-detection
        if verbose:
            print("  Attempting to auto-detect ploidy JSON path...")
        detected_path = detect_ploidy_json_path(input_file, verbose)
        if detected_path:
            ploidy_data = load_ploidy_json(detected_path, verbose)
        else:
            if verbose:
                print("  WARNING: Could not find ploidy JSON. Using fallback removal method.")
                print("           Specify --ploidy-json to use ploidy-based removal.")

    # Remove suffixes
    if verbose:
        print("  Removing copy number suffixes...")
    tree = remove_copy_suffixes(tree, ploidy_data, verbose)
    
    # Remove annotations and branch lengths (BEAST/TreeAnnotator metadata)
    if verbose:
        print("  Removing annotations and branch lengths (keeping only topology)...")
    tree = remove_annotations(tree)

    if verbose:
        cleaned_tips = [leaf.name for leaf in tree.get_terminals()]
        print(f"  Cleaned tree has {len(cleaned_tips)} tips")

        # Show some examples of changes
        changes = []
        for orig, clean in zip(original_tips[:5], cleaned_tips[:5]):
            if orig != clean:
                changes.append(f"    {orig} → {clean}")

        if changes:
            print("  Example changes:")
            for change in changes:
                print(change)

    # Write output
    if verbose:
        print(f"Writing cleaned tree to: {output_file}")

    try:
        Phylo.write(tree, output_file, 'newick')
    except Exception as e:
        print(f"ERROR: Failed to write tree file: {e}")
        sys.exit(1)

    if verbose:
        print("  ✓ Tree processing complete")


def main():
    parser = argparse.ArgumentParser(
        description='Remove copy number suffixes from AlloppNET consensus trees',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python remove_copy_numbers.py alloppnet_consensus.tre alloppnet_final.tre

  # With verbose output
  python remove_copy_numbers.py -v alloppnet_consensus.tre alloppnet_final.tre

Notes:
  - Removes 0 and 1 suffixes added by AlloppNET for polyploid (ploidy=4) homeologs
  - Uses ploidy_level.json to identify polyploid species (only removes suffixes from polyploids)
  - Input: Newick or NEXUS format (TreeAnnotator outputs NEXUS)
  - Output: Newick format
  - Creates final MUL-tree for comparison with other methods
  - If ploidy JSON is not provided, will attempt auto-detection from input path
"""
    )

    parser.add_argument('input_tree',
                       help='Input tree file (Newick or NEXUS format with copy suffixes)')
    parser.add_argument('output_tree',
                       help='Output tree file (Newick format without copy suffixes)')
    parser.add_argument('--ploidy-json', type=str, default=None,
                       help='Path to ploidy_level.json file. If not provided, will attempt auto-detection.')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print progress messages')

    args = parser.parse_args()

    # Verify input file exists
    import os
    if not os.path.exists(args.input_tree):
        print(f"ERROR: Input file not found: {args.input_tree}")
        sys.exit(1)

    # Process tree
    print("=" * 80)
    print("AlloppNET Tree Post-Processing")
    print("=" * 80)
    print(f"Input:  {args.input_tree}")
    print(f"Output: {args.output_tree}")
    print("=" * 80)
    print()

    process_tree_file(args.input_tree, args.output_tree, args.ploidy_json, args.verbose)

    print()
    print("=" * 80)
    print("Post-processing complete!")
    print("=" * 80)
    print(f"Final tree: {args.output_tree}")
    print("=" * 80)


if __name__ == '__main__':
    main()
