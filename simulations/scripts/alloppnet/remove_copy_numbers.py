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
from Bio import Phylo


def remove_copy_suffixes(tree):
    """
    Remove _0 and _1 suffixes from all terminal nodes.

    AlloppNET adds these to track homeologs in tetraploids:
    - TaxonName_0 (genome A)
    - TaxonName_1 (genome B)

    We strip these to get back to original taxon names.

    Args:
        tree: Bio.Phylo tree object

    Returns:
        Modified tree with suffixes removed
    """
    for leaf in tree.get_terminals():
        if leaf.name:
            # Remove _0 or _1 suffix (only at end of name)
            leaf.name = re.sub(r'_[01]$', '', leaf.name)

    return tree


def process_tree_file(input_file, output_file, verbose=False):
    """
    Process tree file to remove copy number suffixes.

    Args:
        input_file (str): Input tree file (Newick format)
        output_file (str): Output tree file (Newick format)
        verbose (bool): Print progress messages
    """
    if verbose:
        print(f"Reading tree from: {input_file}")

    try:
        tree = Phylo.read(input_file, 'newick')
    except Exception as e:
        print(f"ERROR: Failed to read tree file: {e}")
        sys.exit(1)

    if verbose:
        original_tips = [leaf.name for leaf in tree.get_terminals()]
        print(f"  Original tree has {len(original_tips)} tips")

    # Remove suffixes
    tree = remove_copy_suffixes(tree)

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
  - Removes _0 and _1 suffixes added by AlloppNET for tetraploid homeologs
  - Input and output are Newick format
  - Creates final MUL-tree for comparison with other methods
"""
    )

    parser.add_argument('input_tree',
                       help='Input tree file (Newick format with copy suffixes)')
    parser.add_argument('output_tree',
                       help='Output tree file (Newick format without copy suffixes)')
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

    process_tree_file(args.input_tree, args.output_tree, args.verbose)

    print()
    print("=" * 80)
    print("Post-processing complete!")
    print("=" * 80)
    print(f"Final tree: {args.output_tree}")
    print("=" * 80)


if __name__ == '__main__':
    main()
