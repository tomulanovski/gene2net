#!/usr/bin/env python3
"""
Apply taxa replacement map to a phylogenetic tree.

This script reads a taxa map file and applies the replacements to a tree file in place.
Used to fix species trees after processing gene trees with substring issues.

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/fix_substrings_to_astral.py" <tree_file> <taxa_map_file>
"""

import sys
import re
import argparse
from pathlib import Path


def read_taxa_map(map_file):
    """Read taxa replacement map from file."""
    replacements = {}
    
    with open(map_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse tab-delimited format
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"WARNING: Skipping malformed line: {line}")
                continue
            
            original, replacement = parts
            replacements[original] = replacement
    
    return replacements


def apply_replacements_to_tree(tree_string, replacements):
    """Apply taxa replacements to a tree string."""
    if not replacements:
        return tree_string
    
    # Sort by length (longest first) to avoid partial replacements
    sorted_replacements = sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True)
    
    result = tree_string
    for original, replacement in sorted_replacements:
        # Replace only complete taxon names
        # Match: word boundary before, non-alphanumeric after (handles hyphens, colons, etc.)
        pattern = rf'\b{re.escape(original)}(?=[^a-zA-Z0-9_]|$)'
        result = re.sub(pattern, replacement, result)
    
    return result


def apply_taxa_map(tree_file, taxa_map_file, verbose=False):
    """Apply taxa map to tree file (in place)."""
    tree_path = Path(tree_file)
    map_path = Path(taxa_map_file)
    
    # Check if files exist
    if not tree_path.exists():
        print(f"ERROR: Tree file not found: {tree_file}")
        return False
    
    if not map_path.exists():
        print(f"ERROR: Taxa map file not found: {taxa_map_file}")
        return False
    
    # Read the taxa map
    if verbose:
        print(f"Reading taxa map from: {taxa_map_file}")
    
    replacements = read_taxa_map(taxa_map_file)
    
    if not replacements:
        print("WARNING: No replacements found in taxa map")
        return False
    
    if verbose:
        print(f"Found {len(replacements)} replacements:")
        for orig, repl in replacements.items():
            print(f"  {orig} -> {repl}")
    
    # Read the tree
    if verbose:
        print(f"Reading tree from: {tree_file}")
    
    with open(tree_file, 'r') as f:
        tree_string = f.read().strip()
    
    # Apply replacements
    if verbose:
        print("Applying replacements...")
    
    modified_tree = apply_replacements_to_tree(tree_string, replacements)
    
    # Check if any changes were made
    if modified_tree == tree_string:
        print("WARNING: No changes were made to the tree")
        if verbose:
            print("  (This might mean the taxa names don't match the map)")
        return False
    
    # Write back to the same file (in place)
    if verbose:
        print(f"Writing modified tree back to: {tree_file}")
    
    with open(tree_file, 'w') as f:
        f.write(modified_tree)
        if not modified_tree.endswith('\n'):
            f.write('\n')
    
    print(f"Successfully applied taxa map to: {tree_file}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Apply taxa replacement map to a phylogenetic tree",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python apply_taxa_map_to_tree.py species_tree.tre taxa_map.txt
        """
    )
    
    parser.add_argument('tree_file', help='Tree file to modify (will be modified in place)')
    parser.add_argument('taxa_map_file', help='Taxa replacement map file (tab-delimited)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print detailed progress information')
    
    args = parser.parse_args()
    
    success = apply_taxa_map(args.tree_file, args.taxa_map_file, args.verbose)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()