#!/usr/bin/env python3
"""
A general-purpose script to remove accession numbers from node names in Newick tree files.
This script handles multiple trees in a single file and preserves all other parts of the node names.

Accession patterns recognized:
- Standard GenBank/ENA/DDBJ: 1-2 letters followed by digits (e.g., AB123456, DQ274092)
- Extended patterns: up to 4 letters + digits (e.g., ABCD123456)
- Can handle accessions regardless of their position in the name
- Preserves taxonomy information, location data, and variant notations

Dependencies:
- ete3 (install with: pip install ete3)

Usage:
    python remove_accessions_multi.py input_trees.newick -o output_trees.newick
    python remove_accessions_multi.py *.tree --batch
"""

import re
import argparse
from ete3 import Tree

def read_multiple_trees(input_file):
    """
    Read multiple trees from a file.
    Handles both single trees and files with multiple trees.
    """
    trees = []
    
    try:
        with open(input_file, 'r') as f:
            content = f.read().strip()
        
        # Split by semicolons and filter empty strings
        tree_strings = [t.strip() + ';' for t in content.split(';') if t.strip()]
        
        # Parse each tree string
        for i, tree_str in enumerate(tree_strings):
            try:
                tree = Tree(tree_str, format=1)  # Try with branch lengths
                trees.append(tree)
            except:
                try:
                    tree = Tree(tree_str, format=0)  # Try without branch lengths
                    trees.append(tree)
                except:
                    print(f"Warning: Could not parse tree {i+1} in file {input_file}")
                    continue
        
        if trees:
            return trees
            
    except Exception as e:
        # If that fails, try reading as single tree
        try:
            tree = Tree(input_file, format=1)
            return [tree]
        except:
            try:
                tree = Tree(input_file, format=0)
                return [tree]
            except Exception as e2:
                print(f"Error reading tree file {input_file}: {e2}")
                return []
    
    return trees

def remove_accessions_from_name(name, patterns):
    """
    Remove accession numbers from a single node name using multiple patterns.
    
    Args:
        name: The node name to clean
        patterns: List of compiled regex patterns to match accessions
    
    Returns:
        Cleaned name with accessions removed
    """
    if not name:
        return name
    
    original_name = name
    
    # Apply each pattern
    for pattern in patterns:
        name = pattern.sub('', name)
    
    # Clean up any potential double underscores resulting from removal
    name = re.sub(r'__+', '_', name)
    
    # Remove leading underscore if present
    if name.startswith('_'):
        name = name[1:]
    
    # Remove trailing underscore if present
    if name.endswith('_'):
        name = name[:-1]
    
    return name

def get_accession_patterns(custom_patterns=None):
    """
    Define regex patterns to match different accession number formats.
    
    Returns:
        List of compiled regex patterns
    """
    default_patterns = [
        # Standard GenBank/ENA/DDBJ: 1-2 letters + digits (e.g., AB123456, DQ274092)
        r'_([A-Z]{1,2}\d+)(?=_|$)',
        
        # Extended: 3-4 letters + digits (e.g., ABCD123456)
        r'_([A-Z]{3,4}\d+)(?=_|$)',
        
        # Mixed case letters + digits
        r'_([A-Za-z]{1,4}\d+)(?=_|$)',
        
        # Numbers only (if they appear to be accessions)
        r'_(\d{6,})(?=_|$)',  # 6+ digits likely accessions
        
        # Handle accessions at the beginning
        r'^([A-Z]{1,4}\d+)_',
        
        # Handle cases where accession might have dots
        r'_([A-Z]{1,4}\d+\.\d+)(?=_|$)',
    ]
    
    if custom_patterns:
        patterns = custom_patterns
    else:
        patterns = default_patterns
    
    return [re.compile(pattern) for pattern in patterns]

def process_tree(tree, accession_patterns, verbose=False):
    """
    Remove accession numbers from all node names in a tree.
    
    Args:
        tree: ETE3 Tree object
        accession_patterns: List of compiled regex patterns
        verbose: Whether to show detailed output
    
    Returns:
        Modified tree with accessions removed
    """
    nodes_modified = 0
    
    # Process each node
    for node in tree.traverse():
        if node.name:
            original_name = node.name
            cleaned_name = remove_accessions_from_name(node.name, accession_patterns)
            
            if cleaned_name != original_name:
                node.name = cleaned_name
                nodes_modified += 1
                if verbose:
                    print(f"    {original_name} ? {cleaned_name}")
    
    if verbose:
        print(f"    Modified {nodes_modified} node names")
    
    return tree

def process_tree_file(input_file, output_file=None, custom_patterns=None, verbose=False):
    """
    Process a tree file that may contain single or multiple trees.
    """
    trees = read_multiple_trees(input_file)
    
    if not trees:
        print(f"No valid trees found in {input_file}")
        return None
    
    print(f"Found {len(trees)} tree(s) in {input_file}")
    
    # Get accession patterns
    accession_patterns = get_accession_patterns(custom_patterns)
    
    # Process each tree
    processed_trees = []
    
    for i, tree in enumerate(trees):
        tree_label = f"Tree {i+1}" if len(trees) > 1 else "Tree"
        
        if verbose:
            print(f"\nProcessing {tree_label}:")
            print(f"  Original nodes: {len(list(tree.traverse()))}")
            print(f"  Sample names before:")
            for leaf in list(tree.get_leaves())[:3]:
                if leaf.name:
                    print(f"    {leaf.name}")
        
        # Process the tree
        processed_tree = process_tree(tree, accession_patterns, verbose)
        processed_trees.append(processed_tree)
        
        if verbose:
            print(f"  Sample names after:")
            for leaf in list(processed_tree.get_leaves())[:3]:
                if leaf.name:
                    print(f"    {leaf.name}")
        else:
            print(f"{tree_label}: processed {len(list(tree.traverse()))} nodes")
    
    # Write output
    if output_file:
        with open(output_file, 'w') as f:
            for tree in processed_trees:
                f.write(tree.write(format=1) + '\n')
        print(f"\nSaved {len(processed_trees)} processed tree(s) to: {output_file}")
    else:
        print("\nProcessed trees in Newick format:")
        for i, tree in enumerate(processed_trees):
            if len(processed_trees) > 1:
                print(f"Tree {i+1}:")
            print(tree.write(format=1))
        print()
    
    return processed_trees

def main():
    parser = argparse.ArgumentParser(
        description="Remove accession numbers from phylogenetic tree node names",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python remove_accessions_multi.py trees.newick -o clean_trees.newick
  
  # Batch process multiple files
  python remove_accessions_multi.py *.tree --batch
  
  # Verbose output showing changes
  python remove_accessions_multi.py trees.newick -o clean_trees.newick -v
  
  # Test what would be changed without saving
  python remove_accessions_multi.py trees.newick --dry-run -v

Common accession patterns handled:
  - AB123456, DQ274092 (standard GenBank)
  - ABCD123456 (extended)
  - 123456789 (numeric only, 6+ digits)
  - AB123456.1 (with version numbers)
        """
    )
    
    parser.add_argument("input", nargs="+", help="Input tree file(s)")
    parser.add_argument("-o", "--output", help="Output file (for single input)")
    parser.add_argument("--batch", action="store_true",
                       help="Process multiple files in batch mode")
    parser.add_argument("--suffix", default="_clean",
                       help="Suffix for output files in batch mode (default: _clean)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Show detailed output with name changes")
    parser.add_argument("--dry-run", action="store_true",
                       help="Show what would be changed without saving")
    parser.add_argument("--custom-patterns", nargs="+",
                       help="Custom regex patterns for accession matching")
    
    args = parser.parse_args()
    
    # Convert custom patterns to list if provided
    custom_patterns = args.custom_patterns if args.custom_patterns else None
    
    if len(args.input) == 1 and not args.batch:
        # Single file mode
        input_file = args.input[0]
        
        if args.dry_run:
            output_file = None
        else:
            output_file = args.output if args.output else f"{input_file.rsplit('.', 1)[0]}{args.suffix}.tree"
        
        process_tree_file(input_file, output_file, custom_patterns, args.verbose)
    
    else:
        # Batch mode
        for input_file in args.input:
            print(f"\nProcessing {input_file}...")
            
            if args.dry_run:
                output_file = None
            else:
                base_name = input_file.rsplit('.', 1)[0]
                output_file = f"{base_name}{args.suffix}.tree"
            
            process_tree_file(input_file, output_file, custom_patterns, args.verbose)

if __name__ == "__main__":
    main()