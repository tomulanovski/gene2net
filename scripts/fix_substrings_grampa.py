#!/usr/bin/env python3
"""
GRAMPA Substring Label Fixer

This script identifies and fixes tip labels in Newick tree files where one label
is a substring of another, which causes GRAMPA to crash during MUL-tree building.

Usage:
    python fix_grampa_labels.py input_file.nwk [output_file.nwk]
    
If no output file is specified, the script will create one with "_fixed" suffix.
"""

import re
import sys
import argparse
from collections import defaultdict

def extract_labels_from_newick(tree_string):
    """Extract all tip labels from a Newick tree string."""
    # Pattern to match tip labels (alphanumeric, underscore, dash)
    # This matches labels that appear before colons or commas/parentheses
    pattern = r'([A-Za-z0-9_-]+)(?=[:),])'
    labels = re.findall(pattern, tree_string)
    # Remove duplicates while preserving order
    seen = set()
    unique_labels = []
    for label in labels:
        if label not in seen:
            seen.add(label)
            unique_labels.append(label)
    return unique_labels

def find_substring_problems(labels):
    """Find labels that are substrings of other labels."""
    problematic = {}  # {short_label: [longer_labels_containing_it]}
    
    # Sort by length to check shorter labels first
    sorted_labels = sorted(labels, key=len)
    
    for i, short_label in enumerate(sorted_labels):
        containing_labels = []
        for j in range(i + 1, len(sorted_labels)):
            long_label = sorted_labels[j]
            if short_label in long_label and short_label != long_label:
                containing_labels.append(long_label)
        
        if containing_labels:
            problematic[short_label] = containing_labels
    
    return problematic

def generate_safe_replacement(original_label, existing_labels, used_replacements, suffix_style="X"):
    """Generate a safe replacement label that won't create new substring issues."""
    # Different suffix styles to avoid problematic characters
    if suffix_style == "X":
        base_replacement = f"{original_label}X"
    elif suffix_style == "Z":
        base_replacement = f"{original_label}Z"
    elif suffix_style == "number":
        base_replacement = f"{original_label}999"
    else:
        base_replacement = f"{original_label}X"  # default
    
    counter = 0
    
    while True:
        if counter == 0:
            candidate = base_replacement
        else:
            if suffix_style == "number":
                candidate = f"{original_label}{999 + counter}"
            else:
                candidate = f"{base_replacement}{counter}"
        
        # Check if this candidate would create new problems
        is_safe = True
        for existing in existing_labels:
            if existing != original_label:  # Don't compare with the original
                if candidate in existing or existing in candidate:
                    is_safe = False
                    break
        
        # Check against other replacements we've made
        for replacement in used_replacements.values():
            if candidate in replacement or replacement in candidate:
                is_safe = False
                break
        
        if is_safe:
            return candidate
        
        counter += 1
        if counter > 1000:  # Safety valve
            raise Exception(f"Could not generate safe replacement for {original_label}")

def fix_tree_labels(tree_string, label_replacements):
    """Apply label replacements to a tree string."""
    fixed_tree = tree_string
    
    # Sort replacements by length of original label (longest first)
    # to avoid partial replacements
    sorted_replacements = sorted(label_replacements.items(), 
                               key=lambda x: len(x[0]), reverse=True)
    
    for original, replacement in sorted_replacements:
        # Use word boundaries to ensure we're replacing complete labels
        # Pattern matches the label when it's followed by : or ) or ,
        pattern = f'\\b{re.escape(original)}(?=[:),])'
        fixed_tree = re.sub(pattern, replacement, fixed_tree)
    
    return fixed_tree

def process_newick_file(input_file, output_file=None, suffix_style="X"):
    """Process a Newick file and fix substring label issues."""
    
    if output_file is None:
        if input_file.endswith('.nwk'):
            output_file = input_file.replace('.nwk', '_fixed.nwk')
        else:
            output_file = input_file + '_fixed'
    
    print(f"Reading trees from: {input_file}")
    
    try:
        with open(input_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: Could not find file {input_file}")
        return False
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return False
    
    # Split into individual trees (each line should be a tree)
    trees = [line.strip() for line in content.strip().split('\n') if line.strip()]
    
    if not trees:
        print("No trees found in file")
        return False
    
    print(f"Found {len(trees)} tree(s)")
    
    # Collect all labels from all trees
    all_labels = set()
    for tree in trees:
        labels = extract_labels_from_newick(tree)
        all_labels.update(labels)
    
    print(f"Found {len(all_labels)} unique tip labels")
    
    # Find substring problems
    problematic = find_substring_problems(list(all_labels))
    
    if not problematic:
        print("No substring label issues found!")
        print("Your trees should work fine with GRAMPA.")
        return True
    
    print(f"\nFound {len(problematic)} problematic labels:")
    for short_label, containing_labels in problematic.items():
        print(f"  '{short_label}' is contained in: {containing_labels}")
    
    # Generate replacements
    label_replacements = {}
    used_replacements = {}
    
    for problematic_label in problematic.keys():
        safe_replacement = generate_safe_replacement(
            problematic_label, all_labels, used_replacements, suffix_style
        )
        label_replacements[problematic_label] = safe_replacement
        used_replacements[problematic_label] = safe_replacement
    
    print(f"\nProposed replacements:")
    for original, replacement in label_replacements.items():
        print(f"  '{original}' -> '{replacement}'")
    
    # Apply fixes to all trees
    fixed_trees = []
    for tree in trees:
        fixed_tree = fix_tree_labels(tree, label_replacements)
        fixed_trees.append(fixed_tree)
    
    # Write fixed trees
    try:
        with open(output_file, 'w') as f:
            for tree in fixed_trees:
                f.write(tree + '\n')
        print(f"\nFixed trees written to: {output_file}")
        print("You can now use this file with GRAMPA!")
        return True
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Fix substring label issues in Newick tree files for GRAMPA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python fix_grampa_labels.py gene_trees.nwk
    python fix_grampa_labels.py species_tree.nwk species_tree_fixed.nwk
    python fix_grampa_labels.py gene_trees.nwk gene_trees_fixed.nwk
        """
    )
    parser.add_argument('input_file', help='Input Newick tree file')
    parser.add_argument('output_file', nargs='?', 
                       help='Output file (optional, will add _fixed suffix if not provided)')
    parser.add_argument('--suffix', choices=['X', 'Z', 'number'], default='X',
                       help='Suffix style for replacements: X (default), Z, or number')
    parser.add_argument('--dry-run', action='store_true',
                       help='Only check for problems, don\'t fix them')
    
    args = parser.parse_args()
    
    if args.dry_run:
        # Just check for problems
        try:
            with open(args.input_file, 'r') as f:
                content = f.read()
        except FileNotFoundError:
            print(f"Error: Could not find file {args.input_file}")
            return 1
        
        trees = [line.strip() for line in content.strip().split('\n') if line.strip()]
        all_labels = set()
        for tree in trees:
            labels = extract_labels_from_newick(tree)
            all_labels.update(labels)
        
        problematic = find_substring_problems(list(all_labels))
        
        if problematic:
            print(f"Found {len(problematic)} problematic labels:")
            for short_label, containing_labels in problematic.items():
                print(f"  '{short_label}' is contained in: {containing_labels}")
            return 1
        else:
            print("No substring label issues found!")
            return 0
    else:
        # Fix the problems
        success = process_newick_file(args.input_file, args.output_file, args.suffix)
        return 0 if success else 1

if __name__ == '__main__':
    sys.exit(main())