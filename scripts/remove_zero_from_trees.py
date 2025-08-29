#!/usr/bin/env python3
"""
Script to remove ':0;' at the end of Newick tree strings and replace with just ';'
This fixes a common issue with tree conversion programs.

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/remove_zero_from_trees.py" input_file output_file
"""

import sys
import re
import os

def is_balanced(tree_str):
    """Check if parentheses are balanced in a tree string"""
    count = 0
    for char in tree_str:
        if char == '(':
            count += 1
        elif char == ')':
            count -= 1
        if count < 0:
            return False
    return count == 0

def clean_tree(tree_str):
    """Clean up a single tree string"""
    # Skip if empty
    if not tree_str.strip():
        return ""
    
    # Remove any control characters except newlines
    clean = ''.join(c for c in tree_str if c.isprintable() or c == '\n')
    
    # Fix the :0; at the end
    clean = re.sub(r'\)(\d*):0;', r')\1;', clean)
    clean = re.sub(r'\):0;', r');', clean)
    
    # Fix scientific notation at the root
    clean = re.sub(r'\)(\d*):[\d.]+e[-+][\d]+;', r')\1;', clean)
    
    # Ensure the tree ends with a semicolon
    if not clean.strip().endswith(';'):
        clean = clean.strip() + ';'
    
    return clean

def validate_tree(tree_str):
    """Simple validation of a Newick tree string"""
    if not tree_str.strip():
        return False
    
    # Check balanced parentheses
    if not is_balanced(tree_str):
        return False
    
    # Check for valid pattern: should start with ( and end with ;
    if not tree_str.strip().startswith('(') or not tree_str.strip().endswith(';'):
        return False
    
    return True

def fix_tree_file(input_file, output_file):
    """Process the whole file of trees"""
    try:
        # Read the file
        with open(input_file, 'rb') as f:
            content = f.read()
        
        # Detect and normalize line endings
        if b'\r\n' in content:
            print("Windows-style line endings detected, normalizing...")
            content = content.replace(b'\r\n', b'\n')
        
        # Convert to string
        content_str = content.decode('utf-8', errors='replace')
        
        # Split by semicolon to get individual trees
        tree_parts = content_str.split(';')
        
        # Process each tree
        valid_trees = []
        skipped = 0
        fixed = 0
        
        for i, tree_part in enumerate(tree_parts):
            if i < len(tree_parts) - 1:  # Add back the semicolon except for last empty part
                tree = tree_part + ';'
                clean_tree_str = clean_tree(tree)
                
                if validate_tree(clean_tree_str):
                    valid_trees.append(clean_tree_str)
                    if clean_tree_str != tree:
                        fixed += 1
                else:
                    print("Skipping invalid tree (part {}): {}...".format(i, tree[:50]))
                    skipped += 1
        
        # Write the fixed content
        with open(output_file, 'w', newline='\n') as f:
            for tree in valid_trees:
                f.write(tree + '\n')
        
        print("Processing complete:")
        print("- Total tree parts found: {}".format(len(tree_parts) - 1))  # -1 because last split is empty
        print("- Valid trees written: {}".format(len(valid_trees)))
        print("- Trees fixed: {}".format(fixed))
        print("- Trees skipped (invalid): {}".format(skipped))
        
        # Try to verify final file has proper structure
        with open(output_file, 'r') as f:
            final_content = f.read()
            final_trees = final_content.count(';')
            print("- Trees in final file: {}".format(final_trees))
        
        return True
    
    except Exception as e:
        print("Error: {}".format(str(e)))
        return False

def main():
    if len(sys.argv) != 3:
        print("Usage: python fix_tree_file_comprehensive.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print("Processing tree file: {}".format(input_file))
    print("Writing output to: {}".format(output_file))
    
    if fix_tree_file(input_file, output_file):
        print("Success! Now try running your TreeConversion.jar with the fixed file.")
    else:
        print("Failed to process the file.")

if __name__ == "__main__":
    main()