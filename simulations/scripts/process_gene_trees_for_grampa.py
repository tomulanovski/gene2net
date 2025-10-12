#!/usr/bin/env python3
"""
Process gene trees for GRAMPA: Clean, fix substrings, and reformat in one pass.

This script:
1. Removes everything from first underscore onwards in taxa names
2. Fixes substring label issues (where one label is contained in another)
3. Reformats for GRAMPA by adding copy numbers (e.g., 1_species, 2_species)
4. Outputs a taxa replacement map for use with species trees

Usage:
    python "/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/process_gene_trees_for_grampa.py" <input_dir> <output_file> [--taxa-map <map_file>]
"""

import sys
import re
import argparse
from pathlib import Path
from collections import defaultdict


def clean_taxa_names(tree_string):
    """Remove everything from first underscore onwards in taxa names."""
    return re.sub(r'_[^,):]*([,):])', r'\1', tree_string)


def extract_taxa(tree_string):
    """Extract all unique taxa names from a tree string."""
    pattern = r'([a-zA-Z][a-zA-Z0-9_]*)'
    taxa = re.findall(pattern, tree_string)
    return list(set(taxa))


def find_substring_issues(taxa):
    """Find taxa that are substrings of other taxa."""
    problematic = {}
    sorted_taxa = sorted(taxa, key=len)
    
    for i, short_taxon in enumerate(sorted_taxa):
        containing = []
        for j in range(i + 1, len(sorted_taxa)):
            long_taxon = sorted_taxa[j]
            if short_taxon in long_taxon and short_taxon != long_taxon:
                containing.append(long_taxon)
        
        if containing:
            problematic[short_taxon] = containing
    
    return problematic


def generate_replacements(problematic, all_taxa):
    """Generate safe replacement names that won't create new substring issues."""
    replacements = {}
    
    for taxon in problematic.keys():
        # Try adding 'X' suffix
        candidate = f"{taxon}X"
        
        # Check if this creates new issues
        safe = True
        for other_taxon in all_taxa:
            if other_taxon != taxon:
                if candidate in other_taxon or other_taxon in candidate:
                    safe = False
                    break
        
        # Check against other replacements
        for other_replacement in replacements.values():
            if candidate in other_replacement or other_replacement in candidate:
                safe = False
                break
        
        if safe:
            replacements[taxon] = candidate
        else:
            # Fallback: try 'Z' or add numbers
            counter = 1
            while not safe and counter < 1000:
                candidate = f"{taxon}X{counter}"
                safe = True
                for other_taxon in all_taxa:
                    if other_taxon != taxon:
                        if candidate in other_taxon or other_taxon in candidate:
                            safe = False
                            break
                if safe:
                    replacements[taxon] = candidate
                    break
                counter += 1
    
    return replacements


def fix_substring_issues(tree_string, replacements):
    """Apply replacements to fix substring issues."""
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


def reformat_for_grampa(tree_string):
    """Add copy numbers to taxa for GRAMPA format (e.g., 1_species, 2_species)."""
    taxa_counts = defaultdict(int)
    result = []
    i = 0
    
    while i < len(tree_string):
        char = tree_string[i]
        
        # Check if we're at the start of a taxon name (letter or digit)
        if char.isalnum():
            # Check if this position could be a taxon (after '(' or ',')
            if i == 0 or tree_string[i-1] in '(,':
                # Extract the full taxon name
                taxon = ''
                j = i
                while j < len(tree_string) and (tree_string[j].isalnum() or tree_string[j] == '_'):
                    taxon += tree_string[j]
                    j += 1
                
                # Increment count and add reformatted name
                taxa_counts[taxon] += 1
                result.append(f"{taxa_counts[taxon]}_{taxon}")
                i = j
                continue
        
        result.append(char)
        i += 1
    
    return ''.join(result)


def process_single_tree(tree_string, replacements):
    """Process a single tree through all steps."""
    # Step 1: Clean taxa names (remove underscores)
    tree = clean_taxa_names(tree_string)
    
    # Step 2: Fix substring issues
    tree = fix_substring_issues(tree, replacements)
    
    # Step 3: Reformat for GRAMPA
    tree = reformat_for_grampa(tree)
    
    return tree


def write_taxa_map(replacements, output_file, verbose=False):
    """Write the taxa replacement map to a file."""
    if not replacements:
        if verbose:
            print(f"No taxa replacements needed, skipping map file: {output_file}")
        return
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Taxa Replacement Map\n")
        f.write("# Generated by process_gene_trees_for_grampa.py\n")
        f.write("# Format: ORIGINAL_TAXON\\tREPLACEMENT_TAXON\n")
        f.write("#\n")
        
        # Write replacements (sorted for readability)
        for original in sorted(replacements.keys()):
            replacement = replacements[original]
            f.write(f"{original}\t{replacement}\n")
    
    if verbose:
        print(f"Taxa map saved to: {output_file}")


def process_gene_trees(input_dir, output_file, taxa_map_file=None, verbose=False):
    """Process all gene trees in a directory."""
    input_path = Path(input_dir)
    
    if not input_path.exists():
        print(f"ERROR: Input directory not found: {input_dir}")
        return False
    
    # Find all gene tree files (starting with 'g_')
    gene_tree_files = sorted(input_path.glob('g_*'))
    
    if not gene_tree_files:
        print(f"WARNING: No gene tree files found in {input_dir}")
        return False
    
    print(f"Found {len(gene_tree_files)} gene trees")
    
    # Step 1: Collect all trees and clean them
    if verbose:
        print("Step 1: Cleaning taxa names...")
    
    cleaned_trees = []
    for tree_file in gene_tree_files:
        with open(tree_file, 'r') as f:
            tree = f.read().strip()
            cleaned_tree = clean_taxa_names(tree)
            cleaned_trees.append(cleaned_tree)
    
    # Step 2: Analyze all taxa across all trees to find substring issues
    if verbose:
        print("Step 2: Analyzing for substring issues...")
    
    all_taxa = set()
    for tree in cleaned_trees:
        taxa = extract_taxa(tree)
        all_taxa.update(taxa)
    
    print(f"Found {len(all_taxa)} unique taxa")
    
    problematic = find_substring_issues(list(all_taxa))
    
    if problematic:
        print(f"Found {len(problematic)} taxa with substring issues:")
        for taxon, containing in problematic.items():
            print(f"  '{taxon}' is contained in: {containing}")
        
        replacements = generate_replacements(problematic, all_taxa)
        print("Generated replacements:")
        for original, replacement in replacements.items():
            print(f"  '{original}' -> '{replacement}'")
    else:
        print("No substring issues found")
        replacements = {}
    
    # Step 2.5: Write taxa map if requested
    if taxa_map_file:
        write_taxa_map(replacements, taxa_map_file, verbose)
    
    # Step 3: Process all trees and write to output
    if verbose:
        print("Step 3: Reformatting for GRAMPA and writing output...")
    
    with open(output_file, 'w') as out:
        for i, cleaned_tree in enumerate(cleaned_trees, 1):
            # Fix substrings and reformat
            processed_tree = fix_substring_issues(cleaned_tree, replacements)
            final_tree = reformat_for_grampa(processed_tree)
            out.write(final_tree + '\n')
            
            if verbose and i % 100 == 0:
                print(f"  Processed {i}/{len(cleaned_trees)} trees...")
    
    print(f"Successfully processed {len(cleaned_trees)} trees")
    print(f"Output saved to: {output_file}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Process gene trees for GRAMPA: clean, fix substrings, and reformat",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python process_trees_for_grampa.py /path/to/replicate_dir output_trees.tre
    python process_trees_for_grampa.py /path/to/replicate_dir output_trees.tre --taxa-map taxa_map.txt
        """
    )
    
    parser.add_argument('input_dir', help='Directory containing gene tree files (g_*)')
    parser.add_argument('output_file', help='Output file for processed trees')
    parser.add_argument('--taxa-map', '-m', dest='taxa_map_file',
                       help='Output file for taxa replacement map (optional)')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Print detailed progress information')
    
    args = parser.parse_args()
    
    success = process_gene_trees(args.input_dir, args.output_file, 
                                 args.taxa_map_file, args.verbose)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()