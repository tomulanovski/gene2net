#!/usr/bin/env python3
"""
Process SimPhy gene trees for MP-SUGAR input

This script:
1. Reads gene trees from SimPhy output
2. Strips everything from first underscore in taxon names
3. Adds copy numbers to ALL taxa (e.g., RS04_1, RS04_2), even single occurrences
4. Saves all gene trees in NEXUS format for MP-SUGAR
5. Generates taxon mapping dictionary for MP-SUGAR as JSON

Usage:
    python process_simphy_to_mpsugar.py -i <data_dir> -o <output_file> -m <map_file>
"""

import argparse
import os
import json
from ete3 import Tree
from collections import defaultdict


def strip_underscores(tree_string):
    """
    Strip everything from first underscore onwards in taxa names.
    
    Args:
        tree_string (str): Newick tree string
        
    Returns:
        str: Modified tree string
    """
    # This regex replaces taxon_anything with taxon (keeping the delimiter)
    import re
    return re.sub(r'_[^,):]*([,):])', r'\1', tree_string)


def add_copy_numbers_to_tree(tree, global_taxa_info):
    """
    Modifies a tree by adding sequential copy numbers to ALL taxa.
    Also tracks all taxa and their copies globally.
    
    Args:
        tree (ete3.Tree): The tree to modify
        global_taxa_info (dict): Dictionary tracking all taxa and their max copy numbers
        
    Returns:
        ete3.Tree: The modified tree
    """
    # Count occurrences of each taxon in this tree
    taxa_counts = defaultdict(int)
    taxa_current = defaultdict(int)
    
    # First pass: count occurrences
    for leaf in tree.get_leaves():
        if leaf.name:
            taxa_counts[leaf.name] += 1
    
    # Second pass: rename ALL leaves with copy numbers and track globally
    for leaf in tree.get_leaves():
        if leaf.name:
            original_name = leaf.name
            
            # Track this taxon globally
            if original_name not in global_taxa_info:
                global_taxa_info[original_name] = 0
            
            # Always add copy number (starting from 1)
            taxa_current[original_name] += 1
            leaf.name = f"{original_name}_{taxa_current[original_name]}"
            
            # Update max copy number seen globally
            if taxa_current[original_name] > global_taxa_info[original_name]:
                global_taxa_info[original_name] = taxa_current[original_name]
    
    return tree


def create_taxon_map(global_taxa_info):
    """
    Create MP-SUGAR taxon mapping from global taxa information.
    Maps each taxon to all its copy-numbered versions that appear in trees.
    
    Args:
        global_taxa_info (dict): Dictionary with taxa as keys and max copy numbers as values
        
    Returns:
        dict: Taxon mapping for MP-SUGAR
    """
    taxon_map = {}
    
    for taxon, max_copies in sorted(global_taxa_info.items()):
        # Always create list of copies from _1 to _max
        copies = [f"{taxon}_{i}" for i in range(1, max_copies + 1)]
        taxon_map[taxon] = copies
    
    return taxon_map


def write_taxon_map(taxon_map, output_file):
    """
    Write taxon map to JSON file.
    
    Args:
        taxon_map (dict): Taxon mapping dictionary
        output_file (str): Path to output JSON file
    """
    with open(output_file, 'w') as f:
        json.dump(taxon_map, f, indent=4, sort_keys=True)


def process_gene_trees(data_dir, output_file, map_file):
    """
    Process all gene trees from SimPhy output and create NEXUS file and taxon map.
    
    Args:
        data_dir (str): Path to SimPhy data directory (contains g_* files)
        output_file (str): Path to output NEXUS file
        map_file (str): Path to output taxon map JSON file
        
    Returns:
        int: Number of trees processed
    """
    gene_trees = []
    tree_count = 0
    global_taxa_info = {}  # Track all taxa and their max copy numbers across all trees
    
    # Find all gene tree files (starting with 'g_')
    gene_files = sorted([f for f in os.listdir(data_dir) if f.startswith('g_')])
    
    if not gene_files:
        print(f"WARNING: No gene tree files found in {data_dir}")
        return 0
    
    print(f"Found {len(gene_files)} gene tree files")
    
    for gene_file in gene_files:
        gene_path = os.path.join(data_dir, gene_file)
        
        try:
            # Read the tree
            with open(gene_path, 'r') as f:
                tree_string = f.read().strip()
            
            # Step 1: Strip underscores
            stripped_tree_string = strip_underscores(tree_string)
            
            # Step 2: Parse tree
            tree = Tree(stripped_tree_string, format=1)
            
            # Step 3: Add copy numbers and track globally
            modified_tree = add_copy_numbers_to_tree(tree, global_taxa_info)
            
            # Step 4: Get Newick string (format=9 gives only topology with leaf names)
            newick_string = modified_tree.write(format=9)
            
            gene_trees.append(newick_string)
            tree_count += 1
            
        except Exception as e:
            print(f"ERROR processing {gene_file}: {e}")
            continue
    
    # Write NEXUS file
    if gene_trees:
        with open(output_file, 'w') as out:
            out.write("#NEXUS\n")
            out.write("BEGIN TREES;\n")
            
            for i, newick in enumerate(gene_trees, 1):
                out.write(f"Tree geneTree{i} = {newick}\n")
            
            out.write("END;\n")
        
        print(f"Successfully wrote {tree_count} trees to {output_file}")
    
    # Create and write taxon map
    if global_taxa_info:
        taxon_map = create_taxon_map(global_taxa_info)
        write_taxon_map(taxon_map, map_file)
        print(f"Successfully wrote taxon map to {map_file}")
        print(f"Found {len(taxon_map)} unique taxa")
        
        # Print summary
        single_copy_taxa = {k: v for k, v in taxon_map.items() if len(v) == 1}
        multi_copy_taxa = {k: v for k, v in taxon_map.items() if len(v) > 1}
        
        if single_copy_taxa:
            print(f"Taxa with single copy ({len(single_copy_taxa)}): always appear once per tree")
        if multi_copy_taxa:
            print(f"Taxa with multiple copies ({len(multi_copy_taxa)}):")
            for taxon, copies in sorted(multi_copy_taxa.items()):
                print(f"  {taxon}: max {len(copies)} copies per tree")
    
    return tree_count


def main():
    """
    Main function to parse arguments and process gene trees.
    """
    parser = argparse.ArgumentParser(
        description="Process SimPhy gene trees for MP-SUGAR input"
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to SimPhy replicate data directory (contains g_* files)"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output NEXUS file"
    )
    
    parser.add_argument(
        "-m", "--map",
        required=True,
        help="Path to output taxon map JSON file"
    )
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.isdir(args.input):
        print(f"ERROR: Input directory not found: {args.input}")
        return 1
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    map_dir = os.path.dirname(args.map)
    if map_dir:
        os.makedirs(map_dir, exist_ok=True)
    
    # Process the gene trees
    tree_count = process_gene_trees(args.input, args.output, args.map)
    
    if tree_count == 0:
        print("ERROR: No trees were processed")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())