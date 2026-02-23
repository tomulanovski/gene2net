#!/usr/bin/env python3
import re
import json
from ete3 import Tree
import argparse

'''
Diversity-based pruning of phylogenetic trees.

Keeps the k most diverse copies per species, where k is specified via --taxa-json.

Usage:
  # With a JSON file specifying copy numbers per taxon
  python reduce_diverse_k.py input.tre -o output.tre --taxa-json taxa_copies.json

  # With default copies per taxon (e.g., keep 2 for all)
  python reduce_diverse_k.py input.tre -o output.tre --keep-per-taxon 2

The taxa JSON file should map taxon names to copy numbers:
  {"Ephedra_sinica": 2, "Ephedra_foeminea": 1, ...}

Taxon name extraction: by default, everything before the last underscore
(e.g., Ephedra_sinica_KT033298 -> Ephedra_sinica).
Override with --taxon-separator and --accession-position if needed.
'''

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

def extract_taxon_name(leaf_name, taxa_dict=None):
    """
    Extract the taxon name from a leaf name.

    Strategy:
    1. If taxa_dict is provided, try matching leaf_name against known taxon names
       (longest match first). This handles cases like 'Genus_species_var_subsp'.
    2. Fallback: everything before the last underscore.
       e.g., Ephedra_sinica_KT033298 -> Ephedra_sinica

    Examples:
      Ephedra_sinica_KT033298 -> Ephedra_sinica
      Ephedra_saxatilis_var_mairei_KT033359 -> Ephedra_saxatilis_var_mairei (if in taxa_dict)
      Brachypodium_distachyon_acc456 -> Brachypodium_distachyon
    """
    # If we have a taxa dictionary, try to match known taxon names (longest first)
    if taxa_dict:
        # Sort by length descending so longer names match first
        for taxon in sorted(taxa_dict.keys(), key=len, reverse=True):
            if leaf_name.startswith(taxon + '_') or leaf_name == taxon:
                return taxon

    # If no underscore, return as-is
    if '_' not in leaf_name:
        return leaf_name

    # Fallback: everything before the last underscore
    parts = leaf_name.rsplit('_', 1)
    return parts[0]

def get_pairwise_distances(nodes):
    """Calculate pairwise distances between all nodes in the list"""
    distances = {}
    for i, node1 in enumerate(nodes):
        for node2 in nodes[i+1:]:
            try:
                ancestor = node1.get_common_ancestor(node2)
                dist = node1.get_distance(ancestor) + node2.get_distance(ancestor)
                distances[(node1, node2)] = dist
            except:
                # Fallback if distance calculation fails
                distances[(node1, node2)] = 1.0
    return distances

def select_most_diverse_leaves(tree, taxon, nodes, num_to_keep):
    """
    Select the most diverse leaves for a given taxon using a greedy approach.
    """
    # If we have fewer leaves than requested, keep all of them
    if len(nodes) <= num_to_keep:
        return nodes
    
    # If only one leaf needed, find the one with maximum distance to the root
    if num_to_keep == 1:
        try:
            return [max(nodes, key=lambda leaf: leaf.get_distance(tree))]
        except:
            return [nodes[0]]  # Fallback
    
    # Calculate all pairwise distances
    all_distances = get_pairwise_distances(nodes)
    
    # If no distances (only one node), return that node
    if not all_distances:
        return nodes[:num_to_keep]
    
    # Start with the two most distant leaves
    try:
        most_distant_pair = max(all_distances.items(), key=lambda x: x[1])[0]
        selected_leaves = list(most_distant_pair)
    except:
        selected_leaves = nodes[:2]  # Fallback
    
    # Iteratively add the remaining leaves
    while len(selected_leaves) < num_to_keep:
        remaining_leaves = [leaf for leaf in nodes if leaf not in selected_leaves]
        if not remaining_leaves:
            break
            
        # Find the leaf with maximum minimum distance to any already selected leaf
        max_min_distance = -1
        leaf_to_add = None
        
        for leaf in remaining_leaves:
            min_distance = float('inf')
            for selected_leaf in selected_leaves:
                # Get the key for the distance dictionary (order doesn't matter)
                key = (leaf, selected_leaf) if id(leaf) < id(selected_leaf) else (selected_leaf, leaf)
                if key in all_distances:
                    dist = all_distances[key]
                    min_distance = min(min_distance, dist)
            
            if min_distance > max_min_distance:
                max_min_distance = min_distance
                leaf_to_add = leaf
        
        if leaf_to_add:
            selected_leaves.append(leaf_to_add)
        else:
            # Fallback: just add any remaining leaf
            selected_leaves.extend(remaining_leaves[:num_to_keep - len(selected_leaves)])
            break
    
    return selected_leaves

def process_single_tree(tree, taxa_dict, tree_num=None, verbose=False):
    """
    Process a single tree with diversity-based pruning.
    """
    tree_label = f"Tree {tree_num}" if tree_num is not None else "Tree"
    
    if verbose:
        print(f"\nProcessing {tree_label}:")
        print(f"  Original leaves: {len(tree.get_leaves())}")
    
    # Collect all leaves and group by taxon
    taxa_to_leaves = {}
    all_leaves = list(tree.get_leaves())
    
    for leaf in all_leaves:
        taxon = extract_taxon_name(leaf.name, taxa_dict)
        if taxon not in taxa_to_leaves:
            taxa_to_leaves[taxon] = []
        taxa_to_leaves[taxon].append(leaf)
    
    if verbose:
        print(f"  Found taxa: {list(taxa_to_leaves.keys())}")
    
    # Determine leaves to keep
    leaves_to_keep = []
    
    for taxon, leaves in taxa_to_leaves.items():
        # Default to keeping 1 copy if taxon not in dictionary
        num_to_keep = taxa_dict.get(taxon, 1)
        
        if len(leaves) > num_to_keep:
            # Select the most diverse leaves to keep
            selected_leaves = select_most_diverse_leaves(tree, taxon, leaves, num_to_keep)
            leaves_to_keep.extend(selected_leaves)
            
            if verbose:
                print(f"    {taxon}: keeping {len(selected_leaves)} of {len(leaves)} leaves")
        else:
            # Keep all leaves if we have fewer than or equal to the required number
            leaves_to_keep.extend(leaves)
            if verbose:
                print(f"    {taxon}: keeping all {len(leaves)} leaves")
    
    # Get the names of leaves to keep
    names_to_keep = [leaf.name for leaf in leaves_to_keep]
    
    # Prune the tree to keep only the selected leaves
    try:
        tree.prune(names_to_keep)
        if verbose:
            print(f"  Pruned to: {len(tree.get_leaves())} leaves")
    except Exception as e:
        print(f"  Error pruning {tree_label}: {e}")
        return None
    
    return tree

def main():
    parser = argparse.ArgumentParser(
        description="Diversity-based pruning of phylogenetic trees",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using a JSON taxa dictionary (recommended)
  python reduce_diverse_k.py input.tre -o output.tre --taxa-json copies.json -v

  # Using default copies per taxon
  python reduce_diverse_k.py input.tre -o output.tre --keep-per-taxon 2

Taxa JSON format:
  {"Ephedra_sinica": 2, "Ephedra_foeminea": 1, ...}
        """
    )

    parser.add_argument("input", help="Input tree file (single or multiple trees)")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument("--taxa-json", help="JSON file mapping taxon names to copy numbers")
    parser.add_argument("--keep-per-taxon", type=int, default=1,
                       help="Default number of copies to keep per taxon (default: 1)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output")

    args = parser.parse_args()

    # Load taxa dictionary from JSON or use empty dict
    if args.taxa_json:
        with open(args.taxa_json, 'r') as f:
            TAXA_DICT = json.load(f)
        print(f"Loaded copy numbers for {len(TAXA_DICT)} taxa from {args.taxa_json}")
    else:
        TAXA_DICT = {}
        print(f"No --taxa-json provided, using --keep-per-taxon={args.keep_per_taxon} for all taxa")

    # Read trees
    trees = read_multiple_trees(args.input)

    if not trees:
        print(f"No valid trees found in {args.input}")
        return

    print(f"Found {len(trees)} tree(s) in {args.input}")

    # Process each tree
    processed_trees = []
    for i, tree in enumerate(trees):
        tree_num = i + 1 if len(trees) > 1 else None

        # For taxa not in the dictionary, use default keep-per-taxon value
        temp_taxa_dict = TAXA_DICT.copy()

        # Add any new taxa found in tree with default value
        for leaf in tree.get_leaves():
            taxon = extract_taxon_name(leaf.name, temp_taxa_dict)
            if taxon not in temp_taxa_dict:
                temp_taxa_dict[taxon] = args.keep_per_taxon
                if args.verbose:
                    print(f"  Note: {taxon} not in taxa JSON, using default k={args.keep_per_taxon}")

        processed_tree = process_single_tree(tree, temp_taxa_dict, tree_num, args.verbose)
        if processed_tree:
            processed_trees.append(processed_tree)

    # Write output
    if processed_trees:
        with open(args.output, 'w') as f:
            for tree in processed_trees:
                f.write(tree.write(format=1) + '\n')

        print(f"\nSaved {len(processed_trees)} processed tree(s) to {args.output}")
    else:
        print("No trees were successfully processed")

if __name__ == "__main__":
    main()