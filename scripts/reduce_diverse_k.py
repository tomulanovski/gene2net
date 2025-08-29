#!/usr/bin/env python3
import re
from ete3 import Tree
import argparse

'''
python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/reduce_diverse_k.py" input_trees.tre -o output_trees.tre

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

def extract_taxon_name(leaf_name):
    """
    Extract the taxon name from a leaf name with Brachypodium naming convention.
    
    Examples:
    Brachypodium_arbuscula_accession123 -> Brachypodium_arbuscula
    Brachypodium_cf_pinnatum_DLA-2015_accession123 -> Brachypodium_cf_pinnatum_DLA-2015
    Brachypodium_distachyon_accession456 -> Brachypodium_distachyon
    Brachypodium_sylvaticum_subsp_glaucovirens_acc789 -> Brachypodium_sylvaticum_subsp_glaucovirens
    """
    
    # Handle special case: Brachypodium_cf_pinnatum_DLA-2015
    cf_pinnatum_match = re.match(r'(Brachypodium_cf_pinnatum_DLA-2015)', leaf_name)
    if cf_pinnatum_match:
        return cf_pinnatum_match.group(1)
    
    # Handle special case: Brachypodium_sylvaticum_subsp_glaucovirens
    subsp_match = re.match(r'(Brachypodium_sylvaticum_subsp_glaucovirens)', leaf_name)
    if subsp_match:
        return subsp_match.group(1)
    
    # Handle other subspecies patterns (in case you have more)
    general_subsp_match = re.match(r'(Brachypodium_[a-z]+_subsp_[a-z]+)', leaf_name)
    if general_subsp_match:
        return general_subsp_match.group(1)
    
    # Handle other "cf" patterns (in case you have more)
    cf_match = re.match(r'(Brachypodium_cf_[a-z]+_[A-Z0-9-]+)', leaf_name)
    if cf_match:
        return cf_match.group(1)
    
    # Handle standard species: Brachypodium_species_accession -> Brachypodium_species
    standard_match = re.match(r'(Brachypodium_[a-z]+)_', leaf_name)
    if standard_match:
        return standard_match.group(1)
    
    # If no underscore (just species name), return as-is
    if '_' not in leaf_name:
        return leaf_name
    
    # Fallback: return everything before the last underscore (assumes last part is accession)
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
        taxon = extract_taxon_name(leaf.name)
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
  # Single tree or multi-tree file with default (keep 1 per taxon)
  python diversity_prune.py input.newick -o output.newick
  
  # Specify number to keep per taxon
  python diversity_prune.py input.newick -o output.newick --keep-per-taxon 2
  
  # Use custom taxa dictionary (modify script)
  python diversity_prune.py input.newick -o output.newick --verbose
        """
    )
    
    parser.add_argument("input", help="Input tree file (single or multiple trees)")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    parser.add_argument("--keep-per-taxon", type=int, default=1,
                       help="Default number of copies to keep per taxon (default: 1)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # Taxa dictionary with desired copy numbers for Brachypodium species
    TAXA_DICT = {
        'Brachypodium_arbuscula': 1,
        'Brachypodium_boissieri': 3, 
        'Brachypodium_bolusii': 2,
        'Brachypodium_distachyon': 1,
        'Brachypodium_flexum': 3,
        'Brachypodium_genuense': 1,
        'Brachypodium_sylvaticum_subsp_glaucovirens': 1,
        'Brachypodium_hybridum': 2,
        'Brachypodium_kawakamii': 3,
        'Brachypodium_madagascariense': 2,
        'Brachypodium_mexicanum': 2,
        'Brachypodium_phoenicoides': 2,
        'Brachypodium_pinnatum': 1,
        'Brachypodium_retusum': 3,
        'Brachypodium_rupestre': 1,
        'Brachypodium_stacei': 1,
        'Brachypodium_sylvaticum': 1,
        'Brachypodium_cf_pinnatum_DLA-2015': 1
    }
    
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
        
        # For taxa not in the predefined dictionary, use default keep-per-taxon value
        temp_taxa_dict = TAXA_DICT.copy()
        
        # Add any new taxa found in tree with default value
        for leaf in tree.get_leaves():
            taxon = extract_taxon_name(leaf.name)
            if taxon not in temp_taxa_dict:
                temp_taxa_dict[taxon] = args.keep_per_taxon
        
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