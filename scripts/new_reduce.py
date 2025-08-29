#!/usr/bin/env python3
"""
Leaf-Centered Tree Reducer
This script starts from each leaf and checks if it's part of a monophyletic group
or a polytomy that should be reduced, with special handling for names ending in 'a' or 'b'.
"""

import dendropy
import itertools
from collections import defaultdict

def parse_taxon_name(taxon_label):
    """
    Parse taxon label with space separators
    Returns (species, variant_type, variant_name, full_name)
    """
    parts = taxon_label.split()
    if len(parts) < 2:
        print(f"WARNING: Couldn't parse '{taxon_label}' - too few parts")
        return taxon_label, None, None, taxon_label
    
    genus = parts[0]
    species = parts[1]
    base_species = f"{genus} {species}"
    
    # Look for subspecies or variety designations
    variant_type = None
    variant_name = None
    i = 2
    while i < len(parts) - 1:  # Leave at least one part for accession
        if parts[i] == 'subsp' and i + 1 < len(parts):
            variant_type = 'subsp'
            variant_name = parts[i+1]
            i += 2
        elif parts[i] == 'var' and i + 1 < len(parts):
            variant_type = 'var'
            variant_name = parts[i+1]
            i += 2
        else:
            # Assume this is part of the accession
            break
    
    return base_species, variant_type, variant_name, taxon_label

def clean_taxon_labels(tree):
    """
    Clean taxon labels by removing accession numbers.
    Leaves only genus, species, and variant information.
    """
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label
        parts = label.split()
        
        if len(parts) < 2:
            continue  # Skip if can't parse properly
            
        genus = parts[0]
        species = parts[1]
        clean_label = f"{genus} {species}"
        
        # Add subspecies or variety if present
        i = 2
        while i < len(parts):
            if parts[i] == 'subsp' and i + 1 < len(parts):
                clean_label += f" subsp {parts[i+1]}"
                i += 2
            elif parts[i] == 'var' and i + 1 < len(parts):
                clean_label += f" var {parts[i+1]}"
                i += 2
            else:
                i += 1
        
        # Update the taxon label
        leaf.taxon.label = clean_label

def get_variant_key(leaf):
    """Get a key representing the species+variant of a leaf"""
    species, var_type, var_name, _ = parse_taxon_name(leaf.taxon.label)
    if var_type and var_name:
        return f"{species} {var_type} {var_name}"
    return species

def is_monophyletic(tree, taxa_subset):
    """Check if a subset of taxa forms a monophyletic group"""
    if len(taxa_subset) <= 1:
        return True
        
    subset_mrca = tree.mrca(taxa=taxa_subset)
    subset_descendants = set(leaf.taxon for leaf in subset_mrca.leaf_nodes())
    return set(taxa_subset) == subset_descendants

def ends_with_a_or_b(taxon_label):
    """Check if a taxon label ends with 'a' or 'b'"""
    return taxon_label[-1] in ['a', 'b']

def select_nodes_to_remove(candidates, leaf_to_keep, verbose=True):
    """
    Select which nodes to remove from candidates, considering:
    - If a name ends with 'a' or 'b', don't remove it unless all candidates end with 'a' or 'b'
    - Always keep the specified leaf
    """
    # Always keep the specified leaf
    to_remove = [node for node in candidates if node != leaf_to_keep]
    
    # Check if any candidate ends with 'a' or 'b'
    special_candidates = [node for node in to_remove if ends_with_a_or_b(node.taxon.label)]
    
    # If some (but not all) candidates end with 'a' or 'b', keep those
    if 0 < len(special_candidates) < len(to_remove):
        if verbose:
            print(f"  Special case: Keeping nodes ending with 'a' or 'b'")
            print(f"  Special nodes: {[node.taxon.label for node in special_candidates]}")
        
        # Remove special nodes from the to_remove list
        return [node for node in to_remove if not ends_with_a_or_b(node.taxon.label)]
    
    # Otherwise, remove all candidates
    return to_remove

def reduce_tree(input_path, output_path, clean_labels=True, verbose=True):
    """
    Reduce the tree using a leaf-centered approach.
    For each leaf, check if it belongs to a monophyletic group or polytomy
    that should be reduced. Special handling for names ending with 'a' or 'b'.
    """
    print("\n==== STARTING LEAF-CENTERED TREE REDUCTION ====\n")
    
    # Read the tree
    print(f"Reading tree from {input_path}")
    tree = dendropy.Tree.get(path=input_path, schema="newick")
    original_taxa_count = len(tree.leaf_nodes())
    
    print(f"Original tree has {original_taxa_count} taxa")
    
    # Group all leaves by variant key (species + variant)
    variant_groups = defaultdict(list)
    all_leaves = list(tree.leaf_node_iter())
    
    for leaf in all_leaves:
        variant_key = get_variant_key(leaf)
        variant_groups[variant_key].append(leaf)
    
    # Report variant grouping
    print(f"\nFound {len(variant_groups)} distinct species/variant combinations")
    for variant, leaves in variant_groups.items():
        if len(leaves) > 1:
            print(f"  {variant}: {len(leaves)} accessions")
    
    # Nodes to remove
    nodes_to_remove = set()
    processed_leaves = set()
    
    print("\n==== ANALYZING INDIVIDUAL LEAVES ====\n")
    
    # Iterate through all leaves
    for leaf in all_leaves:
        # Skip if already processed
        if leaf in processed_leaves or leaf in nodes_to_remove:
            continue
            
        variant_key = get_variant_key(leaf)
        same_variant_leaves = variant_groups[variant_key]
        
        # Skip if this is the only leaf with this variant
        if len(same_variant_leaves) == 1:
            if verbose:
                print(f"Leaf '{leaf.taxon.label}' is the only one with variant '{variant_key}' - keeping it")
            continue
        
        # Check if leaf is part of a monophyletic group
        found_monophyletic_group = False
        
        # Try to find largest monophyletic group containing this leaf
        for size in range(len(same_variant_leaves), 1, -1):
            if found_monophyletic_group:
                break
                
            # Skip combinations already being processed
            valid_leaves = [l for l in same_variant_leaves if l not in processed_leaves and l not in nodes_to_remove]
            if len(valid_leaves) < size:
                continue
                
            # Generate combinations containing our leaf
            for subset in itertools.combinations(valid_leaves, size):
                if leaf not in subset:
                    continue
                    
                subset_taxa = [node.taxon for node in subset]
                if is_monophyletic(tree, subset_taxa):
                    # Found monophyletic group!
                    found_monophyletic_group = True
                    
                    # Select which nodes to remove considering 'a'/'b' rule
                    to_remove = select_nodes_to_remove(subset, leaf, verbose)
                    
                    if verbose:
                        print(f"Leaf '{leaf.taxon.label}' is part of monophyletic group of size {size}")
                        print(f"  Keeping: {leaf.taxon.label}")
                        if to_remove:
                            print(f"  Removing: {[node.taxon.label for node in to_remove]}")
                        else:
                            print(f"  No nodes to remove due to special rule")
                    
                    nodes_to_remove.update(to_remove)
                    processed_leaves.update(subset)
                    break
        
        if found_monophyletic_group:
            continue
            
        # Check if leaf is part of a polytomy
        parent = leaf.parent_node
        if parent:
            # Get direct children of parent that are leaves
            sibling_leaves = [child for child in parent.child_nodes() if child.is_leaf()]
            
            # Check if any siblings have same variant
            same_variant_siblings = [sib for sib in sibling_leaves 
                                   if sib != leaf and get_variant_key(sib) == variant_key]
            
            if same_variant_siblings:
                # Found polytomy with same variant!
                if verbose:
                    print(f"Leaf '{leaf.taxon.label}' is part of polytomy with {len(same_variant_siblings)} siblings of same variant")
                
                # Consider all leaves in this polytomy (including our leaf)
                polytomy_group = [leaf] + same_variant_siblings
                
                # Select which nodes to remove considering 'a'/'b' rule
                to_remove = select_nodes_to_remove(polytomy_group, leaf, verbose)
                
                if verbose:
                    print(f"  Keeping: {leaf.taxon.label}")
                    if to_remove:
                        print(f"  Removing: {[node.taxon.label for node in to_remove]}")
                    else:
                        print(f"  No nodes to remove due to special rule")
                
                nodes_to_remove.update(to_remove)
                processed_leaves.add(leaf)
                processed_leaves.update(same_variant_siblings)
    
    # Summarize and remove nodes
    removal_count = len(nodes_to_remove)
    
    if removal_count == 0:
        print("\n==== NO NODES TO REMOVE ====")
        print("The tree will remain unchanged.")
        
        # Clean taxon labels if requested
        if clean_labels:
            print("\n==== CLEANING TAXON LABELS ====")
            clean_taxon_labels(tree)
            print("Removed accession numbers from taxon labels")
            
        tree.write(path=output_path, schema="newick")
        print(f"Original tree had {original_taxa_count} taxa")
        print(f"Reduced tree has {original_taxa_count} taxa (unchanged)")
        return tree
    
    # Remove identified nodes
    print(f"\n==== REMOVING {removal_count} NODES ====")
    taxa_to_prune = [node.taxon for node in nodes_to_remove]
    
    for node in nodes_to_remove:
        print(f"Removing: {node.taxon.label}")
    
    tree.prune_taxa(taxa_to_prune)
    
    # Clean taxon labels if requested
    if clean_labels:
        print("\n==== CLEANING TAXON LABELS ====")
        clean_taxon_labels(tree)
        print("Removed accession numbers from taxon labels")
    
    # Write the reduced tree
    tree.write(path=output_path, schema="newick")
    final_count = len(tree.leaf_nodes())
    
    print(f"\n==== REDUCTION COMPLETE ====")
    print(f"Original tree had {original_taxa_count} taxa")
    print(f"Reduced tree has {final_count} taxa")
    print(f"Removed {original_taxa_count - final_count} taxa")
    print(f"Final tree written to {output_path}")
    
    return tree

if __name__ == "__main__":
    # Define paths
    input_path = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Bendiksby_2011/gene_trees/extended_tree.tre"
    output_path = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Bendiksby_2011/gene_trees/reduced_tree.tre"
    
    # Run the reduction with clean labels
    reduce_tree(input_path, output_path, clean_labels=True, verbose=True)