#!/usr/bin/env python3

from ete3 import Tree
import re
from itertools import combinations

tree = Tree("/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/gene_trees/reduced_nia.tre")

def extract_species_id(name):
    """Extract species and ID from taxon name"""
    # The pattern is genus_species_ID_letter
    match = re.match(r'(\w+)_(\w+)_([^_]+)_\w+', name)
    if match:
        genus, species, identifier = match.groups()
        return f"{genus}_{species}_{identifier}"
    
    # Handle special cases like 4x in the identifier
    match = re.match(r'(\w+)_(\w+)_(\w+)_([^_]+)_\w+', name)
    if match:
        genus, species, modifier, identifier = match.groups()
        return f"{genus}_{species}_{modifier}_{identifier}"
    
    return name  # Return the original name if pattern doesn't match

def get_taxa_by_species_id():
    """Group all taxa by species+ID (without the letter suffix)"""
    taxa_by_species_id = {}
    
    for leaf in tree.get_leaves():
        name = leaf.name
        species_id = extract_species_id(name)
        
        if species_id not in taxa_by_species_id:
            taxa_by_species_id[species_id] = []
            
        taxa_by_species_id[species_id].append(leaf)
        
    return taxa_by_species_id

def is_monophyletic(tree, taxa):
    """Check if a set of taxa forms a monophyletic group"""
    if len(taxa) <= 1:
        return False  # Single taxon can't form a monophyletic group for reduction
        
    # Get the most recent common ancestor
    mrca = tree.get_common_ancestor(taxa)
    
    # Get all leaves under the MRCA
    mrca_leaves = set(mrca.get_leaves())
    taxa_set = set(taxa)
    
    # If all leaves under MRCA are in our taxa list, it's monophyletic
    return mrca_leaves.issubset(taxa_set) or mrca_leaves == taxa_set

def find_all_monophyletic_subgroups():
    """Find all monophyletic subgroups for each species group"""
    taxa_by_species_id = get_taxa_by_species_id()
    monophyletic_groups = []
    
    for species_id, taxa_list in taxa_by_species_id.items():
        # Skip if there's only one taxon
        if len(taxa_list) <= 1:
            continue
        
        # Try all possible subgroups, starting from the largest
        for size in range(len(taxa_list), 1, -1):
            for subset in combinations(taxa_list, size):
                if is_monophyletic(tree, subset):
                    # Convert to list for easier handling
                    subset_list = list(subset)
                    monophyletic_groups.append((species_id, subset_list))
                    
                    # Remove these taxa from further consideration
                    for taxon in subset_list:
                        if taxon in taxa_list:
                            taxa_list.remove(taxon)
                    
                    # Break the size loop and move on to the next set of remaining taxa
                    break
        
    return monophyletic_groups

def reduce_tree_based_on_monophyletic_groups(monophyletic_groups):
    """Create a reduced tree keeping only one representative from each monophyletic group"""
    nodes_to_keep = set()
    
    # Add all leaves to nodes_to_keep initially
    for leaf in tree.get_leaves():
        nodes_to_keep.add(leaf.name)
    
    # For each monophyletic group, keep only one representative
    removed_nodes = []
    monophyletic_reductions = []
    
    for species_id, taxa_list in monophyletic_groups:
        # Keep only the first taxon
        kept_node = taxa_list[0].name
        to_remove = [taxon.name for taxon in taxa_list[1:]]
        
        # Update the set of nodes to keep
        for node_name in to_remove:
            if node_name in nodes_to_keep:
                nodes_to_keep.remove(node_name)
                removed_nodes.append(node_name)
        
        monophyletic_reductions.append((species_id, kept_node, to_remove))
    
    # Create a pruned tree with only the nodes to keep
    pruned_tree = tree.copy()
    pruned_tree.prune(list(nodes_to_keep), preserve_branch_length=True)
    
    return pruned_tree, monophyletic_reductions, removed_nodes

# Find all monophyletic groups
monophyletic_groups = find_all_monophyletic_subgroups()

# Reduce the tree
reduced_tree, reductions, removed_nodes = reduce_tree_based_on_monophyletic_groups(monophyletic_groups)

# Print results
print("=== MONOPHYLETIC GROUPS REDUCED ===")
for i, (species_id, kept, removed) in enumerate(reductions, 1):
    print(f"Group {i}: {species_id}")
    print(f"  Kept: {kept}")
    print(f"  Removed: {', '.join(removed)}")
    print()

print("=== NODES REMOVED ===")
print(", ".join(removed_nodes))
print(f"Total nodes removed: {len(removed_nodes)}")

print("\n=== REDUCED TREE (NEWICK FORMAT) ===")
print(reduced_tree.write(format=1))

# Save to the specified path
output_path = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Lawrence_2016/gene_trees/reduced_nia.tre"
with open(output_path, "w") as f:
    f.write(reduced_tree.write(format=1))

print(f"\nReduced tree saved to '{output_path}'")