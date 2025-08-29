import dendropy
import itertools

def get_species_name(taxon_label):
    """Extract species name from the full taxon label."""
    parts = taxon_label.split(' ')
    if len(parts) < 2:
        return taxon_label  # Can't split, return as is
    genus = parts[0]
    species = parts[1]
    
    # Handle subspecies and varieties
    if len(parts) > 2 and parts[2] in ['subsp', 'var']:
        return f"{genus} {species}"  # Return just genus+species
    return f"{genus} {species}"

def reduce_tree(input_tree_path, output_tree_path):
    # Read the tree
    tree = dendropy.Tree.get(path=input_tree_path, schema="newick")
    
    # Group taxa by species
    species_groups = {}
    for leaf in tree.leaf_node_iter():
        species = get_species_name(leaf.taxon.label)
        if species not in species_groups:
            species_groups[species] = []
        species_groups[species].append(leaf)
    
    # Print species groups
    print(f"Found {len(species_groups)} species")
    for species, nodes in species_groups.items():
        print(f"Species: {species}, Accessions: {len(nodes)}")
        if len(nodes) > 1:
            print(f"  Accessions: {[node.taxon.label for node in nodes]}")
    
    # For each species group, find all monophyletic subgroups
    nodes_to_remove = set()
    subgroups_found = 0
    
    for species, nodes in species_groups.items():
        if len(nodes) > 1:  # Multiple accessions exist
            # Get all possible subgroups with at least 2 members
            species_nodes_to_remove = set()
            
            # Check all subgroups, starting with largest
            for size in range(len(nodes), 1, -1):
                for subgroup in itertools.combinations(nodes, size):
                    # Skip if any node in this subgroup is already marked for removal
                    if any(node in species_nodes_to_remove for node in subgroup):
                        continue
                        
                    # Get the mrca of this subgroup
                    mrca = tree.mrca(taxa=[node.taxon for node in subgroup])
                    # Check if subgroup is monophyletic
                    descendants = set([leaf.taxon for leaf in mrca.leaf_nodes()])
                    subgroup_taxa = set([node.taxon for node in subgroup])
                    
                    if descendants == subgroup_taxa:  # Subgroup is monophyletic
                        subgroups_found += 1
                        print(f"Found monophyletic subgroup for {species} with {len(subgroup)} accessions")
                        print(f"  Keeping: {subgroup[0].taxon.label}")
                        print(f"  Removing: {[node.taxon.label for node in subgroup[1:]]}")
                        # Mark for removal (keep first one)
                        species_nodes_to_remove.update(subgroup[1:])
            
            # Add this species' nodes to the global removal set
            nodes_to_remove.update(species_nodes_to_remove)
            
            if not species_nodes_to_remove:
                print(f"No monophyletic subgroups found for {species}")
    
    if subgroups_found == 0:
        print("No monophyletic subgroups found! The tree will remain unchanged.")
    else:
        print(f"Found {subgroups_found} monophyletic subgroups to reduce")
    
    # Prune unwanted nodes
    print(f"Removing {len(nodes_to_remove)} nodes")
    for node in nodes_to_remove:
        tree.prune_taxa([node.taxon])
    
    # Write reduced tree
    tree.write(path=output_tree_path, schema="newick")
    print(f"Original tree had {len(tree.leaf_nodes()) + len(nodes_to_remove)} taxa")
    print(f"Reduced tree has {len(tree.leaf_nodes())} taxa")

# Define paths
input_path = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Bendiksby_2011/gene_trees/extended_tree.tre"
output_path = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Bendiksby_2011/gene_trees/reduced_tree.tre"

# Run the reduction
reduce_tree(input_path, output_path)