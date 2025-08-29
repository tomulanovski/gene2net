import re
from ete3 import Tree

# File paths
INPUT_FILE = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/gene_trees/all_trees.tre"
OUTPUT_FILE = "/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/gene_trees/pruned_trees.tre"
# Taxa to keep (all copies will be kept)
TAXA_TO_KEEP = [
    'RSC04',
    'RSC03',
    'RSC05',
    'RSC02',
    'RSC06',
    'RSC07',
    'RS244',
    'RS421',
    'RSC01',
    'Senec',
    'RS247',
    'RS416',
    'RS265',  # Will keep RS265 but not RS265.1 or RS265.main
    'RS418',
    'RS306',
    'Helia'
]

# Instead of specific exclusions, we'll use a pattern to exclude any taxon.something
def should_keep_leaf(leaf_name):
    """
    Determine if a leaf should be kept based on its name
    Returns: (should_keep, taxon_matched)
    """
    # Check if the leaf name has a dot pattern (taxon.something)
    for taxon in TAXA_TO_KEEP:
        if leaf_name.startswith(taxon + '.'):
            print(f"    Excluding dot-notation leaf: {leaf_name}")
            return False, None
    
    # Then check if it matches any of our taxa to keep
    for taxon in TAXA_TO_KEEP:
        # Check if the leaf name is exactly the taxon or starts with the taxon followed by non-alphabetic character
        # This prevents matching "RS265" with "RS2650" but allows "RS265_something"
        if leaf_name == taxon or leaf_name.startswith(taxon + '_') or leaf_name.startswith(taxon + '.'):
            return True, taxon
            
    # Default case: don't keep
    return False, None

def process_tree(tree_str, tree_idx):
    """Process a single tree string"""
    try:
        # Parse the tree string
        tree = Tree(tree_str, format=1)
        
        # Collect all leaves
        all_leaves = list(tree.get_leaves())
        
        print(f"Tree {tree_idx}: Total leaves: {len(all_leaves)}")
        
        # Identify leaves to keep
        leaves_to_keep = []
        taxa_counts = {taxon: 0 for taxon in TAXA_TO_KEEP}  # To count how many of each taxon we keep
        
        for leaf in all_leaves:
            keep_leaf, matched_taxon = should_keep_leaf(leaf.name)
            
            if keep_leaf:
                leaves_to_keep.append(leaf)
                if matched_taxon:
                    taxa_counts[matched_taxon] += 1
        
        # Report on taxa found
        for taxon, count in taxa_counts.items():
            if count > 0:
                print(f"  {taxon}: keeping all {count} copies")
            else:
                print(f"  Warning: Taxon {taxon} not found in tree {tree_idx}")
        
        # Get the names of leaves to keep
        names_to_keep = [leaf.name for leaf in leaves_to_keep]
        
        # Prune the tree to keep only the selected leaves
        print(f"  Pruning tree {tree_idx} to keep {len(names_to_keep)} leaves...")
        
        # If we have any leaves to keep, prune and return the tree
        if names_to_keep:
            tree.prune(names_to_keep)
            return tree.write(format=1)
        else:
            print(f"  Warning: No leaves to keep in tree {tree_idx}")
            return None
            
    except Exception as e:
        print(f"Error processing tree {tree_idx}: {e}")
        return None

def main():
    # Read all trees from input file
    try:
        with open(INPUT_FILE, 'r') as f:
            content = f.read().strip()
            
            # Check if the file contains multiple trees or a single tree
            if content.count(';') > 1:
                # Multiple trees - split on semicolons
                tree_strings = content.split(';')
                # Remove empty strings and add back the semicolon
                tree_strings = [t.strip() + ';' for t in tree_strings if t.strip()]
            else:
                # Single tree
                tree_strings = [content]
                
    except Exception as e:
        print(f"Error reading input file: {e}")
        return
    
    print(f"Found {len(tree_strings)} trees in input file")
    
    # Process each tree and collect the pruned versions
    pruned_trees = []
    for i, tree_str in enumerate(tree_strings):
        print(f"Processing tree {i+1}/{len(tree_strings)}")
        pruned_tree = process_tree(tree_str, i+1)
        if pruned_tree:
            pruned_trees.append(pruned_tree)
    
    # Write all pruned trees to the output file
    with open(OUTPUT_FILE, 'w') as f:
        for tree in pruned_trees:
            f.write(tree + '\n')
    
    print(f"Wrote {len(pruned_trees)} pruned trees to {OUTPUT_FILE}")
    
    # Print instructions for verification
    print("\nTo verify specific exclusions, check if RS265.1 or RS265.main appear in the output file:")
    print(f"  grep -E 'RS265.1|RS265.main' {OUTPUT_FILE}")

if __name__ == "__main__":
    main()