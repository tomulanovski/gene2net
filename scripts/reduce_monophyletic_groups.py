#!/usr/bin/env python3
"""
Script to reduce identical tip names in phylogenetic trees.
Handles cases where the same species name appears multiple times as tips.
Example: 16 copies of "Brachypodiumarbuscula" -> 1 copy

python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/reduce_monophyletic_groups.py"
"""

from ete3 import Tree
import argparse
from collections import defaultdict, Counter
import sys

def get_species_counts_in_clade(node):
    """
    Get counts of each species name in a clade.
    Returns dict: {species_name: count}
    """
    if node.is_leaf():
        return {node.name: 1}
    
    species_counts = Counter()
    for leaf in node.get_leaves():
        species_counts[leaf.name] += 1
    
    return dict(species_counts)

def is_single_species_clade(node):
    """
    Check if all tips in this clade have the same name.
    """
    if node.is_leaf():
        return True
    
    leaf_names = {leaf.name for leaf in node.get_leaves()}
    return len(leaf_names) == 1

def should_reduce_clade(node, min_copies=2):
    """
    Check if this clade should be reduced:
    - All tips have same name (single species)
    - Has >= min_copies of that species
    """
    if node.is_leaf():
        return False
    
    if not is_single_species_clade(node):
        return False
    
    num_tips = len(node.get_leaves())
    return num_tips >= min_copies

def should_reduce_polytomy(node):
    """
    Check if a polytomy has multiple copies of the same species that should be reduced.
    Example: (Brachy, Brachy, Oryza, Oryza) -> should reduce to (Brachy, Oryza)
    """
    if node.is_leaf() or len(node.children) < 2:
        return False
    
    # Get the species name for each direct child
    child_species = []
    for child in node.children:
        if child.is_leaf():
            child_species.append(child.name)
        elif is_single_species_clade(child):
            # Get the single species name from this monophyletic subclade
            first_leaf = next(iter(child.get_leaves()))
            child_species.append(first_leaf.name)
        else:
            # Mixed clade - can't reduce this polytomy simply
            return False
    
    # Check if any species appears more than once
    species_counts = Counter(child_species)
    return any(count > 1 for count in species_counts.values())

def reduce_single_species_clade(node):
    """
    Replace a clade containing multiple copies of the same species 
    with a single representative tip.
    """
    if node.is_leaf():
        return node
    
    # Get the species name (all tips should have same name)
    first_leaf = next(iter(node.get_leaves()))
    species_name = first_leaf.name
    
    # Create new leaf with the same name
    new_leaf = Tree(name=species_name)
    
    # Preserve branch length if it exists
    if hasattr(node, 'dist') and node.dist is not None:
        new_leaf.dist = node.dist
    
    return new_leaf

def reduce_polytomy(node):
    """
    Reduce a polytomy by keeping only one representative per species.
    Example: (Brachy, Brachy, Oryza, Oryza) -> (Brachy, Oryza)
    """
    if not should_reduce_polytomy(node):
        return node
    
    # Group children by species
    species_to_child = {}
    
    for child in node.children:
        if child.is_leaf():
            species_name = child.name
        elif is_single_species_clade(child):
            first_leaf = next(iter(child.get_leaves()))
            species_name = first_leaf.name
        else:
            # This shouldn't happen given our should_reduce_polytomy check
            continue
            
        # Keep first occurrence of each species
        if species_name not in species_to_child:
            species_to_child[species_name] = child
    
    # Replace children with representatives
    for child in node.children[:]:  # Copy list to avoid modification issues
        node.remove_child(child)
    
    for representative in species_to_child.values():
        node.add_child(representative)
    
    return node

def reduce_tree(tree, min_copies=2):
    """
    Main function to reduce duplicate species in the tree.
    """
    nodes_to_process = []
    
    # First pass: identify nodes that need reduction
    for node in tree.traverse("postorder"):
        if node.is_root():
            continue
            
        if should_reduce_clade(node, min_copies):
            nodes_to_process.append(('clade', node))
        elif should_reduce_polytomy(node):
            nodes_to_process.append(('polytomy', node))
    
    # Second pass: apply reductions
    for reduction_type, node in nodes_to_process:
        parent = node.up
        
        if reduction_type == 'clade' and parent:
            # Replace entire clade with single tip
            reduced_node = reduce_single_species_clade(node)
            parent.remove_child(node)
            parent.add_child(reduced_node)
            
        elif reduction_type == 'polytomy':
            # Reduce polytomy in place
            reduce_polytomy(node)
    
    return tree

def print_tree_stats(tree, label="Tree"):
    """Print basic statistics about the tree."""
    leaf_names = [leaf.name for leaf in tree.get_leaves()]
    species_counts = Counter(leaf_names)
    
    print(f"{label} statistics:")
    print(f"  Total tips: {len(leaf_names)}")
    print(f"  Unique species: {len(species_counts)}")
    print(f"  Species with multiple copies:")
    for species, count in species_counts.items():
        if count > 1:
            print(f"    {species}: {count} copies")
    print()

def read_multiple_trees(input_file):
    """
    Read multiple trees from a file.
    Handles both single trees and files with multiple trees.
    """
    trees = []
    
    # First try to read as multiple trees
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

def process_tree_file(input_file, output_file=None, min_copies=2, verbose=False):
    """
    Process a tree file that may contain single or multiple trees.
    """
    trees = read_multiple_trees(input_file)
    
    if not trees:
        print(f"No valid trees found in {input_file}")
        return None
    
    print(f"Found {len(trees)} tree(s) in {input_file}")
    
    reduced_trees = []
    
    for i, tree in enumerate(trees):
        tree_label = f"Tree {i+1}" if len(trees) > 1 else "Tree"
        
        if verbose:
            print(f"\n{tree_label}:")
            print_tree_stats(tree, f"  Original {tree_label.lower()}")
        else:
            print(f"{tree_label}: {len(tree.get_leaves())} tips", end="")
        
        # Reduce the tree
        reduced_tree = reduce_tree(tree, min_copies)
        reduced_trees.append(reduced_tree)
        
        if verbose:
            print_tree_stats(reduced_tree, f"  Reduced {tree_label.lower()}")
        else:
            print(f" -> {len(reduced_tree.get_leaves())} tips")
    
    # Output
    if output_file:
        with open(output_file, 'w') as f:
            for reduced_tree in reduced_trees:
                f.write(reduced_tree.write(format=1) + '\n')
        print(f"Saved {len(reduced_trees)} reduced tree(s) to: {output_file}")
    else:
        print("\nReduced tree(s) in Newick format:")
        for i, reduced_tree in enumerate(reduced_trees):
            if len(reduced_trees) > 1:
                print(f"Tree {i+1}:")
            print(reduced_tree.write(format=1))
        print()
    
    return reduced_trees

def main():
    parser = argparse.ArgumentParser(
        description="Reduce duplicate species names in phylogenetic trees",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python reduce_duplicates.py gene_tree.newick -o reduced_tree.newick
  
  # Batch process all tree files
  python reduce_duplicates.py *.tree --batch
  
  # Only reduce if ?3 copies of same species
  python reduce_duplicates.py tree.newick --min-copies 3
  
  # Verbose output showing species counts
  python reduce_duplicates.py tree.newick -v
        """
    )
    
    parser.add_argument("input", nargs="+", help="Input tree file(s)")
    parser.add_argument("-o", "--output", help="Output file (for single input)")
    parser.add_argument("--min-copies", type=int, default=2,
                       help="Minimum copies needed to trigger reduction (default: 2)")
    parser.add_argument("--batch", action="store_true",
                       help="Process multiple files in batch mode")
    parser.add_argument("--suffix", default="_reduced",
                       help="Suffix for output files in batch mode (default: _reduced)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Show detailed statistics")
    
    args = parser.parse_args()
    
    if len(args.input) == 1 and not args.batch:
        # Single file mode
        input_file = args.input[0]
        output_file = args.output if args.output else f"{input_file.rsplit('.', 1)[0]}{args.suffix}.tree"
        process_tree_file(input_file, output_file, args.min_copies, args.verbose)
    
    else:
        # Batch mode
        for input_file in args.input:
            base_name = input_file.rsplit('.', 1)[0]
            output_file = f"{base_name}{args.suffix}.tree"
            print(f"Processing {input_file}...")
            process_tree_file(input_file, output_file, args.min_copies, args.verbose)

if __name__ == "__main__":
    main()