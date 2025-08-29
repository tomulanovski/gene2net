#!/usr/bin/env python3
import random
import argparse
from ete3 import Tree
'''
usage: python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/select_random_trees.py" -i input_trees.tre -o output_trees.tre -n number of trees
'''
def read_trees(tree_file):
    """Read trees from a file, assuming one tree per line in Newick format."""
    trees = []
    with open(tree_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                trees.append(line)
    return trees

def count_species(tree_str):
    """Count total number of leaves and unique species in a tree."""
    try:
        tree = Tree(tree_str)
        leaf_names = tree.get_leaf_names()
        total_leaves = len(leaf_names)
        unique_species = len(set(leaf_names))
        return total_leaves, unique_species
    except Exception as e:
        print(f"Warning: Could not parse tree: {str(e)}")
        return 0, 0

def main():
    parser = argparse.ArgumentParser(description='Select trees and show species counts')
    parser.add_argument('-i', '--input', required=True, help='Input file with trees in Newick format')
    parser.add_argument('-n', '--number', type=int, required=True, help='Number of trees to select')
    parser.add_argument('-o', '--output', default='selected_trees.txt', help='Output file for selected trees')
    parser.add_argument('-s', '--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    
    # Read all trees
    all_trees = read_trees(args.input)
    print(f"Read {len(all_trees)} trees from {args.input}")
    
    # Check if requested number is valid
    if args.number > len(all_trees):
        print(f"Warning: Requested number ({args.number}) is larger than available trees ({len(all_trees)})")
        num_trees = len(all_trees)
    else:
        num_trees = args.number
    
    # Select random subset
    selected_trees = random.sample(all_trees, num_trees)
    
    # Count species in each selected tree
    total_counts = []
    unique_counts = []
    
    for tree in selected_trees:
        total, unique = count_species(tree)
        if total > 0:  # Only include valid trees
            total_counts.append(total)
            unique_counts.append(unique)
    
    # Calculate and display averages
    if total_counts:
        avg_total = sum(total_counts) / len(total_counts)
        avg_unique = sum(unique_counts) / len(unique_counts)
        
        print(f"Average number of total leaves per tree: {avg_total:.2f}")
        print(f"Min leaves: {min(total_counts)}, Max leaves: {max(total_counts)}")
        
        print(f"Average number of unique species per tree: {avg_unique:.2f}")
        print(f"Min unique species: {min(unique_counts)}, Max unique species: {max(unique_counts)}")
    else:
        print("No valid trees found in selection")
    
    # Write selected trees to output file
    with open(args.output, 'w') as f:
        for tree in selected_trees:
            f.write(f"{tree}\n")
    
    print(f"Wrote {len(selected_trees)} trees to {args.output}")

if __name__ == "__main__":
    main()