"""
Gene Tree Species Distribution Analyzer

This script analyzes multiple gene trees from a single file to determine 
the distribution of copy numbers for each species across all trees, and 
identifies the most representative copy number for each species using a kernel
smoothing approach.

usage:

python "/groups/itay_mayrose/tomulanovski/gene2net/scripts/copies_smoothing.py" -i input_trees.tre -o output.tsv -m output_multi_set.txt
"""

import argparse
from collections import defaultdict, Counter
from ete3 import Tree

def extract_species_from_leaf_name(leaf_name, extract_mode="full"):
    """
    Extract species identifier from leaf name.
    
    Args:
        leaf_name (str): The name of the leaf node
        extract_mode (str): How to extract the species name:
            - "full": Use the full leaf name
            - "before": Extract part before first underscore
            - "after": Extract part after first underscore
        
    Returns:
        str: The species identifier
    """
    if extract_mode == "after" and "_" in leaf_name:
        return leaf_name.split("_", 1)[1]
    elif extract_mode == "before" and "_" in leaf_name:
        return leaf_name.split("_", 1)[0]
    else:
        return leaf_name

def count_species_copies_in_tree(tree, extract_mode="full"):
    """
    Count the number of copies for each species in a tree.
    
    Args:
        tree (ete3.Tree): An ETE3 Tree object
        extract_mode (str): How to extract the species name ("full", "before", or "after")
        
    Returns:
        dict: Dictionary with species as keys and their copy counts as values
    """
    species_counts = defaultdict(int)
    
    # Count occurrences of each species
    for leaf in tree.get_leaves():
        if leaf.name:
            species = extract_species_from_leaf_name(leaf.name, extract_mode)
            species_counts[species] += 1
            
    return dict(species_counts)

def analyze_trees_from_file(file_path, extract_mode="full"):
    """
    Analyze all trees from a single file containing multiple Newick trees.
    
    Args:
        file_path (str): Path to the file containing multiple trees
        extract_mode (str): How to extract the species name ("full", "before", or "after")
        
    Returns:
        dict: Dictionary with tree index as keys and species counts as values
    """
    tree_species_counts = {}
    tree_index = 0
    
    try:
        # Read the entire file content
        with open(file_path, 'r') as file:
            content = file.read()
        
        # Split by ';' to get individual trees, but add ';' back as it's needed for Newick format
        tree_strings = [ts + ';' for ts in content.split(';') if ts.strip()]
        
        # Process each tree
        for tree_str in tree_strings:
            try:
                tree = Tree(tree_str)
                species_counts = count_species_copies_in_tree(tree, extract_mode)
                tree_species_counts[f"Tree_{tree_index}"] = species_counts
                tree_index += 1
            except Exception as e:
                print(f"Error processing tree #{tree_index}: {e}")
        
        print(f"Successfully processed {tree_index} trees")
        return tree_species_counts
    
    except Exception as e:
        print(f"Error opening or reading file {file_path}: {e}")
        return {}

def get_species_copy_distributions(tree_species_counts):
    """
    Calculate the distribution of copy numbers for each species across all trees.
    
    Args:
        tree_species_counts (dict): Dictionary with tree indices as keys and species counts as values
        
    Returns:
        dict: Dictionary with species as keys and counters of copy numbers as values
    """
    species_copy_distributions = defaultdict(Counter)
    
    for tree_id, species_counts in tree_species_counts.items():
        for species, count in species_counts.items():
            species_copy_distributions[species][count] += 1
            
    return species_copy_distributions

def get_representative_copy_numbers(species_copy_distributions, kernel_width=2):
    """
    Determine the most representative copy number for each species using 
    a kernel smoothing approach to find the peak of the distribution.
    
    Args:
        species_copy_distributions (dict): Dictionary with species as keys and 
                                          counters of copy numbers as values
        kernel_width (int): How far the kernel extends in each direction
        
    Returns:
        dict: Dictionary with species as keys and representative copy numbers as values
    """
    representative_copies = {}
    
    for species, count_distribution in species_copy_distributions.items():
        # Get sorted list of copy numbers and their frequencies
        sorted_counts = sorted(count_distribution.items())
        
        # If there's only one copy number, use it
        if len(sorted_counts) == 1:
            copy_number, _ = sorted_counts[0]
            representative_copies[species] = copy_number
            continue
            
        # If there are just two copy numbers, use the more frequent one
        if len(sorted_counts) == 2:
            if sorted_counts[0][1] >= sorted_counts[1][1]:
                representative_copies[species] = sorted_counts[0][0]
            else:
                representative_copies[species] = sorted_counts[1][0]
            continue
        
        # For 3 or more distinct copy numbers, use kernel smoothing approach
        copy_numbers = [num for num, _ in sorted_counts]
        frequencies = [freq for _, freq in sorted_counts]
        
        # Calculate smoothed distribution - use a triangular kernel
        smoothed_values = {}
        
        for i, center in enumerate(copy_numbers):
            # Apply kernel centered at each observed copy number
            for j, copy_num in enumerate(copy_numbers):
                # Calculate kernel weight based on distance
                distance = abs(center - copy_num)
                if distance <= kernel_width:
                    weight = 1 - (distance / (kernel_width + 1))
                    
                    # Add weighted frequency to the smoothed value
                    if center not in smoothed_values:
                        smoothed_values[center] = 0
                    smoothed_values[center] += frequencies[j] * weight
        
        # Find the copy number with the highest smoothed value
        max_smoothed = 0
        representative = copy_numbers[0]
        
        for copy_num, smoothed_val in smoothed_values.items():
            if smoothed_val > max_smoothed:
                max_smoothed = smoothed_val
                representative = copy_num
        
        representative_copies[species] = representative
            
    return representative_copies

def write_multi_set_file(representative_copies, output_file):
    """
    Write multi_set file with species names repeated according to their representative copy numbers.
    
    Args:
        representative_copies (dict): Dictionary with species as keys and representative copy numbers as values
        output_file (str): Path to the output multi_set file
    """
    with open(output_file, "w") as f:
        # Sort species for consistent output
        for species in sorted(representative_copies.keys()):
            copy_number = representative_copies[species]
            # Write species name on separate lines, repeated by copy number
            for _ in range(copy_number):
                f.write(f"{species}\n")
    
    print(f"Multi-set file written to {output_file}")

def parse_arguments():
    """
    Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Analyze gene trees to determine species copy number distributions."
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Path to the input file containing Newick trees"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=False,
        help="Path to the output TSV file for results (optional)"
    )
    
    parser.add_argument(
        "-m", "--multiset",
        required=False,
        help="Path to the output multi_set.txt file (optional)"
    )
    
    parser.add_argument(
        "-e", "--extract",
        choices=["full", "before", "after"],
        default="full",
        help="How to extract species names from leaf labels: 'full' for complete name, "
             "'before' for text before underscore, 'after' for text after underscore (default: full)"
    )
    
    parser.add_argument(
        "-k", "--kernel-width",
        type=int,
        default=2,
        help="Width of the smoothing kernel (default: 2)"
    )
    
    return parser.parse_args()

def main():
    """
    Main function that analyzes trees and outputs results.
    Uses command line arguments for configuration.
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Check that at least one output is specified
    if not args.output and not args.multiset:
        print("ERROR: At least one output file must be specified (--output or --multiset)")
        return
    
    # Analyze trees
    print(f"Analyzing trees from {args.input}...")
    tree_species_counts = analyze_trees_from_file(
        args.input, 
        extract_mode=args.extract
    )
    
    if not tree_species_counts:
        print(f"No valid trees found in {args.input}")
        return
    
    # Get copy number distributions
    print("Calculating species copy number distributions...")
    species_copy_distributions = get_species_copy_distributions(tree_species_counts)
    
    # Get representative copy numbers using smoothing
    print(f"Finding representative copy numbers using kernel smoothing (width={args.kernel_width})...")
    representative_copies = get_representative_copy_numbers(
        species_copy_distributions,
        kernel_width=args.kernel_width
    )
    
    # Print results
    print("\nResults:")
    for species, copy_number in sorted(representative_copies.items()):
        print(f"{species}: {copy_number} copies")
    
    # Write TSV results to file if specified
    if args.output:
        with open(args.output, "w") as f:
            f.write("Species\tRepresentativeCopyNumber\tDistribution\n")
            for species, copy_number in sorted(representative_copies.items()):
                distribution_str = ", ".join([f"{count}:{freq}" for count, freq in 
                                             species_copy_distributions[species].most_common()])
                f.write(f"{species}\t{copy_number}\t{distribution_str}\n")
        print(f"\nTSV results written to {args.output}")
    
    # Write multi_set file if specified
    if args.multiset:
        write_multi_set_file(representative_copies, args.multiset)
    
    print(f"\nTotal species analyzed: {len(representative_copies)}")
    print(f"Total trees analyzed: {len(tree_species_counts)}")

if __name__ == "__main__":
    main()