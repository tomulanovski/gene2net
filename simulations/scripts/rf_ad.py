#!/usr/bin/env python3
"""
Calculate RF distance for each gene tree vs species tree and compute Average Distance (AD).
AD = Average of all RF distances
Gene tree labels: "speciesname_x_x" -> strips to "speciesname"

Usage: python script.py [directory] [output_file]
  directory: Directory containing tree files (default: current directory)
  output_file: Path for output file (default: rf_distance_results.txt in tree directory)

  python script.py /path/to/trees /path/to/output.txt
"""

from ete3 import Tree
import numpy as np
import glob
import os
import sys

def strip_gene_suffix(label):
    """
    Remove the trailing _x_x pattern from gene tree labels.
    E.g., 'RSC01_subgenome_6_0_0' -> 'RSC01_subgenome_6'
          'Helia_0_0' -> 'Helia'
    All gene trees are assumed to have this _x_x suffix.
    """
    parts = label.rsplit('_', 2)  # Split from the right, max 2 splits
    return parts[0]  # Return everything before the last _x_x

def prepare_gene_tree(gene_tree, species_names):
    """Rename gene tree leaves to match species tree."""
    for leaf in gene_tree.iter_leaves():
        leaf.name = strip_gene_suffix(leaf.name)

def main():
    # Get directory from command line argument or use current directory
    if len(sys.argv) > 1:
        tree_dir = sys.argv[1]
    else:
        tree_dir = "."
    
    # Get output file path from command line or use default
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = os.path.join(tree_dir, "rf_distance_results.txt")
    
    if not os.path.isdir(tree_dir):
        print(f"Error: Directory '{tree_dir}' does not exist!")
        return
    
    print(f"Looking for tree files in: {os.path.abspath(tree_dir)}")
    print(f"Output will be saved to: {os.path.abspath(output_file)}\n")
    
    # File paths
    species_tree_path = os.path.join(tree_dir, "s_tree.trees")
    locus_tree_path = os.path.join(tree_dir, "l_trees.trees")
    
    # Check if species tree exists
    if not os.path.exists(species_tree_path):
        print(f"Error: Species tree '{species_tree_path}' not found!")
        return
    
    # Find gene tree files
    all_tree_files = glob.glob(os.path.join(tree_dir, "*.trees"))
    gene_tree_files = [f for f in all_tree_files 
                       if os.path.basename(f) not in ["s_tree.trees", "l_trees.trees"]]
    
    if not gene_tree_files:
        print("No gene tree files found!")
        return
    
    print(f"Found {len(gene_tree_files)} gene tree files\n")
    
    # Load species tree
    species_tree = Tree(species_tree_path)
    species_names = set([leaf.name for leaf in species_tree.iter_leaves()])
    print(f"Species tree has {len(species_names)} leaves\n")
    
    # Calculate RF distance for each gene tree
    rf_distances = []
    rf_results = []  # Store results for each gene tree
    
    print("Processing gene trees...")
    print("-" * 60)
    
    for i, gene_tree_path in enumerate(gene_tree_files, 1):
        try:
            gene_tree = Tree(gene_tree_path)
            prepare_gene_tree(gene_tree, species_names)
            
            # Calculate RF distance
            rf, max_rf = gene_tree.robinson_foulds(species_tree, unrooted_trees=True)[:2]
            norm_rf = rf / max_rf if max_rf > 0 else 0
            norm_rf_percent = norm_rf * 100
            
            rf_distances.append(rf)
            rf_results.append({
                'filename': os.path.basename(gene_tree_path),
                'rf': rf,
                'max_rf': max_rf,
                'norm_rf': norm_rf,
                'norm_rf_percent': norm_rf_percent
            })
            
            # Show progress
            if i <= 10 or len(gene_tree_files) <= 20:
                print(f"{i}. {os.path.basename(gene_tree_path)}: RF = {rf}/{max_rf} ({norm_rf_percent:.1f}%)")
            elif i % 50 == 0:
                print(f"Processed {i}/{len(gene_tree_files)} trees...")
                
        except Exception as e:
            print(f"{i}. {os.path.basename(gene_tree_path)}: ERROR - {e}")
    
    # Calculate Average Distance (AD)
    if rf_distances:
        ad = np.mean(rf_distances)
        norm_rf_values = [r['norm_rf'] for r in rf_results]
        ad_percent = np.mean(norm_rf_values) * 100
        
        print("\n" + "=" * 60)
        print("RESULTS")
        print("=" * 60)
        print(f"Gene trees analyzed: {len(rf_distances)}/{len(gene_tree_files)}")
        print(f"\nAverage Distance (AD): {ad:.2f} ({ad_percent:.2f}%)")
        print(f"\nRF Distance Statistics:")
        print(f"  Mean (AD):  {ad:.2f} ({ad_percent:.2f}%)")
        print(f"  Std Dev:    {np.std(rf_distances):.2f}")
        print(f"  Min:        {np.min(rf_distances):.0f}")
        print(f"  Max:        {np.max(rf_distances):.0f}")
        print(f"  Median:     {np.median(rf_distances):.2f}")
        
        # Save results
        try:
            with open(output_file, 'w') as f:
                f.write("RF Distance Analysis\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Species tree: {species_tree_path}\n")
                f.write(f"Gene trees analyzed: {len(rf_distances)}/{len(gene_tree_files)}\n\n")
                f.write(f"Average Distance (AD): {ad:.2f} ({ad_percent:.2f}%)\n\n")
                f.write("RF Distance Statistics:\n")
                f.write(f"  Mean (AD):     {ad:.2f} ({ad_percent:.2f}%)\n")
                f.write(f"  Std Dev:       {np.std(rf_distances):.2f}\n")
                f.write(f"  Min:           {np.min(rf_distances):.0f}\n")
                f.write(f"  Max:           {np.max(rf_distances):.0f}\n")
                f.write(f"  Median:        {np.median(rf_distances):.2f}\n\n")
                
                f.write("=" * 60 + "\n")
                f.write("Individual Gene Tree RF Distances:\n")
                f.write("=" * 60 + "\n\n")
                for i, result in enumerate(rf_results, 1):
                    f.write(f"{i:3d}. {result['filename']:<40s} RF: {result['rf']:3.0f}/{result['max_rf']:3.0f}  ({result['norm_rf_percent']:5.1f}%)\n")
            
            print(f"\nResults saved to: {output_file}")
        except Exception as e:
            print(f"\nError saving results: {e}")
    else:
        print("\nNo gene trees successfully analyzed!")

if __name__ == "__main__":
    main()