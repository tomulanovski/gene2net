#!/usr/bin/env python3
"""
Extract tree heights (average root-to-leaf distances) from gene trees.
These will be used to calculate substitution rates once divergence times are available.
Outputs: CSV statistics, PNG visualization, and pickle file for sampling.
"""

import os
import pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from io import StringIO

# Dataset paths
DATASETS = {
    'Zhao_2021': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Zhao_et_al_2021/gene_trees',
    'Ren_2024': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Ren_et_al_2024/gene_trees',
    'Morales_Briones_2021': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Morales-Briones_et_al_2021/gene_trees'
}

OUTPUT_DIR = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions')


def calculate_tree_height(tree_file):
    """
    Calculate average root-to-leaf distance (tree height) from a phylogenetic tree.
    Returns the mean path length from root to all leaves.
    """
    try:
        tree = Phylo.read(tree_file, 'newick')
        
        # Get all terminal (leaf) nodes
        terminals = tree.get_terminals()
        
        if not terminals:
            return None
        
        # Calculate distance from root to each leaf
        root_to_leaf_distances = []
        for terminal in terminals:
            distance = tree.distance(tree.root, terminal)
            root_to_leaf_distances.append(distance)
        
        # Return average root-to-leaf distance
        return np.mean(root_to_leaf_distances)
        
    except Exception as e:
        print(f"Error processing {tree_file}: {e}")
        return None


def collect_tree_heights(base_dir):
    """Collect tree heights from all .treefile files in subdirectories."""
    tree_heights = []
    
    base_path = Path(base_dir)
    
    # Find all .treefile files in subdirectories
    tree_files = list(base_path.rglob('*.treefile'))
    
    print(f"  Found {len(tree_files)} .treefile files")
    
    for tree_file in tree_files:
        height = calculate_tree_height(tree_file)
        if height is not None:
            tree_heights.append({
                'height': height,
                'source_file': str(tree_file)
            })
    
    return tree_heights


def calculate_statistics(heights_data, dataset_name):
    """Calculate summary statistics for tree heights."""
    if not heights_data:
        return None
    
    heights = [d['height'] for d in heights_data]
    
    return {
        'Dataset': dataset_name,
        'N_genes': len(heights),
        'Min_height': np.min(heights),
        'Max_height': np.max(heights),
        'Mean_height': np.mean(heights),
        'Median_height': np.median(heights),
        'Std_height': np.std(heights),
        'Q1_height': np.percentile(heights, 25),
        'Q3_height': np.percentile(heights, 75)
    }


def create_visualizations(all_heights, output_path):
    """Create visualization of tree height distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Prepare data for plotting
    datasets = []
    heights = []
    for dataset_name, heights_data in all_heights.items():
        dataset_heights = [d['height'] for d in heights_data]
        datasets.extend([dataset_name] * len(dataset_heights))
        heights.extend(dataset_heights)
    
    df = pd.DataFrame({'Dataset': datasets, 'Height': heights})
    
    # 1. Histograms for each dataset
    ax1 = axes[0, 0]
    for dataset_name, heights_data in all_heights.items():
        dataset_heights = [d['height'] for d in heights_data]
        ax1.hist(dataset_heights, alpha=0.6, label=dataset_name, bins=30)
    ax1.set_xlabel('Tree Height (substitutions/site)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Tree Height Distribution by Dataset', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Box plots
    ax2 = axes[0, 1]
    df.boxplot(column='Height', by='Dataset', ax=ax2)
    ax2.set_xlabel('Dataset', fontsize=12)
    ax2.set_ylabel('Tree Height (substitutions/site)', fontsize=12)
    ax2.set_title('Tree Height Distributions', fontsize=14, fontweight='bold')
    plt.sca(ax2)
    plt.xticks(rotation=45, ha='right')
    
    # 3. Combined histogram
    ax3 = axes[1, 0]
    all_vals = heights
    ax3.hist(all_vals, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax3.axvline(np.mean(all_vals), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {np.mean(all_vals):.4f}')
    ax3.axvline(np.median(all_vals), color='green', linestyle='--', linewidth=2, 
                label=f'Median: {np.median(all_vals):.4f}')
    ax3.set_xlabel('Tree Height (substitutions/site)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Combined Tree Height Distribution', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Violin plot
    ax4 = axes[1, 1]
    dataset_heights = [
        [d['height'] for d in all_heights[ds]] 
        for ds in all_heights.keys()
    ]
    parts = ax4.violinplot(dataset_heights,
                           positions=range(len(all_heights)),
                           showmeans=True, showmedians=True)
    ax4.set_xticks(range(len(all_heights)))
    ax4.set_xticklabels(all_heights.keys(), rotation=45, ha='right')
    ax4.set_ylabel('Tree Height (substitutions/site)', fontsize=12)
    ax4.set_title('Tree Height Distributions (Violin Plot)', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_path}")
    plt.close()


def main():
    """Main execution function."""
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Collect tree heights from all datasets
    print("Collecting tree heights from gene trees...")
    all_heights = {}
    stats_list = []
    
    for dataset_name, base_dir in DATASETS.items():
        print(f"\nProcessing {dataset_name}...")
        heights_data = collect_tree_heights(base_dir)
        
        if heights_data:
            all_heights[dataset_name] = heights_data
            stats = calculate_statistics(heights_data, dataset_name)
            stats_list.append(stats)
            heights = [d['height'] for d in heights_data]
            print(f"  Processed {len(heights)} trees")
            print(f"  Height range: {min(heights):.6f} - {max(heights):.6f} substitutions/site")
            print(f"  Mean: {np.mean(heights):.6f}, Median: {np.median(heights):.6f}")
        else:
            print(f"  No valid tree files found in {base_dir}")
    
    # Calculate combined statistics
    combined_heights = [d for heights_data in all_heights.values() for d in heights_data]
    combined_stats = calculate_statistics(combined_heights, 'Combined')
    stats_list.append(combined_stats)
    
    # Save statistics to CSV
    csv_path = OUTPUT_DIR / 'tree_heights.csv'
    df_stats = pd.DataFrame(stats_list)
    df_stats.to_csv(csv_path, index=False)
    print(f"\nStatistics saved to {csv_path}")
    
    # Save raw data for sampling
    pkl_path = OUTPUT_DIR / 'tree_heights.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_heights, f)
    print(f"Raw data saved to {pkl_path}")
    
    # Create visualizations
    png_path = OUTPUT_DIR / 'tree_heights.png'
    create_visualizations(all_heights, png_path)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(df_stats.to_string(index=False))
    
    print("\n" + "="*60)
    print("NEXT STEPS: Calculating Substitution Rates")
    print("="*60)
    print("Once you have divergence times (in generations) for each dataset:")
    print("  substitution_rate = tree_height / divergence_time")
    print("\nTo use the tree heights:")
    print("  import pickle")
    print(f"  with open('{pkl_path}', 'rb') as f:")
    print("      heights = pickle.load(f)")
    print("  # For a specific dataset:")
    print("  zhao_heights = [d['height'] for d in heights['Zhao_2021']]")
    print("  # Calculate rates (example with divergence_time = 1000000):")
    print("  rates = [h / 1000000 for h in zhao_heights]")


if __name__ == '__main__':
    main()
