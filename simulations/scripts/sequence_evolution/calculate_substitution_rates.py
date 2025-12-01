#!/usr/bin/env python3
"""
Calculate substitution rates from gene tree heights and TimeTree species trees.
Uses species tree heights (in millions of years) to normalize gene tree heights.
Outputs: CSV statistics, PNG visualization, and pickle file for sampling.
"""

import pickle
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from Bio import Phylo

# Input/Output paths
INPUT_DIR = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions')
OUTPUT_DIR = INPUT_DIR

# Species tree paths (TimeTree dated trees)
SPECIES_TREES = {
    'Zhao_2021': '/groups/itay_mayrose/tomulanovski/gene2net/papers/Zhao_2021/timetree_zhao.nwk',
    'Ren_2024': '/groups/itay_mayrose/tomulanovski/gene2net/papers/Ren_2024/timetree_ren.nwk',
    'Morales_Briones_2021': '/groups/itay_mayrose/tomulanovski/gene2net/papers/Morales-Briones_2021/timetree_morales.nwk'
}

# Generation times (in years) - all set to 1 year
GENERATION_TIMES = {
    'Zhao_2021': 1.0,
    'Ren_2024': 1.0,
    'Morales_Briones_2021': 1.0
}


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate substitution rates from gene trees and species trees',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--gen-time-zhao',
        type=float,
        default=GENERATION_TIMES['Zhao_2021'],
        help='Generation time in years for Zhao_2021 dataset'
    )
    
    parser.add_argument(
        '--gen-time-ren',
        type=float,
        default=GENERATION_TIMES['Ren_2024'],
        help='Generation time in years for Ren_2024 dataset'
    )
    
    parser.add_argument(
        '--gen-time-morales',
        type=float,
        default=GENERATION_TIMES['Morales_Briones_2021'],
        help='Generation time in years for Morales_Briones_2021 dataset'
    )
    
    return parser.parse_args()


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
        print(f"Error calculating tree height for {tree_file}: {e}")
        return None


def calculate_species_tree_heights():
    """
    Calculate species tree heights for each dataset.
    TimeTree branch lengths are in millions of years (MY).
    Returns dict with species tree heights in MY.
    """
    species_heights = {}
    
    print("Calculating species tree heights from TimeTree...")
    print("="*60)
    
    for dataset_name, tree_path in SPECIES_TREES.items():
        tree_path_obj = Path(tree_path)
        
        if not tree_path_obj.exists():
            print(f"ERROR: Species tree not found for {dataset_name}: {tree_path}")
            print(f"Please ensure the TimeTree file exists at this location.")
            return None
        
        height = calculate_tree_height(tree_path)
        
        if height is not None:
            species_heights[dataset_name] = height
            print(f"{dataset_name}:")
            print(f"  Species tree: {tree_path}")
            print(f"  Crown age (tree height): {height:.2f} million years (MY)")
        else:
            print(f"ERROR: Could not calculate height for {dataset_name} species tree")
            return None
    
    print("="*60 + "\n")
    return species_heights


def calculate_substitution_rates(tree_heights_dict, species_heights, generation_times):
    """
    Calculate substitution rates from gene tree heights, species tree heights, and generation times.
    
    Formula: rate = gene_tree_height / (species_tree_height_MY × 1,000,000 / generation_time_years)
              rate = gene_tree_height × generation_time / (species_tree_height × 1,000,000)
    
    Units: substitutions per site per generation
    """
    rates_dict = {}
    
    print("Calculating substitution rates per generation...")
    print("="*60)
    
    for dataset_name, gene_heights_data in tree_heights_dict.items():
        species_height_my = species_heights.get(dataset_name)
        gen_time = generation_times.get(dataset_name)
        
        if species_height_my is None or gen_time is None:
            print(f"Warning: Missing data for {dataset_name}, skipping...")
            continue
        
        # Convert species tree height from MY to generations
        species_height_generations = (species_height_my * 1_000_000) / gen_time
        
        rates = []
        for entry in gene_heights_data:
            gene_height = entry['height']
            
            # Calculate rate per generation
            rate = gene_height / species_height_generations
            
            rates.append({
                'rate': rate,
                'gene_tree_height': gene_height,
                'species_tree_height_MY': species_height_my,
                'species_tree_height_generations': species_height_generations,
                'generation_time_years': gen_time,
                'source_file': entry['source_file']
            })
        
        rates_dict[dataset_name] = rates
        
        rate_values = [r['rate'] for r in rates]
        print(f"{dataset_name}: Calculated {len(rates)} substitution rates")
        print(f"  Species crown age: {species_height_my:.2f} MY")
        print(f"  Generation time: {gen_time} years")
        print(f"  Divergence time: {species_height_generations:.2e} generations")
        print(f"  Rate range: {min(rate_values):.2e} - {max(rate_values):.2e} per site per generation")
        print(f"  Rate mean: {np.mean(rate_values):.2e}, median: {np.median(rate_values):.2e}")
    
    print("="*60 + "\n")
    return rates_dict


def calculate_rate_statistics(rates_data, dataset_name):
    """Calculate summary statistics for substitution rates."""
    if not rates_data:
        return None
    
    rates = [d['rate'] for d in rates_data]
    
    stats = {
        'Dataset': dataset_name,
        'N_genes': len(rates),
        'Min_rate': np.min(rates),
        'Max_rate': np.max(rates),
        'Mean_rate': np.mean(rates),
        'Median_rate': np.median(rates),
        'Std_rate': np.std(rates),
        'Q1_rate': np.percentile(rates, 25),
        'Q3_rate': np.percentile(rates, 75)
    }
    
    # Only add dataset-specific info if not Combined
    if dataset_name != 'Combined':
        stats['Species_crown_age_MY'] = rates_data[0]['species_tree_height_MY'] if rates_data else None
        stats['Generation_time_years'] = rates_data[0]['generation_time_years'] if rates_data else None
        stats['Divergence_time_generations'] = rates_data[0]['species_tree_height_generations'] if rates_data else None
    else:
        # For Combined, these don't make sense so leave as None or indicate mixed
        stats['Species_crown_age_MY'] = 'Mixed'
        stats['Generation_time_years'] = 'Mixed'
        stats['Divergence_time_generations'] = 'Mixed'
    
    return stats


def create_rate_visualizations(all_rates, output_path):
    """Create visualization of substitution rate distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Prepare data for plotting
    datasets = []
    rates = []
    for dataset_name, rates_data in all_rates.items():
        dataset_rates = [d['rate'] for d in rates_data]
        datasets.extend([dataset_name] * len(dataset_rates))
        rates.extend(dataset_rates)
    
    df = pd.DataFrame({'Dataset': datasets, 'Rate': rates})
    
    # 1. Histograms for each dataset
    ax1 = axes[0, 0]
    for dataset_name, rates_data in all_rates.items():
        dataset_rates = [d['rate'] for d in rates_data]
        ax1.hist(dataset_rates, alpha=0.6, label=dataset_name, bins=30)
    ax1.set_xlabel('Substitution Rate (per site per generation)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Substitution Rate Distribution by Dataset', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    
    # 2. Box plots
    ax2 = axes[0, 1]
    df.boxplot(column='Rate', by='Dataset', ax=ax2)
    ax2.set_xlabel('Dataset', fontsize=12)
    ax2.set_ylabel('Substitution Rate (per site per generation)', fontsize=12)
    ax2.set_title('Substitution Rate Distributions', fontsize=14, fontweight='bold')
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    plt.sca(ax2)
    plt.xticks(rotation=45, ha='right')
    
    # 3. Combined histogram
    ax3 = axes[1, 0]
    all_vals = rates
    ax3.hist(all_vals, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax3.axvline(np.mean(all_vals), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {np.mean(all_vals):.2e}')
    ax3.axvline(np.median(all_vals), color='green', linestyle='--', linewidth=2, 
                label=f'Median: {np.median(all_vals):.2e}')
    ax3.set_xlabel('Substitution Rate (per site per generation)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_title('Combined Substitution Rate Distribution', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    
    # 4. Violin plot
    ax4 = axes[1, 1]
    dataset_rates = [
        [d['rate'] for d in all_rates[ds]] 
        for ds in all_rates.keys()
    ]
    parts = ax4.violinplot(dataset_rates,
                           positions=range(len(all_rates)),
                           showmeans=True, showmedians=True)
    ax4.set_xticks(range(len(all_rates)))
    ax4.set_xticklabels(all_rates.keys(), rotation=45, ha='right')
    ax4.set_ylabel('Substitution Rate (per site per generation)', fontsize=12)
    ax4.set_title('Substitution Rate Distributions (Violin Plot)', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_path}")
    plt.close()


def main():
    """Main execution function."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Update generation times from arguments
    generation_times = {
        'Zhao_2021': args.gen_time_zhao,
        'Ren_2024': args.gen_time_ren,
        'Morales_Briones_2021': args.gen_time_morales
    }
    
    print("\n" + "="*60)
    print("Substitution Rate Calculation from TimeTree")
    print("="*60)
    print("Generation times (years):")
    for dataset, gen_time in generation_times.items():
        print(f"  {dataset}: {gen_time}")
    print("="*60 + "\n")
    
    # Step 1: Calculate species tree heights
    species_heights = calculate_species_tree_heights()
    if species_heights is None:
        print("ERROR: Could not calculate species tree heights. Exiting.")
        return 1
    
    # Step 2: Load gene tree heights
    heights_file = INPUT_DIR / 'tree_heights.pkl'
    if not heights_file.exists():
        print(f"ERROR: Gene tree heights file not found: {heights_file}")
        print("Please run extract_tree_heights.py first!")
        return 1
    
    print("Loading gene tree heights...")
    with open(heights_file, 'rb') as f:
        tree_heights = pickle.load(f)
    print(f"Loaded tree heights for {sum(len(v) for v in tree_heights.values())} gene trees\n")
    
    # Step 3: Calculate substitution rates
    all_rates = calculate_substitution_rates(tree_heights, species_heights, generation_times)
    
    if not all_rates:
        print("ERROR: No rates calculated. Please check inputs.")
        return 1
    
    # Step 4: Calculate statistics
    stats_list = []
    for dataset_name, rates_data in all_rates.items():
        stats = calculate_rate_statistics(rates_data, dataset_name)
        stats_list.append(stats)
    
    # Combined statistics
    combined_rates = [d for rates in all_rates.values() for d in rates]
    combined_stats = calculate_rate_statistics(combined_rates, 'Combined')
    stats_list.append(combined_stats)
    
    # Step 5: Save statistics to CSV
    csv_path = OUTPUT_DIR / 'substitution_rates.csv'
    df_stats = pd.DataFrame(stats_list)
    df_stats.to_csv(csv_path, index=False)
    print(f"Statistics saved to {csv_path}")
    
    # Step 6: Save raw data for sampling
    pkl_path = OUTPUT_DIR / 'substitution_rates.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_rates, f)
    print(f"Raw data saved to {pkl_path}")
    
    # Step 7: Create visualizations
    png_path = OUTPUT_DIR / 'substitution_rates.png'
    create_rate_visualizations(all_rates, png_path)
    
    # Print summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    
    # Print detailed table
    print("\nPer-Dataset Statistics:")
    print(df_stats.to_string(index=False))
    
    print("\n" + "="*60)
    print("Usage Instructions")
    print("="*60)
    print("\nTo sample substitution rates in your simulation code:")
    print("  import pickle, random")
    print(f"  with open('{pkl_path}', 'rb') as f:")
    print("      rates = pickle.load(f)")
    print("  # Sample from combined distribution:")
    print("  combined = [r['rate'] for ds in rates.values() for r in ds]")
    print("  sampled_rate = random.choice(combined)")
    print("  # Use sampled_rate as -su parameter in SimPhy")
    print("\nNote: Rates are in substitutions per site per generation")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
