#!/usr/bin/env python3
"""
Infer Tree Height Distribution from Real Datasets

This script:
1. Reads IQ-TREE output files (.treefile) from multiple datasets
2. Calculates average root-to-leaf distance (tree height) for each gene tree
3. Fits candidate distributions (Gamma, Log-normal, Normal)
4. Selects best distribution using AIC and KS test
5. Saves parameters for sampling in simulations

Tree height = average evolutionary distance from root to leaves (substitutions/site)
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
from ete3 import Tree

# Configuration
BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/papers"
DATASETS = {
    "Zhao_2021": f"{BASE_DIR}/Zhao_2021/gene_trees",
    "Diaz-Perez_2018": f"{BASE_DIR}/Diaz-Perez_2018/gene_trees",
    "Sessa_2012b": f"{BASE_DIR}/Sessa_2012b/gene_trees",
    "Hori_2014": f"{BASE_DIR}/Hori_2014/gene_trees"
}

OUTPUT_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/distribution_inferences"
TREE_EXTENSION = ".treefile"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)


def calculate_tree_height(tree_file):
    """
    Calculate average root-to-leaf distance (tree height) from a Newick tree file.
    
    Args:
        tree_file: Path to .treefile (Newick format with branch lengths)
    
    Returns:
        float: Average distance from root to all leaves (substitutions/site)
    """
    try:
        # Load tree with branch lengths (format=1 for Newick with branch lengths)
        tree = Tree(tree_file, format=1)
        
        # Get all leaves
        leaves = list(tree.iter_leaves())
        
        if len(leaves) == 0:
            return None
        
        # Calculate distance from root to each leaf
        distances = []
        for leaf in leaves:
            # get_distance returns branch length sum from node to leaf
            distance = tree.get_distance(tree, leaf)
            distances.append(distance)
        
        # Return average root-to-leaf distance
        return np.mean(distances)
    
    except Exception as e:
        print(f"Warning: Failed to parse {tree_file}: {e}")
        return None


def collect_tree_heights(datasets):
    """
    Collect tree heights from all datasets.
    
    Args:
        datasets: Dictionary mapping dataset names to tree directories
    
    Returns:
        DataFrame with columns: dataset, gene, tree_height
    """
    all_data = []
    
    for dataset_name, tree_dir in datasets.items():
        print(f"\nProcessing {dataset_name}...")
        
        if not os.path.exists(tree_dir):
            print(f"  Warning: Directory not found: {tree_dir}")
            continue
        
        # Find all .treefile files
        tree_files = list(Path(tree_dir).glob(f"*{TREE_EXTENSION}"))
        print(f"  Found {len(tree_files)} tree files")
        
        for tree_file in tree_files:
            height = calculate_tree_height(str(tree_file))
            
            if height is not None:
                all_data.append({
                    'dataset': dataset_name,
                    'gene': tree_file.stem.replace(TREE_EXTENSION, ''),
                    'tree_height': height
                })
    
    return pd.DataFrame(all_data)


def fit_distribution(data, dist_name):
    """
    Fit a distribution to the data using Maximum Likelihood Estimation.
    
    Args:
        data: Array of tree heights
        dist_name: Name of scipy.stats distribution
    
    Returns:
        dict: Contains parameters, AIC, KS test results
    """
    dist = getattr(stats, dist_name)
    
    # Fit distribution parameters using MLE
    params = dist.fit(data)
    
    # Calculate log-likelihood
    log_likelihood = np.sum(dist.logpdf(data, *params))
    
    # Calculate AIC: 2k - 2*log(L)
    k = len(params)  # number of parameters
    aic = 2 * k - 2 * log_likelihood
    
    # Kolmogorov-Smirnov test
    ks_statistic, ks_pvalue = stats.kstest(data, lambda x: dist.cdf(x, *params))
    
    return {
        'distribution': dist_name,
        'parameters': params,
        'aic': aic,
        'ks_statistic': ks_statistic,
        'ks_pvalue': ks_pvalue,
        'log_likelihood': log_likelihood
    }


def plot_results(df, fits, output_file):
    """
    Create visualization of tree height distribution and fits.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Tree Height Distribution Analysis', fontsize=16, fontweight='bold')
    
    heights = df['tree_height'].values
    
    # Panel 1: Histogram with fitted distributions
    ax = axes[0, 0]
    ax.hist(heights, bins=50, density=True, alpha=0.6, color='skyblue', edgecolor='black')
    
    x_range = np.linspace(heights.min(), heights.max(), 1000)
    colors = {'gamma': 'red', 'lognorm': 'green', 'norm': 'blue'}
    
    for fit in fits:
        dist = getattr(stats, fit['distribution'])
        y = dist.pdf(x_range, *fit['parameters'])
        label = f"{fit['distribution'].capitalize()} (AIC={fit['aic']:.1f})"
        ax.plot(x_range, y, label=label, linewidth=2, color=colors.get(fit['distribution'], 'black'))
    
    ax.set_xlabel('Tree Height (substitutions/site)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('Distribution Fits', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # Panel 2: Q-Q plots for best fit
    best_fit = min(fits, key=lambda x: x['aic'])
    ax = axes[0, 1]
    dist = getattr(stats, best_fit['distribution'])
    stats.probplot(heights, dist=dist, sparams=best_fit['parameters'], plot=ax)
    ax.set_title(f'Q-Q Plot: {best_fit["distribution"].capitalize()} (Best Fit)', 
                 fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Panel 3: Distribution by dataset
    ax = axes[1, 0]
    datasets = df['dataset'].unique()
    positions = range(len(datasets))
    
    box_data = [df[df['dataset'] == ds]['tree_height'].values for ds in datasets]
    bp = ax.boxplot(box_data, positions=positions, labels=datasets, patch_artist=True)
    
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    
    ax.set_xlabel('Dataset', fontsize=11)
    ax.set_ylabel('Tree Height (substitutions/site)', fontsize=11)
    ax.set_title('Tree Heights by Dataset', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Panel 4: Summary statistics table
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_stats = df.groupby('dataset')['tree_height'].agg(['count', 'mean', 'std', 'min', 'max'])
    summary_stats.loc['Total'] = [
        df.shape[0],
        df['tree_height'].mean(),
        df['tree_height'].std(),
        df['tree_height'].min(),
        df['tree_height'].max()
    ]
    
    table_data = []
    table_data.append(['Dataset', 'N', 'Mean', 'Std', 'Min', 'Max'])
    for idx, row in summary_stats.iterrows():
        table_data.append([
            idx,
            f"{int(row['count'])}",
            f"{row['mean']:.3f}",
            f"{row['std']:.3f}",
            f"{row['min']:.3f}",
            f"{row['max']:.3f}"
        ])
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.25, 0.12, 0.15, 0.15, 0.15, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Style header row
    for i in range(6):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style total row
    for i in range(6):
        table[(len(table_data)-1, i)].set_facecolor('#E8F5E9')
        table[(len(table_data)-1, i)].set_text_props(weight='bold')
    
    ax.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")


def main():
    print("="*70)
    print("TREE HEIGHT DISTRIBUTION INFERENCE")
    print("="*70)
    
    # Collect tree heights from all datasets
    print("\n1. Collecting tree heights from datasets...")
    df = collect_tree_heights(DATASETS)
    
    if df.empty:
        print("Error: No tree files found!")
        return
    
    print(f"\nTotal trees collected: {len(df)}")
    print(f"Datasets: {df['dataset'].unique()}")
    print(f"\nTree height range: {df['tree_height'].min():.4f} - {df['tree_height'].max():.4f}")
    print(f"Mean tree height: {df['tree_height'].mean():.4f} ± {df['tree_height'].std():.4f}")
    
    # Fit candidate distributions
    print("\n2. Fitting candidate distributions...")
    heights = df['tree_height'].values
    
    distributions = ['gamma', 'lognorm', 'norm']
    fits = []
    
    for dist_name in distributions:
        print(f"\n   Fitting {dist_name}...")
        fit_result = fit_distribution(heights, dist_name)
        fits.append(fit_result)
        print(f"     AIC: {fit_result['aic']:.2f}")
        print(f"     KS p-value: {fit_result['ks_pvalue']:.4f}")
    
    # Select best fit (lowest AIC)
    best_fit = min(fits, key=lambda x: x['aic'])
    print(f"\n3. Best distribution: {best_fit['distribution'].upper()}")
    print(f"   AIC: {best_fit['aic']:.2f}")
    print(f"   KS p-value: {best_fit['ks_pvalue']:.4f}")
    print(f"   Parameters: {best_fit['parameters']}")
    
    # Prepare output data
    output_data = {
        'best_distribution': {
            'name': best_fit['distribution'],
            'parameters': list(best_fit['parameters']),
            'aic': best_fit['aic'],
            'ks_pvalue': best_fit['ks_pvalue']
        },
        'all_distributions': [
            {
                'name': f['distribution'],
                'parameters': list(f['parameters']),
                'aic': f['aic'],
                'ks_pvalue': f['ks_pvalue']
            }
            for f in fits
        ],
        'summary_statistics': {
            'n_trees': int(len(df)),
            'mean': float(df['tree_height'].mean()),
            'std': float(df['tree_height'].std()),
            'min': float(df['tree_height'].min()),
            'max': float(df['tree_height'].max()),
            'median': float(df['tree_height'].median())
        },
        'per_dataset_stats': df.groupby('dataset')['tree_height'].agg(['count', 'mean', 'std', 'min', 'max']).to_dict('index')
    }
    
    # Save outputs
    print("\n4. Saving outputs...")
    
    # JSON parameters
    json_file = os.path.join(OUTPUT_DIR, 'tree_height_distribution_params.json')
    with open(json_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"   Parameters saved to: {json_file}")
    
    # CSV summary - per dataset statistics like alignment length format
    dataset_stats = df.groupby('dataset')['tree_height'].agg([
        ('N_trees', 'count'),
        ('Mean_height', 'mean'),
        ('Std_height', 'std'),
        ('Median_height', 'median'),
        ('Min_height', 'min'),
        ('Max_height', 'max')
    ]).reset_index()
    dataset_stats.columns = ['Dataset', 'N_trees', 'Mean_height', 'Std_height', 'Median_height', 'Min_height', 'Max_height']
    
    # Add overall row
    overall_row = pd.DataFrame([{
        'Dataset': 'OVERALL',
        'N_trees': len(df),
        'Mean_height': df['tree_height'].mean(),
        'Std_height': df['tree_height'].std(),
        'Median_height': df['tree_height'].median(),
        'Min_height': df['tree_height'].min(),
        'Max_height': df['tree_height'].max()
    }])
    
    summary_df = pd.concat([dataset_stats, overall_row], ignore_index=True)
    
    csv_file = os.path.join(OUTPUT_DIR, 'tree_height_distribution_summary.csv')
    summary_df.to_csv(csv_file, index=False)
    print(f"   Summary saved to: {csv_file}")
    
    # Visualization
    plot_file = os.path.join(OUTPUT_DIR, 'tree_height_distribution_plot.png')
    plot_results(df, fits, plot_file)
    
    # Usage example
    print("\n" + "="*70)
    print("HOW TO SAMPLE TREE HEIGHTS IN YOUR SIMULATIONS:")
    print("="*70)
    print("""
import json
from scipy import stats

# Load parameters
params = json.load(open('/groups/itay_mayrose/tomulanovski/gene2net/distribution_inferences/tree_height_distribution_params.json'))

# Get best distribution
dist_name = params['best_distribution']['name']
dist_params = params['best_distribution']['parameters']

# Create distribution object
dist = getattr(stats, dist_name)

# Sample a tree height
tree_height = dist.rvs(*dist_params)
print(f"Sampled tree height: {tree_height:.4f} substitutions/site")
""")
    
    print("="*70)
    print("DONE!")
    print("="*70)


if __name__ == "__main__":
    main()
