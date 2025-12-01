#!/usr/bin/env python3
"""
Print statistics of genes with extreme alpha values.
"""

import pickle
import re
import pandas as pd
import numpy as np
from pathlib import Path

# Paths
GTR_PARAMS_PATH = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/gtr_parameters_unfiltered.pkl')
TREE_HEIGHTS_PATH = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/tree_heights.pkl')


def extract_alignment_info_from_iqtree(iqtree_file):
    """
    Extract alignment metadata from IQ-TREE file header.
    """
    try:
        with open(iqtree_file, 'r') as f:
            content = f.read(2000)
        
        info = {}
        
        # Number of sequences and sites
        pattern = r'Input data:\s+(\d+)\s+sequences?\s+with\s+(\d+)\s+nucleotide sites'
        match = re.search(pattern, content)
        if match:
            info['n_sequences'] = int(match.group(1))
            info['alignment_length'] = int(match.group(2))
        
        # Number of constant sites
        pattern = r'Number of constant sites:\s+(\d+)'
        match = re.search(pattern, content)
        if match:
            info['constant_sites'] = int(match.group(1))
        
        # Number of parsimony informative sites
        pattern = r'Number of parsimony informative sites:\s+(\d+)'
        match = re.search(pattern, content)
        if match:
            info['informative_sites'] = int(match.group(1))
        
        # Calculate percentages and gaps
        if 'alignment_length' in info:
            if 'constant_sites' in info:
                info['constant_percentage'] = (info['constant_sites'] / info['alignment_length']) * 100
            if 'informative_sites' in info:
                info['informative_percentage'] = (info['informative_sites'] / info['alignment_length']) * 100
            
            # Calculate gaps percentage
            if 'constant_sites' in info and 'informative_sites' in info:
                variable_sites = info['alignment_length'] - info['constant_sites']
                non_informative_variable = variable_sites - info['informative_sites']
                # Approximation: non-informative variable sites often contain gaps
                info['gaps_percentage'] = (non_informative_variable / info['alignment_length']) * 100
        
        return info
        
    except Exception as e:
        return {}


def print_statistics(df, label, alpha_threshold=None):
    """Print detailed statistics for a dataframe."""
    print(f"\n{'='*70}")
    print(f"{label}")
    if alpha_threshold:
        print(f"(alpha > {alpha_threshold})")
    print(f"{'='*70}")
    print(f"Number of genes: {len(df)}\n")
    
    features = {
        'Alpha': 'alpha',
        'Alignment Length (bp)': 'alignment_length',
        'Number of Sequences': 'n_sequences',
        'Constant Sites (%)': 'constant_percentage',
        'Informative Sites (%)': 'informative_percentage',
        'Gaps/Variable Non-informative (%)': 'gaps_percentage',
        'Tree Height': 'tree_height'
    }
    
    for display_name, col_name in features.items():
        if col_name in df.columns:
            data = df[col_name].dropna()
            if len(data) > 0:
                print(f"{display_name}:")
                print(f"  Mean:   {data.mean():>10.2f}")
                print(f"  Median: {data.median():>10.2f}")
                print(f"  Std:    {data.std():>10.2f}")
                print(f"  Min:    {data.min():>10.2f}")
                print(f"  Max:    {data.max():>10.2f}")
                print(f"  Q1:     {data.quantile(0.25):>10.2f}")
                print(f"  Q3:     {data.quantile(0.75):>10.2f}")
                print()


def main():
    """Main analysis function."""
    print("="*70)
    print("Statistics of Genes with Extreme Alpha Values")
    print("="*70 + "\n")
    
    # Load GTR parameters
    print("Loading GTR parameters...")
    with open(GTR_PARAMS_PATH, 'rb') as f:
        gtr_params = pickle.load(f)
    
    # Load tree heights
    print("Loading tree heights...")
    with open(TREE_HEIGHTS_PATH, 'rb') as f:
        tree_heights = pickle.load(f)
    
    # Create mapping of source_file to tree height
    height_map = {}
    for dataset, heights_list in tree_heights.items():
        for entry in heights_list:
            source = entry['source_file']
            height_map[source] = entry['height']
    
    # Collect all genes with their features
    print("Extracting features from IQ-TREE files...")
    all_genes = []
    
    for dataset_name, params_list in gtr_params.items():
        print(f"  Processing {dataset_name}...")
        for params in params_list:
            if 'alpha' in params and params['alpha'] is not None:
                source_file = params.get('source_file', '')
                
                # Extract info from IQ-TREE file
                iqtree_info = extract_alignment_info_from_iqtree(source_file)
                
                # Get tree height
                tree_height = height_map.get(source_file)
                
                gene_data = {
                    'dataset': dataset_name,
                    'alpha': params['alpha'],
                    'tree_height': tree_height,
                    'source_file': source_file,
                    **iqtree_info,
                }
                all_genes.append(gene_data)
    
    df_genes = pd.DataFrame(all_genes)
    
    # Define extreme alpha threshold
    extreme_threshold = 3.0
    df_extreme = df_genes[df_genes['alpha'] > extreme_threshold]
    df_normal = df_genes[df_genes['alpha'] <= extreme_threshold]
    
    # Print statistics for extreme alpha genes
    print_statistics(df_extreme, "EXTREME ALPHA GENES", extreme_threshold)
    
    # Print statistics for normal alpha genes (for comparison)
    print_statistics(df_normal, "NORMAL ALPHA GENES (for comparison)", None)
    
    # Print comparison summary
    print(f"\n{'='*70}")
    print("COMPARISON SUMMARY")
    print(f"{'='*70}\n")
    
    print(f"Total genes: {len(df_genes)}")
    print(f"Extreme alpha genes (α > {extreme_threshold}): {len(df_extreme)} ({100*len(df_extreme)/len(df_genes):.1f}%)")
    print(f"Normal alpha genes (α ≤ {extreme_threshold}): {len(df_normal)} ({100*len(df_normal)/len(df_genes):.1f}%)")
    
    print("\n" + "="*70)
    print("Done!")
    print("="*70)


if __name__ == '__main__':
    main()
