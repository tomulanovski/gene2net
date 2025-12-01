#!/usr/bin/env python3
"""
Report number of sequences per dataset from GTR parameters pickle file.
"""

import pickle
import numpy as np
from pathlib import Path

pkl_path = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/gtr_parameters_all.pkl')

print("="*70)
print("NUMBER OF SEQUENCES PER DATASET")
print("="*70)

with open(pkl_path, 'rb') as f:
    all_params = pickle.load(f)

for dataset_name, params_list in all_params.items():
    n_sequences_list = [p['n_sequences'] for p in params_list if 'n_sequences' in p]
    
    if n_sequences_list:
        print(f"\n{dataset_name}:")
        print(f"  N genes: {len(n_sequences_list)}")
        print(f"  Mean sequences per gene: {np.mean(n_sequences_list):.1f}")
        print(f"  Median sequences per gene: {np.median(n_sequences_list):.1f}")
        print(f"  Range: {int(np.min(n_sequences_list))} - {int(np.max(n_sequences_list))}")
        print(f"  Std: {np.std(n_sequences_list):.1f}")
    else:
        print(f"\n{dataset_name}: No sequence count data available")

# Combined statistics
all_n_sequences = [p['n_sequences'] for params in all_params.values() 
                   for p in params if 'n_sequences' in p]

if all_n_sequences:
    print(f"\n{'='*70}")
    print("COMBINED (All Datasets):")
    print(f"  Total genes: {len(all_n_sequences)}")
    print(f"  Mean sequences per gene: {np.mean(all_n_sequences):.1f}")
    print(f"  Median sequences per gene: {np.median(all_n_sequences):.1f}")
    print(f"  Range: {int(np.min(all_n_sequences))} - {int(np.max(all_n_sequences))}")
    print(f"  Std: {np.std(all_n_sequences):.1f}")
    print("="*70)
