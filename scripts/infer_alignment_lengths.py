"""
================================================================================
ALIGNMENT LENGTH DISTRIBUTION INFERENCE
================================================================================

PURPOSE:
    This script infers the distribution of alignment lengths from multiple 
    sequence alignment (MSA) datasets. The inferred distribution can later be 
    used to generate realistic alignment lengths for simulations.

WHAT IT DOES:
    1. Reads all FASTA alignment files from multiple dataset directories
    2. Extracts the alignment length from each file (all sequences in one 
       MSA have the same length due to gaps)
    3. Collects lengths across all datasets (treating each gene independently)
    4. Computes statistics per dataset and overall
    5. Fits three candidate distributions (Gamma, Log-normal, Normal) to the 
       pooled data
    6. Selects the best-fitting distribution using AIC (Akaike Information 
       Criterion - lower is better)
    7. Generates visualizations and saves parameters for future use

INPUT:
    - Multiple directories containing FASTA alignment files (.fasta or .fa)
    - Each file = one gene alignment
    - All sequences in one file = same length (MSA property)

OUTPUT:
    1. Console: Detailed statistics per dataset and summary table
    2. PNG file: 4-panel visualization
       - Histogram with fitted distribution curve
       - Q-Q plot (quality check)
       - Box plot comparing datasets
       - Overlaid histograms per dataset
    3. JSON file: Complete data including:
       - Best distribution name and parameters
       - ALL distributions tested with AIC values
       - Per-dataset statistics (mean, std, median, min, max)
       - All individual lengths
    4. CSV file: Simple summary table for Excel

HOW TO USE THE RESULTS (SAMPLING):
    After running this script, you can sample new realistic alignment lengths:
    
    import json
    import numpy as np
    from scipy import stats
    
    # Load saved parameters
    with open('alignment_length_params.json', 'r') as f:
        params = json.load(f)
    
    # Get the best distribution
    dist_name = params['best_distribution']['name']
    dist_params = params['best_distribution']['parameters']
    dist = getattr(stats, dist_name)
    
    # Sample ONE new length
    new_length = int(dist.rvs(*dist_params))
    
    # Sample MANY new lengths
    new_lengths = [int(dist.rvs(*dist_params)) for _ in range(100)]

METHODOLOGY:
    - Distribution fitting uses Maximum Likelihood Estimation (MLE)
    - MLE finds parameters that maximize the probability of observing your data
    - AIC balances fit quality with model complexity (penalizes extra parameters)
    - All genes are treated as independent observations (pooled across datasets)

AUTHOR: Tom Ulanovski
DATE: 2025
================================================================================
"""

from Bio import AlignIO
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import glob
import os
import json

# ===================================================================
# CONFIGURATION
# ===================================================================

output_dir = "/groups/itay_mayrose/tomulanovski/gene2net/distribution_inferences/"
os.makedirs(output_dir, exist_ok=True)

datasets = [
    "/groups/itay_mayrose/tomulanovski/gene2net/papers/Zhao_2021/genes/",
    "/groups/itay_mayrose/tomulanovski/gene2net/papers/Diaz-Perez_2018/genes/",
    "/groups/itay_mayrose/tomulanovski/gene2net/papers/Sessa_2012b/genes/",
    "/groups/itay_mayrose/tomulanovski/gene2net/papers/Hori_2014/genes/",
]

# ===================================================================
# STEP 1: Collect alignment lengths
# ===================================================================

all_data = {}  # Store per-dataset information
alignment_lengths = []

print("Collecting alignment lengths from all FASTA files...\n")

for dataset_path in datasets:
    # Extract dataset name from the paper directory (2 levels up)
    # e.g., /path/to/papers/Zhao_2021/genes/ -> Zhao_2021
    path_parts = dataset_path.rstrip('/').split('/')
    dataset_name = path_parts[-2]  # Get the paper name, not "genes"
    
    # Find ALL fasta files (any pattern)
    fasta_files = glob.glob(os.path.join(dataset_path, "*.fasta")) + \
                  glob.glob(os.path.join(dataset_path, "*.fa"))
    
    print(f"{'='*70}")
    print(f"Dataset: {dataset_name}")
    print(f"Path: {dataset_path}")
    print(f"Found {len(fasta_files)} alignment files\n")
    
    dataset_lengths = []
    
    # Read each alignment
    for fasta_file in sorted(fasta_files):
        try:
            aln = AlignIO.read(fasta_file, "fasta")
            length = aln.get_alignment_length()
            n_sequences = len(aln)
            alignment_lengths.append(length)
            dataset_lengths.append(length)
            print(f"  {os.path.basename(fasta_file):45s} {length:5d} bp  ({n_sequences} seqs)")
        except Exception as e:
            print(f"  ERROR - {os.path.basename(fasta_file)}: {e}")
    
    # Store dataset stats
    if dataset_lengths:
        all_data[dataset_name] = {
            'lengths': dataset_lengths,
            'count': len(dataset_lengths),
            'mean': np.mean(dataset_lengths),
            'std': np.std(dataset_lengths) if len(dataset_lengths) > 1 else 0.0,
            'median': np.median(dataset_lengths),
            'min': min(dataset_lengths),
            'max': max(dataset_lengths)
        }
        print(f"\n  Statistics for {dataset_name}:")
        print(f"    Number of alignments: {len(dataset_lengths)}")
        print(f"    Average length:       {np.mean(dataset_lengths):.2f} bp")
        print(f"    Std deviation:        {all_data[dataset_name]['std']:.2f} bp")
        print(f"    Median length:        {np.median(dataset_lengths):.2f} bp")
        print(f"    Range:                [{min(dataset_lengths)} - {max(dataset_lengths)}] bp")
    else:
        print(f"  WARNING: No files processed!")
    
    print()

alignment_lengths = np.array(alignment_lengths)

# ===================================================================
# SUMMARY TABLE - ALL DATASETS
# ===================================================================

print("\n" + "="*85)
print("SUMMARY TABLE - ALL DATASETS:")
print("="*85)
print(f"{'Dataset':<25} {'N':>6} {'Mean':>10} {'Std':>10} {'Median':>10} {'Min':>8} {'Max':>8}")
print("-"*85)

for dataset_name in sorted(all_data.keys()):
    data = all_data[dataset_name]
    print(f"{dataset_name:<25} {data['count']:>6} {data['mean']:>10.1f} "
          f"{data['std']:>10.1f} {data['median']:>10.1f} {data['min']:>8} {data['max']:>8}")

print("-"*85)
print(f"{'OVERALL (pooled)':<25} {len(alignment_lengths):>6} {alignment_lengths.mean():>10.1f} "
      f"{alignment_lengths.std():>10.1f} {np.median(alignment_lengths):>10.1f} "
      f"{alignment_lengths.min():>8} {alignment_lengths.max():>8}")
print("="*85)

if len(alignment_lengths) == 0:
    print("\nERROR: No alignments collected! Check your paths.")
    exit(1)

# ===================================================================
# STEP 2: Fit distributions
# ===================================================================

print("\n\nFitting distributions to pooled data...\n")

distributions = {
    'gamma': stats.gamma,
    'lognorm': stats.lognorm,
    'norm': stats.norm,
}

best_dist_name = None
best_dist_obj = None
best_params = None
best_aic = np.inf

# Store all distribution results
all_distribution_results = {}

print(f"{'Distribution':<12} {'AIC':<15} {'KS p-value':<12}")
print("-"*40)

for name, dist in distributions.items():
    params = dist.fit(alignment_lengths)
    log_lik = np.sum(dist.logpdf(alignment_lengths, *params))
    k = len(params)
    aic = 2*k - 2*log_lik
    ks_stat, ks_pval = stats.kstest(alignment_lengths, 
                                     lambda x: dist.cdf(x, *params))
    
    print(f"{name:<12} {aic:<15.2f} {ks_pval:<12.6f}")
    
    # Store results for this distribution
    all_distribution_results[name] = {
        'parameters': [float(p) for p in params],
        'aic': float(aic),
        'log_likelihood': float(log_lik),
        'ks_statistic': float(ks_stat),
        'ks_pvalue': float(ks_pval),
        'n_parameters': k
    }
    
    if aic < best_aic:
        best_aic = aic
        best_dist_name = name
        best_dist_obj = dist
        best_params = params

print("-"*40)
print(f"\nBest distribution: {best_dist_name.upper()}")
print(f"Parameters: {best_params}")

# ===================================================================
# STEP 3: Visualize
# ===================================================================

fig = plt.figure(figsize=(16, 10))

# Main histogram with fitted distribution
ax1 = plt.subplot(2, 2, 1)
ax1.hist(alignment_lengths, bins=50, density=True, alpha=0.7, 
         edgecolor='black', color='skyblue', label='Observed')

x = np.linspace(alignment_lengths.min(), alignment_lengths.max(), 1000)
pdf = best_dist_obj.pdf(x, *best_params)
ax1.plot(x, pdf, 'r-', linewidth=2.5, label=f'Fitted {best_dist_name}')

ax1.set_xlabel('Alignment Length (bp)', fontsize=11)
ax1.set_ylabel('Density', fontsize=11)
ax1.set_title(f'Pooled Alignment Length Distribution (n={len(alignment_lengths)})', 
              fontsize=12, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Q-Q plot
ax2 = plt.subplot(2, 2, 2)
stats.probplot(alignment_lengths, dist=best_dist_obj, sparams=best_params, 
               plot=ax2)
ax2.set_title(f'Q-Q Plot ({best_dist_name})', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)

# Box plot by dataset
ax3 = plt.subplot(2, 2, 3)
dataset_names = sorted(all_data.keys())
dataset_lengths_list = [all_data[name]['lengths'] for name in dataset_names]

bp = ax3.boxplot(dataset_lengths_list, tick_labels=dataset_names, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightblue')

ax3.set_ylabel('Alignment Length (bp)', fontsize=11)
ax3.set_title('Alignment Lengths by Dataset', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha='right')

# Histogram by dataset (overlaid)
ax4 = plt.subplot(2, 2, 4)
colors = plt.cm.Set3(np.linspace(0, 1, len(dataset_names)))
for i, (dataset_name, color) in enumerate(zip(dataset_names, colors)):
    lengths = all_data[dataset_name]['lengths']
    ax4.hist(lengths, bins=30, alpha=0.5, label=f"{dataset_name} (n={len(lengths)})", 
             color=color, edgecolor='black')

ax4.set_xlabel('Alignment Length (bp)', fontsize=11)
ax4.set_ylabel('Frequency', fontsize=11)
ax4.set_title('Alignment Length Distribution by Dataset', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
output_file = os.path.join(output_dir, 'alignment_length_distribution.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nSaved figure: {output_file}")
plt.close()

# ===================================================================
# STEP 4: Save parameters and dataset statistics
# ===================================================================

output = {
    'best_distribution': {
        'name': best_dist_name,
        'parameters': [float(p) for p in best_params],
        'aic': all_distribution_results[best_dist_name]['aic'],
        'ks_pvalue': all_distribution_results[best_dist_name]['ks_pvalue']
    },
    
    # All distribution comparison results
    'all_distributions_tested': all_distribution_results,
    
    'n_alignments_total': int(len(alignment_lengths)),
    
    'per_dataset_statistics': {
        name: {
            'n_alignments': data['count'],
            'mean_length': float(data['mean']),
            'std_length': float(data['std']),
            'median_length': float(data['median']),
            'min_length': int(data['min']),
            'max_length': int(data['max']),
            'all_lengths': [int(x) for x in data['lengths']]
        }
        for name, data in all_data.items()
    },
    
    'overall_statistics': {
        'mean': float(np.mean(alignment_lengths)),
        'std': float(np.std(alignment_lengths)),
        'median': float(np.median(alignment_lengths)),
        'min': int(alignment_lengths.min()),
        'max': int(alignment_lengths.max()),
        'q25': float(np.percentile(alignment_lengths, 25)),
        'q75': float(np.percentile(alignment_lengths, 75))
    }
}

output_json = os.path.join(output_dir, 'alignment_length_params.json')
with open(output_json, 'w') as f:
    json.dump(output, f, indent=2)

print(f"Saved parameters: {output_json}")

# Also save a simple CSV summary
csv_file = os.path.join(output_dir, 'alignment_length_summary.csv')
with open(csv_file, 'w') as f:
    # Header row
    f.write("Dataset,N_alignments,Mean_length,Std_length,Median_length,Min_length,Max_length\n")
    
    # Data rows - one per dataset (sorted alphabetically)
    for name in sorted(all_data.keys()):
        data = all_data[name]
        f.write(f"{name},{data['count']},{data['mean']:.2f},{data['std']:.2f},"
                f"{data['median']:.2f},{data['min']},{data['max']}\n")
    
    # Overall row
    f.write(f"OVERALL,{len(alignment_lengths)},{alignment_lengths.mean():.2f},"
            f"{alignment_lengths.std():.2f},{np.median(alignment_lengths):.2f},"
            f"{alignment_lengths.min()},{alignment_lengths.max()}\n")

print(f"Saved CSV summary: {csv_file}")

# ===================================================================
# STEP 5: Test sampling
# ===================================================================

def sample_alignment_length():
    """Sample a new alignment length from the fitted distribution"""
    return int(best_dist_obj.rvs(*best_params))

print(f"\n{'='*70}")
print("Sample alignment lengths (10 random draws):")
print(f"{'='*70}")
samples = [sample_alignment_length() for _ in range(10)]
for i, s in enumerate(samples, 1):
    print(f"  Sample {i:2d}: {s:5d} bp")

print(f"\n{'='*70}")
print("DONE! All outputs saved to:")
print(f"  {output_dir}")
print(f"{'='*70}")
