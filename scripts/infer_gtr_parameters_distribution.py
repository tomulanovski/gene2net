#!/usr/bin/env python3
"""
Infer GTR Parameters Distribution from Real Datasets (IMPROVED VERSION)

This script:
1. Reads IQ-TREE output files (.iqtree) from multiple datasets
2. Extracts substitution model parameters INDEPENDENTLY for each gene
3. Handles different models flexibly (GTR, HKY, TPM3u, etc.)
4. Fits distributions to parameters with maximum available data
5. Saves parameters for sampling in simulations

KEY IMPROVEMENT:
- Extracts parameters independently (not requiring all parameters present)
- More genes for each parameter type → better distributions
- Handles models with/without +F, +G4, +I, +R2, etc.
"""

import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
from collections import Counter

# Configuration
BASE_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/papers"
DATASETS = {
    "Zhao_2021": f"{BASE_DIR}/Zhao_2021/iqtree",
    "Diaz-Perez_2018": f"{BASE_DIR}/Diaz-Perez_2018/iqtree",
    "Sessa_2012b": f"{BASE_DIR}/Sessa_2012b/iqtree",
    "Hori_2014": f"{BASE_DIR}/Hori_2014/iqtree"
}

OUTPUT_DIR = "/groups/itay_mayrose/tomulanovski/gene2net/distribution_inferences"
IQTREE_EXTENSION = ".iqtree"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)


def parse_iqtree_file(iqtree_file):
    """
    Parse an IQ-TREE output file to extract model parameters.
    NOW: Extracts parameters independently - doesn't require all to be present.
    
    Args:
        iqtree_file: Path to .iqtree file
    
    Returns:
        dict: Contains available parameters (model name, rates, base frequencies, alpha)
    """
    try:
        with open(iqtree_file, 'r') as f:
            content = f.read()
        
        result = {
            'model': None,
            'rates': {},
            'base_freqs': {},
            'alpha': None
        }
        
        # Extract model name (e.g., "GTR+F+G4", "HKY+F+G4", "K2P+I")
        model_match = re.search(r'Model of substitution:\s+(\S+)', content)
        if model_match:
            result['model'] = model_match.group(1)
        
        # Extract rate parameters (if present)
        rate_section = re.search(r'Rate parameter R:(.*?)(?:State frequencies:|Rate matrix Q:)', content, re.DOTALL)
        if rate_section:
            rate_text = rate_section.group(1)
            
            # Parse each rate
            rate_patterns = {
                'A-C': r'A-C:\s+([\d.]+)',
                'A-G': r'A-G:\s+([\d.]+)',
                'A-T': r'A-T:\s+([\d.]+)',
                'C-G': r'C-G:\s+([\d.]+)',
                'C-T': r'C-T:\s+([\d.]+)',
                'G-T': r'G-T:\s+([\d.]+)'
            }
            
            for rate_name, pattern in rate_patterns.items():
                match = re.search(pattern, rate_text)
                if match:
                    result['rates'][rate_name] = float(match.group(1))
        
        # Extract base frequencies (if model has +F)
        # Try multiple patterns to catch different formats
        freq_patterns = [
            # Pattern 1: Standard format with "State frequencies:" header
            r'State frequencies:.*?\n\n(.*?)(?:\n\n|Rate matrix Q:)',
            # Pattern 2: Direct empirical counts format
            r'State frequencies: \(empirical counts from alignment\)(.*?)(?:\n\n|Rate matrix Q:)',
            # Pattern 3: Equal frequencies format
            r'State frequencies: \(equal frequencies\)'
        ]
        
        freq_found = False
        for pattern in freq_patterns:
            freq_section = re.search(pattern, content, re.DOTALL)
            if freq_section:
                if 'equal frequencies' in pattern:
                    # Equal frequencies: 0.25 for each
                    result['base_freqs'] = {
                        'pi_A': 0.25,
                        'pi_C': 0.25,
                        'pi_G': 0.25,
                        'pi_T': 0.25
                    }
                    freq_found = True
                    break
                else:
                    freq_text = freq_section.group(1)
                    
                    # Extract individual frequencies
                    pi_patterns = {
                        'pi_A': r'pi\(A\)\s*=\s*([\d.]+)',
                        'pi_C': r'pi\(C\)\s*=\s*([\d.]+)',
                        'pi_G': r'pi\(G\)\s*=\s*([\d.]+)',
                        'pi_T': r'pi\(T\)\s*=\s*([\d.]+)'
                    }
                    
                    temp_freqs = {}
                    for freq_name, pi_pattern in pi_patterns.items():
                        match = re.search(pi_pattern, freq_text)
                        if match:
                            temp_freqs[freq_name] = float(match.group(1))
                    
                    # Only accept if we got all 4 frequencies
                    if len(temp_freqs) == 4:
                        result['base_freqs'] = temp_freqs
                        freq_found = True
                        break
        
        # Extract Gamma alpha parameter (if model has +G or +G4)
        alpha_match = re.search(r'Gamma shape alpha:\s+([\d.]+)', content)
        if alpha_match:
            result['alpha'] = float(alpha_match.group(1))
        
        # Return even if not all parameters are present
        return result
    
    except Exception as e:
        print(f"Warning: Failed to parse {iqtree_file}: {e}")
        return None


def collect_parameters(datasets):
    """
    Collect parameters from all .iqtree files.
    NOW: Creates separate records for each parameter type.
    
    Args:
        datasets: Dictionary mapping dataset names to tree directories
    
    Returns:
        tuple: (all_data DataFrame, counts dict)
    """
    all_data = []
    counts = {
        'total_files': 0,
        'with_model': 0,
        'with_base_freqs': 0,
        'with_alpha': 0,
        'with_rates': 0
    }
    
    for dataset_name, tree_dir in datasets.items():
        print(f"\nProcessing {dataset_name}...")
        
        if not os.path.exists(tree_dir):
            print(f"  Warning: Directory not found: {tree_dir}")
            continue
        
        # Find all .iqtree files
        iqtree_files = list(Path(tree_dir).glob(f"*{IQTREE_EXTENSION}"))
        print(f"  Found {len(iqtree_files)} .iqtree files")
        counts['total_files'] += len(iqtree_files)
        
        for iqtree_file in iqtree_files:
            params = parse_iqtree_file(str(iqtree_file))
            
            if params:
                row = {
                    'dataset': dataset_name,
                    'gene': iqtree_file.stem.replace(IQTREE_EXTENSION, ''),
                    'model': params['model']
                }
                
                # Track what we have
                if params['model']:
                    counts['with_model'] += 1
                
                # Add alpha if present
                if params['alpha'] is not None:
                    row['alpha'] = params['alpha']
                    counts['with_alpha'] += 1
                
                # Add base frequencies if present
                if params['base_freqs']:
                    for freq_name, freq_value in params['base_freqs'].items():
                        row[freq_name] = freq_value
                    counts['with_base_freqs'] += 1
                
                # Add rates if present
                if params['rates']:
                    for rate_name, rate_value in params['rates'].items():
                        row[rate_name] = rate_value
                    counts['with_rates'] += 1
                
                all_data.append(row)
    
    return pd.DataFrame(all_data), counts


def is_gtr_model(model_name):
    """
    Check if a model is GTR or GTR-like.
    """
    if not model_name:
        return False
    return 'GTR' in model_name.upper()


def fit_distribution(data, dist_name):
    """
    Fit a distribution to the data using Maximum Likelihood Estimation.
    
    Returns:
        dict: Contains parameters, AIC, KS test results
        None: If fitting fails
    """
    try:
        # Check for constant data or insufficient variation
        if len(np.unique(data)) < 3:
            return None
        
        if np.std(data) < 1e-10:
            return None
        
        dist = getattr(stats, dist_name)
        
        # Fit distribution parameters using MLE
        params = dist.fit(data)
        
        # Check if parameters are valid
        if any(np.isnan(p) or np.isinf(p) for p in params):
            return None
        
        # Calculate log-likelihood
        log_likelihood = np.sum(dist.logpdf(data, *params))
        
        # Check if log-likelihood is valid
        if np.isnan(log_likelihood) or np.isinf(log_likelihood):
            return None
        
        # Calculate AIC: 2k - 2*log(L)
        k = len(params)
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
    
    except Exception as e:
        return None


def fit_parameter_distributions(df):
    """
    Fit distributions to each parameter using ALL available data.
    
    Returns:
        dict: Best fits for each parameter
    """
    results = {}
    
    # Distributions to try
    distributions = ['gamma', 'lognorm', 'norm']
    
    # Parameters to fit
    params_to_fit = {
        'alpha': 'Gamma Alpha',
        'pi_A': 'Base Frequency A',
        'pi_C': 'Base Frequency C',
        'pi_G': 'Base Frequency G',
        'pi_T': 'Base Frequency T'
    }
    
    print("\nFitting distributions to parameters...")
    
    for param_name, param_label in params_to_fit.items():
        if param_name not in df.columns:
            print(f"  {param_label}: Not found in data")
            continue
        
        data = df[param_name].dropna().values
        if len(data) < 10:
            print(f"  {param_label}: Too few observations ({len(data)})")
            continue
        
        print(f"\n  {param_label} ({len(data)} observations):")
        
        fits = []
        for dist_name in distributions:
            fit_result = fit_distribution(data, dist_name)
            if fit_result is not None:
                fits.append(fit_result)
                print(f"    {dist_name}: AIC={fit_result['aic']:.2f}, KS p-value={fit_result['ks_pvalue']:.4f}")
            else:
                print(f"    {dist_name}: Failed - invalid fit")
        
        if fits:
            best_fit = min(fits, key=lambda x: x['aic'])
            results[param_name] = {
                'best_fit': best_fit,
                'all_fits': fits,
                'n_obs': len(data)
            }
    
    return results


def analyze_gtr_rates(df):
    """
    Analyze GTR rate parameters using ALL available data.
    
    Returns:
        dict: GTR rate analysis and fitted distributions
    """
    # Check how many GTR genes we have
    gtr_df = df[df['model'].apply(is_gtr_model)].copy()
    n_gtr = len(gtr_df)
    
    print(f"\n{'='*70}")
    print(f"GTR RATE PARAMETER ANALYSIS")
    print(f"{'='*70}")
    print(f"Total genes with GTR model: {n_gtr}")
    
    rate_params = ['A-C', 'A-G', 'A-T', 'C-G', 'C-T', 'G-T']
    
    results = {
        'n_gtr_genes': n_gtr,
        'use_gtr_only': n_gtr >= 50,
        'rates': {}
    }
    
    if n_gtr >= 50:
        print(f"Using only GTR genes for rate parameter inference (n={n_gtr})")
        rate_df = gtr_df
    else:
        print(f"Too few GTR genes ({n_gtr}). Using all genes with rate parameters.")
        # Use all genes that have ANY rate parameters
        rate_df = df.copy()
    
    # Fit distributions to each rate parameter independently
    distributions = ['gamma', 'lognorm', 'norm']
    
    for rate_name in rate_params:
        if rate_name not in rate_df.columns:
            continue
        
        data = rate_df[rate_name].dropna().values
        if len(data) < 10:
            print(f"\n  {rate_name}: Too few observations ({len(data)})")
            continue
        
        print(f"\n  {rate_name} ({len(data)} observations):")
        
        fits = []
        for dist_name in distributions:
            fit_result = fit_distribution(data, dist_name)
            if fit_result is not None:
                fits.append(fit_result)
                print(f"    {dist_name}: AIC={fit_result['aic']:.2f}, KS p-value={fit_result['ks_pvalue']:.4f}")
            else:
                print(f"    {dist_name}: Failed - invalid fit")
        
        if fits:
            best_fit = min(fits, key=lambda x: x['aic'])
            results['rates'][rate_name] = {
                'best_fit': best_fit,
                'all_fits': fits,
                'n_obs': len(data)
            }
    
    return results


def plot_results(df, param_fits, rate_analysis, output_file):
    """
    Create visualization of parameter distributions.
    """
    # Count models
    model_counts = Counter(df['model'].dropna())
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle('GTR Parameter Distribution Analysis (IMPROVED)', fontsize=16, fontweight='bold')
    
    # Plot 1: Model distribution (top-left, larger)
    ax1 = fig.add_subplot(gs[0, :2])
    models = list(model_counts.keys())[:10]  # Top 10 models
    counts = [model_counts[m] for m in models]
    
    bars = ax1.barh(range(len(models)), counts, color='skyblue', edgecolor='black')
    ax1.set_yticks(range(len(models)))
    ax1.set_yticklabels(models)
    ax1.set_xlabel('Number of Genes', fontsize=11)
    ax1.set_title('Top 10 Substitution Models Used', fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Highlight GTR
    for i, model in enumerate(models):
        if 'GTR' in model.upper():
            bars[i].set_color('lightcoral')
    
    # Add text
    gtr_count = sum(1 for m in df['model'] if is_gtr_model(m))
    ax1.text(0.98, 0.02, f'GTR genes: {gtr_count}/{len(df)}', 
             transform=ax1.transAxes, fontsize=10, va='bottom', ha='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot 2: Alpha distribution
    ax2 = fig.add_subplot(gs[0, 2])
    if 'alpha' in param_fits:
        alpha_data = df['alpha'].dropna()
        ax2.hist(alpha_data, bins=30, density=True, alpha=0.6, color='skyblue', edgecolor='black')
        
        best_fit = param_fits['alpha']['best_fit']
        x_range = np.linspace(alpha_data.min(), alpha_data.max(), 100)
        dist = getattr(stats, best_fit['distribution'])
        y = dist.pdf(x_range, *best_fit['parameters'])
        ax2.plot(x_range, y, 'r-', linewidth=2, label=best_fit['distribution'])
        ax2.set_xlabel('Alpha', fontsize=10)
        ax2.set_ylabel('Density', fontsize=10)
        ax2.set_title(f'Gamma Alpha (n={len(alpha_data)})', fontsize=11, fontweight='bold')
        ax2.legend()
        ax2.grid(alpha=0.3)
    
    # Plots 3-6: Base frequencies
    base_freq_params = ['pi_A', 'pi_C', 'pi_G', 'pi_T']
    base_freq_labels = ['π(A)', 'π(C)', 'π(G)', 'π(T)']
    
    for idx, (param, label) in enumerate(zip(base_freq_params, base_freq_labels)):
        ax = fig.add_subplot(gs[1, idx % 3] if idx < 3 else gs[2, 0])
        
        if param in param_fits:
            data = df[param].dropna()
            ax.hist(data, bins=30, density=True, alpha=0.6, color='lightgreen', edgecolor='black')
            
            best_fit = param_fits[param]['best_fit']
            x_range = np.linspace(data.min(), data.max(), 100)
            dist = getattr(stats, best_fit['distribution'])
            y = dist.pdf(x_range, *best_fit['parameters'])
            ax.plot(x_range, y, 'r-', linewidth=2, label=best_fit['distribution'])
            ax.set_xlabel(label, fontsize=10)
            ax.set_ylabel('Density', fontsize=10)
            ax.set_title(f'{label} (n={len(data)})', fontsize=11, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
    
    # Plot 7: Rate parameters summary (if available)
    ax7 = fig.add_subplot(gs[2, 1:])
    if rate_analysis['rates']:
        rate_names = list(rate_analysis['rates'].keys())
        rate_means = [df[rate].dropna().mean() for rate in rate_names]
        rate_stds = [df[rate].dropna().std() for rate in rate_names]
        
        x_pos = np.arange(len(rate_names))
        ax7.bar(x_pos, rate_means, yerr=rate_stds, alpha=0.6, color='coral', 
                edgecolor='black', capsize=5)
        ax7.set_xticks(x_pos)
        ax7.set_xticklabels(rate_names)
        ax7.set_ylabel('Rate Value', fontsize=11)
        ax7.set_title('GTR Rate Parameters (Mean ± Std)', fontsize=12, fontweight='bold')
        ax7.grid(axis='y', alpha=0.3)
        
        # Add sample size info
        n_obs = rate_analysis['rates'][rate_names[0]]['n_obs']
        ax7.text(0.98, 0.98, f'n={n_obs}', transform=ax7.transAxes, 
                fontsize=10, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        ax7.text(0.5, 0.5, 'Insufficient GTR rate data', 
                transform=ax7.transAxes, ha='center', va='center', fontsize=12)
        ax7.axis('off')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()  # Close figure to prevent hanging
    print(f"\nPlot saved to: {output_file}")


def main():
    print("="*70)
    print("GTR PARAMETER DISTRIBUTION INFERENCE (IMPROVED)")
    print("="*70)
    
    # Collect parameters from all datasets
    print("\n1. Collecting parameters from .iqtree files...")
    df, counts = collect_parameters(DATASETS)
    
    if df.empty:
        print("Error: No .iqtree files found or parsed successfully!")
        return
    
    print(f"\n{'='*70}")
    print("DATA COLLECTION SUMMARY")
    print(f"{'='*70}")
    print(f"Total .iqtree files processed: {counts['total_files']}")
    print(f"Genes with model information: {counts['with_model']}")
    print(f"Genes with base frequencies: {counts['with_base_freqs']}")
    print(f"Genes with alpha parameter: {counts['with_alpha']}")
    print(f"Genes with rate parameters: {counts['with_rates']}")
    print(f"Total genes in dataset: {len(df)}")
    print(f"Datasets: {df['dataset'].unique()}")
    
    # Analyze models used
    print("\n2. Analyzing substitution models used...")
    model_counts = Counter(df['model'].dropna())
    print("\nTop 10 models:")
    for model, count in model_counts.most_common(10):
        print(f"  {model}: {count}")
    
    # Fit distributions to base frequencies and alpha
    print("\n3. Fitting distributions to base frequencies and alpha...")
    param_fits = fit_parameter_distributions(df)
    
    # Analyze GTR rates
    print("\n4. Analyzing GTR rate parameters...")
    rate_analysis = analyze_gtr_rates(df)
    
    # Prepare output data
    output_data = {
        'summary': {
            'n_total_genes': len(df),
            'n_total_files': counts['total_files'],
            'n_with_base_freqs': counts['with_base_freqs'],
            'n_with_alpha': counts['with_alpha'],
            'n_with_rates': counts['with_rates'],
            'n_gtr_genes': rate_analysis['n_gtr_genes'],
            'using_gtr_only_for_rates': rate_analysis['use_gtr_only'],
            'models_used': dict(model_counts)
        },
        'base_frequencies': {},
        'alpha': {},
        'gtr_rates': {}
    }
    
    # Add base frequency fits
    for param in ['pi_A', 'pi_C', 'pi_G', 'pi_T']:
        if param in param_fits:
            best_fit = param_fits[param]['best_fit']
            output_data['base_frequencies'][param] = {
                'distribution': best_fit['distribution'],
                'parameters': list(best_fit['parameters']),
                'aic': best_fit['aic'],
                'ks_pvalue': best_fit['ks_pvalue'],
                'n_obs': param_fits[param]['n_obs']
            }
    
    # Add alpha fit
    if 'alpha' in param_fits:
        best_fit = param_fits['alpha']['best_fit']
        output_data['alpha'] = {
            'distribution': best_fit['distribution'],
            'parameters': list(best_fit['parameters']),
            'aic': best_fit['aic'],
            'ks_pvalue': best_fit['ks_pvalue'],
            'n_obs': param_fits['alpha']['n_obs']
        }
    
    # Add rate fits
    if rate_analysis['rates']:
        for rate_name, rate_data in rate_analysis['rates'].items():
            best_fit = rate_data['best_fit']
            output_data['gtr_rates'][rate_name] = {
                'distribution': best_fit['distribution'],
                'parameters': list(best_fit['parameters']),
                'aic': best_fit['aic'],
                'ks_pvalue': best_fit['ks_pvalue'],
                'n_obs': rate_data['n_obs']
            }
    
    # Save outputs
    print("\n5. Saving outputs...")
    
    # JSON parameters
    json_file = os.path.join(OUTPUT_DIR, 'gtr_parameters_distribution_params.json')
    with open(json_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"   Parameters saved to: {json_file}")
    
    # CSV summary - per dataset statistics
    summary_rows = []
    for dataset in df['dataset'].unique():
        dataset_df = df[df['dataset'] == dataset]
        row = {
            'Dataset': dataset,
            'N_genes': len(dataset_df),
            'N_with_alpha': dataset_df['alpha'].notna().sum(),
            'N_with_base_freqs': dataset_df['pi_A'].notna().sum(),
            'Mean_alpha': dataset_df['alpha'].mean(),
            'Std_alpha': dataset_df['alpha'].std(),
            'Mean_pi_A': dataset_df['pi_A'].mean(),
            'Mean_pi_C': dataset_df['pi_C'].mean(),
            'Mean_pi_G': dataset_df['pi_G'].mean(),
            'Mean_pi_T': dataset_df['pi_T'].mean()
        }
        summary_rows.append(row)
    
    # Add overall row
    summary_rows.append({
        'Dataset': 'OVERALL',
        'N_genes': len(df),
        'N_with_alpha': df['alpha'].notna().sum(),
        'N_with_base_freqs': df['pi_A'].notna().sum(),
        'Mean_alpha': df['alpha'].mean(),
        'Std_alpha': df['alpha'].std(),
        'Mean_pi_A': df['pi_A'].mean(),
        'Mean_pi_C': df['pi_C'].mean(),
        'Mean_pi_G': df['pi_G'].mean(),
        'Mean_pi_T': df['pi_T'].mean()
    })
    
    summary_df = pd.DataFrame(summary_rows)
    csv_file = os.path.join(OUTPUT_DIR, 'gtr_parameters_distribution_summary.csv')
    summary_df.to_csv(csv_file, index=False)
    print(f"   Summary saved to: {csv_file}")
    
    # Visualization
    plot_file = os.path.join(OUTPUT_DIR, 'gtr_parameters_distribution_plot.png')
    plot_results(df, param_fits, rate_analysis, plot_file)
    
    # Usage example
    print("\n" + "="*70)
    print("HOW TO SAMPLE GTR PARAMETERS IN YOUR SIMULATIONS:")
    print("="*70)
    print("""
import json
from scipy import stats

# Load parameters
params = json.load(open('/groups/itay_mayrose/tomulanovski/gene2net/distribution_inferences/gtr_parameters_distribution_params.json'))

# Sample alpha
alpha_dist = getattr(stats, params['alpha']['distribution'])
alpha = alpha_dist.rvs(*params['alpha']['parameters'])

# Sample base frequencies
base_freqs = {}
for base in ['pi_A', 'pi_C', 'pi_G', 'pi_T']:
    dist = getattr(stats, params['base_frequencies'][base]['distribution'])
    base_freqs[base] = dist.rvs(*params['base_frequencies'][base]['parameters'])

# Normalize base frequencies to sum to 1
total = sum(base_freqs.values())
base_freqs = {k: v/total for k, v in base_freqs.items()}

# Sample GTR rates
gtr_rates = {}
for rate in ['A-C', 'A-G', 'A-T', 'C-G', 'C-T', 'G-T']:
    if rate in params['gtr_rates']:
        dist = getattr(stats, params['gtr_rates'][rate]['distribution'])
        gtr_rates[rate] = dist.rvs(*params['gtr_rates'][rate]['parameters'])
    else:
        # If G-T failed, use 1.0 as reference rate
        gtr_rates[rate] = 1.0

print(f"Sampled alpha: {alpha:.4f}")
print(f"Sampled base frequencies: {base_freqs}")
print(f"Sampled GTR rates: {gtr_rates}")
""")
    
    print("="*70)
    print("DONE!")
    print("="*70)


if __name__ == "__main__":
    main()
