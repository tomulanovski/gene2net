#!/usr/bin/env python3
"""
Filter alignments by quality criteria, then extract GTR parameters only from filtered genes.
No need to re-run IQ-TREE - we just select which existing IQ-TREE files to use.
"""

import argparse
from pathlib import Path
from Bio import SeqIO
import re
import pickle
import pandas as pd
import numpy as np

# Dataset configurations
DATASETS = {
    'Zhao_2021': {
        'alignment_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Zhao_et_al_2021',
        'iqtree_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Zhao_et_al_2021/gene_trees',
        'extensions': ['.fasta', '.fa', '.fna', '.fas']
    },
    'Ren_2024': {
        'alignment_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Ren_et_al_2024',
        'iqtree_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Ren_et_al_2024/gene_trees',
        'extensions': ['.fas', '.fasta', '.fa', '.fna']
    },
    'Morales_Briones_2021': {
        'alignment_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Morales-Briones_et_al_2021',
        'iqtree_dir': '/groups/itay_mayrose/ronenshtein/gene2network/papers/Morales-Briones_et_al_2021/gene_trees',
        'extensions': ['.aln-cln']
    }
}

OUTPUT_DIR = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions')


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter alignments by quality and extract GTR parameters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--min-length', type=int, default=None,
                       help='Minimum alignment length (bp). Default: no filtering')
    parser.add_argument('--min-sequences', type=int, default=None,
                       help='Minimum number of sequences. Default: no filtering')
    parser.add_argument('--max-gaps', type=float, default=None,
                       help='Maximum gap percentage. Default: no filtering')
    parser.add_argument('--min-informative', type=float, default=None,
                       help='Minimum informative sites percentage. Default: no filtering')
    parser.add_argument('--max-constant', type=float, default=None,
                       help='Maximum constant sites percentage. Default: no filtering')
    
    # Convenience flag to enable all filters with default values
    parser.add_argument('--enable-filters', action='store_true',
                       help='Enable all quality filters with default values (length≥350, sequences≥15, gaps≤30%%, informative≥20%%, constant≤60%%)')
    
    return parser.parse_args()


def extract_gene_id(alignment_file, dataset_name):
    """Extract gene ID from alignment filename - just use the stem (filename without extension)."""
    return alignment_file.stem


def find_iqtree_file(gene_id, iqtree_dir, dataset_name):
    """
    Find IQ-TREE file for a gene ID using known patterns.
    Much faster than recursive search.
    """
    iqtree_dir = Path(iqtree_dir)
    
    # Use dataset-specific patterns based on examples
    if dataset_name == 'Zhao_2021':
        # Three possible patterns:
        # 1. aligned_10005_old.fasta → iqt_10005_old/run_10005_old.iqt.iqtree
        # 2. aligned_10210.fasta → iqt_10210/run_10210.iqt.iqtree
        # 3. aligned_R044790.fasta → iqt_R044790/run_R044790.iqt.iqtree
        # Extract the ID part (remove "aligned_" prefix if present)
        
        # Remove "aligned_" prefix if it exists
        if gene_id.startswith('aligned_'):
            id_part = gene_id.replace('aligned_', '', 1)
        else:
            id_part = gene_id
        
        # Check if it contains "_old"
        if '_old' in id_part:
            # Use pattern with "_old"
            pattern = iqtree_dir / f'iqt_{id_part}' / f'run_{id_part}.iqt.iqtree'
        else:
            # Use pattern without "_old"
            pattern = iqtree_dir / f'iqt_{id_part}' / f'run_{id_part}.iqt.iqtree'
        
        if pattern.exists():
            return pattern
    
    elif dataset_name == 'Ren_2024':
        # Pattern: gene_trees/iqt_2000088/run_2000088.iqt.iqtree
        # Gene ID from 2000088.fas is "2000088"
        pattern = iqtree_dir / f'iqt_{gene_id}' / f'run_{gene_id}.iqt.iqtree'
        if pattern.exists():
            return pattern
    
    elif dataset_name == 'Morales_Briones_2021':
        # Pattern: gene_trees/iqt_gene00141_158_14.inclade1.ortho1.rt.mo.out.NT.fs/run_gene00141_158_14.inclade1.ortho1.rt.mo.out.NT.fs.iqt.iqtree
        # Gene ID is the full filename stem
        pattern = iqtree_dir / f'iqt_{gene_id}' / f'run_{gene_id}.iqt.iqtree'
        if pattern.exists():
            return pattern
    
    return None


def extract_alignment_info_from_iqtree(iqtree_file):
    """
    Extract alignment metadata from IQ-TREE file header.
    This gives us constant sites and informative sites percentages.
    """
    try:
        with open(iqtree_file, 'r') as f:
            content = f.read(2000)  # Just read first 2000 chars (header section)
        
        info = {}
        
        # Extract number of sequences and sites
        pattern = r'Input data:\s+(\d+)\s+sequences?\s+with\s+(\d+)\s+nucleotide sites'
        match = re.search(pattern, content)
        if match:
            info['n_sequences'] = int(match.group(1))
            info['alignment_length'] = int(match.group(2))
        
        # Extract number of constant sites
        pattern = r'Number of constant sites:\s+(\d+)'
        match = re.search(pattern, content)
        if match:
            info['constant_sites'] = int(match.group(1))
        
        # Extract number of parsimony informative sites
        pattern = r'Number of parsimony informative sites:\s+(\d+)'
        match = re.search(pattern, content)
        if match:
            info['informative_sites'] = int(match.group(1))
        
        # Calculate percentages if we have the data
        if 'alignment_length' in info:
            if 'constant_sites' in info:
                info['constant_percentage'] = (info['constant_sites'] / info['alignment_length']) * 100
            if 'informative_sites' in info:
                info['informative_percentage'] = (info['informative_sites'] / info['alignment_length']) * 100
        
        return info
        
    except Exception as e:
        return {}


def check_alignment_quality(alignment_file, iqtree_file, min_length, min_sequences, max_gaps, 
                           min_informative, max_constant):
    """
    Check if alignment meets quality criteria.
    Returns (passes, stats_dict).
    
    If all filters are None, always returns True (no filtering).
    """
    try:
        sequences = list(SeqIO.parse(alignment_file, 'fasta'))
        
        if not sequences:
            return False, None
        
        n_sequences = len(sequences)
        alignment_length = len(sequences[0].seq)
        
        # Calculate gap percentage from alignment
        total_positions = n_sequences * alignment_length
        gap_count = sum(str(seq.seq).count('-') + str(seq.seq).count('N') 
                       for seq in sequences)
        gap_percentage = (gap_count / total_positions) * 100
        
        # Get constant and informative sites from IQ-TREE file
        iqtree_info = extract_alignment_info_from_iqtree(iqtree_file)
        
        stats = {
            'n_sequences': n_sequences,
            'alignment_length': alignment_length,
            'gap_percentage': gap_percentage,
            'constant_percentage': iqtree_info.get('constant_percentage'),
            'informative_percentage': iqtree_info.get('informative_percentage'),
            'constant_sites': iqtree_info.get('constant_sites'),
            'informative_sites': iqtree_info.get('informative_sites')
        }
        
        # Check if any filtering is enabled
        filtering_enabled = any([
            min_length is not None,
            min_sequences is not None,
            max_gaps is not None,
            min_informative is not None,
            max_constant is not None
        ])
        
        if not filtering_enabled:
            # No filtering - accept all alignments
            return True, stats
        
        # Apply only the filters that are set
        passes = True
        if min_length is not None:
            passes = passes and (alignment_length >= min_length)
        if min_sequences is not None:
            passes = passes and (n_sequences >= min_sequences)
        if max_gaps is not None:
            passes = passes and (gap_percentage <= max_gaps)
        if min_informative is not None and stats['informative_percentage'] is not None:
            passes = passes and (stats['informative_percentage'] >= min_informative)
        if max_constant is not None and stats['constant_percentage'] is not None:
            passes = passes and (stats['constant_percentage'] <= max_constant)
        
        return passes, stats
        
    except Exception as e:
        print(f"Error reading {alignment_file}: {e}")
        return False, None


def parse_iqtree_file(iqtree_file):
    """Parse IQ-TREE output file to extract GTR+Gamma parameters."""
    try:
        with open(iqtree_file, 'r') as f:
            content = f.read()
        
        params = {}
        
        # Extract GTR rate parameters
        rate_pattern = r'Rate parameter R:\s+A-C:\s+([\d.]+)\s+A-G:\s+([\d.]+)\s+A-T:\s+([\d.]+)\s+C-G:\s+([\d.]+)\s+C-T:\s+([\d.]+)\s+G-T:\s+([\d.]+)'
        rate_match = re.search(rate_pattern, content)
        if rate_match:
            params['AC'] = float(rate_match.group(1))
            params['AG'] = float(rate_match.group(2))
            params['AT'] = float(rate_match.group(3))
            params['CG'] = float(rate_match.group(4))
            params['CT'] = float(rate_match.group(5))
            params['GT'] = float(rate_match.group(6))
        else:
            return None
        
        # Extract base frequencies
        freq_pattern = r'pi\(A\) = ([\d.]+)\s+pi\(C\) = ([\d.]+)\s+pi\(G\) = ([\d.]+)\s+pi\(T\) = ([\d.]+)'
        freq_match = re.search(freq_pattern, content)
        if freq_match:
            params['pi_A'] = float(freq_match.group(1))
            params['pi_C'] = float(freq_match.group(2))
            params['pi_G'] = float(freq_match.group(3))
            params['pi_T'] = float(freq_match.group(4))
        else:
            return None
        
        # Extract Gamma alpha parameter
        alpha = None
        alpha_pattern = r'Gamma shape alpha:\s+([\d.]+)'
        alpha_match = re.search(alpha_pattern, content)
        if alpha_match:
            alpha = float(alpha_match.group(1))
            params['alpha'] = alpha
            params['rate_model'] = 'Gamma'
        
        return params
        
    except Exception as e:
        print(f"Error parsing {iqtree_file}: {e}")
        return None


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Handle --enable-filters flag
    if args.enable_filters:
        if args.min_length is None:
            args.min_length = 350
        if args.min_sequences is None:
            args.min_sequences = 15
        if args.max_gaps is None:
            args.max_gaps = 30.0
        if args.min_informative is None:
            args.min_informative = 20.0
        if args.max_constant is None:
            args.max_constant = 60.0
    
    # Determine if filtering is active
    filtering_active = any([
        args.min_length is not None,
        args.min_sequences is not None,
        args.max_gaps is not None,
        args.min_informative is not None,
        args.max_constant is not None
    ])
    
    print("="*70)
    print("GTR Parameter Extraction")
    print("="*70)
    
    if filtering_active:
        print("\n⚠ Quality filtering ENABLED:")
        if args.min_length is not None:
            print(f"  Minimum alignment length: {args.min_length} bp")
        if args.min_sequences is not None:
            print(f"  Minimum sequences: {args.min_sequences}")
        if args.max_gaps is not None:
            print(f"  Maximum gap percentage: {args.max_gaps}%")
        if args.min_informative is not None:
            print(f"  Minimum informative sites: {args.min_informative}%")
        if args.max_constant is not None:
            print(f"  Maximum constant sites: {args.max_constant}%")
    else:
        print("\n✓ No filtering - processing ALL genes")
    
    print("="*70 + "\n")
    
    # Process each dataset
    all_parameters = {}
    filter_stats = []
    
    for dataset_name, dataset_config in DATASETS.items():
        print(f"Processing {dataset_name}...")
        
        alignment_dir = Path(dataset_config['alignment_dir'])
        extensions = dataset_config['extensions']
        
        # Find all alignment files with any of the extensions
        alignment_files = []
        for ext in extensions:
            alignment_files.extend(list(alignment_dir.glob(f'*{ext}')))
        
        print(f"  Found {len(alignment_files)} alignment files")
        
        passed_alignments = []
        failed_alignments = []
        parameters_list = []
        
        for alignment_file in alignment_files:
            # Extract gene ID
            gene_id = extract_gene_id(alignment_file, dataset_name)
            if not gene_id:
                continue
            
            # Get corresponding IQ-TREE file first (we need it for filtering now)
            iqtree_path = find_iqtree_file(gene_id, dataset_config['iqtree_dir'], dataset_name)
            
            if not iqtree_path or not iqtree_path.exists():
                failed_alignments.append((gene_id, 'IQ-TREE file not found'))
                continue
            
            # Check alignment quality (now including informative/constant sites)
            passes, stats = check_alignment_quality(
                alignment_file,
                iqtree_path,
                args.min_length, 
                args.min_sequences, 
                args.max_gaps,
                args.min_informative,
                args.max_constant
            )
            
            if passes and stats:
                # Parse GTR parameters
                params = parse_iqtree_file(iqtree_path)
                
                if params:
                    params['source_file'] = str(iqtree_path)
                    params['alignment_file'] = str(alignment_file)
                    params['gene_id'] = gene_id
                    params.update(stats)
                    parameters_list.append(params)
                    passed_alignments.append((gene_id, stats))
            else:
                if stats and filtering_active:
                    reason = []
                    if args.min_length is not None and stats['alignment_length'] < args.min_length:
                        reason.append(f"length={stats['alignment_length']}")
                    if args.min_sequences is not None and stats['n_sequences'] < args.min_sequences:
                        reason.append(f"n_seq={stats['n_sequences']}")
                    if args.max_gaps is not None and stats['gap_percentage'] > args.max_gaps:
                        reason.append(f"gaps={stats['gap_percentage']:.1f}%")
                    if args.min_informative is not None and stats['informative_percentage'] is not None and stats['informative_percentage'] < args.min_informative:
                        reason.append(f"informative={stats['informative_percentage']:.1f}%")
                    if args.max_constant is not None and stats['constant_percentage'] is not None and stats['constant_percentage'] > args.max_constant:
                        reason.append(f"constant={stats['constant_percentage']:.1f}%")
                    failed_alignments.append((gene_id, ', '.join(reason)))
        
        all_parameters[dataset_name] = parameters_list
        
        # Report statistics
        n_passed = len(passed_alignments)
        n_failed = len(failed_alignments)
        n_total = n_passed + n_failed
        
        if filtering_active:
            print(f"  ✓ Passed filter: {n_passed}/{n_total} ({100*n_passed/n_total:.1f}%)")
            print(f"  ✗ Failed filter: {n_failed}/{n_total} ({100*n_failed/n_total:.1f}%)")
        else:
            print(f"  ✓ Processed: {n_passed} genes")
            if n_failed > 0:
                print(f"  ⚠ Skipped (missing IQ-TREE): {n_failed}")
        
        if parameters_list:
            df = pd.DataFrame(parameters_list)
            print(f"  Alpha range: [{df['alpha'].min():.2f}, {df['alpha'].max():.2f}]")
            print(f"  Alpha mean: {df['alpha'].mean():.2f}, median: {df['alpha'].median():.2f}")
        
        filter_stats.append({
            'Dataset': dataset_name,
            'Total': n_total,
            'Passed': n_passed,
            'Failed': n_failed,
            'Pass_Rate': n_passed/n_total if n_total > 0 else 0
        })
        
        print()
    
    # Overall statistics
    total_genes = sum(len(params) for params in all_parameters.values())
    print("="*70)
    print(f"TOTAL: Extracted GTR parameters from {total_genes} genes")
    print("="*70 + "\n")
    
    # Save filter statistics
    df_filter_stats = pd.DataFrame(filter_stats)
    suffix = '_filtered' if filtering_active else '_all'
    filter_csv = OUTPUT_DIR / f'alignment_stats{suffix}.csv'
    df_filter_stats.to_csv(filter_csv, index=False)
    print(f"✓ Statistics saved to {filter_csv}")
    
    # Save GTR parameters
    pkl_path = OUTPUT_DIR / f'gtr_parameters{suffix}.pkl'
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_parameters, f)
    print(f"✓ GTR parameters saved to {pkl_path}")
    
    # Compare alpha distributions
    print("\n" + "="*70)
    print("ALPHA DISTRIBUTION COMPARISON")
    print("="*70 + "\n")
    
    for dataset_name, params_list in all_parameters.items():
        if params_list:
            alphas = [p['alpha'] for p in params_list if 'alpha' in p]
            if alphas:
                print(f"{dataset_name}:")
                print(f"  N genes: {len(alphas)}")
                print(f"  Alpha: mean={np.mean(alphas):.3f}, median={np.median(alphas):.3f}")
                print(f"  Alpha: min={np.min(alphas):.3f}, max={np.max(alphas):.3f}")
                print(f"  Alpha > 3.0: {sum(1 for a in alphas if a > 3.0)} genes")
                print()
    
    # Summary
    all_alphas = [p['alpha'] for params in all_parameters.values() 
                  for p in params if 'alpha' in p]
    
    if all_alphas:
        print(f"Combined ({len(all_alphas)} genes):")
        print(f"  Alpha: mean={np.mean(all_alphas):.3f}, median={np.median(all_alphas):.3f}")
        print(f"  Alpha: min={np.min(all_alphas):.3f}, max={np.max(all_alphas):.3f}")
        print(f"  Alpha > 3.0: {sum(1 for a in all_alphas if a > 3.0)} genes ({100*sum(1 for a in all_alphas if a > 3.0)/len(all_alphas):.1f}%)")
    
    suffix_text = suffix.replace('_', ' ')
    print(f"\n✓ Done! Use gtr_parameters{suffix}.pkl for your simulations.")


if __name__ == '__main__':
    main()
