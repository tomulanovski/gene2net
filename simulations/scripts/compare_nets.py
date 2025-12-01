from pathlib import Path
from reticulate_tree import ReticulateTree
from compare_reticulations import one_stop_compare
import pandas as pd

def detect_tree_format(newick_str):
    """
    Detect if a tree string is extended Newick (with #H markers) or standard Newick.
    
    Returns:
        'enewick' if extended Newick with reticulations
        'newick' if standard Newick (potentially MUL-tree)
    """
    newick_str = newick_str.strip()
    if '#H' in newick_str or '#UID' in newick_str or '//' in newick_str:
        return 'enewick'
    return 'newick'

def collect_mul_trees(main_directory, ground_truth_path, output_filename='output.tre'):
    """
    Collect trees from subdirectories and ground truth.
    Automatically detects format (standard Newick vs extended Newick).
    
    Args:
        main_directory: Path to main directory containing method subdirectories
        ground_truth_path: Path to ground truth tree file
        output_filename: Name of output files to look for (default: 'output.tre')
    
    Returns:
        Dictionary with tree names as keys and tree strings as values,
        plus a format dictionary indicating each tree's format
    """
    main_dir = Path(main_directory)
    trees = {}
    formats = {}
    
    # Load ground truth
    try:
        with open(ground_truth_path, 'r') as f:
            tree_str = f.read().strip()
            trees['ground_truth'] = tree_str
            formats['ground_truth'] = detect_tree_format(tree_str)
        print(f"Loaded ground truth from {ground_truth_path}")
        print(f"  Format detected: {formats['ground_truth']}")
    except Exception as e:
        print(f"Error loading ground truth: {e}")
        return None, None
    
    # Iterate through subdirectories
    for subdir in main_dir.iterdir():
        if subdir.is_dir():
            output_file = subdir / output_filename
            
            # Check if output.tre exists
            if output_file.exists():
                try:
                    with open(output_file, 'r') as f:
                        tree_str = f.read().strip()
                        trees[subdir.name] = tree_str
                        formats[subdir.name] = detect_tree_format(tree_str)
                    print(f"Loaded {subdir.name} from {output_file}")
                    print(f"  Format detected: {formats[subdir.name]}")
                except Exception as e:
                    print(f"Error loading {subdir.name}: {e}")
            else:
                print(f"Skipping {subdir.name} - no {output_filename} found")
    
    return trees, formats

def run_comparison_analysis(main_directory, ground_truth_path, 
                           output_dir=None, 
                           use_polyphest=False,
                           threshold=0.2,
                           normalize=True,
                           visualize_trees=False):
    """
    Complete comparison pipeline for phylogenetic trees/networks.
    Automatically handles both MUL-trees (standard Newick) and 
    networks (extended Newick with #H markers).
    
    Args:
        main_directory: Path to directory with method subdirectories
        ground_truth_path: Path to ground truth tree
        output_dir: Where to save results (default: main_directory/comparison_results)
        use_polyphest: Use Polyphest fuzzy matching for MUL-trees (ignored for extended Newick)
        threshold: Threshold for Polyphest mode (only used if use_polyphest=True)
        normalize: Normalize edit distance in Polyphest mode
        visualize_trees: Generate visualizations for each tree
    """
    
    # Setup output directory
    if output_dir is None:
        output_dir = Path(main_directory) / 'comparison_results'
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print("="*60)
    print("PHYLOGENETIC NETWORK COMPARISON PIPELINE")
    print("="*60)
    print(f"Main directory: {main_directory}")
    print(f"Ground truth: {ground_truth_path}")
    print(f"Output directory: {output_dir}")
    print(f"MUL-tree conversion mode: {'POLYPHEST' if use_polyphest else 'STRICT'}")
    if use_polyphest:
        print(f"  - Threshold: {threshold}")
        print(f"  - Normalize: {normalize}")
    print("="*60 + "\n")
    
    # Step 1: Collect all trees
    print("STEP 1: Collecting trees...")
    trees, formats = collect_mul_trees(main_directory, ground_truth_path)
    
    if trees is None or len(trees) <= 1:
        print("\nNot enough trees to compare. Exiting.")
        return None
    
    print(f"\nCollected {len(trees)} trees: {list(trees.keys())}")
    print(f"  Formats: {formats}\n")
    
    # Step 2: Convert and analyze
    print("STEP 2: Converting to networks and computing measurements...")
    
    results = []
    logs = []
    
    for name, tree_str in trees.items():
        try:
            print(f"  Processing {name} ({formats[name]})...")
            
            # Choose initialization parameters based on format
            if formats[name] == 'enewick':
                # Extended Newick - direct network representation
                # No MUL-tree conversion needed
                rt = ReticulateTree(tree_str)
                print(f"    Loaded as extended Newick (network)")
                
            else:  # 'newick'
                # Standard Newick - might be MUL-tree
                rt_temp = ReticulateTree(tree_str)
                
                if rt_temp.check_duplicated():
                    # It's a MUL-tree - apply conversion
                    conversion_params = {'is_multree': True}
                    if use_polyphest:
                        conversion_params['threshold'] = threshold
                        conversion_params['normalize'] = normalize
                    rt = ReticulateTree(tree_str, **conversion_params)
                    print(f"    Detected MUL-tree, applied {('POLYPHEST' if use_polyphest else 'STRICT')} conversion")
                else:
                    # Simple tree - no conversion needed
                    rt = rt_temp
                    print(f"    Loaded as simple tree (no duplicates)")
            
            measurements = rt.measure(printout=False)
            results.append({'name': name, 'object': rt, **measurements})
            
            # Optional visualization
            if visualize_trees:
                viz_file = output_dir / f"{name}_network.png"
                rt.visualize(filename=str(viz_file))
                print(f"    Saved visualization to {viz_file}")
                
        except Exception as e:
            error_msg = f"Error processing {name}: {e}"
            print(f"    {error_msg}")
            logs.append(error_msg)
    
    if not results:
        print("\nNo trees successfully processed. Exiting.")
        return None
    
    # Step 3: Create comparison matrices
    print("\nSTEP 3: Computing pairwise comparisons...")
    try:
        df = pd.DataFrame(results).set_index('name')
        from compare_reticulations import pairwise_comparison
        comparison_matrices = pairwise_comparison(df, debug=False)
        print("Pairwise comparisons complete\n")
    except Exception as e:
        print(f"Error in pairwise comparison: {e}")
        return None
    
    # Step 4: Save results
    print("STEP 4: Saving results...")
    
    # Save measurements table
    measurements_file = output_dir / 'measurements.csv'
    df_export = df.drop(columns=['object'])  # Don't try to save the objects
    
    # Convert complex columns to string representation for CSV
    for col in df_export.columns:
        if df_export[col].dtype == 'object':
            df_export[col] = df_export[col].astype(str)
    
    df_export.to_csv(measurements_file)
    print(f"  Saved measurements to {measurements_file}")
    
    # Save comparison matrices
    for metric_name, matrix in comparison_matrices.items():
        matrix_file = output_dir / f'{metric_name}_matrix.csv'
        matrix.to_csv(matrix_file)
        print(f"  Saved {metric_name} matrix to {matrix_file}")
    
    # Step 5: Generate summary focusing on ground truth
    print("\nSTEP 5: Generating summary report...")
    summary_file = output_dir / 'summary_report.txt'
    
    with open(summary_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("PHYLOGENETIC NETWORK COMPARISON SUMMARY REPORT\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"Ground truth: {ground_truth_path}\n")
        f.write(f"Ground truth format: {formats.get('ground_truth', 'unknown')}\n")
        f.write(f"MUL-tree conversion mode: {'POLYPHEST' if use_polyphest else 'STRICT'}\n")
        if use_polyphest:
            f.write(f"  Threshold: {threshold}, Normalize: {normalize}\n")
        f.write(f"\nMethods compared: {[name for name in df.index if name != 'ground_truth']}\n")
        f.write(f"Formats: {formats}\n\n")
        
        f.write("="*60 + "\n")
        f.write("DISTANCE FROM GROUND TRUTH (lower is better)\n")
        f.write("="*60 + "\n\n")
        
        # Extract ground truth comparisons
        for metric_name, matrix in comparison_matrices.items():
            f.write(f"\n{metric_name.upper().replace('_', ' ')}:\n")
            f.write("-" * 40 + "\n")
            
            if 'ground_truth' in matrix.index:
                gt_row = matrix.loc['ground_truth']
                # Sort by distance (excluding self-comparison)
                sorted_methods = gt_row[gt_row.index != 'ground_truth'].sort_values()
                
                for method, distance in sorted_methods.items():
                    f.write(f"  {method:20s}: {distance:.4f}\n")
            else:
                f.write("  Ground truth not found in comparison\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write("BASIC STATISTICS\n")
        f.write("="*60 + "\n\n")
        
        for name in df.index:
            f.write(f"\n{name}:\n")
            f.write(f"  Format: {formats.get(name, 'unknown')}\n")
            f.write(f"  Reticulation count: {df.loc[name, 'reticulation_count']}\n")
            f.write(f"  Leaf counts: {df.loc[name, 'leaf_counts']}\n")
    
    print(f"  Saved summary report to {summary_file}")
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)
    print(f"\nAll results saved to: {output_dir}")
    print(f"\nKey files:")
    print(f"  - summary_report.txt: Quick overview of results")
    print(f"  - measurements.csv: Detailed measurements for each tree")
    print(f"  - *_matrix.csv: Full pairwise comparison matrices")
    
    return {
        'data': df,
        'comparisons': comparison_matrices,
        'formats': formats,
        'logs': logs
    }


if __name__ == '__main__':
    
    # Configuration
    MAIN_DIR = '/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/Shahrestani_2015/results/ILS_low_dup_low/'
    GROUND_TRUTH = '/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks/Shahrestani_2015.tre'
    
    # Run analysis with STRICT mode (exact matching for MUL-trees)
    results = run_comparison_analysis(
        main_directory=MAIN_DIR,
        ground_truth_path=GROUND_TRUTH,
        use_polyphest=False,  # Change to True for fuzzy matching on MUL-trees
        visualize_trees=False  # Change to True to generate network visualizations
    )
    
    # Optional: Run with POLYPHEST mode for comparison
    # results_polyphest = run_comparison_analysis(
    #     main_directory=MAIN_DIR,
    #     ground_truth_path=GROUND_TRUTH,
    #     output_dir=Path(MAIN_DIR) / 'comparison_results_polyphest',
    #     use_polyphest=True,
    #     threshold=0.2,
    #     normalize=True,
    #     visualize_trees=False
    # )