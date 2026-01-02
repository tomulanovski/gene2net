# Real Data Analysis Pipeline Guide

Complete guide for analyzing real phylogenetic network data from published papers, comparing different inference methods pairwise.

## Overview

This pipeline analyzes real phylogenetic network data where:
- **No ground truth**: We compare methods against each other (not against a known truth)
- **No replicates**: Each network has a single result per method
- **Optional published networks**: Some papers include published network files for comparison

The pipeline performs **pairwise comparisons** between all available methods for each network, generating summary tables and publication-quality figures.

## Directory Structure

```
gene2net/
├── papers/                    # Real data directory (on cluster)
│   └── {network_name}/        # e.g., Bendiksby_2011
│       └── networks/
│           ├── grampa/
│           │   └── network.tre
│           ├── polyphest/
│           │   └── network.tre
│           ├── padre/
│           │   └── network.tre
│           ├── mpallop/
│           │   └── network.tre
│           ├── alloppnet/
│           │   └── network.tre
│           └── paper/          # Optional: published network
│               └── network.tre
├── scripts/
│   ├── papers_config.yaml      # Configuration file
│   ├── collect_results.py       # Phase 1: Data collection
│   ├── compute_comparisons.py  # Phase 2: Pairwise comparisons
│   ├── summarize_results.py    # Phase 3: Summary tables
│   ├── create_analysis_figures.py  # Phase 4: Visualizations
│   ├── run_analysis.py         # Main orchestrator
│   └── REAL_DATA_ANALYSIS_GUIDE.md  # This file
└── analysis/                   # Output directory
    └── {timestamp}/
        ├── inventory.csv
        ├── comparisons_raw.csv
        ├── method_availability.csv
        ├── pairwise_summary.csv
        ├── per_network_comparisons.csv
        ├── method_rankings.csv
        ├── plots/
        │   ├── 01_method_availability.pdf/png
        │   ├── 02_availability_heatmap.pdf/png
        │   ├── 03_pairwise_*.pdf/png
        │   └── ...
        └── cache/               # Cached comparison results
```

## Prerequisites

1. **Python packages**: pandas, numpy, matplotlib, seaborn, pyyaml
2. **Dependencies**: The pipeline uses modules from `simulations/scripts/`:
   - `reticulate_tree.py`
   - `compare_reticulations.py`
3. **Data**: Papers directory with network files in the expected structure

## Configuration

Edit `scripts/papers_config.yaml` to configure:

```yaml
# Methods to analyze
methods:
  grampa:
    directory: "grampa"
    output_file: "network.tre"
  polyphest:
    directory: "polyphest"
    output_file: "network.tre"
  # ... etc

# Networks to analyze (must match directory names in papers/)
networks:
  - Bendiksby_2011
  - Koenen_2020
  # ... etc

# Paths (relative to gene2net root)
papers_dir: "papers"
output_dir: "analysis"
```

## Running the Pipeline

### Option 1: One-Command Execution (Recommended)

Run the complete pipeline with a single command:

```bash
# From gene2net root directory
python scripts/run_analysis.py
```

**Options:**
- `--force-recompute`: Force recompute all comparisons (ignore cache)
- `--output DIR`: Custom output directory (default: `analysis/`)
- `--dry-run`: Show what would be processed without computing
- `--no-figures`: Skip figure generation (faster for testing)
- `--config PATH`: Use custom config file (default: `scripts/papers_config.yaml`)

**Examples:**
```bash
# Force recompute everything
python scripts/run_analysis.py --force-recompute

# Custom output location
python scripts/run_analysis.py --output analysis/latest/

# Dry run to see what will be processed
python scripts/run_analysis.py --dry-run

# Skip figures (faster)
python scripts/run_analysis.py --no-figures
```

### Option 2: Step-by-Step Execution

Run each phase individually for more control:

#### Phase 1: Data Collection
```bash
python scripts/collect_results.py --export analysis/inventory.csv
```

**Output**: `inventory.csv` with columns:
- `network`: Network name
- `method`: Method name
- `network_path`: Path to network file
- `exists`: Boolean (file exists)
- `file_size`: File size in bytes

#### Phase 2: Compute Comparisons
```bash
python scripts/compute_comparisons.py \
    analysis/inventory.csv \
    analysis/cache/ \
    --export analysis/comparisons_raw.csv \
    --force-recompute
```

**Output**: `comparisons_raw.csv` with columns:
- `network`: Network name
- `method1`: First method in comparison
- `method2`: Second method in comparison
- `metric`: Metric name (e.g., `edit_distance_multree`, `rf_distance`)
- `value`: Metric value
- `status`: SUCCESS or FAILED

**Caching**: Results are cached in `cache/` directory. Use `--force-recompute` to ignore cache.

#### Phase 3: Generate Summary Tables
```bash
python scripts/summarize_results.py \
    analysis/comparisons_raw.csv \
    analysis/inventory.csv \
    analysis/
```

**Output**: Multiple CSV files in `analysis/`:
- `method_availability.csv`: Completion rates per method
- `pairwise_summary.csv`: Aggregated metrics per method pair
- `per_network_comparisons.csv`: Detailed per-network results
- `method_rankings.csv`: Overall method performance rankings

#### Phase 4: Generate Figures
```bash
python scripts/create_analysis_figures.py \
    --comparisons analysis/comparisons_raw.csv \
    --inventory analysis/inventory.csv \
    --output analysis/
```

**Output**: PDF and PNG figures in `analysis/plots/`:
- `01_method_availability.pdf/png`: Bar chart of completion rates
- `02_availability_heatmap.pdf/png`: Method × network availability matrix
- `03_pairwise_*_heatmap.pdf/png`: Pairwise distance heatmaps
- `04_pairwise_*_boxplot.pdf/png`: Distribution of distances by method pair
- `05_per_network_*.pdf/png`: Per-network comparison plots
- `06_method_rankings.pdf/png`: Method ranking chart

## Understanding the Output

### Summary Tables

#### `method_availability.csv`
Shows which methods succeeded for which networks:
- `method`: Method name
- `total_networks`: Total number of networks
- `available`: Number of networks with successful output
- `missing`: Number of networks with missing output
- `completion_rate`: Percentage of networks with available output

#### `pairwise_summary.csv`
Aggregated metrics across all networks for each method pair:
- `method1`, `method2`: Methods being compared
- `metric`: Metric name
- `mean`, `std`, `min`, `max`: Statistics across networks
- `n_networks`: Number of networks with this comparison

#### `per_network_comparisons.csv`
Wide-format table with one row per network:
- `network`: Network name
- `{method1}_vs_{method2}_{metric}`: Comparison values

#### `method_rankings.csv`
Overall method performance:
- `method`: Method name
- `avg_edit_distance_multree`: Average distance to other methods
- `rank_edit_distance_multree`: Ranking (lower = better)

### Key Metrics

The pipeline computes the following metrics (same as simulation analysis):

1. **`edit_distance_multree`**: Edit distance on MUL-trees (PRIMARY METRIC)
2. **`rf_distance`**: Robinson-Foulds distance on MUL-trees (PRIMARY METRIC)
3. **`edit_distance`**: Edit distance on folded networks (legacy)
4. **`num_rets_diff`**: Absolute difference in reticulation counts
5. **`num_rets_bias`**: Signed difference in reticulation counts (bias)
6. **`ploidy_diff.*`**: Polyploid identification metrics
7. **`ret_leaf_jaccard.*`**: Reticulation leaf set Jaccard similarity
8. **`ret_sisters_jaccard.*`**: Sister relationship Jaccard similarity

### Figures

- **Method Availability**: Shows which methods work on which networks
- **Pairwise Heatmaps**: Visualize average distances between methods
- **Boxplots**: Show distribution of distances for each method pair
- **Per-Network Plots**: Compare methods within each network
- **Rankings**: Overall method performance comparison

## Workflow on Cluster

Since the papers directory is on the cluster, follow this workflow:

1. **Local Development**: Test scripts locally (they'll report missing files)
2. **Push to GitHub**: Commit and push your changes
3. **Cluster Setup**: SSH to cluster and pull latest changes
4. **Run Analysis**: Execute pipeline on cluster where papers directory exists
5. **Download Results**: Transfer results back to local machine if needed

```bash
# On cluster
cd /path/to/gene2net
git pull
python scripts/run_analysis.py --output analysis/
```

## Troubleshooting

### Common Issues

#### 1. "Configuration file not found"
**Solution**: Make sure you're running from the gene2net root directory, or use `--config` to specify the full path.

#### 2. "No valid comparisons to summarize"
**Cause**: No method pairs have both outputs available for any network.
**Solution**: Check that at least 2 methods have outputs for some networks.

#### 3. "Failed to load network"
**Cause**: Network file is corrupted or in wrong format.
**Solution**: Check the network file manually, ensure it's valid Newick format.

#### 4. Import errors for `reticulate_tree` or `compare_reticulations`
**Solution**: The pipeline automatically adds `simulations/scripts/` to the path. If this fails, ensure you're running from the gene2net root directory.

#### 5. Empty output directory
**Cause**: All comparisons failed or no methods have outputs.
**Solution**: Check the completion report from Phase 1 to see what's available.

### Debugging Tips

1. **Start with dry run**: `python scripts/run_analysis.py --dry-run`
2. **Check inventory first**: Run Phase 1 and examine `inventory.csv`
3. **Test single comparison**: Manually test loading networks with Python
4. **Check cache**: Look in `cache/` directory to see what's been computed
5. **Read reports**: Check `comparison_report.txt` for detailed error messages

### Performance

- **Caching**: Comparisons are cached by default. First run is slow, subsequent runs are fast.
- **Parallelization**: Currently single-threaded. For large datasets, consider running in parallel manually.
- **Memory**: Each comparison loads two networks into memory. For many large networks, monitor memory usage.

## Comparison with Simulation Analysis

| Aspect | Simulation | Real Data |
|--------|-----------|-----------|
| Ground truth | Yes (true network) | No (optional published) |
| Replicates | Multiple per network | Single per network |
| Comparisons | GT vs Inferred | Method vs Method |
| Aggregation | Across replicates | Across networks |
| Focus | Accuracy metrics | Method differences |
| Missing data | Rare (controlled) | Common (methods fail) |

## Next Steps

After running the analysis:

1. **Review summary tables**: Check `method_availability.csv` and `pairwise_summary.csv`
2. **Examine figures**: Look at heatmaps and boxplots to understand method differences
3. **Identify patterns**: Which methods are most similar? Which networks are problematic?
4. **Compare to published**: If `paper` method exists, see how methods compare to published networks
5. **Generate reports**: Use the summary tables for publication

## Example Workflow

```bash
# 1. Check what data is available
python scripts/collect_results.py --verbose

# 2. Run complete analysis
python scripts/run_analysis.py

# 3. Check results
ls -lh analysis/*/plots/
cat analysis/*/method_availability.csv

# 4. If needed, recompute with different settings
python scripts/run_analysis.py --force-recompute --output analysis/latest/
```

## Support

For issues or questions:
1. Check the completion report from Phase 1
2. Review `comparison_report.txt` for errors
3. Examine the inventory CSV to see what's available
4. Test individual network loading manually

---

**Last Updated**: 2024
**Pipeline Version**: 1.0



