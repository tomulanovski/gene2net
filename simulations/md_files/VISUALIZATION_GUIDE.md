# Visualization and Analysis Guide

Complete guide for generating publication-quality figures and tables from phylogenetic network inference results.

## Overview

This guide covers the complete workflow from simulation results to publication-ready visualizations. It includes:

1. **Ground Truth Statistics**: Network characterization with WGD event detection
2. **Analysis Figures**: Comprehensive per-configuration visualizations
3. **Method Performance Evaluation**: Success rates, accuracy metrics, and difficulty analysis

## Quick Start

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Generate analysis figures for one configuration
python create_analysis_figures.py --config conf_ils_low_10M

# Generate for multiple configurations
python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M
```

**Output location**: `simulations/analysis/summary/{config}/plots/` and `tables/`

## Prerequisites

Before generating visualizations, ensure you have:

1. **Completed simulations** for all ILS levels (low, medium, high)
2. **Run postprocessing** for each configuration:
   ```bash
   python postprocess_results.py conf_ils_low_10M
   python postprocess_results.py conf_ils_medium_10M
   python postprocess_results.py conf_ils_high_10M
   ```
3. **Run summary pipeline** for each configuration:
   ```bash
   python run_full_summary.py conf_ils_low_10M
   python run_full_summary.py conf_ils_medium_10M
   python run_full_summary.py conf_ils_high_10M
   ```
4. **Updated network statistics** (automatically loaded from default location)

## Network Ground Truth Statistics

### What Are They?

Ground truth statistics characterize each of the 21 networks used in simulations. They distinguish between different types of polyploidy events.

**File**: `simulations/networks/mul_tree_final_stats.csv`

### Key Columns

- **Filename**: Network file name (e.g., `Ding_2023.tre`)
- **Num_Species**: Total unique species in the network
- **Num_Polyploids**: Number of species appearing >1 time in MUL-tree
- **Max_Copies**: Maximum number of copies any species has
- **H_Strict**: Number of reticulations (allopolyploidization/hybridization events)
- **H_Relaxed**: Reticulations with relaxed folding algorithm
- **H_Diff**: Difference between strict and relaxed
- **Num_Autopolyploidization_Events**: Number of autopolyploidization events
- **Total_WGD**: Total whole genome duplication events (auto + allo)
- **Polyploid_Names**: Comma-separated list of polyploid species

### Critical Distinctions

#### Autopolyploidization EVENT
- **Definition**: A single WGD that duplicates an entire clade
- **Detection**: Identical sibling subtrees in MUL-tree
- **Example**: Ding_2023 has ONE event (Rch clade duplication) → creates 7 autopolyploid species
- **Network folding**: These are simplified away (do NOT create reticulations)

#### Allopolyploidization/Hybridization
- **Definition**: Reticulation events where different lineages merge
- **Detection**: Identical non-sibling duplications persist as reticulations after network folding
- **Measure**: H_Strict count
- **Network folding**: These persist as reticulation nodes (2+ parents)

#### Total WGD Events
```
Total_WGD = Num_Autopolyploidization_Events + H_Strict
```

**Example (Ding_2023.tre)**:
- Num_Polyploids: 7 (Aro, CiP, Jma, Jmi, Jre, Pla, Rch all appear twice)
- H_Strict: 0 (no reticulations/hybridization)
- Num_Autopolyploidization_Events: 1 (ONE event duplicated the Rch clade)
- Total_WGD: 1 (1 auto + 0 allo)

**Example (Liu_2023.tre)**:
- Num_Polyploids: 8
- H_Strict: 6 (6 hybridization events)
- Num_Autopolyploidization_Events: 2 (2 WGD events)
- Total_WGD: 8 (2 auto + 6 allo)

### Regenerating Statistics

If you need to update the network statistics (e.g., after modifying network files):

```bash
cd simulations/scripts
python run_reticulation_stats.py ../networks/ ../networks/mul_tree_final_stats.csv
```

### Testing Autopolyploid Detection

Verify the autopolyploidization event detection is working correctly:

```bash
cd simulations/scripts
python test_autopolyploid_detection.py
```

Expected output for Ding_2023:
- ✓ Num_Autopolyploidization_Events: 1
- ✓ H_Strict: 0
- ✓ Total_WGD: 1

## Generating Analysis Figures

### The Analysis Script

**File**: `simulations/scripts/create_analysis_figures.py`

**Key Features**:
- Per-configuration output organization
- Non-interactive matplotlib backend (no X11 required - works on cluster)
- Automatic path resolution relative to script location
- Default network stats path (no need to specify every time)

### Usage

```bash
cd simulations/scripts

# Single configuration
python create_analysis_figures.py --config conf_ils_low_10M

# Multiple configurations (recommended)
python create_analysis_figures.py --config conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

# Custom network stats path (optional - has sensible default)
python create_analysis_figures.py --config conf_ils_low_10M --network-stats /custom/path/stats.csv
```

### Output Structure

Each configuration gets its own organized directory:

```
simulations/analysis/summary/{config}/
├── plots/
│   ├── 01_success_vs_reticulations.pdf      # Success rate vs H_Strict
│   ├── 01_success_vs_reticulations.png
│   ├── 02_success_vs_polyploids.pdf         # Success rate vs Num_Polyploids
│   ├── 02_success_vs_polyploids.png
│   ├── 03_success_vs_wgd.pdf                # Success rate vs Total_WGD
│   ├── 03_success_vs_wgd.png
│   ├── 04_method_performance_overview.pdf   # Edit distance & success panels
│   ├── 04_method_performance_overview.png
│   ├── 05_method_network_heatmap.pdf        # Edit distance heatmap
│   ├── 05_method_network_heatmap.png
│   ├── 06_reticulation_accuracy.pdf         # Scatter: true vs inferred H
│   ├── 06_reticulation_accuracy.png
│   ├── 07_difficulty_correlations.pdf       # Network properties correlation
│   └── 07_difficulty_correlations.png
└── tables/
    ├── summary_table.csv                     # Detailed per-method stats
    └── summary_simple.csv                    # Simplified summary
```

## Figure Descriptions

### Figure 1: Success Rate vs Reticulations
**File**: `01_success_vs_reticulations.pdf/png`

**Purpose**: Shows how the number of reticulation events affects each method's success rate.

**X-axis**: H_Strict (number of reticulations/hybridization events)
**Y-axis**: Success rate (%)

**Interpretation**:
- Methods with flat lines handle reticulations well
- Steep downward slopes indicate methods struggle with complex networks
- Compare methods to identify which are robust to high reticulation counts

### Figure 2: Success Rate vs Polyploids
**File**: `02_success_vs_polyploids.pdf/png`

**Purpose**: Shows how the number of polyploid species affects success rate.

**X-axis**: Num_Polyploids (species appearing >1 time)
**Y-axis**: Success rate (%)

**Interpretation**:
- Identifies methods robust to polyploidy complexity
- High polyploid counts indicate complex MUL-tree structures

### Figure 3: Success Rate vs Total WGD
**File**: `03_success_vs_wgd.pdf/png`

**Purpose**: Shows combined effect of all genome duplication events.

**X-axis**: Total_WGD (autopolyploidization + allopolyploidization events)
**Y-axis**: Success rate (%)

**Interpretation**:
- Combined view of genome duplication impact
- Shows overall WGD tolerance regardless of event type

### Figure 4: Method Performance Overview
**File**: `04_method_performance_overview.pdf/png`

**Layout**: Two panel figure
- **Left panel**: Mean edit distance per method (bar plot)
- **Right panel**: Success rate per method (bar plot)

**Purpose**: Quick comparison of method accuracy and reliability.

**Interpretation**:
- Lower edit distance = more accurate network inference
- Higher success rate = more robust/reliable method
- Ideal methods: Low edit distance + high success rate

### Figure 5: Method-Network Heatmap
**File**: `05_method_network_heatmap.pdf/png`

**Purpose**: Shows edit distance for each method on each network.

**Rows**: Methods (grampa, polyphest, padre, mpsugar)
**Columns**: Networks (21 networks sorted by complexity)
**Color**: Edit distance (blue = good, red = poor)

**Interpretation**:
- Identify problematic networks (red columns)
- Identify struggling methods (red rows)
- White cells = method failed to produce output

### Figure 6: Reticulation Accuracy
**File**: `06_reticulation_accuracy.pdf/png`

**Purpose**: How accurately does each method recover the correct number of reticulations?

**X-axis**: True H_Strict (ground truth)
**Y-axis**: Inferred H (from method output)
**Diagonal line**: Perfect accuracy

**Interpretation**:
- Points on diagonal = perfect inference
- Points above diagonal = overestimating reticulations
- Points below diagonal = underestimating reticulations

### Figure 7: Network Difficulty Correlations
**File**: `07_difficulty_correlations.pdf/png`

**Purpose**: Which network properties make inference difficult?

**Layout**: Correlation matrix heatmap
**Rows/Columns**: Network properties (H_Strict, Num_Polyploids, Total_WGD, etc.)
**Color**: Correlation strength (red = positive, blue = negative)

**Interpretation**:
- Positive correlation with edit distance = property makes inference harder
- Identify which network characteristics are most challenging

## Tables

### Summary Table
**File**: `summary_table.csv`

**Columns**:
- `method`: Method name
- `total_networks`: Total number of networks tested
- `successful_runs`: Number of successful inferences
- `success_rate`: Percentage of successful runs
- `mean_edit_distance`: Average edit distance (accuracy)
- `std_edit_distance`: Standard deviation of edit distance
- `min_edit_distance`: Best edit distance achieved
- `max_edit_distance`: Worst edit distance achieved

**Use**: Comprehensive method comparison for paper tables

### Summary Simple Table
**File**: `summary_simple.csv`

Simplified version with just essential metrics for quick overview.

## Data Aggregation

### Success Rate Calculation

For plots showing success rate vs network characteristics:

```python
# For each characteristic value (e.g., H_Strict = 5)
networks_with_value = networks where H_Strict == 5
total_runs = sum of all runs across these networks
successful_runs = sum of successful runs across these networks
success_rate = (successful_runs / total_runs) * 100
```

**Example**: If 3 networks all have H_Strict=5, and each has 5 replicates:
- Total runs = 3 networks × 5 replicates = 15
- If 12 succeed: success_rate = 12/15 × 100 = 80%

### Edit Distance Calculation

Edit distances are calculated per network:
1. For each network, average edit distance across successful replicates
2. Plot distribution across all networks for each method
3. Heatmap shows per-network mean edit distance

## Customization

### Color Schemes

Edit `create_analysis_figures.py` to customize method colors:

```python
METHOD_COLORS = {
    'grampa': '#1f77b4',    # Blue
    'polyphest': '#ff7f0e', # Orange
    'padre': '#2ca02c',     # Green
    'mpsugar': '#d62728'    # Red
}
```

### Figure Sizes

Modify publication settings at top of script:

```python
plt.rcParams['figure.dpi'] = 300        # Resolution
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11          # Base font size
plt.rcParams['axes.labelsize'] = 12     # Axis labels
plt.rcParams['axes.titlesize'] = 13     # Titles
plt.rcParams['legend.fontsize'] = 10    # Legend
```

## Troubleshooting

### Error: "FileNotFoundError: inventory.csv"

**Problem**: Summary pipeline not run for specified configuration.

**Solution**:
```bash
cd simulations/scripts
python run_full_summary.py conf_ils_low_10M
```

### Error: "KeyError: 'inferred_exists'"

**Problem**: Old inventory format without boolean success column.

**Solution**: Re-run postprocessing pipeline:
```bash
python postprocess_results.py conf_ils_low_10M
python run_full_summary.py conf_ils_low_10M
```

### Empty plots or missing data

**Problem**: Network names don't match between inventory and stats files.

**Solution**:
- Check that network names in `mul_tree_final_stats.csv` match those in `inventory.csv`
- Names should be without `.tre` extension
- Case-sensitive matching required

### X11 display errors on cluster

**Problem**: Matplotlib trying to open display window.

**Solution**: Script already uses non-interactive backend (`matplotlib.use('Agg')`). If still having issues:
```bash
export DISPLAY=
python create_analysis_figures.py --config conf_ils_low_10M
```

## Workflow Checklist

Use this checklist to ensure all steps are completed:

- [ ] Simulations completed for all ILS levels (low, medium, high)
- [ ] Postprocessing run for each configuration
- [ ] Summary pipeline run for each configuration
- [ ] Network statistics file exists: `simulations/networks/mul_tree_final_stats.csv`
- [ ] Generated analysis figures for each configuration
- [ ] Reviewed all 7 plots per configuration
- [ ] Checked tables for method comparison
- [ ] Identified key figures for main paper vs supplementary

## Next Steps

After generating visualizations:

1. **Review all figures**: Check each plot makes sense and shows expected patterns

2. **Select key figures for paper**: Typically 3-4 main figures + rest in supplementary

3. **Extract insights**:
   - Which methods perform best overall?
   - How does ILS level affect performance?
   - What network properties correlate with difficulty?
   - Are autopolyploids harder to infer than allopolyploids?

4. **Statistical analysis** (if needed):
   - Pairwise method comparisons (Wilcoxon signed-rank test)
   - ANOVA for ILS level effects
   - Regression analysis for network property effects

5. **Customize for publication**:
   - Adjust colors to journal requirements
   - Modify labels/legends as needed
   - Export at required resolution

## File Locations Reference

```
gene2net/
└── simulations/
    ├── networks/
    │   ├── *.tre                           # 21 input MUL-trees
    │   └── mul_tree_final_stats.csv        # Ground truth statistics
    ├── scripts/
    │   ├── run_reticulation_stats.py       # Generate network stats
    │   ├── test_autopolyploid_detection.py # Test autopolyploid detection
    │   ├── postprocess_results.py          # Required preprocessing
    │   ├── run_full_summary.py             # Required preprocessing
    │   └── create_analysis_figures.py      # Main visualization script
    ├── analysis/
    │   └── summary/
    │       ├── conf_ils_low_10M/
    │       │   ├── inventory.csv           # Input: run inventory
    │       │   ├── aggregated_metrics.csv  # Input: per-network metrics
    │       │   ├── comparisons/            # Input: network comparisons
    │       │   ├── plots/                  # Output: 7 figures (PDF + PNG)
    │       │   └── tables/                 # Output: summary tables
    │       ├── conf_ils_medium_10M/        # Same structure
    │       └── conf_ils_high_10M/          # Same structure
    └── md_files/
        └── VISUALIZATION_GUIDE.md          # This file
```

## Related Documentation

- **COMPLETE_PIPELINE_GUIDE.md**: Full simulation workflow from start to finish
- **SUMMARY_PIPELINE_GUIDE.md**: Details on summary pipeline (prerequisite)
- **POSTPROCESSING_GUIDE.md**: Details on postprocessing (prerequisite)
- **CLAUDE.md**: Repository-wide guidance

## Support

For issues or questions:
1. Check this guide first
2. Review related documentation (see above)
3. Examine example outputs in `simulations/analysis/summary/`
4. Check `CLAUDE.md` for repository context
