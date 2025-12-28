# Comprehensive Visualization and Analysis Guide

Complete guide for generating publication-quality figures and tables from phylogenetic network inference results.

## Overview

This guide covers the complete workflow from simulation results to publication-ready visualizations. After completing simulations for ILS low, medium, and high levels, this pipeline generates **42 comprehensive plots** organized into 4 categories:

1. **Completion Rate Analysis (12 plots)**: How network characteristics affect method success
2. **Accuracy vs Network Properties (12 plots)**: How network properties affect reconstruction accuracy
3. **Advanced Performance Metrics (7 plots)**: Fine-grained accuracy analysis (Jaccard distances, F1 scores, distributions)
4. **Distributions & Summaries (11 plots)**: Overall performance, correlations, and comparisons

Plus **2 comprehensive tables** for publication.

## Quick Start

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Generate analysis figures for one configuration
python create_analysis_figures.py --config conf_ils_low_10M

# Generate for multiple configurations (recommended)
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
- **H_Strict**: Number of reticulations using Holm folding (strict)
- **H_Relaxed**: Reticulations with Polyphest relaxed folding (θ=0.2)
- **H_Diff**: Difference between strict and relaxed
- **Num_Autopolyploidization_Events**: Number of autopolyploidization events
- **Total_WGD**: Total whole genome duplication events (auto + allo)
- **Polyploid_Names**: Comma-separated list of polyploid species

### Derived Metrics

The pipeline automatically calculates:
- **Polyploid_Ratio**: Num_Polyploids / Num_Species (proportion of polyploid species)
- **Ret_Density**: H_Strict / Num_Species (reticulations per species)

## Generating Analysis Figures

### The Analysis Script

**File**: `simulations/scripts/create_analysis_figures.py`

**Key Features**:
- Per-configuration output organization
- Non-interactive matplotlib backend (no X11 required - works on cluster)
- Automatic path resolution relative to script location
- Default network stats path (no need to specify every time)
- Both combined plots (all methods) AND faceted plots (per-method subplots)
- Error bars showing variability across networks
- Scatter plots for discrete data (completion rates)
- Automatic cleanup of old plots and tables before generation
- Organized into 4 clear categories with progress tracking

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
│   # CATEGORY 1: Completion Rate vs Network Characteristics (12 plots)
│   ├── 01_combined_completion_vs_h_strict.pdf/png
│   ├── 02_combined_completion_vs_h_relaxed.pdf/png
│   ├── 03_combined_completion_vs_polyploids.pdf/png
│   ├── 04_combined_completion_vs_total_wgd.pdf/png
│   ├── 05_combined_completion_vs_num_species.pdf/png
│   ├── 06_combined_completion_vs_max_copies.pdf/png
│   │
│   # CATEGORY 2: Accuracy vs Network Characteristics (12 plots)
│   ├── 11_combined_editdist_multree_vs_num_species.pdf/png
│   ├── 12_combined_editdist_multree_vs_h_strict.pdf/png
│   ├── 13_combined_editdist_multree_vs_polyploids.pdf/png
│   ├── 14_combined_editdist_multree_vs_max_copies.pdf/png
│   ├── 15_combined_rf_vs_num_species.pdf/png
│   ├── 16_combined_rf_vs_h_strict.pdf/png
│   │
│   # CATEGORY 3: Advanced Metrics (7 plots)
│   ├── 21_combined_ret_leaf_jaccard_vs_h_strict.pdf/png
│   ├── 22_combined_ret_sisters_jaccard_vs_h_strict.pdf/png
│   ├── 23_polyploid_f1_performance.pdf/png
│   ├── 08d_ret_leaf_jaccard_distribution.pdf/png
│   ├── 08e_ret_sisters_jaccard_distribution.pdf/png
│   │
│   # CATEGORY 4: Distributions, Comparisons, Summaries (11 plots)
│   ├── 05_folding_comparison.pdf/png
│   ├── 06_folding_accuracy_comparison.pdf/png
│   ├── 07_reticulation_bias_boxplot.pdf/png
│   ├── 08_edit_distance_multree_boxplot.pdf/png
│   ├── 08b_edit_distance_multree_distribution.pdf/png
│   ├── 08c_rf_distance_distribution.pdf/png
│   ├── 09_per_network_breakdown.pdf/png
│   ├── 10_method_summary.pdf/png
│   ├── 31_comprehensive_correlation_heatmap.pdf/png
│   ├── 32_per_method_correlation_heatmap.pdf/png
│   │
│   └── individual_methods/                                    # Faceted subplot versions
│       ├── 01_faceted_completion_vs_h_strict.pdf/png
│       ├── 02_faceted_completion_vs_h_relaxed.pdf/png
│       ├── 03_faceted_completion_vs_polyploids.pdf/png
│       ├── 04_faceted_completion_vs_total_wgd.pdf/png
│       ├── 05_faceted_completion_vs_num_species.pdf/png
│       ├── 06_faceted_completion_vs_max_copies.pdf/png
│       ├── 11_faceted_editdist_multree_vs_num_species.pdf/png
│       ├── 12_faceted_editdist_multree_vs_h_strict.pdf/png
│       ├── 13_faceted_editdist_multree_vs_polyploids.pdf/png
│       ├── 14_faceted_editdist_multree_vs_max_copies.pdf/png
│       ├── 15_faceted_rf_vs_num_species.pdf/png
│       ├── 16_faceted_rf_vs_h_strict.pdf/png
│       ├── 21_faceted_ret_leaf_jaccard_vs_h_strict.pdf/png
│       └── 22_faceted_ret_sisters_jaccard_vs_h_strict.pdf/png
│
└── tables/
    ├── 01_method_performance_summary.csv
    └── 02_per_network_performance.csv
```

**Total output per configuration**: 42 plots (22 combined + 20 faceted) + 2 tables

**Note**: The script automatically cleans old plots and tables from the output directories before generating new ones, ensuring you always have the latest visualizations. Summary files from `run_full_summary` (like `inventory.csv`, `aggregated_metrics.csv`) are preserved.

## Figure Descriptions

### CATEGORY 1: Completion Rate vs Network Characteristics

These plots show how often each method successfully completes inference for networks with different characteristics. **Data is discrete** (scatter plots with error bars, no connecting lines).

#### Combined vs Faceted Views

For each characteristic, TWO versions are generated:

- **Combined plots** (`plots/XX_combined_*.pdf`): All methods on one graph
  - Best for: Direct method comparison
  - Includes: Error bars, legend, scatter points (no lines)

- **Faceted plots** (`plots/individual_methods/XX_faceted_*.pdf`): Grid of subplots, one per method
  - Best for: Detailed per-method trends without overlapping points
  - Cleaner for complex patterns

#### Figure 01-02: Completion Rate vs Reticulations (Holm Fold)
**Files**:
- `01_combined_completion_vs_h_strict.pdf/png`
- `individual_methods/01_faceted_completion_vs_h_strict.pdf/png`

**X-axis**: Number of Reticulations (Holm Fold)
**Y-axis**: Completion Rate (%)

**Interpretation**:
- Flat patterns = method robust to reticulations
- Downward trends = method struggles with network complexity
- Error bars show variability across networks with same H_Strict value
- Scatter points (no lines) reflect discrete data nature

#### Figure 03-04: Completion Rate vs Reticulations (Polyphest Fold)
**Files**:
- `02_combined_completion_vs_h_relaxed.pdf/png`
- `individual_methods/02_faceted_completion_vs_h_relaxed.pdf/png`

**X-axis**: Number of Reticulations (Polyphest Fold)
**Y-axis**: Completion Rate (%)

**Interpretation**:
- Compare with Holm Fold to see if folding algorithm affects difficulty
- H_Relaxed may differ from H_Strict for the same network
- Identifies folding-sensitive methods

#### Figure 05-06: Completion Rate vs Polyploids
**Files**:
- `03_combined_completion_vs_polyploids.pdf/png`
- `individual_methods/03_faceted_completion_vs_polyploids.pdf/png`

**X-axis**: Number of Polyploid Species
**Y-axis**: Completion Rate (%)

**Interpretation**:
- Identifies methods robust to polyploidy
- High polyploid counts = complex MUL-tree structure
- Different from reticulation count (species-level vs event-level)

#### Figure 07-08: Completion Rate vs Total WGD
**Files**:
- `04_combined_completion_vs_total_wgd.pdf/png`
- `individual_methods/04_faceted_completion_vs_total_wgd.pdf/png`

**X-axis**: Total WGD Events (autopolyploidization + allopolyploidization)
**Y-axis**: Completion Rate (%)

**Interpretation**:
- Combined view of all genome duplication events
- Shows overall WGD tolerance regardless of event type

#### Figure 09-10: Completion Rate vs Num Species
**Files**:
- `05_combined_completion_vs_num_species.pdf/png`
- `individual_methods/05_faceted_completion_vs_num_species.pdf/png`

**X-axis**: Number of Species
**Y-axis**: Completion Rate (%)

**Why important**:
- Phylogeny size is a fundamental complexity factor
- Identifies methods that scale well vs those that struggle with large datasets
- Independent of polyploidy/reticulation complexity

**Interpretation**:
- Upward trend = method benefits from more data
- Downward trend = method struggles with computational complexity
- Flat pattern = method unaffected by phylogeny size

#### Figure 11-12: Completion Rate vs Max Copies
**Files**:
- `06_combined_completion_vs_max_copies.pdf/png`
- `individual_methods/06_faceted_completion_vs_max_copies.pdf/png`

**X-axis**: Maximum Copies per Species
**Y-axis**: Completion Rate (%)

**Why important**:
- Max_Copies indicates MUL-tree structural complexity
- Species with 9 copies create highly complex tree structures
- Different from Num_Polyploids (one highly duplicated species vs many 2× species)

**Interpretation**:
- High Max_Copies = extreme duplication of single lineage
- Tests method robustness to local vs distributed complexity

---

### CATEGORY 2: Accuracy vs Network Characteristics

These plots show how accurately methods reconstruct networks, and how accuracy varies with network properties.

**Key concept**: Edit distance and RF distance measure topological accuracy (0 = perfect, higher = worse). These plots reveal which network properties make accurate reconstruction harder.

#### Figure 13-14: Edit Distance vs Num Species
**Files**:
- `11_combined_editdist_multree_vs_num_species.pdf/png`
- `individual_methods/11_faceted_editdist_multree_vs_num_species.pdf/png`

**X-axis**: Number of Species
**Y-axis**: Normalized Edit Distance (MUL-tree)

**Why important**: Does phylogeny size affect accuracy?

**Interpretation**:
- Upward trend = larger phylogenies are harder to reconstruct accurately
- Flat line = method maintains accuracy regardless of size
- Compare slopes between methods to identify which scale better

#### Figure 15-16: Edit Distance vs H_Strict
**Files**:
- `12_combined_editdist_multree_vs_h_strict.pdf/png`
- `individual_methods/12_faceted_editdist_multree_vs_h_strict.pdf/png`

**X-axis**: Number of Reticulations (Holm Fold)
**Y-axis**: Normalized Edit Distance (MUL-tree)

**Why important**: Do more reticulations make accurate inference harder?

**Interpretation**:
- Upward trend = reticulation complexity degrades accuracy
- Method comparison reveals which methods handle reticulations better
- Critical for understanding reticulation limits

#### Figure 17-18: Edit Distance vs Num Polyploids
**Files**:
- `13_combined_editdist_multree_vs_polyploids.pdf/png`
- `individual_methods/13_faceted_editdist_multree_vs_polyploids.pdf/png`

**X-axis**: Number of Polyploid Species
**Y-axis**: Normalized Edit Distance (MUL-tree)

**Why important**: Does polyploid count affect reconstruction quality?

**Interpretation**:
- Tests if methods can maintain accuracy with many polyploid species
- Identifies polyploidy-robust methods

#### Figure 19-20: Edit Distance vs Max Copies
**Files**:
- `14_combined_editdist_multree_vs_max_copies.pdf/png`
- `individual_methods/14_faceted_editdist_multree_vs_max_copies.pdf/png`

**X-axis**: Maximum Copies per Species
**Y-axis**: Normalized Edit Distance (MUL-tree)

**Why important**: Does extreme duplication affect accuracy?

**Interpretation**:
- High Max_Copies = extreme local complexity
- Tests if methods struggle with highly duplicated lineages

#### Figure 21-22: RF Distance vs Num Species
**Files**:
- `15_combined_rf_vs_num_species.pdf/png`
- `individual_methods/15_faceted_rf_vs_num_species.pdf/png`

**X-axis**: Number of Species
**Y-axis**: Robinson-Foulds Distance (MUL-tree)

**Interpretation**: Similar to edit distance, but using RF metric. RF distance measures splits/topology differences.

#### Figure 23-24: RF Distance vs H_Strict
**Files**:
- `16_combined_rf_vs_h_strict.pdf/png`
- `individual_methods/16_faceted_rf_vs_h_strict.pdf/png`

**X-axis**: Number of Reticulations (Holm Fold)
**Y-axis**: Robinson-Foulds Distance (MUL-tree)

**Interpretation**: How reticulation complexity affects RF-based accuracy measurement.

---

### CATEGORY 3: Advanced Performance Metrics

Beyond overall topology accuracy, these plots assess specific aspects of network inference quality.

#### Figure 25-26: Reticulation Leaf Jaccard Distance vs H_Strict
**Files**:
- `21_combined_ret_leaf_jaccard_vs_h_strict.pdf/png`
- `individual_methods/21_faceted_ret_leaf_jaccard_vs_h_strict.pdf/png`

**X-axis**: Number of Reticulations (Holm Fold)
**Y-axis**: Reticulation Leaf Set Distance (0 = perfect match, 1 = no match)

**What it measures**: How well do inferred reticulation leaf sets match true reticulation leaf sets?

**Why important**:
- Edit distance is global; this is reticulation-specific
- Distance = 1 - Jaccard similarity
- 0.0 = perfect reticulation leaf recovery
- 1.0 = no overlap

**Interpretation**:
- Lower values = method correctly identifies which leaves are under reticulations
- Upward trend = more reticulations make identification harder
- Method comparison reveals reticulation detection accuracy

#### Figure 27-28: Sister Relationship Jaccard Distance vs H_Strict
**Files**:
- `22_combined_ret_sisters_jaccard_vs_h_strict.pdf/png`
- `individual_methods/22_faceted_ret_sisters_jaccard_vs_h_strict.pdf/png`

**X-axis**: Number of Reticulations (Holm Fold)
**Y-axis**: Sister Relationship Distance (0 = perfect match, 1 = no match)

**What it measures**: How well do inferred sister relationships match true sister relationships?

**Why important**:
- Sister relationships define network topology
- More stringent than leaf sets (requires correct pairing)
- Critical for phylogenetic accuracy

**Interpretation**:
- Lower values = correct topological relationships
- Higher than leaf Jaccard = method gets leaves right but relationships wrong

#### Figure 29: Polyploid Identification F1 Score
**File**: `23_polyploid_f1_performance.pdf/png`

**Layout**: Two-panel figure
- **Left panel**: F1 score per method (bar chart)
- **Right panel**: Precision vs Recall per method (grouped bars)

**What it measures**: How accurately do methods identify which species are polyploid?

**Metrics**:
- **Precision**: TP / (TP + FP) - of identified polyploids, how many are correct?
- **Recall**: TP / (TP + FN) - of true polyploids, how many did we find?
- **F1**: Harmonic mean of precision and recall (overall accuracy)

**Why important**:
- Polyploid identification is a key biological question
- High precision = few false positives (reliable calls)
- High recall = few false negatives (don't miss polyploids)
- F1 balances both

**Interpretation**:
- F1 = 1.0: Perfect identification
- F1 < 0.5: Poor identification
- Compare precision vs recall to understand method behavior

#### Figure 30: Reticulation Leaf Jaccard Distance Distribution
**File**: `08d_ret_leaf_jaccard_distribution.pdf/png`

**Layout**: Boxplot showing distribution of reticulation leaf set distances per method

**Interpretation**:
- Lower values = better performance
- Box shows quartiles, whiskers show range
- Outliers indicate networks where method struggled
- Compare medians between methods

#### Figure 31: Sister Relationship Jaccard Distance Distribution
**File**: `08e_ret_sisters_jaccard_distribution.pdf/png`

**Layout**: Boxplot showing distribution of sister relationship distances per method

**Interpretation**:
- Lower values = better performance
- More stringent metric than leaf Jaccard
- Shows consistency of topological relationship recovery

---

### CATEGORY 4: Distributions, Comparisons, and Summaries

#### Figure 32: Folding Method Comparison
**File**: `05_folding_comparison.pdf/png`

**Purpose**: Direct comparison of Holm Fold vs Polyphest Fold completion rates

**Interpretation**: Which folding algorithm is more permissive for method completion?

#### Figure 33: Folding Accuracy Comparison
**File**: `06_folding_accuracy_comparison.pdf/png`

**Purpose**: Reticulation count errors for both folding methods

**Interpretation**: Which folding method provides better ground truth for comparison?

#### Figure 34: Reticulation Error Distribution (Percentage Bias)
**File**: `07_reticulation_bias_boxplot.pdf/png`

**Layout**: Box plots showing distribution of percentage bias: (Inferred - True) / True × 100

**Key Features**:
- **Percentage bias**: Signed error as percentage of true reticulation count
- **Positive values**: Over-estimation (inferred > true)
- **Negative values**: Under-estimation (inferred < true)
- **Can exceed ±100%**: If true=2 and inferred=4, bias = +100%; if inferred=0, bias = -100%
- **Mean percentage bias** shown above each box

**Interpretation**:
- Centered at 0% = accurate
- Positive = systematic overestimation
- Negative = systematic underestimation
- Box width shows consistency across networks
- Outliers indicate networks with extreme errors

**Example**: If true reticulations = 2:
- Inferred = 4 → bias = +100% (doubled)
- Inferred = 0 → bias = -100% (missed all)
- Inferred = 1 → bias = -50% (underestimated by half)

#### Figure 35: Edit Distance Distribution (MUL-tree)
**File**: `08_edit_distance_multree_boxplot.pdf/png`

**Layout**: Box plots showing edit distance distribution per method

**Interpretation**:
- Lower = better
- Box width shows consistency
- See EDIT_DISTANCE_EXPLAINED.md for scale interpretation

#### Figure 36: 3-Way Distance Metric Comparison
**File**: `08a_distance_metrics_comparison.pdf/png`

**Layout**: Three side-by-side boxplots comparing:
- Network Edit Distance
- MUL-tree Edit Distance
- RF Distance (MUL-tree)

**Purpose**: Compare different distance metrics to understand which aspects of topology are captured

#### Figure 37: Edit Distance MUL-tree Distribution (Detailed)
**File**: `08b_edit_distance_multree_distribution.pdf/png`

**Layout**: Boxplot distribution of MUL-tree edit distance

**Interpretation**: Detailed view of MUL-tree accuracy distribution

#### Figure 38: RF Distance Distribution
**File**: `08c_rf_distance_distribution.pdf/png`

**Layout**: Boxplot distribution of Robinson-Foulds distance

**Interpretation**: RF-based accuracy distribution across methods

#### Figure 39: Per-Network Completion Breakdown
**File**: `09_per_network_breakdown.pdf/png`

**Layout**: Bar plot with networks on x-axis, grouped by method

**Interpretation**:
- Identifies problematic networks (all methods fail)
- Identifies method-specific failures
- Shows network-by-network performance

#### Figure 40: Method Performance Summary
**File**: `10_method_summary.pdf/png`

**Layout**: Four-panel bar chart
1. **Completion rate**: Overall success rate
2. **Mean edit distance (MUL-tree)**: Topological accuracy
3. **Mean absolute error (MAE)**: Reticulation count absolute error
4. **Mean percentage bias**: Reticulation count signed error as percentage

**Purpose**: High-level overview for paper main text

**Key Features**:
- Panel 4 shows percentage bias with color coding:
  - **Red bars**: Over-estimation (positive bias)
  - **Blue bars**: Under-estimation (negative bias)
  - **Gray bars**: No data
- Reference line at 0% for perfect accuracy

#### Figure 41: Aggregated Correlation Heatmap
**File**: `31_comprehensive_correlation_heatmap.pdf/png`

**Layout**: Heatmap showing correlations between:
- **Rows**: All network properties (Num_Species, H_Strict, H_Relaxed, Num_Polyploids, Max_Copies, Total_WGD, Polyploid_Ratio)
- **Columns**: All performance metrics (completion_rate, edit_distance_multree, rf_distance, num_rets_bias_pct)

**Color scale**: Red = positive correlation, Blue = negative correlation

**Why important**:
- Identifies which network properties most strongly predict difficulty
- Shows relationships between properties and multiple performance aspects
- Reveals unexpected correlations
- **Aggregated across all methods** (see Figure 42 for per-method view)

**Interpretation**:
- Strong positive correlation with edit_distance = property makes accurate inference harder
- Strong negative correlation with completion_rate = property causes failures
- Example: If H_Strict strongly correlates with edit_distance, reticulations are a major accuracy challenge

#### Figure 42: Per-Method Correlation Heatmaps
**File**: `32_per_method_correlation_heatmap.pdf/png`

**Layout**: Grid of heatmaps, one subplot per method

**Purpose**: Shows correlation patterns **separately for each method**

**Why important**:
- Different methods may have different difficulty factors
- Reveals method-specific strengths and weaknesses
- Complements aggregated view (Figure 41)

**Interpretation**:
- Compare heatmaps across methods to see if difficulty factors are consistent
- Method-specific patterns may indicate algorithmic differences
- Strong correlations in one method but not others = method-specific challenge

---

## Tables

### Table 1: Method Performance Summary
**File**: `01_method_performance_summary.csv`

**Columns**:
- `Method`: Method name
- `Total_Runs`: Total network attempts
- `Completed_Runs`: Successful completions
- `Completion_Rate_%`: Success percentage
- `Mean_Edit_Distance`: Average accuracy (MUL-tree)
- `Median_Edit_Distance`: Median accuracy
- `Std_Edit_Distance`: Accuracy variability
- `Mean_Reticulation_Error`: Average reticulation count error (absolute)
- `Median_Reticulation_Error`: Median reticulation count error
- `Mean_Reticulation_Bias_%`: Average percentage bias (signed)

**Use**: Main results table for paper

### Table 2: Per-Network Performance
**File**: `02_per_network_performance.csv`

**Columns**:
- Network properties (H_Strict, H_Relaxed, Num_Polyploids, Total_WGD)
- Per-method completion rates
- Per-method edit distances
- Per-method reticulation errors and biases

**Use**: Supplementary material - detailed breakdown

## Understanding the Metrics

### Edit Distance
- **Range**: 0 (identical) to ~1.0 (completely different)
- **Calculation**: NetworkX graph edit distance, normalized by max(nodes + edges)
- **Why not Robinson-Foulds**: RF requires unique leaf labels (not MUL-trees) and doesn't handle reticulations
- **See**: `EDIT_DISTANCE_EXPLAINED.md` for details

### Jaccard Distance
- **Range**: 0 (perfect match) to 1 (no overlap)
- **Calculation**: Distance = 1 - Jaccard similarity
- **Jaccard similarity**: |True ∩ Inferred| / |True ∪ Inferred|
- **Variants**:
  - `ret_leaf_jaccard.dist`: Reticulation leaf set distance
  - `ret_sisters_jaccard.dist`: Sister relationship distance
- **Interpretation**: Lower values indicate better performance

### Reticulation Error Metrics

#### Absolute Error (`num_rets_diff`)
- **Definition**: |Inferred - True|
- **Range**: 0 to infinity
- **Interpretation**: Magnitude of error, regardless of direction

#### Percentage Bias (`num_rets_bias_pct`)
- **Definition**: (Inferred - True) / True × 100
- **Range**: Can exceed ±100%
- **Interpretation**:
  - **Positive**: Over-estimation (inferred > true)
  - **Negative**: Under-estimation (inferred < true)
  - **0%**: Perfect accuracy
- **Examples**:
  - True=2, Inferred=4 → +100% (doubled)
  - True=2, Inferred=0 → -100% (missed all)
  - True=10, Inferred=12 → +20% (slight over-estimation)
- **Edge case**: If True=0, percentage cannot be calculated (uses absolute error instead)

### F1 Score (Polyploid Identification)
- **Range**: 0 (worst) to 1 (perfect)
- **Calculation**: 2 × (Precision × Recall) / (Precision + Recall)
- **Based on**: True Positives, False Positives, False Negatives
- **Precision**: TP / (TP + FP) - accuracy of positive predictions
- **Recall**: TP / (TP + FN) - coverage of true positives

## Data Aggregation

### Completion Rate Calculation
```python
# For each characteristic value
networks_with_value = networks where characteristic == value
total_runs = sum of all runs across these networks
completed_runs = sum of completed runs
completion_rate = (completed_runs / total_runs) × 100
```

### Error Bars
Error bars represent standard error across networks:
```python
# Calculate per-network rates first
network_rates = [rate for each network with characteristic value]
mean_rate = mean(network_rates)
std_error = std(network_rates) / sqrt(number_of_networks)
```

### Percentage Bias Calculation
```python
# For each network-method pair
bias = inferred_reticulations - true_reticulations
if true_reticulations > 0:
    bias_pct = (bias / true_reticulations) × 100
else:
    bias_pct = bias  # Use absolute value if true=0
```

## Customization

### Color Schemes
Edit `create_analysis_figures.py`:
```python
METHOD_COLORS = {
    'grampa': '#0173B2',      # Blue
    'polyphest_p60': '#DE8F05',   # Orange
    'padre': '#ECE133',       # Yellow
    'mpsugar': '#CA9161'      # Tan
}
```

### Adding New Characteristics
To add analysis for a new network property:

1. Ensure column exists in `mul_tree_final_stats.csv`
2. Add completion plot calls in `generate_all_figures()`:
   ```python
   self.plot_completion_vs_characteristic_combined(
       'YourColumn', 'Your Label', 'prefix')
   ```
3. Add accuracy plot calls if desired
4. Update documentation

## Troubleshooting

### Error: "FileNotFoundError: inventory.csv"
**Solution**: Run summary pipeline:
```bash
python run_full_summary.py conf_ils_low_10M
```

### Error: "KeyError: 'inferred_exists'"
**Solution**: Re-run postprocessing:
```bash
python postprocess_results.py conf_ils_low_10M
python run_full_summary.py conf_ils_low_10M
```

### Missing Jaccard/F1 plots
**Problem**: Metrics not calculated in comparison pipeline

**Solution**: Ensure `compare_nets.py` calculates all metrics:
- `ret_leaf_jaccard` and `ret_leaf_jaccard.dist`
- `ret_sisters_jaccard` and `ret_sisters_jaccard.dist`
- `ploidy_diff` (TP/FP/FN)

### Plots look crowded
**Solution**: Use faceted versions in `individual_methods/` subdirectory

### Old plots still present
**Solution**: The script automatically cleans old plots before generation. If you see old files, they may be from a different configuration or the cleanup may have failed. Re-run the script.

## Workflow Checklist

- [ ] Simulations completed (low, medium, high)
- [ ] Postprocessing run for each configuration
- [ ] Summary pipeline run for each configuration
- [ ] Network statistics file exists
- [ ] Generated all 42 plots successfully
- [ ] Reviewed combined plots for overview
- [ ] Reviewed faceted plots for details
- [ ] Checked correlation heatmaps (aggregated and per-method) for insights
- [ ] Reviewed tables
- [ ] Selected figures for main text vs supplementary

## Figure Selection Guide for Publications

### Main Text (4-6 key figures)

**Essential**:
1. **Fig 40 (method summary)**: Comprehensive overview - completion, accuracy, reticulation error
2. **Fig 35 (edit distance distribution)**: Overall accuracy comparison
3. **Fig 34 (reticulation bias)**: Percentage bias shows systematic over/under-estimation
4. **Fig 41 (aggregated correlation)**: Identifies difficulty predictors

**Optional** (choose based on story):
5. **Fig 09 (combined)**: Completion vs Num_Species - if size matters
6. **Fig 13 (combined)**: Edit distance vs Num_Species - if accuracy degrades with size
7. **Fig 29 (polyploid F1)**: If polyploid identification is key
8. **Fig 42 (per-method correlation)**: If method-specific patterns are important

### Supplementary Material

**Category 1 - Completion Analysis**:
- All faceted completion plots (detailed per-method trends)
- All combined completion plots (Figs 01-12)

**Category 2 - Accuracy Analysis**:
- All edit distance vs characteristic plots (Figs 13-20)
- All RF distance plots (Figs 21-24)
- Shows how accuracy varies with network properties

**Category 3 - Advanced Metrics**:
- Jaccard distance plots (Figs 25-28)
- Jaccard distribution boxplots (Figs 30-31)
- Polyploid F1 if not in main text

**Category 4 - Distributions**:
- Folding comparisons (Figs 32-33)
- Error distributions (Figs 34-38)
- Per-network breakdown (Fig 39)

**Tables**:
- Table 1: Main supplementary table
- Table 2: Detailed per-network results

## Analysis Workflow

After generating all figures:

### 1. Overall Performance
- Review Fig 40 (method summary) - which methods perform best overall?
- Review Fig 35 (edit distance distribution) - how accurate are methods?
- Review Table 1 - quantitative comparison

### 2. Difficulty Factors
- Review Fig 41 (aggregated correlation) - which properties predict difficulty?
- Review Fig 42 (per-method correlation) - are difficulty factors consistent across methods?
- Review Category 1 combined plots - how does each property affect completion?
- Review Category 2 plots - how does each property affect accuracy?

### 3. Method-Specific Analysis
- Use faceted plots to see per-method trends clearly
- Identify method strengths/weaknesses from slopes
- Compare methods on specific network types (Fig 39)
- Compare per-method correlation patterns (Fig 42)

### 4. Advanced Insights
- Review Jaccard distance plots (Figs 25-28) - are topological relationships correct?
- Review Jaccard distributions (Figs 30-31) - consistency of performance
- Review F1 plot (Fig 29) - how well are polyploids identified?
- Review percentage bias (Fig 34) - systematic over/under-estimation patterns
- Cross-reference with biological expectations

### 5. ILS Level Comparison
- Generate plots for all three ILS levels
- Compare how ILS affects each metric
- Identify ILS-robust vs ILS-sensitive methods

### 6. Statistical Analysis (if needed)
- Wilcoxon signed-rank test for pairwise method comparison
- ANOVA for ILS level effects
- Regression for network property effects

## File Locations Reference

```
gene2net/
└── simulations/
    ├── networks/
    │   ├── *.tre                           # 21 input MUL-trees
    │   └── mul_tree_final_stats.csv        # Ground truth statistics
    ├── scripts/
    │   ├── create_analysis_figures.py      # Main visualization script (42 plots)
    │   ├── run_full_summary.py             # Prerequisite
    │   └── postprocess_results.py          # Prerequisite
    ├── analysis/
    │   └── summary/
    │       └── conf_ils_low_10M/
    │           ├── inventory.csv           # Input
    │           ├── aggregated_metrics.csv  # Input
    │           ├── comparisons_raw.csv     # Input
    │           ├── plots/                  # 22 combined plots (PDF + PNG)
    │           │   └── individual_methods/ # 20 faceted plots (PDF + PNG)
    │           └── tables/                 # 2 summary tables (CSV)
    └── md_files/
        ├── VISUALIZATION_GUIDE.md          # This file
        └── EDIT_DISTANCE_EXPLAINED.md      # Metric details
```

## Related Documentation

- **COMPLETE_PIPELINE_GUIDE.md**: Full simulation workflow
- **SUMMARY_PIPELINE_GUIDE.md**: Summary pipeline details
- **POSTPROCESSING_GUIDE.md**: Postprocessing details
- **EDIT_DISTANCE_EXPLAINED.md**: Edit distance metric explanation
- **CLAUDE.md**: Repository-wide guidance

## Quick Reference: All 42 Plots

### Category 1: Completion Rate (12 plots)
1-2. vs H_Strict (combined + faceted)
3-4. vs H_Relaxed (combined + faceted)
5-6. vs Num_Polyploids (combined + faceted)
7-8. vs Total_WGD (combined + faceted)
9-10. vs Num_Species (combined + faceted)
11-12. vs Max_Copies (combined + faceted)

### Category 2: Accuracy Metrics (12 plots)
13-14. Edit Distance vs Num_Species (combined + faceted)
15-16. Edit Distance vs H_Strict (combined + faceted)
17-18. Edit Distance vs Num_Polyploids (combined + faceted)
19-20. Edit Distance vs Max_Copies (combined + faceted)
21-22. RF Distance vs Num_Species (combined + faceted)
23-24. RF Distance vs H_Strict (combined + faceted)

### Category 3: Advanced Metrics (7 plots)
25-26. Reticulation Leaf Jaccard Distance vs H_Strict (combined + faceted)
27-28. Sister Relationship Jaccard Distance vs H_Strict (combined + faceted)
29. Polyploid F1 Score
30. Reticulation Leaf Jaccard Distance Distribution
31. Sister Relationship Jaccard Distance Distribution

### Category 4: Summaries (11 plots)
32. Folding comparison
33. Folding accuracy
34. Reticulation error distribution (percentage bias)
35. Edit distance distribution (MUL-tree)
36. 3-Way distance metric comparison
37. Edit distance MUL-tree distribution (detailed)
38. RF distance distribution
39. Per-network breakdown
40. Method summary
41. Aggregated correlation heatmap
42. Per-method correlation heatmaps

## Support

For issues or questions:
1. Check this guide
2. Review EDIT_DISTANCE_EXPLAINED.md for metric details
3. Check CLAUDE.md for repository context
4. Verify prerequisites completed
5. Review example outputs in `simulations/analysis/summary/`
