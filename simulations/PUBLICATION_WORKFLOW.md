# Publication Workflow Guide

Complete guide for generating publication-quality figures and tables from phylogenetic network inference results.

## Overview

This workflow takes you from raw simulation results to publication-ready figures and tables. It includes:

1. **Ground Truth Statistics**: Enhanced network characterization including autopolyploid detection
2. **Method Performance Analysis**: Success rates, edit distances, and rankings
3. **Publication Figures**: Comprehensive visualizations for papers/conferences
4. **Statistical Analysis**: Correlations between network properties and method performance

## Quick Start

```bash
# 1. Re-run reticulation stats with autopolyploid detection (on cluster)
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts
python run_reticulation_stats.py

# 2. Generate publication figures (locally or on cluster)
python create_publication_figures.py \
  --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
  --network-stats ../networks/mul_tree_final_stats.csv \
  --output ../analysis/publication_figures
```

## Step-by-Step Workflow

### Step 1: Update Network Ground Truth Statistics

**What**: Enhanced reticulation statistics that distinguish between autopolyploid and allopolyploid events.

**File**: `simulations/scripts/run_reticulation_stats.py`

**Key Features**:
- **Autopolyploids**: WGD events creating identical subtrees (e.g., Rch+Rch in Ding_2023)
- **Allopolyploids**: Hybridization events between different lineages (H_Strict)
- **Total WGD**: Sum of autopolyploid and allopolyploid events

**Run on cluster**:
```bash
ssh cluster
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts
python run_reticulation_stats.py
```

**Output**: `simulations/networks/mul_tree_final_stats.csv`

**Columns**:
- `Filename`: Network file name
- `Num_Species`: Total unique species
- `Num_Polyploids`: Species with >1 copy
- `Max_Copies`: Maximum copies of any species
- `H_Strict`: Reticulations (allopolyploid events)
- `H_Relaxed`: Reticulations with relaxed folding
- `H_Diff`: Difference between strict and relaxed
- `Num_Autopolyploids`: Autopolyploid WGD events
- `Num_Allopolyploids`: Allopolyploid events (= H_Strict)
- `Total_WGD`: Total whole genome duplications
- `Polyploid_Names`: List of polyploid species

**Expected Output Example** (Ding_2023.tre):
```
Num_Polyploids: 7
H_Strict: 0 (no hybridization)
Num_Autopolyploids: 1 (the duplicated Rch clade)
Total_WGD: 1
```

### Step 2: Verify Autopolyploid Detection (Optional)

**Test locally** (if Python/ete3 environment is compatible):
```bash
cd simulations/scripts
python test_autopolyploid_detection.py
```

This tests that Ding_2023.tre correctly identifies the duplicated Rch clade as an autopolyploid event.

### Step 3: Generate Publication Figures

**What**: Comprehensive publication-quality visualizations.

**File**: `simulations/scripts/create_publication_figures.py`

**Prerequisites**:
- Updated network stats (from Step 1)
- Completed summary pipeline for all ILS levels:
  ```bash
  python simulations/scripts/run_full_summary.py conf_ils_low_10M
  python simulations/scripts/run_full_summary.py conf_ils_medium_10M
  python simulations/scripts/run_full_summary.py conf_ils_high_10M
  ```

**Run**:
```bash
cd simulations/scripts

# Generate all figures for all ILS levels
python create_publication_figures.py \
  --configs conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
  --network-stats ../networks/mul_tree_final_stats.csv \
  --output ../analysis/publication_figures
```

**Generated Figures**:

1. **fig1_success_vs_reticulations.pdf/png**
   - Success rate vs number of reticulations
   - Separate panels for each ILS level (low, medium, high)
   - Shows which methods handle high-reticulation networks better

2. **fig2_success_vs_polyploids.pdf/png**
   - Success rate vs number of polyploid species
   - Identifies methods robust to polyploidy complexity

3. **fig3_success_vs_total_wgd.pdf/png**
   - Success rate vs total WGD events (autopolyploids + allopolyploids)
   - Combined view of whole genome duplication impact

4. **fig4_success_vs_autopolyploids.pdf/png**
   - Success rate vs number of autopolyploid events
   - Specifically examines WGD (non-hybridization) impact

5. **fig5_edit_distance_boxplot.pdf/png**
   - Distribution of edit distances for each method
   - Boxplots with individual data points overlay
   - Compares accuracy across methods and ILS levels

6. **fig6_correlation_heatmap_{config}.pdf/png** (one per ILS level)
   - Correlation matrix: network properties vs edit distance
   - Identifies which network characteristics make inference difficult
   - Separate heatmap for each ILS level

**Generated Tables**:

1. **table1_method_summary.csv** (and .tex for LaTeX)
   - Comprehensive summary table
   - Columns: ILS_Level, Method, Num_Networks, Num_Success, Success_Rate_%, Mean_Edit_Distance, Std_Edit_Distance
   - Ready for publication insertion

## Figure Descriptions for Paper

### Success Rate Figures (1-4)

**Purpose**: Demonstrate how network complexity affects method performance.

**Key Insights**:
- **Reticulations**: Shows which methods scale to high-reticulation networks
- **Polyploids**: Identifies methods robust to polyploidy
- **Total WGD**: Combined effect of all genome duplications
- **Autopolyploids**: Specific impact of WGD (non-hybrid) events

**Aggregation Strategy**:
When multiple networks have the same characteristic value (e.g., 5 reticulations), the success rate is calculated as:
```
success_rate = (total successful runs) / (total runs) * 100
```
across all networks with that value.

### Edit Distance Comparison (Figure 5)

**Purpose**: Compare method accuracy across ILS levels.

**Key Features**:
- Boxplot shows distribution (median, quartiles, outliers)
- Individual points overlaid for transparency
- Separate panel per ILS level for direct comparison

### Correlation Heatmap (Figure 6)

**Purpose**: Identify which network properties correlate with method difficulty.

**Interpretation**:
- **Positive correlation** (red): Higher values → Higher edit distance (worse performance)
- **Negative correlation** (green): Higher values → Lower edit distance (better performance)
- **No correlation** (white): Network property doesn't affect this method

## Data Aggregation Strategy

### For Success Rate Plots

When multiple networks share the same characteristic (e.g., 3 networks all have 5 reticulations):

**Approach**: Aggregate across all runs
```python
# For a given characteristic value and method:
num_success = count(successful runs across all matching networks)
num_total = count(total runs across all matching networks)
success_rate = num_success / num_total * 100
```

**Why**: Provides robust estimate of success probability for that characteristic value.

### For Edit Distance Plots

**Approach**: Show distribution across networks
- Boxplot displays median ± quartiles
- Individual points show per-network means
- Standard deviation shows variability within network replicates

## Additional Analysis Ideas

### Suggested Extensions

1. **Runtime Analysis**
   - Plot runtime vs network size/complexity
   - Identify computational bottlenecks

2. **Accuracy vs Runtime Trade-off**
   - Scatter plot: edit distance vs runtime
   - Pareto frontier of method performance

3. **Network Difficulty Classification**
   - Define "easy", "medium", "hard" based on success rates
   - Analyze what makes networks hard

4. **Method Failure Analysis**
   - Investigate specific networks where methods fail
   - Common failure patterns

5. **Statistical Tests**
   - Pairwise method comparisons (Wilcoxon signed-rank)
   - ANOVA for ILS level effects

## File Locations

```
gene2net/
├── simulations/
│   ├── networks/
│   │   ├── *.tre                                    # Input MUL-trees
│   │   └── mul_tree_final_stats.csv                 # Ground truth (output of Step 1)
│   ├── scripts/
│   │   ├── run_reticulation_stats.py                # Step 1: Enhanced stats
│   │   ├── test_autopolyploid_detection.py          # Optional testing
│   │   ├── run_full_summary.py                      # Prerequisites
│   │   └── create_publication_figures.py            # Step 3: Figure generation
│   └── analysis/
│       ├── summary/                                 # Summary pipeline outputs
│       │   ├── conf_ils_low_10M/
│       │   ├── conf_ils_medium_10M/
│       │   └── conf_ils_high_10M/
│       └── publication_figures/                     # Final outputs (Step 3)
│           ├── fig1_success_vs_reticulations.pdf
│           ├── fig2_success_vs_polyploids.pdf
│           ├── fig3_success_vs_total_wgd.pdf
│           ├── fig4_success_vs_autopolyploids.pdf
│           ├── fig5_edit_distance_boxplot.pdf
│           ├── fig6_correlation_heatmap_*.pdf
│           └── table1_method_summary.csv
```

## Troubleshooting

### "No valid data found" when running run_reticulation_stats.py

**Problem**: Network files not found or inaccessible.

**Solution**:
- Verify path in script: `DIR_PATH = "/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks/"`
- Check file extensions: Should be `.tre`, `.nwk`, or `.tree`

### "Summary directory not found" when generating figures

**Problem**: Summary pipeline not run for specified configs.

**Solution**:
```bash
# Run summary pipeline first for each config
python simulations/scripts/run_full_summary.py conf_ils_low_10M
python simulations/scripts/run_full_summary.py conf_ils_medium_10M
python simulations/scripts/run_full_summary.py conf_ils_high_10M
```

### Empty plots or missing data points

**Problem**: Network stats not merged correctly.

**Solution**:
- Ensure network names match between:
  - `mul_tree_final_stats.csv` (network column)
  - `inventory.csv` (network column)
- Names should be without `.tre` extension

### Python environment issues (ete3, matplotlib)

**Problem**: Missing dependencies or Python version incompatibility.

**Solution**:
```bash
# On cluster, activate gene2net environment
conda activate gene2net

# If missing packages
pip install pandas numpy matplotlib seaborn scipy
```

## Next Steps After Figure Generation

1. **Review figures** in `simulations/analysis/publication_figures/`

2. **Select key figures** for paper:
   - Typically 3-4 main figures
   - Rest go to supplementary material

3. **Customize as needed**:
   - Edit Python script color schemes, sizes
   - Adjust labels for journal style

4. **Extract insights**:
   - Which methods perform best overall?
   - How does ILS level affect performance?
   - What network properties correlate with difficulty?
   - Are autopolyploids harder to infer than allopolyploids?

5. **Statistical tests** (if needed):
   - Run additional analyses based on figures
   - Add significance markers to plots

## Key Concepts

### Autopolyploid vs Allopolyploid

**Autopolyploid (WGD)**:
- Genome duplication within a species (A + A → AA)
- Creates identical subtrees in MUL-tree
- Example: Ding_2023.tre has duplicated Rch clade

**Allopolyploid (Hybridization)**:
- Genome combination between species (A + B → AB)
- Creates reticulations in network (H_Strict)
- Different parent lineages merge

**Total WGD**:
```
Total_WGD = Num_Autopolyploids + Num_Allopolyploids
```

### Success Rate Calculation

For a given network characteristic value and method:
```
Success Rate = (# successful runs) / (# total runs) × 100
```

Where:
- **Successful run**: Method produced output file
- **Total runs**: All attempts for that characteristic × method combination

### Edit Distance (Normalized)

Measures similarity between inferred and true network:
```
Edit Distance = graph_edit_distance(inferred, true) / normalization_factor
```

Where:
- `graph_edit_distance`: Minimum edits (add/remove nodes/edges) to transform one network to another
- Normalization: Typically max(|nodes| + |edges|) of the two networks
- **Range**: 0 (identical) to 1+ (very different)

## Citation

If using this workflow, cite:
- **GRAMPA**: [Upcoming publication]
- **Polyphest**: Ogilvie et al. (2016)
- **MPSUGAR**: Heydari et al. (2020)
- **PADRE**: [Publication details]
- **ReticulateTree**: [This repository]

## Contact

For questions about this workflow:
- Check existing documentation in `simulations/` directory
- Review `CLAUDE.md` for general repository guidance
- Examine example outputs in `simulations/analysis/`
