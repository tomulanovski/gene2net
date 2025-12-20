# Summary Pipeline Usage Guide

This pipeline evaluates phylogenetic network reconstruction methods by comparing inferred networks to ground truth across multiple configurations, methods, and replicates.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Expected Input Files](#expected-input-files)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Metrics Explained](#metrics-explained)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Post-Processing Required

**IMPORTANT**: Before running the summary pipeline, you must post-process raw program outputs to extract clean MUL-trees.

```bash
# Step 1: Post-process raw outputs (REQUIRED)
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# Step 2: Run summary pipeline
python simulations/scripts/run_full_summary.py conf_ils_low_10M
```

**Why?** Each program outputs results in different formats:
- **Polyphest**: Text file with "multree:" and "network:" lines
- **GRAMPA**: TSV with multiple ranked trees (need best one)
- **MPSUGAR**: Detailed report with embedded Newick tree
- **PADRE**: Already clean (just needs copying)

Post-processing extracts/copies these to standardized `{program}_result.tre` files.

**See**: `simulations/md_files/POSTPROCESSING_GUIDE.md` for detailed instructions

---

## Expected Input Files

### Required Directory Structure

```
/groups/itay_mayrose/tomulanovski/gene2net/simulations/
├── networks/                                    # Ground truth networks (21 networks)
│   ├── Bendiksby_2011.tre
│   ├── Koenen_2020.tre
│   ├── Brysting_2007.tre
│   ├── ... (18 more networks)
│   └── mul_tree_final_stats.csv                 # Network characteristics (optional for correlations)
│
└── {NetworkName}/                               # For each of 21 networks
    └── results/
        └── {CONFIG}/                            # e.g., conf_ils_low_10M
            ├── grampa/
            │   ├── replicate_1/
            │   │   └── grampa_result.tre    # Post-processed result
            │   │     
            │   ├── replicate_2/
            │   │   └ grampa_result.tre
            │   ├── ... (replicates 3-5)
            │
            ├── polyphest_p50/
            │   ├── replicate_1/
            │   │   └── polyphest_result.tre
            │   ├── replicate_2/
            │   │   └── polyphest_result.tre
            │   ├── ... (replicates 3-5)
            │
            ├── polyphest_p70/
            │   └── ... (5 replicates with polyphest_result.tre)
            │
            ├── polyphest_p90/
            │   └── ... (5 replicates with polyphest_result.tre)
            │
            ├── mpsugar/
            │   ├── replicate_1/
            │   │   └── mpsugar_result.tre
            │   └── ... (replicates 2-5)
            │
            └── padre/
                ├── replicate_1/
                │   └── padre_result.tre
                └── ... (replicates 2-5)
```

### Ground Truth Networks

**Location**: `/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks/`

**Format**: Newick format (`.tre` files)

**21 Networks**:
- Bendiksby_2011, Koenen_2020, Brysting_2007, Lawrence_2016
- Diaz-Perez_2018, Wisecaver_2023, Ding_2023, Liang_2019
- Popp_2005, Wu_2015, Liu_2023, Ren_2024
- Marcussen_2011, Marcussen_2012, Sessa_2012b, Zhao_2021
- Hori_2014, Marcussen_2015, Shahrestani_2015
- Morales-Briones_2021, Soza_2014

### Inferred Network Files

**Naming Convention**: `{program}_result.tre` (standardized post-processed MUL-trees)

**Required Files per Configuration** (created by `postprocess_results.py`):
- GRAMPA: `grampa/replicate_{1-5}/grampa_result.tre`
- Polyphest (p50): `polyphest_p50/replicate_{1-5}/polyphest_result.tre`
- Polyphest (p70): `polyphest_p70/replicate_{1-5}/polyphest_result.tre`
- Polyphest (p90): `polyphest_p90/replicate_{1-5}/polyphest_result.tre`
- MPSUGAR: `mpsugar/replicate_{1-5}/mpsugar_result.tre`
- PADRE: `padre/replicate_{1-5}/padre_result.tre`

**Total Files per Configuration**: 21 networks × 6 methods × 5 replicates = **630 files**

**Note**: These files are created by post-processing raw program outputs (see [Prerequisites](#prerequisites))

### Network Characteristics File (Optional)

**Location**: `simulations/simulations/networks/mul_tree_final_stats.csv`

**Purpose**: Used for Level 3 correlation analysis

**Required Columns**:
- `Filename` - Network filename (e.g., "Bendiksby_2011.tre")
- `Num_Species` - Number of species in network
- `Num_Polyploids` - Number of polyploid species
- `Max_Copies` - Maximum copy number
- `H_Strict` - Reticulation count (strict conversion)

**If missing**: Pipeline skips Level 3 (correlation analysis) but completes other levels

---

## Quick Start

```bash
# Navigate to gene2net root directory
cd /groups/itay_mayrose/tomulanovski/gene2net

# STEP 1: Post-process raw program outputs (REQUIRED FIRST)
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# STEP 2: Run summary pipeline
python simulations/scripts/run_full_summary.py conf_ils_low_10M
```

**Expected Runtime**:
- Post-processing: 1-2 minutes
- Summary pipeline: 5-10 minutes (first run), <30 seconds (with cache)

**Expected Output**: ~8 files + Level 1 directory with 15 CSV files

---

## Pipeline Overview

The pipeline consists of 4 modules that run sequentially:

### 1. Data Collection (`collect_results.py`)
- Scans directory structure for all expected files
- Creates inventory of available data
- Reports completion rates per method/network

### 2. Comparison Engine (`compute_comparisons.py`)
- Loads ground truth and inferred networks using `ReticulateTree`
- Computes 5 metrics using `pairwise_compare()` from `compare_reticulations.py`
- Caches results with file hash validation
- Handles missing files gracefully

### 3. Aggregation (`aggregate_and_summarize.py`)
- Aggregates across 5 replicates (mean ± std)
- Generates 3 summary levels
- Handles partial data (< 5 replicates)

### 4. Orchestrator (`run_full_summary.py`)
- Runs all 3 modules with one command
- Manages output directories
- Provides progress reporting

---

## Usage Examples

### Basic Usage

#### Process Single Configuration

```bash
# Process conf_ils_low_10M (all methods, all networks)
python simulations/scripts/run_full_summary.py conf_ils_low_10M

# Expected output:
#   - 630 comparisons (21 networks × 6 methods × 5 replicates)
#   - 3 summary levels
#   - Output in: simulations/analysis/summary/conf_ils_low_10M/
```

#### Process All Configurations

```bash
# Run each configuration separately
python simulations/scripts/run_full_summary.py conf_ils_low_10M
python simulations/scripts/run_full_summary.py conf_ils_medium_10M
python simulations/scripts/run_full_summary.py conf_ils_high_10M

# Or use a loop
for config in conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M; do
    echo "Processing $config..."
    python simulations/scripts/run_full_summary.py $config
done
```

### Method Filtering

#### Process Specific Methods

```bash
# Only GRAMPA
python simulations/scripts/run_full_summary.py conf_ils_low_10M --methods grampa

# GRAMPA and PADRE
python simulations/scripts/run_full_summary.py conf_ils_low_10M --methods grampa padre

# All Polyphest variants
python simulations/scripts/run_full_summary.py conf_ils_low_10M --methods polyphest_p50 polyphest_p70 polyphest_p90

# GRAMPA, MPSUGAR, and one Polyphest
python simulations/scripts/run_full_summary.py conf_ils_low_10M --methods grampa mpsugar polyphest_p50
```

**Available Methods**:
- `grampa`
- `polyphest_p50`
- `polyphest_p70`
- `polyphest_p90`
- `mpsugar`
- `padre`

### Configuration Options

#### Custom Configuration File

```bash
# Use custom config file (default: simulations/summary_config.yaml)
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --config /path/to/custom_config.yaml
```

#### Custom Output Directory

```bash
# Change output location (default: from config file)
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --output /custom/output/path/

# Results will be in: /custom/output/path/conf_ils_low_10M/
```

#### Include Network Characteristics

```bash
# Specify network stats file for correlation analysis
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --network-stats simulations/simulations/networks/mul_tree_final_stats.csv

# If not specified, pipeline looks for default location automatically
```

### Cache Management

#### Force Recompute

```bash
# Ignore cache and recompute all comparisons
python simulations/scripts/run_full_summary.py conf_ils_low_10M --force-recompute

# Use when:
#   - Result files were updated
#   - You suspect cached comparisons are incorrect
#   - Debugging comparison issues
```

#### Using Cache (Default)

```bash
# Automatically uses cache if files haven't changed
python simulations/scripts/run_full_summary.py conf_ils_low_10M

# Cache validation:
#   - Checks file hashes (SHA256)
#   - Recomputes if ground truth or inferred network changed
#   - Skips recomputation if files unchanged
```

### Dry Run

```bash
# Preview what would be processed without computing
python simulations/scripts/run_full_summary.py conf_ils_low_10M --dry-run

# Shows:
#   - Which files exist
#   - Completion rates
#   - What would be processed
#   - Exits before running comparisons
```

### Complete Examples

#### Example 1: Quick Test on Subset

```bash
# Test pipeline on just GRAMPA and PADRE for one config
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --methods grampa padre \
    --output test_results/

# Expected: ~210 comparisons (21 networks × 2 methods × 5 replicates)
```

#### Example 2: Full Analysis with All Options

```bash
# Complete analysis with custom paths and stats
python simulations/scripts/run_full_summary.py conf_ils_medium_10M \
    --config simulations/summary_config.yaml \
    --network-stats simulations/simulations/networks/mul_tree_final_stats.csv \
    --output simulations/analysis/summary/ \
    --force-recompute

# This will:
#   1. Process all 6 methods
#   2. Use custom config
#   3. Include correlation analysis
#   4. Force recompute all comparisons
#   5. Output to specified directory
```

#### Example 3: Compare Different Polyphest Percentiles

```bash
# Process only Polyphest variants to compare percentiles
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --methods polyphest_p50 polyphest_p70 polyphest_p90

# Then check level2_method_rankings.csv to see which percentile performs best
```

#### Example 4: Process Multiple Configs in Parallel

```bash
# Terminal 1
python simulations/scripts/run_full_summary.py conf_ils_low_10M &

# Terminal 2
python simulations/scripts/run_full_summary.py conf_ils_medium_10M &

# Terminal 3
python simulations/scripts/run_full_summary.py conf_ils_high_10M &

# Wait for all to complete
wait

echo "All configurations processed!"
```

#### Example 5: Check Data Before Running

```bash
# First, check what raw outputs are available
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run

# Post-process available outputs
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# Review completion rates, then decide which methods to process
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --methods grampa padre  # Only process complete methods
```

---

## Output Files

### Output Directory Structure

```
simulations/analysis/summary/
└── {CONFIG}/                            # e.g., conf_ils_low_10M
    ├── cache/                           # Cached comparisons (auto-managed)
    │   ├── Bendiksby_2011_grampa_rep1.pkl
    │   ├── Bendiksby_2011_grampa_rep2.pkl
    │   └── ... (630 cache files)
    │
    ├── inventory.csv                    # Full inventory of expected files
    ├── comparison_report.txt            # Detailed success/failure report
    ├── comparisons_raw.csv              # All replicate comparisons
    ├── aggregated_metrics.csv           # Mean ± std across replicates
    │
    ├── level1_detailed_per_network/     # Per-network pivot tables
    │   ├── edit_distance.csv            # Combined (mean ± std)
    │   ├── edit_distance_mean.csv       # Mean only
    │   ├── edit_distance_std.csv        # Std only
    │   ├── edit_distance_n.csv          # Sample sizes
    │   ├── num_rets_diff.csv
    │   ├── num_rets_diff_mean.csv
    │   ├── num_rets_diff_std.csv
    │   ├── num_rets_diff_n.csv
    │   ├── ploidy_diff.csv
    │   ├── ploidy_diff_mean.csv
    │   ├── ploidy_diff_std.csv
    │   ├── ploidy_diff_n.csv
    │   ├── ret_leaf_jaccard.csv
    │   ├── ret_leaf_jaccard_mean.csv
    │   ├── ret_leaf_jaccard_std.csv
    │   └── ret_leaf_jaccard_n.csv
    │   └── ret_sisters_jaccard.csv
    │   └── ... (3 more _mean/_std/_n files)
    │
    ├── level2_method_rankings.csv       # Overall method rankings
    └── level3_network_correlations.csv  # Network property correlations
```

### Key Output Files Explained

#### inventory.csv
**Purpose**: Lists all expected file combinations and availability

**Columns**:
- `network` - Network name (e.g., "Bendiksby_2011")
- `config` - Configuration name (e.g., "conf_ils_low_10M")
- `method` - Method name (e.g., "grampa")
- `replicate` - Replicate number (1-5)
- `gt_path` - Path to ground truth network
- `inferred_path` - Path to inferred network
- `gt_exists` - Boolean: ground truth file exists
- `inferred_exists` - Boolean: inferred file exists and non-empty
- `file_size` - Size of inferred file in bytes

**Use**: Check data completeness before running comparisons

#### comparisons_raw.csv
**Purpose**: All individual replicate comparisons

**Columns**:
- `network`, `config`, `method`, `replicate` - Identifiers
- `metric` - Metric name (e.g., "edit_distance", "ploidy_diff.dist")
- `value` - Metric value
- `status` - "SUCCESS" or "FAILED"

**Rows**: One per (network × method × replicate × metric)
**Example Count**: 21 networks × 6 methods × 5 reps × 13 metrics = ~8,190 rows

#### aggregated_metrics.csv
**Purpose**: Statistics across 5 replicates

**Columns**:
- `network`, `config`, `method`, `metric` - Identifiers
- `mean` - Mean across replicates
- `std` - Standard deviation
- `min` - Minimum value
- `max` - Maximum value
- `n_valid` - Number of successful replicates (1-5)

**Use**: Quick overview of method performance with uncertainty

**Example Row**:
```
Bendiksby_2011, conf_ils_low_10M, grampa, edit_distance, 0.234, 0.045, 0.189, 0.298, 5
```
Interpretation: GRAMPA had mean edit distance 0.234 ± 0.045 on Bendiksby_2011 (all 5 replicates succeeded)

#### level1_detailed_per_network/*.csv
**Purpose**: Detailed comparison tables for each metric

**Format**: Pivot tables
- **Rows**: 21 networks
- **Columns**: 6 methods
- **Values**: mean ± std (n_valid)

**Files**:
1. `edit_distance.csv` - Overall network structure similarity
2. `num_rets_diff.csv` - Reticulation count accuracy
3. `ploidy_diff.csv` - Polyploid detection accuracy
4. `ret_leaf_jaccard.csv` - Reticulation leaf set accuracy
5. `ret_sisters_jaccard.csv` - Reticulation topology accuracy

**Plus**: Separate `*_mean.csv`, `*_std.csv`, `*_n.csv` for each metric

**Use**: Compare methods per network, identify which networks are hard/easy

**Example Cell**: `0.2340 ± 0.0450 (5)` = mean 0.234, std 0.045, 5 replicates

#### level2_method_rankings.csv
**Purpose**: Overall method performance rankings

**Columns**:
- `method` - Method name
- `metric` - Metric name
- `avg_value` - Average metric value across all networks
- `avg_rank` - Average rank (1=best, lower is better)
- `num_best` - How many networks this method ranked #1
- `num_worst` - How many networks this method ranked last
- `num_networks` - Number of networks with data

**Use**: Identify best overall methods, compare method consistency

**Example**: If GRAMPA has avg_rank=2.3 for edit_distance, it typically ranks 2nd-3rd best

#### level3_network_correlations.csv
**Purpose**: How network characteristics affect reconstruction difficulty

**Columns**:
- `method` - Method name
- `metric` - Metric name (edit_distance, num_rets_diff)
- `network_property` - Property name (num_species, num_reticulations, etc.)
- `correlation` - Pearson correlation coefficient (-1 to +1)
- `p_value` - Statistical significance
- `significant` - Boolean: p < 0.05
- `n_networks` - Number of networks in analysis

**Use**: Understand which properties make networks harder for each method

**Example**: correlation=0.75, p=0.001 for (grampa, edit_distance, num_reticulations)
→ GRAMPA struggles more with networks that have more reticulations

#### comparison_report.txt
**Purpose**: Detailed success/failure breakdown

**Contents**:
- Total comparisons attempted
- Success/failure counts and percentages
- Cache usage statistics
- Detailed list of all failed comparisons with error messages

**Use**: Debugging, identifying problematic networks/methods

---

## Metrics Explained

### Primary Metrics (Lower is Better)

#### edit_distance
- **Range**: 0-1 (normalized)
- **Meaning**: Graph edit distance measuring overall network structure similarity
- **0** = Perfect match (identical networks)
- **1** = Maximum difference
- **Use**: Overall reconstruction quality

#### num_rets_diff
- **Range**: 0-infinity (integer)
- **Meaning**: Absolute difference in reticulation count: |H_truth - H_inferred|
- **0** = Exact reticulation count match
- **Positive** = Method over/under-estimated reticulations
- **Use**: Did method detect right number of reticulations?

#### ploidy_diff.dist
- **Range**: 0-1 (Jaccard distance)
- **Meaning**: How well polyploid species were identified
- **0** = All polyploid copies correctly identified
- **1** = No polyploid copies correctly identified
- **Use**: Polyploidy detection accuracy

#### ret_leaf_jaccard.dist
- **Range**: 0-1 (Jaccard distance)
- **Meaning**: Accuracy of reticulation leaf sets (which species involved)
- **0** = Perfect match (right species in each reticulation)
- **1** = No overlap
- **Use**: Did method identify correct species as reticulated?

#### ret_sisters_jaccard.dist
- **Range**: 0-1 (Jaccard distance)
- **Meaning**: Reticulation topology accuracy including sister clades
- **0** = Perfect topology match
- **1** = No topology match
- **Use**: Most stringent topology comparison

### Secondary Metrics (Per-Metric Breakdown)

For each of `ploidy_diff`, `ret_leaf_jaccard`, `ret_sisters_jaccard`:

#### .FP (False Positive Rate)
- **Range**: 0-1
- **Meaning**: Proportion of incorrectly identified features
- **0** = No false positives (all detections correct)
- **1** = All detections were wrong

#### .FN (False Negative Rate)
- **Range**: 0-1
- **Meaning**: Proportion of missed features
- **0** = No false negatives (found everything)
- **1** = Missed everything

#### .TP (True Positive Rate)
- **Range**: 0-1
- **Meaning**: Proportion correctly identified (approximate as 1 - dist)
- **1** = Perfect detection
- **0** = Nothing correct

---

## Advanced Usage

### Run Individual Modules

For debugging or customization:

```bash
# 1. Data collection only
python simulations/scripts/collect_results.py conf_ils_low_10M \
    --export inventory.csv

# 2. Comparisons only (requires inventory)
python simulations/scripts/compute_comparisons.py inventory.csv cache/ \
    --export comparisons.csv \
    --report report.txt

# 3. Aggregation only (requires comparisons)
python simulations/scripts/aggregate_and_summarize.py comparisons.csv output/ \
    --network-stats simulations/simulations/networks/mul_tree_final_stats.csv \
    --export-aggregated aggregated.csv
```

### Modify Configuration File

Edit `simulations/summary_config.yaml` to customize:

```yaml
# Add new method
methods:
  my_new_method:
    directory: "my_method"
    output_file: "my_method_result.tre"

# Change paths
base_dir: "/custom/path/to/simulations"

# Change number of replicates
num_replicates: 10  # if you have 10 replicates instead of 5
```

### Custom Analysis Scripts

Use outputs programmatically:

```python
import pandas as pd

# Load aggregated metrics
agg = pd.read_csv('summary/conf_ils_low_10M/aggregated_metrics.csv')

# Find best method per network for edit_distance
best_methods = (agg[agg['metric'] == 'edit_distance']
                .loc[agg.groupby('network')['mean'].idxmin()]
                [['network', 'method', 'mean']])

print(best_methods)
```

---

## Troubleshooting

### Common Issues

#### "Configuration file not found"

**Error**: `FileNotFoundError: simulations/summary_config.yaml`

**Solution**:
```bash
# Make sure you're in gene2net root
cd /groups/itay_mayrose/tomulanovski/gene2net

# Or specify full path
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --config /full/path/to/summary_config.yaml
```

#### Low Completion Rate

**Symptom**: `comparison_report.txt` shows many missing files

**Diagnosis**:
```bash
# Check what files exist
python simulations/scripts/collect_results.py conf_ils_low_10M
```

**Solutions**:
- Check if method jobs finished running (check SLURM logs)
- Verify post-processing created `{program}_result.tre` files
- Check file naming matches config exactly
- Re-run incomplete method jobs

#### "Failed to load network" Errors

**Error**: `ReticulateTree failed to load: {path}`

**Possible Causes**:
1. File doesn't exist or is empty
2. Invalid Newick format
3. Parsing errors
4. **Post-processing not run** (most common!)

**Debug**:
```bash
# Check if post-processed files exist
ls -lh results/conf_ils_low_10M/grampa/replicate_1/*_result.tre

# If missing, run post-processing first
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# Check file exists and has content
ls -lh /path/to/network/file.tre

# Try loading manually in Python
python
>>> from reticulate_tree import ReticulateTree
>>> tree = ReticulateTree("/path/to/file.tre")
```

#### High Variance in Replicates

**Symptom**: Large `std` values in aggregated_metrics.csv

**Interpretation**: Method has high variability across replicates

**Check**:
- Is this expected for stochastic methods?
- Are all replicates from same parameter settings?
- Check individual replicate values in comparisons_raw.csv

#### Missing Level 3 (Correlations)

**Symptom**: `level3_network_correlations.csv` not generated

**Cause**: Network stats file not found or not specified

**Solution**:
```bash
# Specify network stats file explicitly
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --network-stats simulations/simulations/networks/mul_tree_final_stats.csv
```

### Performance Issues

#### Slow First Run

**Expected**: 5-10 minutes for 630 comparisons

**If slower**:
- Check disk I/O (cluster storage)
- Verify network files aren't corrupted (causing re-reads)
- Check cache directory is writable

#### Slow Subsequent Runs

**Expected**: <30 seconds with cache

**If slower**:
- Cache might be invalidated (files changed)
- Check `comparison_report.txt` for "Newly computed" count
- Use `--force-recompute` intentionally invalidates cache

---

## Example Workflow

### Complete Analysis Workflow

```bash
# 1. Navigate to project
cd /groups/itay_mayrose/tomulanovski/gene2net

# 2. Check which methods completed
echo "Checking raw outputs for conf_ils_low_10M..."
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run

# 3. Post-process raw outputs to extract clean MUL-trees
echo "Post-processing outputs..."
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# 4. Run summary pipeline
echo "Running summary pipeline..."
python simulations/scripts/run_full_summary.py conf_ils_low_10M

# 4. Quick check of results
cd simulations/analysis/summary/conf_ils_low_10M
echo "Files created:"
ls -lh

# 5. View method rankings
echo "Method rankings:"
column -t -s, level2_method_rankings.csv | head -20

# 6. Copy results for local analysis
scp -r conf_ils_low_10M/ local_machine:/path/to/analysis/

# 7. Repeat for other configurations
cd /groups/itay_mayrose/tomulanovski/gene2net
for config in conf_ils_medium_10M conf_ils_high_10M; do
    echo "Processing $config..."
    python simulations/scripts/run_full_summary.py $config
done
```

### Quick Comparison Workflow

```bash
# Compare just two methods quickly
python simulations/scripts/run_full_summary.py conf_ils_low_10M \
    --methods grampa padre \
    --output quick_test/

# View results
cd quick_test/conf_ils_low_10M
cat level2_method_rankings.csv
```

---

## Getting Help

### Built-in Help

```bash
# Main orchestrator help
python simulations/scripts/run_full_summary.py --help

# Individual module help
python simulations/scripts/collect_results.py --help
python simulations/scripts/compute_comparisons.py --help
python simulations/scripts/aggregate_and_summarize.py --help
```

### Debugging Tips

1. **Start small**: Test on subset with `--methods grampa`
2. **Use dry-run**: Preview with `--dry-run` first
3. **Check logs**: Read `comparison_report.txt` for errors
4. **Verify inputs**: Run `collect_results.py` standalone
5. **Check cache**: Delete cache directory to force fresh start

### Additional Resources

- **Post-processing guide**: `simulations/md_files/POSTPROCESSING_GUIDE.md`
- **Implementation plan**: `.claude/plans/toasty-churning-harbor.md`
- **Pipeline checker**: `simulations/scripts/check_pipeline_status.py`
- **Main documentation**: `CLAUDE.md`
