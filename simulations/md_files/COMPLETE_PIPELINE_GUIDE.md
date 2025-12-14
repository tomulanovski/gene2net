# Complete Pipeline Guide: SimPhy to Summary

This guide covers the **entire workflow** for phylogenetic network inference simulations, from generating gene trees through final analysis.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Step 1: Generate Simulations (SimPhy)](#step-1-generate-simulations-simphy)
3. [Step 2: Run Network Inference Methods](#step-2-run-network-inference-methods)
4. [Step 3: Post-Process Results](#step-3-post-process-results)
5. [Step 4: Analyze and Compare](#step-4-analyze-and-compare)
6. [Optional: Additional Analyses](#optional-additional-analyses)
7. [Complete Examples](#complete-examples)
8. [Troubleshooting](#troubleshooting)

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│ COMPLETE SIMULATION PIPELINE                                        │
└─────────────────────────────────────────────────────────────────────┘

STEP 1: GENERATE SIMULATIONS
   │
   ├─> submit_simphy.sh
   │       ↓
   │   21 networks × 5 replicates × 1000 gene trees
   │       ↓
   └─> check_pipeline_status.py --step simphy
           ✓ Validate 105,000 gene trees generated

STEP 2: RUN NETWORK INFERENCE METHODS
   │
   ├─> submit_all_methods.sh (RECOMMENDED)
   │   OR submit individual methods
   │       ↓
   │   Prep: Format trees for each method (21 × 5 = 105 prep jobs)
   │       ↓
   │   Run: Infer networks (21 × 5 = 105 inference jobs per method)
   │       ↓
   └─> check_pipeline_status.py --step prep/run
           ✓ Validate inputs prepared and outputs generated

STEP 3: POST-PROCESS RESULTS
   │
   └─> postprocess_results.py
           ↓
       Extract clean MUL-trees from method-specific formats
           ✓ 630 standardized result files (21 × 6 methods × 5 reps)

STEP 4: ANALYZE AND COMPARE
   │
   └─> run_full_summary.py
           ↓
       Collect → Compare → Aggregate → Summarize
           ✓ Performance metrics, rankings, correlations

OPTIONAL: ADDITIONAL ANALYSES
   │
   ├─> submit_rf.sh (Robinson-Foulds distances)
   └─> submit_dup_loss_extraction.sh (Duplication/loss events)
```

---

## Step 1: Generate Simulations (SimPhy)

### What This Does

SimPhy generates gene trees under different evolutionary scenarios:
- **ILS (Incomplete Lineage Sorting):** Different effective population sizes
- **Dup/Loss:** Different duplication and loss rates
- Each simulation: 1000 gene trees × 5 replicates × 21 networks

### Running SimPhy

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# ILS scenarios
./submit_simphy.sh conf_ils_low_10M       # Ne = 200,000 (low ILS)
./submit_simphy.sh conf_ils_medium_10M    # Ne = 1,000,000 (medium ILS)
./submit_simphy.sh conf_ils_high_10M      # Ne = 2,000,000 (high ILS)

# Duplication/loss scenarios
./submit_simphy.sh conf_dup_loss_low_10M      # dup=0.001, loss=0.001
./submit_simphy.sh conf_dup_loss_medium_10M   # dup=0.005, loss=0.005
./submit_simphy.sh conf_dup_loss_high_10M     # dup=0.01, loss=0.01
```

**Dry-run option available:** Add `--dry-run` to preview

### Monitor SimPhy Jobs

```bash
# Check job status
squeue -u $USER | grep simphy

# Monitor logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_*.out

# Check for errors
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/simphy_*.err
```

### Validate SimPhy Completion

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Validate simulations completed
python check_pipeline_status.py conf_ils_low_10M --step simphy

# Expected output: "Success: 105 / 105 (100.0%)"
# This means: 21 networks × 5 replicates = 105 successful simulations
```

**What validation checks:**
- Each replicate has exactly 1000 gene trees
- Handles both single-batch and multi-batch SimPhy outputs
- Reports missing or incomplete simulations

### Expected Outputs

```
simulations/simulations/{NETWORK}/data/{CONFIG}/
└── replicate_{1-5}/
    ├── 1/                      # Single batch mode
    │   ├── g_trees.trees       # 1000 gene trees
    │   ├── s_tree.trees        # Species tree
    │   └── l_trees.trees       # Locus trees
    OR
    ├── batch_1/                # Multi-batch mode
    │   └── 1/
    │       └── g_trees.trees   # Subset of trees
    ├── batch_2/
    │   └── 1/
    │       └── g_trees.trees   # More trees
    └── ...                     # Until 1000 total trees
```

### Troubleshooting SimPhy

**Problem:** Simulations timeout

**Cause:** SimPhy may retry with multi-batch mode automatically

**Solution:** Check logs - multi-batch is normal fallback behavior

**Problem:** Some replicates show < 1000 trees

**Solution:**
```bash
# Re-run specific networks by array index
# Network indices are 1-21 for the 21 networks
./submit_simphy.sh conf_ils_low_10M 10M 0 0 200000 "13"  # Just network 13
```

---

## Step 2: Run Network Inference Methods

### What This Does

Runs phylogenetic network inference on simulated gene trees using:
- **GRAMPA:** Species tree + reticulation graph
- **Polyphest:** Maximum parsimony network
- **PADRE:** Parsimony-based network
- **MPSUGAR:** Bayesian network inference

### Option A: Master Script (RECOMMENDED)

Run all methods at once for one or more configurations:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Single configuration
./submit_all_methods.sh conf_ils_low_10M

# Multiple configurations
./submit_all_methods.sh conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

# Preview first (dry-run)
./submit_all_methods.sh conf_ils_low_10M --dry-run

# Specific methods only
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest
```

**Advantages:**
- Single command for all methods
- Parallel submission
- Consistent parameters
- Clear status matrix

### Option B: Two-Step Workflow (Prep then Run)

Prepare inputs, validate, then run:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Step 2a: Prep all method inputs
./submit_all_methods.sh conf_ils_low_10M --prep-only

# Wait for prep jobs, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step prep

# Step 2b: Run all methods
cd ../jobs
./submit_all_methods.sh conf_ils_low_10M --run-only
```

**Advantages:**
- Can validate inputs before running
- Can run methods at different times
- Easier to debug prep issues

### Option C: Individual Methods

Run methods separately with custom parameters:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# GRAMPA
./submit_grampa_pipeline.sh conf_ils_low_10M

# Polyphest with custom percentile
./submit_polyphest_pipeline.sh conf_ils_low_10M --percentile 50

# PADRE with more memory
./submit_padre_pipeline.sh conf_ils_low_10M --java-mem 8g

# MPSUGAR with more iterations
./submit_mpsugar_pipeline.sh conf_ils_low_10M --iterations 1000 --chains 2
```

**See:** `METHODS_GUIDE.md` for detailed method documentation

### Monitor Method Jobs

```bash
# All your jobs
squeue -u $USER

# Prep jobs
squeue -u $USER | grep prep

# Run jobs
squeue -u $USER | grep -E "grampa|polyphest|padre|mpsugar"

# Monitor logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grampa_*.out
```

### Validate Method Completion

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Validate prep inputs ready
python check_pipeline_status.py conf_ils_low_10M --step prep

# Validate run outputs generated
python check_pipeline_status.py conf_ils_low_10M --step run

# Detailed validation (shows all networks)
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# Validate specific method
python check_pipeline_status.py conf_ils_low_10M --step run --method grampa
```

### Expected Outputs

```
simulations/simulations/{NETWORK}/
├── processed/{CONFIG}/          # Prep stage outputs
│   ├── grampa_input/replicate_{1-5}/
│   ├── polyphest_input/replicate_{1-5}/
│   ├── padre_input/replicate_{1-5}/
│   └── mpsugar_input/replicate_{1-5}/
└── results/{CONFIG}/            # Run stage outputs
    ├── grampa/replicate_{1-5}/
    │   └── grampa-scores.txt
    ├── polyphest_p{PERCENTILE}/replicate_{1-5}/
    │   └── polyphest_trees-polyphest.txt
    ├── padre/replicate_{1-5}/
    │   └── padre_trees-result.tre
    └── mpsugar/replicate_{1-5}/
        └── mpsugar_results.txt
```

---

## Step 3: Post-Process Results

### What This Does

Extracts clean MUL-tree strings from program-specific output formats. Each method outputs results in different formats - this step standardizes them all to `{method}_result.tre` files.

### Running Post-Processing

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Post-process all methods for one config
python postprocess_results.py conf_ils_low_10M

# Specific methods only
python postprocess_results.py conf_ils_low_10M --methods grampa polyphest_p60

# Multiple configs (run sequentially or in parallel)
python postprocess_results.py conf_ils_low_10M
python postprocess_results.py conf_ils_medium_10M
python postprocess_results.py conf_ils_high_10M

# Or in parallel (faster)
python postprocess_results.py conf_ils_low_10M &
python postprocess_results.py conf_ils_medium_10M &
python postprocess_results.py conf_ils_high_10M &
wait

# Dry-run to preview
python postprocess_results.py conf_ils_low_10M --dry-run
```

### What Gets Processed

For each method, raw outputs are converted:

| Method | Raw File | Extracted To |
|--------|----------|--------------|
| **GRAMPA** | `grampa-scores.txt` (TSV) | `grampa_result.tre` |
| **Polyphest** | `polyphest_trees-polyphest.txt` | `polyphest_result.tre` |
| **MPSUGAR** | `mpsugar_results.txt` | `mpsugar_result.tre` |
| **PADRE** | `padre_trees-result.tre` (copy) | `padre_result.tre` |

### Expected Output

```
simulations/simulations/{NETWORK}/results/{CONFIG}/
├── grampa/replicate_{1-5}/
│   ├── grampa-scores.txt       # Original
│   └── grampa_result.tre       # Extracted ✓
├── polyphest_p60/replicate_{1-5}/
│   ├── polyphest_trees-polyphest.txt  # Original
│   └── polyphest_result.tre           # Extracted ✓
...
```

**Total files:** 630 result files (21 networks × 6 methods × 5 replicates)

### Validation

The script reports:
- **Total:** Expected files (630 for full pipeline)
- **Success:** Successfully extracted
- **Skipped:** Missing raw outputs or already extracted
- **Failed:** Extraction errors

---

## Step 4: Analyze and Compare

### What This Does

The full summary pipeline performs comprehensive analysis:
1. **Collect:** Inventory all result files
2. **Compare:** Compute metrics between inferred and true networks
3. **Aggregate:** Summarize across replicates
4. **Rank:** Method performance and consistency

### Running Full Summary

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Full analysis for one config
python run_full_summary.py conf_ils_low_10M

# Specific methods only
python run_full_summary.py conf_ils_low_10M --methods grampa polyphest_p60

# Force recompute (ignore cache)
python run_full_summary.py conf_ils_low_10M --force-recompute

# Multiple configs (run separately)
python run_full_summary.py conf_ils_low_10M
python run_full_summary.py conf_ils_medium_10M
python run_full_summary.py conf_ils_high_10M

# Or in parallel (faster, but uses more CPU)
python run_full_summary.py conf_ils_low_10M &
python run_full_summary.py conf_ils_medium_10M &
python run_full_summary.py conf_ils_high_10M &
wait
```

### Metrics Computed

For each inferred network vs. true network:

1. **edit_distance** - Overall structural similarity (0-1)
2. **num_rets_diff** - Difference in reticulation count
3. **ploidy_diff.dist** - Polyploidy detection accuracy
4. **ret_leaf_jaccard.dist** - Reticulation leaf set overlap
5. **ret_sisters_jaccard.dist** - Reticulation topology accuracy

**Plus:** False positives, false negatives, true positives for applicable metrics

### Output Files Generated

```
simulations/analysis/summary/{CONFIG}/
├── inventory.csv                   # File existence inventory
├── comparisons_raw.csv             # All replicate comparisons (~8,190 rows)
├── aggregated_metrics.csv          # Mean ± std across replicates
├── comparison_report.txt           # Success/failure summary
├── cache/                          # Cached comparisons (auto-managed)
│   └── *.pkl
├── level1_detailed_per_network/    # Per-network pivot tables
│   ├── edit_distance.csv
│   ├── edit_distance_mean.csv
│   ├── edit_distance_std.csv
│   ├── edit_distance_n.csv
│   ├── num_rets_diff*.csv
│   ├── ploidy_diff*.csv
│   ├── ret_leaf_jaccard*.csv
│   └── ret_sisters_jaccard*.csv
├── level2_method_rankings.csv      # Method rankings
└── level3_network_correlations.csv # Network property correlations
```

### Understanding the Results

**Level 1 - Per-Network Metrics:**
- Each CSV shows one metric across networks (rows) and methods (columns)
- Format: `mean ± std (n)` where n = number of successful replicates
- Lower values generally better (distances)

**Level 2 - Method Rankings:**
- Overall performance ranking per metric
- Consistency scores (lower std = more consistent)

**Level 3 - Network Correlations:**
- How network properties affect method performance
- Requires `mul_tree_final_stats.csv` network characteristics file

---

## Optional: Additional Analyses

These are **separate tools** not required for the main pipeline:

### Robinson-Foulds Distances

Calculate RF distances between inferred and true species trees:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

./submit_rf.sh conf_ils_low_10M

# Outputs: rf_distance_results.txt in each replicate directory
```

### Duplication/Loss Event Extraction

Extract dup/loss events from SimPhy databases:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

./submit_dup_loss_extraction.sh conf_ils_low_10M

# Outputs: dup_loss_summary.txt per replicate + aggregated summary
```

**Note:** These are typically run **after** SimPhy for characterizing the simulation scenarios, not for evaluating methods.

---

## Complete Examples

### Example 1: Single Config, All Methods

```bash
# Navigate to jobs directory
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# 1. Generate simulations
./submit_simphy.sh conf_ils_low_10M

# Wait for SimPhy jobs, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step simphy

# 2. Run all methods
cd ../jobs
./submit_all_methods.sh conf_ils_low_10M

# Wait for method jobs, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# 3. Post-process
python postprocess_results.py conf_ils_low_10M

# 4. Analyze
python run_full_summary.py conf_ils_low_10M

# 5. View results
cd ../analysis/summary/conf_ils_low_10M/
ls level1_detailed_per_network/edit_distance.csv
```

### Example 2: Multiple Configs, Prep Then Run

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# 1. Simulations (already done)
# Assume conf_dup_loss_{low,medium,high}_10M already simulated

# 2a. Prep all methods for all configs
./submit_all_methods.sh conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M --prep-only

# Wait and validate prep
cd ../scripts
python check_pipeline_status.py conf_dup_loss_low_10M --step prep
python check_pipeline_status.py conf_dup_loss_medium_10M --step prep
python check_pipeline_status.py conf_dup_loss_high_10M --step prep

# 2b. Run all methods
cd ../jobs
./submit_all_methods.sh conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M --run-only

# Wait and validate run
cd ../scripts
python check_pipeline_status.py conf_dup_loss_low_10M --step run
python check_pipeline_status.py conf_dup_loss_medium_10M --step run
python check_pipeline_status.py conf_dup_loss_high_10M --step run

# 3. Post-process all configs
for config in conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M; do
    python postprocess_results.py $config
done

# 4. Analyze all configs (can run in parallel)
python run_full_summary.py conf_dup_loss_low_10M &
python run_full_summary.py conf_dup_loss_medium_10M &
python run_full_summary.py conf_dup_loss_high_10M &
wait

echo "All analyses complete!"
```

### Example 3: Polyphest Parameter Sweep

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Run Polyphest with different percentiles
for percentile in 50 60 70 90; do
    ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile $percentile
done

# Wait for completion, then validate each
cd ../scripts
for percentile in 50 60 70 90; do
    python check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile $percentile
done

# Post-process all percentile variants
python postprocess_results.py conf_ils_low_10M --methods polyphest_p50 polyphest_p60 polyphest_p70 polyphest_p90

# Compare all variants
python run_full_summary.py conf_ils_low_10M --methods polyphest_p50 polyphest_p60 polyphest_p70 polyphest_p90

# Results show performance across percentile choices
```

---

## Troubleshooting

### General Workflow Issues

**Problem:** Jobs pending for long time

**Cause:** Cluster queue congestion

**Solution:**
```bash
# Check queue status
squeue | wc -l  # Total jobs in queue

# Check your jobs
squeue -u $USER

# Consider running smaller batches or during off-peak hours
```

**Problem:** Out of disk space

**Cause:** SimPhy and method outputs are large

**Solution:**
```bash
# Check disk usage
du -sh /groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/

# Clean up old logs if needed
rm /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/*_old_*.{out,err}
```

### Pipeline-Specific Issues

See individual guides:
- **SimPhy issues:** Check SLURM logs in `/simulations/logs/simphy_*.err`
- **Method issues:** See `METHODS_GUIDE.md`
- **Validation issues:** See `PIPELINE_STATUS_GUIDE.md`
- **Analysis issues:** Check summary pipeline logs

### Getting Help

1. **Check logs:** `/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/`
2. **Validate:** Use `check_pipeline_status.py --verbose` for detailed diagnostics
3. **Dry-run:** Use `--dry-run` flags to preview commands before executing

---

## Time Estimates

Approximate wall-clock times (cluster dependent):

| Step | Time | Notes |
|------|------|-------|
| **SimPhy** | 2-12 hours | Depends on batch mode, retries |
| **Methods Prep** | 1-3 hours | Parallel across networks |
| **GRAMPA** | 3-6 hours | ASTRAL + GRAMPA |
| **Polyphest** | 2-4 hours | Fast |
| **PADRE** | 2-5 hours | Medium |
| **MPSUGAR** | 6-24 hours | Slowest, depends on iterations |
| **Post-processing** | 5-10 minutes | Fast Python script |
| **Summary** | 10-30 minutes | Comparison computation |

**Total for full pipeline:** ~1-2 days from SimPhy to final summary

---

## Quick Reference Card

```bash
# ====================
# COMPLETE WORKFLOW
# ====================

# 1. SIMPHY
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs
./submit_simphy.sh CONFIG
cd ../scripts && python check_pipeline_status.py CONFIG --step simphy

# 2. METHODS
cd ../jobs
./submit_all_methods.sh CONFIG
cd ../scripts && python check_pipeline_status.py CONFIG --step run

# 3. POST-PROCESS
python postprocess_results.py CONFIG

# 4. ANALYZE
python run_full_summary.py CONFIG

# 5. VIEW RESULTS
cd ../analysis/summary/CONFIG/
```

---

## Summary

This complete pipeline transforms raw species trees into comprehensive method evaluations:

1. **SimPhy** generates realistic gene trees (105,000 per config)
2. **Methods** infer networks from gene trees (4 methods × 21 networks × 5 replicates)
3. **Post-processing** standardizes outputs for comparison
4. **Analysis** computes performance metrics and rankings

**Result:** Quantitative evaluation of network inference method performance under different evolutionary scenarios.

For detailed information on specific steps:
- **Methods:** See `METHODS_GUIDE.md`
- **Validation:** See `PIPELINE_STATUS_GUIDE.md`
- **Repository structure:** See `../CLAUDE.md`
