# Complete Pipeline Guide: SimPhy to Summary

This guide covers the **entire workflow** for phylogenetic network inference simulations, from generating gene trees through final analysis.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Step 1: Generate Simulations (SimPhy)](#step-1-generate-simulations-simphy)
3. [Step 2: Simulate Sequences (AliSim)](#step-2-simulate-sequences-alisim)
4. [Step 3: Run Network Inference Methods](#step-3-run-network-inference-methods)
5. [Step 4: Post-Process Results](#step-4-post-process-results)
6. [Step 5: Analyze and Compare](#step-5-analyze-and-compare)
7. [Optional: Additional Analyses](#optional-additional-analyses)
8. [Complete Examples](#complete-examples)
9. [Troubleshooting](#troubleshooting)

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│ COMPLETE SIMULATION PIPELINE                                        │
└─────────────────────────────────────────────────────────────────────┘

STEP 1: GENERATE SIMULATIONS (Gene Trees)
   │
   ├─> submit_simphy.sh
   │       ↓
   │   21 networks × 5 replicates × 1000 gene trees
   │       ↓
   └─> check_pipeline_status.py --step simphy
           ✓ Validate 105,000 gene trees generated

STEP 2: SIMULATE SEQUENCES (DNA Alignments)
   │
   ├─> submit_sequences.sh
   │       ↓
   │   Sample GTR+Gamma parameters from 2,709 empirical genes
   │       ↓
   │   21 networks × 5 replicates × 1000 alignments
   │       ↓
   └─> Verification: Check alignments generated
           ✓ 105,000 sequence alignments with realistic parameters

STEP 3: RUN NETWORK INFERENCE METHODS
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

STEP 4: POST-PROCESS RESULTS
   │
   └─> postprocess_results.py
           ↓
       Extract clean MUL-trees from method-specific formats
           ✓ 630 standardized result files (21 × 6 methods × 5 reps)

STEP 5: ANALYZE AND COMPARE
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

## Step 2: Simulate Sequences (AliSim)

### What This Does

After SimPhy generates gene trees, this step simulates realistic DNA sequences along those trees using:
- **GTR+Gamma substitution model** with parameters sampled from empirical data
- **2,709 real genes** from three datasets (Zhao_2021, Ren_2024, Morales_Briones_2021)
- **Alignment lengths** sampled from empirical distribution (147-6,487 bp, median 609 bp)

**Key Design Decision:** Each gene tree gets the same GTR+Gamma parameters across all replicates (just like SimPhy uses the same substitution rate across replicates). This ensures:
- Replicates represent different stochastic realizations of the **same** evolutionary scenario
- Variance comes from random mutations, not parameter variation
- Proper statistical evaluation of method performance under specific conditions

### Prerequisites

```bash
# 1. Ensure GTR parameters are extracted (one-time setup)
ls -lh /groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/gtr_parameters_all.pkl

# If file doesn't exist, extract parameters:
cd /groups/itay_mayrose/tomulanovski/gene2net
python simulations/scripts/sequence_evolution/filter_and_extract_gtr.py

# 2. Verify SimPhy completed successfully
cd simulations/scripts
python check_pipeline_status.py conf_ils_low_10M --step simphy

# That's it! No tree splitting needed - the script reads directly from SimPhy output
```

### Running Sequence Simulation

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Simulate sequences for all 21 networks
./submit_sequences.sh conf_ils_low_10M

# Test on single network first (recommended)
./submit_sequences.sh conf_ils_low_10M Ding_2023

# Different configurations
./submit_sequences.sh conf_ils_med_10M
./submit_sequences.sh conf_ils_high_10M
```

### What Happens During Execution

For each gene tree (e.g., gene tree #42 in network "Ding_2023"):

1. **Detect batch structure** - Script automatically determines:
   - Single batch: Reads `replicate_1/1/g_trees0042.trees`
   - Batches of 10: Calculates batch 5, reads `replicate_1/batch_5/1/g_trees2.trees`
   - Batches of 1: Reads `replicate_1/batch_42/1/g_trees1.trees`

2. **Sample parameters once** from the 2,709 empirical genes:
   ```
   Example sampled parameters:
   - GTR rates: AC=1.52, AG=4.31, AT=1.08, CG=0.67, CT=5.29, GT=1.00
   - Base frequencies: πA=0.289, πC=0.211, πG=0.211, πT=0.289
   - Alpha (Gamma shape): 0.543
   - Alignment length: 609 bp
   ```

3. **Use same parameters across all replicates**, but with **different random seeds**:
   ```
   Replicate 1: Same model, seed=421234 → alignment_0042.phy
   Replicate 2: Same model, seed=422567 → alignment_0042.phy (different sequences!)
   Replicate 3: Same model, seed=423891 → alignment_0042.phy (different sequences!)
   ...
   ```

**Why same parameters?** Just like SimPhy uses the same substitution rate across replicates, sequence simulation uses the same GTR+Gamma parameters. Replicates = different random outcomes of the **same** evolutionary process.

### Monitor Sequence Simulation Jobs

```bash
# Check running jobs
squeue -u $USER | grep alisim

# Check job status
sacct -X --format=JobName,State,ExitCode | grep alisim_conf_ils_low_10M

# Monitor logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_conf_ils_low_10M_Ding_2023_1.out
```

### Validate Sequence Completion

**Recommended: Use the pipeline status checker**

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Check sequence alignments for all networks
python check_pipeline_status.py conf_ils_low_10M --step sequences

# Verbose output to see details
python check_pipeline_status.py conf_ils_low_10M --step sequences --verbose
```

**Manual verification (alternative)**

```bash
# Check alignments were generated for a network
NETWORK="Ding_2023"
CONFIG="conf_ils_low_10M"
BASE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Count alignments per replicate (should be 1000 each)
for rep in {1..5}; do
    count=$(ls ${BASE}/${NETWORK}/data/${CONFIG}/replicate_${rep}/1/alignments/*.phy 2>/dev/null | wc -l)
    echo "Replicate $rep: $count alignments"
done

# Inspect a sample alignment
head -20 ${BASE}/${NETWORK}/data/${CONFIG}/replicate_1/1/alignments/alignment_0001.phy
```

### Expected Outputs

**Single batch mode:**
```
simulations/simulations/{NETWORK}/data/{CONFIG}/
└── replicate_{1-5}/
    └── 1/
        ├── g_trees0001.trees        # Gene trees (from SimPhy)
        ├── g_trees0002.trees
        ├── ...
        └── alignments/              # NEW: Sequence alignments
            ├── alignment_0001.phy   # Sequences for gene tree 0001
            ├── alignment_0002.phy   # Sequences for gene tree 0002
            └── ...
```

**Multi-batch mode:**
```
simulations/simulations/{NETWORK}/data/{CONFIG}/
└── replicate_{1-5}/
    ├── batch_1/1/
    │   ├── g_trees1.trees ... g_treesN.trees   # Gene trees (SimPhy)
    ├── batch_2/1/
    │   ├── g_trees1.trees ... g_treesN.trees
    ├── ...
    └── 1/
        └── alignments/              # NEW: Sequence alignments (all here!)
            ├── alignment_0001.phy   # Sequences for gene tree 0001
            ├── alignment_0002.phy   # Sequences for gene tree 0002
            └── ...
```

**Total output:** 105,000 sequence alignments (21 networks × 5 replicates × 1000 genes)

**Note:** Alignments are always placed in `replicate_N/1/alignments/`, regardless of batch mode.

### Parameter Sampling Details

**Special handling for extreme alpha values:**
- For α ≤ 3.0: Use `GTR+Gamma` model (rate heterogeneity)
- For α > 3.0: Use `GTR` only (uniform rates across sites)
  - Affects 90 genes (3.3% of empirical data)
  - Biologically reasonable: high α ≈ uniform rates

**Empirical parameter statistics:**
- **GTR rates:** Transitions (AG, CT) higher than transversions (AC, AT, CG)
- **Base composition:** Moderately AT-rich (56% AT, 44% GC on average)
- **Alpha distribution:** Median 0.54, mean 18.47 (right-skewed due to outliers)
- **Alignment lengths:** Range 147-6,487 bp, median 609 bp

### Troubleshooting Sequence Simulation

**Problem:** "GTR parameters file not found"

**Solution:**
```bash
cd /groups/itay_mayrose/tomulanovski/gene2net
python simulations/scripts/sequence_evolution/filter_and_extract_gtr.py
```

**Problem:** "No replicate directories found"

**Solution:** Run SimPhy first for that configuration:
```bash
./submit_simphy.sh conf_ils_low_10M
```

**Problem:** Some alignments missing

**Solution:** Check job logs for errors:
```bash
grep -l "ERROR" /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alisim_*.err
```

---

## Step 3: Run Network Inference Methods

### What This Does

Runs phylogenetic network inference on simulated gene trees using:
- **GRAMPA:** Species tree + reticulation graph
- **Polyphest:** Maximum parsimony network
- **PADRE:** Parsimony-based network
- **MPSUGAR:** Bayesian network inference
- **AlloppNET:** Bayesian MCMC network for allopolyploids (8 compatible networks only)

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

# Include AlloppNET (runs on 8 compatible networks only)
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest,alloppnet
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

# AlloppNET (all-in-one pipeline, ~5 days per replicate)
./submit_alloppnet_pipeline.sh conf_ils_low_10M

# AlloppNET - specific replicates only
./submit_alloppnet_pipeline.sh conf_ils_low_10M --replicates 1,3,5
```

**See:**
- `METHODS_GUIDE.md` for GRAMPA, Polyphest, PADRE, MPSUGAR
- `ALLOPPNET_GUIDE.md` for AlloppNET detailed documentation

### Monitor Method Jobs

```bash
# All your jobs
squeue -u $USER

# Prep jobs
squeue -u $USER | grep prep

# Run jobs
squeue -u $USER | grep -E "grampa|polyphest|padre|mpsugar|alloppnet"

# Monitor logs (GRAMPA example)
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grampa_*.out

# Monitor AlloppNET progress (~5 days per job)
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/alloppnet_*.out
grep "STATE_" results/<Network>/<config>/alloppnet/replicate_1/sampledmultrees.txt | tail
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
│   ├── mpsugar_input/replicate_{1-5}/
│   └── alloppnet_input/replicate_{1-5}/  # Only for 8 compatible networks
│       ├── alignment_*.nex (1000 files)
│       ├── taxa_table.txt
│       └── ploidy_level.json
└── results/{CONFIG}/            # Run stage outputs
    ├── grampa/replicate_{1-5}/
    │   └── grampa-scores.txt
    ├── polyphest_p{PERCENTILE}/replicate_{1-5}/
    │   └── polyphest_trees-polyphest.txt
    ├── padre/replicate_{1-5}/
    │   └── padre_trees-result.tre
    ├── mpsugar/replicate_{1-5}/
    │   └── mpsugar_results.txt
    └── alloppnet/replicate_{1-5}/  # Only for 8 compatible networks
        ├── alloppnet_final.tre  # Final consensus tree
        ├── sampledmultrees.txt  # BEAST MCMC output
        └── alloppnet.XML        # BEAST configuration
```

---

## Step 4: Post-Process Results

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

## Step 5: Analyze and Compare

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

# Per-replicate output: rf_distance_results.txt in each replicate directory
#   Location: simulations/{NETWORK}/data/{CONFIG}/replicate_{1-5}/rf_distance_results.txt
#
# Aggregated output (AUTOMATICALLY CREATED): aggregate_AD_{CONFIG}.csv
#   Location: /groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/aggregate_AD_conf_ils_low_10M.csv
#   Contains: Mean ± std AD (Average Distance) across replicates for each network
```

**What it does:**
- Calculates Robinson-Foulds distances for each replicate (per-replicate analysis)
- **Automatically aggregates** results across all networks and replicates (no separate command needed)
- **Handles multi-batch SimPhy outputs:** Automatically detects and merges batches of 1, 10, or 100 trees
- Single SLURM job processes all 21 networks × 5 replicates

### Duplication/Loss Event Extraction

Extract dup/loss events from SimPhy databases:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

./submit_dup_loss_extraction.sh conf_ils_low_10M

# Per-replicate output: dup_loss_summary.txt
#   Location: simulations/{NETWORK}/data/{CONFIG}/replicate_{1-5}/dup_loss_summary.txt
#   Contains: Duplication and loss event counts for that replicate
#
# Aggregated output (AUTOMATICALLY CREATED): dup_loss_summary_{CONFIG}.txt
#   Location: /groups/itay_mayrose/tomulanovski/gene2net/simulations/dup_loss_summary_conf_ils_low_10M.txt
#   Contains: Summary statistics (mean ± std) across all networks and replicates
```

**What it does:**
- Extracts gene trees, locus trees, and species trees from SimPhy SQLite databases
- Counts duplication and loss events per replicate
- **Automatically aggregates** results across all networks and replicates (no separate command needed)
- **Handles multi-batch SimPhy outputs:** Automatically detects and merges all batches (1, 10, or 100 trees per batch)
- Single SLURM job processes all 21 networks × 5 replicates

**Note:** These are typically run **after** SimPhy for characterizing the simulation scenarios, not for evaluating methods.

---

## Complete Examples

### Example 1: Single Config, All Methods

```bash
# Navigate to jobs directory
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# 1. Generate gene trees (SimPhy)
./submit_simphy.sh conf_ils_low_10M

# Wait for SimPhy jobs, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step simphy

# 2. Simulate sequences (AliSim)
cd ../jobs
./submit_sequences.sh conf_ils_low_10M

# Wait for sequence jobs, then validate
NETWORK="Ding_2023"
CONFIG="conf_ils_low_10M"
BASE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"
for rep in {1..5}; do
    count=$(ls ${BASE}/${NETWORK}/data/${CONFIG}/replicate_${rep}/1/alignments/*.phy 2>/dev/null | wc -l)
    echo "Replicate $rep: $count alignments"
done

# 3. Run all methods
cd ../jobs
./submit_all_methods.sh conf_ils_low_10M

# Wait for method jobs, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# 4. Post-process
python postprocess_results.py conf_ils_low_10M

# 5. Analyze
python run_full_summary.py conf_ils_low_10M

# 6. View results
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
| **Sequence Simulation** | 1-3 hours | Parallel array jobs (1000 genes) |
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

# 1. SIMPHY (Gene Trees)
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs
./submit_simphy.sh CONFIG
cd ../scripts && python check_pipeline_status.py CONFIG --step simphy

# 2. SEQUENCES (DNA Alignments)
cd ../jobs
./submit_sequences.sh CONFIG
cd ../scripts && python check_pipeline_status.py CONFIG --step sequences

# 3. METHODS (Network Inference)
cd ../jobs
./submit_all_methods.sh CONFIG
cd ../scripts && python check_pipeline_status.py CONFIG --step run

# 4. POST-PROCESS (Extract Results)
python postprocess_results.py CONFIG

# 5. ANALYZE (Compare & Summarize)
python run_full_summary.py CONFIG

# 6. VIEW RESULTS
cd ../analysis/summary/CONFIG/
```

---

## Summary

This complete pipeline transforms raw species trees into comprehensive method evaluations:

1. **SimPhy** generates realistic gene trees (105,000 per config)
2. **Sequence Simulation** generates DNA alignments using empirically-sampled GTR+Gamma parameters (105,000 alignments per config)
3. **Methods** infer networks from gene trees (4 methods × 21 networks × 5 replicates)
4. **Post-processing** standardizes outputs for comparison
5. **Analysis** computes performance metrics and rankings

**Result:** Quantitative evaluation of network inference method performance under different evolutionary scenarios, using biologically realistic sequence data.

For detailed information on specific steps:
- **Sequence Simulation:** See `SEQUENCE_SIMULATION_GUIDE.md`
- **Methods:** See `METHODS_GUIDE.md`
- **Validation:** See `PIPELINE_STATUS_GUIDE.md`
- **Repository structure:** See `../CLAUDE.md`
