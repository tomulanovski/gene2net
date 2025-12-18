# Methods Guide: Running Phylogenetic Network Inference

This guide covers running network inference methods (GRAMPA, Polyphest, PADRE, MPSUGAR, AlloppNET) on simulated gene trees.

**Prerequisites:** SimPhy simulations must be completed and validated.

**Note:** AlloppNET has a dedicated guide (`ALLOPPNET_GUIDE.md`) due to its unique all-in-one pipeline and network compatibility requirements.

---

## Table of Contents

1. [Quick Start: Master Orchestration Script](#quick-start-master-orchestration-script)
2. [Running Individual Methods](#running-individual-methods)
3. [Method-Specific Parameters](#method-specific-parameters)
4. [Common Usage Patterns](#common-usage-patterns)
5. [Validation and Monitoring](#validation-and-monitoring)
6. [Troubleshooting](#troubleshooting)

---

## Quick Start: Master Orchestration Script

**RECOMMENDED:** Use the master script to run all or selected methods across multiple configurations.

### Basic Usage

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Run all methods for one config
./submit_all_methods.sh conf_ils_low_10M

# Run all methods for multiple configs
./submit_all_methods.sh conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M

# Run specific methods only
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest

# Dry-run to preview
./submit_all_methods.sh conf_ils_low_10M --dry-run
```

### Two-Step Workflow (Prep then Run)

```bash
# Step 1: Prep inputs for all methods
./submit_all_methods.sh conf_ils_low_10M --prep-only

# Wait for jobs to complete, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step prep

# Step 2: Run all methods
cd ../jobs
./submit_all_methods.sh conf_ils_low_10M --run-only
```

### Command-Line Options

```
Usage: ./submit_all_methods.sh [CONFIG1 CONFIG2 ...] [OPTIONS]

Configuration:
  CONFIGS...                      One or more configuration names
                                  Default: conf_ils_low/medium/high_10M

Method Control:
  --methods METHOD1,METHOD2,...   Methods to run (grampa,polyphest,padre,mpsugar,alloppnet)
                                  Default: all methods
  --prep-only                     Only prepare inputs
  --run-only                      Only run methods (assumes prep done)

Common Parameters:
  --replicates N                  Number of replicates (default: 5)
  --dry-run                       Preview without submitting
  --verbose                       Show detailed output

Method-Specific:
  --polyphest-percentile N        Polyphest percentile (default: 60)
  --polyphest-iso-threshold N     Iso threshold (default: 0.2)
  --padre-java-mem SIZE           PADRE Java heap (default: 4g)
  --mpsugar-iterations N          MPSUGAR iterations (default: 500)
  --mpsugar-chains N              MPSUGAR chains (default: 1)
```

---

## Running Individual Methods

If you need to run methods separately (e.g., different timing or parameters):

### GRAMPA

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Full pipeline (prep + ASTRAL + GRAMPA)
./submit_grampa_pipeline.sh conf_ils_low_10M

# Prep only
./submit_grampa_pipeline.sh conf_ils_low_10M --prep-only

# ASTRAL only (assumes prep done)
./submit_grampa_pipeline.sh conf_ils_low_10M --astral-only

# GRAMPA only (assumes prep + ASTRAL done)
./submit_grampa_pipeline.sh conf_ils_low_10M --run-only

# Skip prep, run ASTRAL + GRAMPA
./submit_grampa_pipeline.sh conf_ils_low_10M --skip-prep
```

**GRAMPA Pipeline Stages:**
1. **Prep** (21 tasks): Fix substring issues, clean trees
2. **ASTRAL** (21 tasks): Infer species trees
3. **GRAMPA** (105 tasks): Network inference (21 networks × 5 replicates)

**Critical:** GRAMPA requires `fix_substrings_grampa.py` preprocessing in prep stage.

### Polyphest

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Full pipeline (default percentile 60)
./submit_polyphest_pipeline.sh conf_ils_low_10M

# Custom percentile
./submit_polyphest_pipeline.sh conf_ils_low_10M --percentile 50

# Custom iso-threshold
./submit_polyphest_pipeline.sh conf_ils_low_10M --iso-threshold 0.3

# Prep only
./submit_polyphest_pipeline.sh conf_ils_low_10M --prep-only

# Run only
./submit_polyphest_pipeline.sh conf_ils_low_10M --run-only
```

**Polyphest Pipeline Stages:**
1. **Prep** (21 tasks): Create multi-set files
2. **Run** (105 tasks): Network inference with specified percentile

**Note:** Output directory includes percentile: `polyphest_p{PERCENTILE}/`

### PADRE

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Full pipeline (default 4g memory)
./submit_padre_pipeline.sh conf_ils_low_10M

# Custom Java heap size
./submit_padre_pipeline.sh conf_ils_low_10M --java-mem 8g

# Prep only
./submit_padre_pipeline.sh conf_ils_low_10M --prep-only

# Run only
./submit_padre_pipeline.sh conf_ils_low_10M --run-only
```

**PADRE Pipeline Stages:**
1. **Prep** (21 tasks): Format conversion
2. **Run** (105 tasks): Network inference

### MPSUGAR

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Full pipeline (default 500 iterations, 1 chain)
./submit_mpsugar_pipeline.sh conf_ils_low_10M

# Custom iterations and chains
./submit_mpsugar_pipeline.sh conf_ils_low_10M --iterations 1000 --chains 2

# Prep only
./submit_mpsugar_pipeline.sh conf_ils_low_10M --prep-only

# Run only
./submit_mpsugar_pipeline.sh conf_ils_low_10M --run-only
```

**MPSUGAR Pipeline Stages:**
1. **Prep** (21 tasks): Convert to NEXUS, create taxon map JSON
2. **Run** (105 tasks): Bayesian network inference

### AlloppNET

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Full pipeline (prep + BEAST run, ~5 days per network)
./submit_alloppnet_pipeline.sh conf_ils_low_10M

# Two-step workflow (RECOMMENDED)
# Step 1: Prep (fast, ~minutes)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --prep-only

# Step 2: Run BEAST (slow, ~5 days)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --run-only

# Specific replicates only
./submit_alloppnet_pipeline.sh conf_ils_low_10M --replicates 1,3,5 --prep-only
```

**AlloppNET Pipeline Stages:**
1. **Prep** (8 tasks, ~minutes each):
   - Convert PHY → NEXUS (1000 alignments per network)
   - Analyze copy number distributions using kernel smoothing
   - Generate `ploidy_level.json` (robust ploidy inference)
   - Create `taxa_table.txt` (homeolog pairing for allotetraploids)
   - Generate BEAST XML using AlloppDT scripts
2. **Run** (8 tasks, ~5 days each):
   - Run BEAST (100M iterations)
   - Summarize with TreeAnnotator (10% burnin)
   - Post-process (remove copy suffixes from tree)

**Important Notes:**
- **Network compatibility:** Only runs on 8 networks with max 2 copies (diploid/tetraploid species):
  - Bendiksby_2011, Ding_2023, Koenen_2020, Liu_2023
  - Shahrestani_2015, Wisecaver_2023, Wu_2015, Zhao_2021
- **Total jobs:** 8 networks × 5 replicates = 40 BEAST runs
- **Runtime:** ~5 days per BEAST run (use `--prep-only` to validate inputs first)
- **Input:** Uses sequence alignments (not just gene trees like other methods)
- **Validation:** Use `--method alloppnet` when checking status

---

## Method-Specific Parameters

### Polyphest Percentile

**What it does:** Controls the threshold for reticulation detection.

**Common values:**
- `50` - Conservative (fewer reticulations)
- `60` - Default (balanced)
- `70` - Moderate (more reticulations)
- `90` - Liberal (most reticulations)

**Example comparing percentiles:**
```bash
./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 50
./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 70
./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile 90
```

### PADRE Memory

**What it does:** Sets Java heap size for PADRE execution.

**Guidelines:**
- `4g` - Default (works for most cases)
- `8g` - Large networks or many gene trees
- `16g` - Very large datasets

**Example:**
```bash
./submit_all_methods.sh conf_ils_low_10M --methods padre --padre-java-mem 8g
```

### MPSUGAR Iterations and Chains

**What they do:**
- **Iterations:** MCMC chain length (convergence quality)
- **Chains:** Number of independent MCMC runs

**Guidelines:**
- Default: `500 iterations, 1 chain` - Fast, good for testing
- Conservative: `1000 iterations, 2 chains` - Better convergence
- Thorough: `5000 iterations, 4 chains` - Publication quality

**Example:**
```bash
./submit_all_methods.sh conf_ils_low_10M --methods mpsugar --mpsugar-iterations 1000 --mpsugar-chains 2
```

---

## Common Usage Patterns

### Run All Methods Across Multiple Configs

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# ILS configs
./submit_all_methods.sh conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M

# Dup/loss configs
./submit_all_methods.sh conf_dup_loss_low_10M conf_dup_loss_medium_10M conf_dup_loss_high_10M
```

### Parameter Sweep for Single Method

```bash
# Test different Polyphest percentiles
for percentile in 50 60 70 90; do
    ./submit_all_methods.sh conf_ils_low_10M --methods polyphest --polyphest-percentile $percentile
done

# Test different MPSUGAR iterations
for iter in 500 1000 2000; do
    ./submit_all_methods.sh conf_ils_low_10M --methods mpsugar --mpsugar-iterations $iter
done
```

### Run Subset of Methods

```bash
# Only GRAMPA and Polyphest
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest

# Only fast methods (exclude MPSUGAR)
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest,padre
```

### Staged Execution

```bash
# Stage 1: Prep all methods at once
./submit_all_methods.sh conf_ils_low_10M --prep-only

# Wait for prep to complete, validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step prep

# Stage 2: Run methods separately (different timing/resources)
cd ../jobs
./submit_grampa_pipeline.sh conf_ils_low_10M --run-only &
./submit_polyphest_pipeline.sh conf_ils_low_10M --run-only &
./submit_padre_pipeline.sh conf_ils_low_10M --run-only &
# Wait for fast methods before submitting slow MPSUGAR
wait
./submit_mpsugar_pipeline.sh conf_ils_low_10M --run-only
```

---

## Validation and Monitoring

### Check Job Status

```bash
# All your jobs
squeue -u $USER

# Specific configuration
squeue -u $USER | grep conf_ils_low

# Specific method
squeue -u $USER | grep grampa
```

### Validate Pipeline Progress

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# Check prep inputs ready
python check_pipeline_status.py conf_ils_low_10M --step prep

# Check run outputs generated
python check_pipeline_status.py conf_ils_low_10M --step run

# Verbose output (shows all networks)
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# Check specific method
python check_pipeline_status.py conf_ils_low_10M --step run --method grampa

# Check Polyphest with specific percentile
python check_pipeline_status.py conf_ils_low_10M --step run --method polyphest --percentile 50

# Check AlloppNET (8 networks only)
python check_pipeline_status.py conf_ils_low_10M --step prep --method alloppnet
python check_pipeline_status.py conf_ils_low_10M --step run --method alloppnet --verbose
```

### Monitor Logs

```bash
# Tail prep logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/prep_grampa_*.out

# Tail run logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_polyphest_*.out

# Check for errors
grep -i error /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/*.err
```

---

## Troubleshooting

### Prep Stage Issues

**Problem:** Prep jobs fail for GRAMPA

**Solution:** Check that gene trees don't have substring issues
```bash
# The prep stage runs fix_substrings_grampa.py automatically
# Check prep logs for details:
cat /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/prep_grampa_*.err
```

**Problem:** MPSUGAR prep fails

**Solution:** Ensure NEXUS format is correct and taxon map is created
```bash
# Check prep output
ls -lh simulations/Bendiksby_2011/processed/conf_ils_low_10M/mpsugar_input/replicate_1/
# Should contain: mpsugar_trees.nex, taxon_map.json
```

### Run Stage Issues

**Problem:** PADRE runs out of memory

**Solution:** Increase Java heap size
```bash
./submit_all_methods.sh conf_ils_low_10M --methods padre --padre-java-mem 8g
```

**Problem:** MPSUGAR taking too long

**Solution:** Reduce iterations for testing
```bash
./submit_all_methods.sh conf_ils_low_10M --methods mpsugar --mpsugar-iterations 100
```

**Problem:** GRAMPA produces no output

**Solution:** Check that ASTRAL step completed successfully
```bash
# Validate ASTRAL outputs exist
ls simulations/Bendiksby_2011/processed/conf_ils_low_10M/grampa_input/replicate_1/species.tre

# If missing, re-run ASTRAL step
./submit_grampa_pipeline.sh conf_ils_low_10M --astral-only
```

### Validation Issues

**Problem:** check_pipeline_status.py shows missing outputs

**Solution:** Check if jobs are still running or failed
```bash
# Check job status
squeue -u $USER

# Check error logs for failed jobs
tail /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/*_error*.err
```

**Problem:** Polyphest output not found with different percentile

**Solution:** Specify the percentile you used
```bash
# If you ran with --percentile 50, validate with:
python check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 50
```

---

## Expected Outputs

### Prep Stage Output Locations

```
simulations/simulations/{NETWORK}/processed/{CONFIG}/
├── grampa_input/replicate_{1-5}/
│   ├── grampa_trees.tre      # Fixed substring trees
│   ├── clean_trees.tre        # Clean gene trees
│   └── species.tre            # ASTRAL species tree (after ASTRAL step)
├── polyphest_input/replicate_{1-5}/
│   ├── polyphest_trees.tre    # Formatted trees
│   └── multi_set.txt          # Multi-set file
├── padre_input/replicate_{1-5}/
│   └── padre_trees.tre        # Formatted trees
├── mpsugar_input/replicate_{1-5}/
│   ├── mpsugar_trees.nex      # NEXUS format
│   └── taxon_map.json         # Taxon mapping
└── alloppnet_input/replicate_{1-5}/
    ├── alignment_*.nex        # NEXUS alignments (1000 files)
    ├── ploidy_level.json      # Ploidy assignments (kernel smoothing)
    └── taxa_table.txt         # Homeolog pairing
```

### Run Stage Output Locations

```
simulations/simulations/{NETWORK}/results/{CONFIG}/
├── grampa/replicate_{1-5}/
│   └── grampa-scores.txt      # Best network with scores
├── polyphest_p{PERCENTILE}/replicate_{1-5}/
│   └── polyphest_trees-polyphest.txt  # Inferred network
├── padre/replicate_{1-5}/
│   └── padre_trees-result.tre  # Inferred network
├── mpsugar/replicate_{1-5}/
│   └── mpsugar_results.txt     # Inferred network with posterior
└── alloppnet/replicate_{1-5}/
    ├── alloppnet_final.tre    # Final consensus network
    ├── sampledmultrees.txt    # BEAST MCMC samples (100K trees)
    └── alloppnet.XML          # BEAST configuration
```

### Summary

- **Standard methods (GRAMPA, Polyphest, PADRE, MPSUGAR):**
  - **21 networks** × **5 replicates** = **105 total runs per method**
  - **4 methods** = **420 total output files**

- **AlloppNET (special case):**
  - **8 compatible networks** × **5 replicates** = **40 total runs**
  - Only runs on networks with max 2 copies (diploid/tetraploid only)
  - All-in-one pipeline: prep → BEAST → summarize (~5 days per run)
  - See `ALLOPPNET_GUIDE.md` for complete details

- All outputs validated with `check_pipeline_status.py --step run`

---

## Next Steps After Methods Complete

After all methods complete and are validated:

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts

# 1. Post-process raw outputs (extract clean trees)
python postprocess_results.py conf_ils_low_10M

# 2. Run full analysis and comparison
python run_full_summary.py conf_ils_low_10M
```

See `COMPLETE_PIPELINE_GUIDE.md` for the full workflow.

---

## Quick Reference

### Most Common Commands

```bash
# Run everything for one config
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs
./submit_all_methods.sh conf_ils_low_10M

# Check status
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# Post-process
python postprocess_results.py conf_ils_low_10M
python run_full_summary.py conf_ils_low_10M
```

### Help Commands

```bash
# Master script help
./submit_all_methods.sh --help

# Individual method help
./submit_grampa_pipeline.sh --help
./submit_polyphest_pipeline.sh --help
./submit_padre_pipeline.sh --help
./submit_mpsugar_pipeline.sh --help
./submit_alloppnet_pipeline.sh --help

# AlloppNET has its own dedicated guide
cat ../md_files/ALLOPPNET_GUIDE.md

# Validation help
python check_pipeline_status.py --help
```
