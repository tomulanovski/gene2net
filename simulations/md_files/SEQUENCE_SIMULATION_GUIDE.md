# Sequence Simulation with Empirical GTR+Gamma Parameters

## Overview

This guide describes how to simulate DNA sequences along gene trees using GTR+Gamma model parameters sampled from empirical distributions. The parameters come from 2,709 genes across three datasets: Zhao_2021 (982 genes), Ren_2024 (727 genes), and Morales_Briones_2021 (1000 genes).

## Key Features

- **Empirically-sampled parameters**: Each gene tree gets a randomly sampled parameter set from 2,709 real genes
- **Realistic sequence evolution**: GTR rate parameters, base frequencies, alpha, and alignment lengths all come from real data
- **Consistent across replicates**: The same parameter set is used for a given gene tree across all replicates
- **Handles extreme alpha values**: When alpha > 3.0, uses uniform rate distribution (no +G) as agreed with supervisors

## Empirical Parameter Distributions

### Source Data
- **Zhao_2021**: 982 alignments, mean length 968 bp, mean 14.6 sequences
- **Ren_2024**: 727 alignments, mean length 888 bp, mean 78.4 sequences
- **Morales_Briones_2021**: 1000 alignments, mean length 376 bp, mean 44.3 sequences

### Combined Statistics
- **Total genes**: 2,709 parameter sets
- **Alignment length**: 147-6,487 bp (median 609 bp, mean 728 bp)
- **Alpha parameter**: median 0.54, mean 18.47
  - 90 genes (3.3%) have α > 3.0 (treated as uniform rate)
  - Maximum α = 999.0

### What Gets Sampled
For each gene tree, one parameter set is randomly sampled containing:
- **GTR rate parameters**: AC, AG, AT, CG, CT, GT (relative to GT=1.0)
- **Base frequencies**: π_A, π_C, π_G, π_T
- **Alpha parameter**: Gamma shape for rate heterogeneity
- **Alignment length**: Number of base pairs to simulate

## Prerequisites

### 1. Ensure GTR Parameters Are Extracted

The empirical parameters must be extracted and saved to a pickle file:

```bash
# Check if the pickle file exists
ls -lh /groups/itay_mayrose/tomulanovski/gene2net/simulations/distributions/gtr_parameters_all.pkl

# If not, run the extraction script
cd /groups/itay_mayrose/tomulanovski/gene2net
python simulations/scripts/sequence_evolution/filter_and_extract_gtr.py
```

### 2. SimPhy Simulations Must Be Complete

Before simulating sequences, you need gene trees from SimPhy:

```bash
# Check if SimPhy completed for your configuration
python simulations/scripts/check_pipeline_status.py ils_low_10M --step simphy
```

That's it! The sequence simulation script automatically reads gene trees from SimPhy output (whether in single-batch or multi-batch mode).

## Usage

### Quick Start: Simulate for All Networks

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# Submit sequence simulation for all 21 networks
./submit_sequences.sh ils_low_10M
```

### Test on a Single Network

```bash
# Submit for just one network to test
./submit_sequences.sh ils_low_10M Ding_2023
```

### Submit for Different Configurations

```bash
# Low ILS
./submit_sequences.sh ils_low_10M

# Medium ILS
./submit_sequences.sh ils_med_10M

# High ILS
./submit_sequences.sh ils_high_10M

# With duplication/loss
./submit_sequences.sh ils_low_dup_loss_low
```

## Script Details

### 1. `sample_gtr_parameters.py`

**Location**: `simulations/scripts/sequence_evolution/sample_gtr_parameters.py`

**Purpose**: Randomly samples one GTR+Gamma parameter set from the empirical distribution.

**Usage**:
```bash
# Get bash variables (default)
python sample_gtr_parameters.py gtr_parameters_all.pkl

# Output:
# IQTREE_MODEL='GTR{1.234/4.567/...}+FU{0.25/0.25/...}+G{0.543}'
# ALIGNMENT_LENGTH=609
# ALPHA=0.543000

# Get just the model string
python sample_gtr_parameters.py gtr_parameters_all.pkl --output-format iqtree

# Get JSON format
python sample_gtr_parameters.py gtr_parameters_all.pkl --output-format json

# Use specific random seed
python sample_gtr_parameters.py gtr_parameters_all.pkl --seed 12345
```

**Model String Format**:
- For α ≤ 3.0: `GTR{rates}+FU{freqs}+G{alpha}`
- For α > 3.0: `GTR{rates}+FU{freqs}` (uniform rate, no +G)

### 2. `simulate_sequences_1_dataset.sh`

**Location**: `simulations/jobs/simulate_sequences_1_dataset.sh`

**Purpose**: SLURM array job that simulates sequences for one network.

**Key Features**:
- Array tasks 1-1000 (one per gene tree)
- Samples parameters once per gene tree
- Uses same parameters across all replicates
- Creates unique seeds for each replicate's simulation

**Direct Usage** (not recommended, use submit_sequences.sh instead):
```bash
sbatch --job-name=alisim_Ding_2023_ils_low_10M \
    simulate_sequences_1_dataset.sh Ding_2023 ils_low_10M
```

**Resource Requirements**:
- Time: 30 minutes per array task
- Memory: 2GB
- CPUs: 1

### 3. `simulate_sequences_all_networks.sh`

**Location**: `simulations/jobs/simulate_sequences_all_networks.sh`

**Purpose**: Wrapper that submits jobs for all 21 networks.

**Usage**:
```bash
# Submit all networks
./simulate_sequences_all_networks.sh ils_low_10M

# Shows progress and summary
```

### 4. `submit_sequences.sh` (Recommended)

**Location**: `simulations/jobs/submit_sequences.sh`

**Purpose**: Master submission script with validation and flexible options.

**Advantages**:
- Validates that GTR parameters exist
- Validates that scripts exist
- Provides detailed summary
- Supports single network or all networks
- Consistent with `submit_simphy.sh` interface

**Usage**:
```bash
# All networks (recommended)
./submit_sequences.sh ils_low_10M

# Single network for testing
./submit_sequences.sh ils_low_10M Ding_2023

# Different configurations
./submit_sequences.sh ils_high_10M
./submit_sequences.sh ils_low_dup_loss_low
```

## Workflow Example

### Complete Pipeline from SimPhy to Sequences

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/jobs

# 1. Run SimPhy simulations
./submit_simphy.sh ils_low_10M 10M 0 0 200000

# 2. Wait for SimPhy to complete, then verify
cd ../scripts
python check_pipeline_status.py ils_low_10M --step simphy

# 3. If SimPhy succeeded, simulate sequences
cd ../jobs
./submit_sequences.sh ils_low_10M

# 4. Monitor sequence simulation jobs
squeue -u $USER | grep alisim_ils_low_10M

# 5. Check logs for any errors
tail -f ../logs/alisim_ils_low_10M_Ding_2023_*.err
```

## Output Structure

After running sequence simulations for a configuration:

```
simulations/
└── <NetworkName>/
    └── data/
        └── <configuration>/           # e.g., ils_low_10M
            ├── replicate_1/
            │   └── 1/
            │       ├── g_trees0001.trees
            │       ├── g_trees0002.trees
            │       ├── ...
            │       └── alignments/
            │           ├── alignment_0001.phy
            │           ├── alignment_0002.phy
            │           └── ...
            ├── replicate_2/
            │   └── 1/
            │       └── alignments/
            │           └── ...
            └── replicate_N/
                └── 1/
                    └── alignments/
                        └── ...
```

Each `alignment_XXXX.phy` file contains sequences simulated using:
- Gene tree from `g_treesXXXX.trees`
- Empirically-sampled GTR+Gamma parameters
- Empirically-sampled alignment length

## Monitoring and Troubleshooting

### Check Job Status

```bash
# See running jobs
squeue -u $USER | grep alisim

# Check completed jobs
sacct -X --format=JobName,State,ExitCode | grep alisim_ils_low_10M

# See detailed job info
sacct -j <JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

### View Logs

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs

# Check for errors
grep -l "ERROR" alisim_ils_low_10M_*.err

# View a specific log
less alisim_ils_low_10M_Ding_2023_1.out
```

### Common Issues

**Issue**: "GTR parameters file not found"
```bash
# Solution: Extract parameters first
cd /groups/itay_mayrose/tomulanovski/gene2net
python simulations/scripts/sequence_evolution/filter_and_extract_gtr.py
```

**Issue**: "No replicate directories found"
```bash
# Solution: Run SimPhy first
cd simulations/jobs
./submit_simphy.sh ils_low_10M
```

**Issue**: "Gene tree file not found"
```bash
# Solution: Check if SimPhy completed successfully
python simulations/scripts/check_pipeline_status.py ils_low_10M --step simphy
```

## Parameter Sampling Details

### Seed Strategy

1. **Parameter sampling seed**: `SLURM_ARRAY_TASK_ID * 1000 + RANDOM`
   - Ensures different genes get different parameter sets
   - Reproducible if RANDOM is seeded

2. **Simulation seed**: `SLURM_ARRAY_TASK_ID * 10000 + replicate * 1000 + RANDOM % 1000`
   - Different for each replicate
   - Ensures different sequence realizations

### Why Same Parameters Across Replicates?

For a given gene tree number (e.g., gene tree 42):
- **Same GTR+Gamma parameters** across all replicates
- **Different simulation seeds** create different sequence realizations

This design:
- Represents biological reality (parameters are gene-specific, not replicate-specific)
- Maintains variance through different sequence realizations
- Simplifies downstream analysis

## Verification

### Check That Sequences Were Generated

**Recommended: Use the pipeline status checker**

```bash
# Check sequence alignments for all networks
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts
python check_pipeline_status.py ils_low_10M --step sequences

# Check only specific network with verbose output
python check_pipeline_status.py ils_low_10M --step sequences --verbose
```

**Manual verification (alternative)**

```bash
# For a specific network and configuration
NETWORK="Ding_2023"
CONFIG="ils_low_10M"
BASE="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

# Count alignments in replicate 1
ls ${BASE}/${NETWORK}/data/${CONFIG}/replicate_1/1/alignments/*.phy | wc -l
# Should be 1000 (or however many gene trees you have)

# Check all replicates
for rep in {1..5}; do
    count=$(ls ${BASE}/${NETWORK}/data/${CONFIG}/replicate_${rep}/1/alignments/*.phy 2>/dev/null | wc -l)
    echo "Replicate $rep: $count alignments"
done
```

### Inspect a Sample Alignment

```bash
head -20 ${BASE}/${NETWORK}/data/${CONFIG}/replicate_1/1/alignments/alignment_0001.phy
```

## Next Steps

After sequences are generated, you can:

1. **Run phylogenetic network inference methods**:
   ```bash
   cd simulations/jobs
   ./submit_all_methods.sh ils_low_10M
   ```

2. **Validate pipeline status**:
   ```bash
   cd simulations/scripts
   python check_pipeline_status.py ils_low_10M --step run
   ```

3. **Analyze results**:
   ```bash
   python postprocess_results.py ils_low_10M
   python run_full_summary.py ils_low_10M
   ```

## Technical Notes

### IQ-TREE Model String Format

The model strings use IQ-TREE's parameter specification format:

```
GTR{rAC/rAG/rAT/rCG/rCT/rGT}+FU{piA/piC/piG/piT}+G{alpha}
```

Where:
- `GTR{...}`: Six rate parameters (AC, AG, AT, CG, CT, GT)
- `+FU{...}`: User-defined base frequencies (A, C, G, T)
- `+G{alpha}`: Gamma distribution with shape parameter alpha
- For α > 3.0: Omit `+G` to use uniform rate distribution

### Decision on Alpha > 3.0

Following discussion with supervisors, genes with α > 3.0 (90 genes, 3.3% of data) are treated as having uniform rate distribution. This is biologically reasonable as high alpha indicates nearly uniform rates across sites.

## Summary

The improved sequence simulation pipeline:

1. ✅ Uses empirically-inferred GTR+Gamma parameters from 2,709 real genes
2. ✅ Samples alignment lengths from empirical distribution (147-6,487 bp)
3. ✅ Maintains biological realism in parameter combinations
4. ✅ Handles extreme alpha values appropriately (α > 3.0 → uniform)
5. ✅ Works across all configurations (ILS, dup/loss, etc.)
6. ✅ Integrates seamlessly with existing pipeline
7. ✅ Provides comprehensive logging and error checking

For questions or issues, consult:
- This guide
- Script comments and usage messages
- Pipeline status checker: `check_pipeline_status.py`
