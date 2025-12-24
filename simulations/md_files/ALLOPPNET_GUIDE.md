# AlloppNET Guide

**Complete guide to running AlloppNET phylogenetic network inference on simulation data**

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Quick Start](#quick-start)
4. [Detailed Workflow](#detailed-workflow)
5. [Network Compatibility](#network-compatibility)
6. [Directory Structure](#directory-structure)
7. [Troubleshooting](#troubleshooting)
8. [Output Files](#output-files)
9. [Post-Processing](#post-processing)
10. [Runtime Considerations](#runtime-considerations)

---

## Overview

AlloppNET is a Bayesian MCMC method for inferring phylogenetic networks in allopolyploid species. It uses BEAST 1.8+ to infer species networks from multi-locus sequence alignments, designed specifically for diploids and allotetraploids.

**Key Features:**
- Infers species networks with hybridization
- Handles diploid (2n) and allotetraploid (4n) species
- Uses BEAST for MCMC sampling
- 100M iterations (~5 days per replicate)
- Only runs on 8 compatible networks (max 2 copies in topology)

**Method Information:**
- **Authors:** Jones et al.
- **Software:** BEAST 1.8+, AlloppDT R scripts
- **Type:** Bayesian MCMC network inference
- **Ploidy Support:** Diploid (2n), Allotetraploid (4n)

---

## Prerequisites

### 1. Completed SimPhy Simulations

AlloppNET requires completed sequence alignments from SimPhy:

```bash
# Verify simulations are complete
python ../scripts/check_pipeline_status.py conf_ils_low_10M --step sequences
```

Expected output structure:
```
simulations/<Network>/data/<config>/replicate_N/1/alignments/
├── alignment_0001.phy ... alignment_1000.phy
```

### 2. AlloppDT R Scripts (REQUIRED)

**CRITICAL:** Before running AlloppNET, you must copy AlloppDT scripts from the cluster:

```bash
# On the cluster
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts/alloppnet/

# Copy AlloppDT R scripts (required for BEAST XML generation)
cp "/groups/itay_mayrose/tomulanovski/gene2net/jones_examples/2013-05-15-manual (2)/AlloppDT_5beastxml_toplevel.r" .
cp "/groups/itay_mayrose/tomulanovski/gene2net/jones_examples/2013-05-15-manual (2)/AlloppDT_6beastxml_bits.r" .

# Verify files exist
ls -l AlloppDT_*.r
```

You should see:
- `AlloppDT_5beastxml_toplevel.r`
- `AlloppDT_6beastxml_bits.r`

### 3. Conda Environment

AlloppNET requires a dedicated conda environment with BEAST:

```bash
# Create environment (if not exists)
conda create -n alloppnet -c bioconda beast=1.8.4

# Activate and verify
conda activate alloppnet
beast -help
treeannotator -help
```

---

## Quick Start

### Two-Step Workflow (Recommended)

```bash
cd simulations/jobs

# Step 1: Prep (fast, ~minutes)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --prep-only

# Wait for prep to complete, then validate
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step prep

# Step 2: Run (slow, ~5 days)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --run-only

# Validate after completion
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step run
```

### Why Two Steps?

**PREP (fast):**
- Converts alignments, analyzes copy distributions (kernel smoothing)
- Generates BEAST XML - validates before expensive BEAST run
- Runtime: ~minutes per network

**RUN (slow):**
- Runs BEAST MCMC (~5 days per network-replicate)
- Summarizes and post-processes results

### Submit Specific Replicates

```bash
# Prep only replicate 1
./submit_alloppnet_pipeline.sh conf_ils_low_10M --replicates 1 --prep-only

# Run replicates 1,3,5 only (after prep validated)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --replicates 1,3,5 --run-only
```

### Monitor Progress

```bash
# Check SLURM queue
squeue -u $USER | grep alloppnet

# Check prep logs
tail -f ../logs/alloppnet_prep_conf_ils_low_10M_rep1_*.out

# Check run logs
tail -f ../logs/alloppnet_run_conf_ils_low_10M_rep1_*.out

# Validate prep outputs
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step prep

# Validate final results
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step run
```

---

## Detailed Workflow

AlloppNET runs as a two-step pipeline:
- **PREP step** (fast, ~minutes): Generate inputs and BEAST XML
- **RUN step** (slow, ~5 days): Run BEAST and post-process

### Step 1: Prepare Input

**Script:** `prepare_alloppnet_input.py`

**What it does:**
1. Analyzes copy number distributions per taxon across all 1000 alignments
2. Uses kernel smoothing to determine representative copy numbers (robust to outliers)
3. Generates `ploidy_level.json` based on representative copy counts:
   - 1 copy → diploid (ploidy=2)
   - 2+ copies → tetraploid (ploidy=4)
4. Converts PHY alignments to NEXUS format
5. Creates `taxa_table.txt` with homeolog pairing

**Inputs:**
- `data/<config>/replicate_N/1/alignments/alignment_*.phy` (1000 files)

**Outputs:**
- `processed/<config>/alloppnet_input/replicate_N/alignment_*.nex` (1000 files)
- `processed/<config>/alloppnet_input/replicate_N/ploidy_level.json`
- `processed/<config>/alloppnet_input/replicate_N/taxa_table.txt`

**Sequence Name Format (SimPhy):**
SimPhy gene tree leaves follow the naming scheme `species_locusid_individualid`:
- `species`: The species/taxon name
- `locusid`: The locus/gene number (corresponds to alignment file number)
- `individualid`: The individual/copy number within that locus
- Examples: `Galeopsisspeciosa_0_0`, `Lamiummoschatum3_1_1`, `Complex_Taxon_Name_5_2`
- Conversion: `species_locusid_individualid` → `species_individualid` for NEXUS
- Handles complex taxon names with underscores correctly

**Ploidy Assignment (Kernel Smoothing):**
- Analyzes copy number distribution across all 1000 genes
- Uses triangular kernel smoothing to find representative copy number
- Robust to outliers from occasional gene duplication or loss
- Example with gene loss:
  - Distribution: {1 copy: 950 genes, 2 copies: 50 genes}
  - Representative: 1 copy → diploid (ploidy=2)
- Example with stable tetraploidy:
  - Distribution: {2 copies: 980 genes, 1 copy: 20 genes}
  - Representative: 2 copies → tetraploid (ploidy=4)

### Step 2: Generate BEAST XML

**Script:** `generate_beast_xml.r`

**What it does:**
1. Sources AlloppDT R scripts
2. Reads NEXUS alignments and taxa_table.txt
3. Generates BEAST XML configuration

**BEAST Configuration:**
- Chain length: 100M iterations
- Screen logging: every 10K iterations
- Parameter sampling: every 1K iterations
- Gene tree sampling: every 1K iterations
- Species network sampling: every 1K iterations

**Output:**
- `results/<config>/alloppnet/replicate_N/alloppnet.XML`

### Step 3: Run BEAST

**Runtime:** ~5 days (120 hours) per network-replicate

**What it does:**
1. Runs BEAST MCMC inference
2. Samples species networks and gene trees
3. Logs parameters and diagnostics

**Outputs:**
- `sampledmultrees.txt` - Species networks (100K trees, 1 per 1K iterations)
- `sampledgtrees1.txt ... sampledgtrees1000.txt` - Gene trees
- `sampledparams.txt` - MCMC parameters
- `DBUGTUNE.txt` - Debug/tuning information

**Note:** Post-processing (TreeAnnotator + copy number removal) is handled separately by `postprocess_results.py` after BEAST completes or times out.

**Monitoring:**
```bash
# Check BEAST progress (look for "STATE_XXXXX")
grep "STATE_" results/<Network>/<config>/alloppnet/replicate_1/sampledmultrees.txt | tail

# Check for errors
tail -100 ../logs/alloppnet_*_rep1_*.err
```

### Step 4: Post-Process Results

**Script:** `submit_alloppnet_postprocess.sh` (sbatch job)

**What it does:**
1. Runs TreeAnnotator on `sampledmultrees.txt` (10% burnin, mean heights)
2. Removes `_0` and `_1` suffixes from tip labels using `remove_copy_numbers.py`
3. Produces final MUL-tree for comparison

**Burnin Calculation:**
```
Total trees: 100,000 (STATE_1000 through STATE_100000000)
Usable trees: 99,999 (exclude STATE_0)
Burnin: 9,999 trees (10%)
```

**Output:**
- `alloppnet_result.tre` - Final consensus tree (ready for comparison)

**Usage:**
```bash
# Option 1: Use postprocess_results.py (automatically submits sbatch jobs on cluster)
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods alloppnet

# Option 2: Manually submit sbatch jobs
cd simulations/jobs
./submit_alloppnet_postprocess.sh conf_ils_low_10M

# Post-process specific replicates (manual only)
./submit_alloppnet_postprocess.sh conf_ils_low_10M --replicates 1,3,5

# Dry run (preview commands)
./submit_alloppnet_postprocess.sh conf_ils_low_10M --dry-run
```

**Why sbatch?**
- TreeAnnotator processing 100K+ trees is CPU/memory intensive
- Requires proper resource allocation (16GB RAM, 2 hours)
- Runs in parallel (8 networks × 5 replicates = 40 jobs)
- Better error handling and logging

**Note:** `postprocess_results.py` automatically detects if you're on the cluster and submits sbatch jobs for AlloppNET. If you're not on the cluster, it will show instructions to run manually.

**Benefits:**
- Can run post-processing even if BEAST job times out
- Proper resource allocation via SLURM
- Parallel processing (faster than sequential)
- Can re-run post-processing without re-running BEAST

---

## Network Compatibility

**AlloppNET only runs on 8 of 21 networks** - those with max 2 copies in the network topology:

| Compatible Networks |
|---------------------|
| Bendiksby_2011      |
| Ding_2023           |
| Koenen_2020         |
| Liu_2023            |
| Shahrestani_2015    |
| Wisecaver_2023      |
| Wu_2015             |
| Zhao_2021           |

**Why only 8 networks?**
- AlloppNET is designed for diploid/allotetraploid species
- Networks with >2 copies (e.g., hexaploids, octoploids) are not supported
- These 8 networks have maximum 2 copies of each locus in the topology

**How to check compatibility:**
```bash
# View network metadata
cat ../networks/mul_tree_final_stats.csv | grep -E "Bendiksby|Ding|Koenen|Liu|Shahrestani|Wisecaver|Wu|Zhao"
```

---

## Directory Structure

```
simulations/<Network>/
├── data/<config>/
│   └── replicate_N/1/alignments/
│       └── alignment_*.phy (1000 files)
│           # Input: PHY alignments from SimPhy
│
├── processed/<config>/alloppnet_input/
│   └── replicate_N/
│       ├── alignment_0001.nex ... alignment_1000.nex
│       ├── ploidy_level.json
│       └── taxa_table.txt
│
└── results/<config>/alloppnet/
    └── replicate_N/
        ├── alloppnet.XML
        ├── sampledmultrees.txt
        ├── sampledgtrees1.txt ... sampledgtrees1000.txt
        ├── sampledparams.txt
        ├── DBUGTUNE.txt
        └── alloppnet_result.tre (created by postprocess_results.py)
```

**Key Files:**

| File | Purpose |
|------|---------|
| `ploidy_level.json` | Ploidy assignment for each taxon (2 or 4) |
| `taxa_table.txt` | Maps taxa to individuals and genomes (A/B) |
| `alloppnet.XML` | BEAST configuration file |
| `sampledmultrees.txt` | Species networks from MCMC (input for post-processing) |
| `alloppnet_result.tre` | **Final output for comparison** (created by postprocess_results.py) |

---

## Troubleshooting

### BEAST Fails to Start

**Symptoms:**
```
ERROR: Could not activate alloppnet environment
ERROR: beast command not found
```

**Solutions:**
```bash
# Recreate conda environment
conda create -n alloppnet -c bioconda beast=1.8.4 -y

# Verify installation
conda activate alloppnet
beast -help
```

### Missing AlloppDT Scripts

**Symptoms:**
```
ERROR: AlloppDT_5beastxml_toplevel.r not found
```

**Solution:**
Follow [Prerequisites](#2-alloppdt-r-scripts-required) to copy scripts from cluster.

### Convergence Issues

**Symptoms:**
- Very low effective sample sizes (ESS < 200)
- Poor mixing in parameters
- High autocorrelation

**Solutions:**
1. Check `DBUGTUNE.txt` for warnings
2. Consider increasing chain length (requires editing `generate_beast_xml.r`)
3. Check for incompatible ploidy assignments in `ploidy_level.json`

### Memory Errors

**Symptoms:**
```
ERROR: Out of memory
slurmstepd: error: Exceeded job memory limit
```

**Solutions:**
```bash
# Edit alloppnet_array.sh to increase memory
#SBATCH --mem=64G   # Increase from 32G to 64G

# Resubmit job
sbatch --array=N alloppnet_array.sh conf_ils_low_10M 1
```

### Empty or Missing Output

**Symptoms:**
- `alloppnet_result.tre` not created
- `sampledmultrees.txt` is empty or has only STATE_0

**Diagnosis:**
```bash
# Check SLURM logs for errors
tail -200 ../logs/alloppnet_*_repN_*.err

# Check if BEAST completed (may have timed out)
grep "BEAST run complete" ../logs/alloppnet_*_repN_*.out

# Count trees in MCMC output
grep -c "tree STATE_" results/<Network>/<config>/alloppnet/replicate_N/sampledmultrees.txt

# If sampledmultrees.txt exists, try post-processing manually
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods alloppnet
```

### Taxa Table Format Errors

**Symptoms:**
```
ERROR: Taxa table format invalid
ERROR: Missing genome assignments
```

**Solution:**
```bash
# Verify taxa_table.txt format
head processed/<Network>/<config>/alloppnet_input/replicate_1/taxa_table.txt

# Expected format:
# ID species individual genome
# Taxon1_0 Taxon1 Taxon1_ind1 A
# Taxon1_1 Taxon1 Taxon1_ind1 B
```

---

## Output Files

### ploidy_level.json

Maps taxon names to ploidy levels:

```json
{
  "Galeopsisspeciosa": 2,
  "Lamiumalbum": 4,
  "Lamiummoschatum3": 4
}
```

- `2` = diploid (1 copy in alignments)
- `4` = tetraploid (2+ copies in alignments)

### taxa_table.txt

Maps sequence IDs to species, individuals, and genomes:

```
ID species individual genome
Galeopsisspeciosa_0 Galeopsisspeciosa Galeopsisspeciosa_ind1 A
Lamiumalbum_0 Lamiumalbum Lamiumalbum_ind1 A
Lamiumalbum_1 Lamiumalbum Lamiumalbum_ind1 B
```

**Format:**
- **Diploid:** Single row per taxon (genome A)
- **Tetraploid:** Pairs of rows (genomes A/B for each individual)
- Homeologs paired: _0 with _1, _2 with _3, etc.

### sampledmultrees.txt

BEAST MCMC output with species networks:

```
#NEXUS

BEGIN TREES;
tree STATE_0 = [...]
tree STATE_1000 = [...]
tree STATE_2000 = [...]
...
tree STATE_100000000 = [...]
END;
```

- 100,000 trees (sampled every 1K iterations)
- Use for convergence diagnostics
- TreeAnnotator summarizes to consensus

### alloppnet_final.tre

**Final consensus tree ready for comparison:**

```
((TaxonA:0.5,TaxonB:0.5):0.3,(TaxonC:0.4,TaxonD:0.4):0.4);
```

- Newick format
- Copy number suffixes removed
- Mean node heights from TreeAnnotator
- Can be compared with other methods (GRAMPA, Polyphest, etc.)

---

## Post-Processing

### Validate Results

```bash
# Check input preparation
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step prep

# Check BEAST outputs
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step run --verbose

# Post-process results (if not done automatically)
python ../scripts/postprocess_results.py conf_ils_low_10M --methods alloppnet
```

**Expected output:**
```
AlloppNET Input Validation (prep):
  Total: 40 (8 networks × 5 replicates)
  Success: 40 (100.0%)
  Missing files: 0

AlloppNET Output Validation (run):
  Total: 40
  Success: 40 (100.0%)
  Required files present:
    - alloppnet_result.tre
    - sampledmultrees.txt
    - alloppnet.XML
```

### Compare with True Networks

```bash
# Run network comparison (if post-processing pipeline is set up)
python ../scripts/compare_nets.py \
    --true-network ../networks/<Network>.tre \
    --inferred results/<Network>/<config>/alloppnet/replicate_1/alloppnet_final.tre \
    --output comparison_results.csv
```

### Extract Statistics

```bash
# Parse BEAST log files for parameter estimates
grep "posterior" results/<Network>/<config>/alloppnet/replicate_1/sampledparams.txt

# Calculate effective sample sizes
# (Requires Tracer or custom script)
```

---

## Runtime Considerations

### Computational Resources

**Per network-replicate:**
- **Time:** ~5 days (120 hours)
- **Memory:** 32GB (increase to 64GB if needed)
- **CPUs:** 1 (BEAST is single-threaded)
- **Disk:** ~10GB (MCMC output files)

**Total for all 8 networks × 5 replicates (40 jobs):**
- **Sequential runtime:** 200 days
- **Parallel runtime (10 simultaneous jobs):** ~20 days
- **Storage:** ~400GB

### SLURM Job Configuration

```bash
#SBATCH --time=120:00:00       # 5 days
#SBATCH --mem=32G              # 32GB memory
#SBATCH --cpus-per-task=1      # 1 CPU
#SBATCH --partition=itaym-pool
```

### Optimization Tips

1. **Parallel Submission:**
   - Submit all replicates simultaneously
   - Each replicate is an array job (8 networks)
   - Max 40 parallel jobs (cluster permitting)

2. **Storage Management:**
   - Gene trees (~1GB each) can be deleted after consensus if needed
   - Keep `sampledmultrees.txt` and `alloppnet_final.tre`

3. **Monitoring:**
   - Check progress every 12-24 hours
   - Look for convergence issues early (ESS, autocorrelation)
   - Be prepared to restart if BEAST crashes

### Restart Strategies

BEAST supports checkpointing, but it's not enabled by default in AlloppDT scripts. If a job times out:

1. **Increase time limit** (if cluster allows >5 days)
2. **Reduce chain length** (edit `generate_beast_xml.r`)
3. **Resume from checkpoint** (requires manual BEAST XML editing)

---

## Integration with Pipeline

### Use with submit_all_methods.sh

AlloppNET is integrated into the master pipeline orchestrator:

```bash
# Prep all methods including AlloppNET
cd simulations/jobs
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest,alloppnet --prep-only

# Run all methods after validation
./submit_all_methods.sh conf_ils_low_10M --methods grampa,polyphest,alloppnet --run-only

# Submit only AlloppNET (both steps)
./submit_all_methods.sh conf_ils_low_10M --methods alloppnet
```

**Note:** AlloppNET now supports `--prep-only` and `--run-only` flags like other methods!

### Validation in Pipeline

```bash
# Check AlloppNET prep outputs
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step prep

# Check AlloppNET final results
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step run --verbose

# Check all methods including AlloppNET
python ../scripts/check_pipeline_status.py conf_ils_low_10M --step run
```

---

## References

- **AlloppNET Paper:** Jones et al. (method description)
- **BEAST Documentation:** http://beast.community/
- **AlloppDT Scripts:** Original implementation by Jones lab

---

## Summary

**Quick Reference:**

```bash
# 1. Copy AlloppDT scripts (ONE TIME SETUP)
cp ".../AlloppDT_5beastxml_toplevel.r" simulations/scripts/alloppnet/
cp ".../AlloppDT_6beastxml_bits.r" simulations/scripts/alloppnet/

# 2. Prep step (fast, ~minutes)
cd simulations/jobs
./submit_alloppnet_pipeline.sh conf_ils_low_10M --prep-only

# 3. Validate prep outputs
python ../scripts/check_pipeline_status.py conf_ils_low_10M --method alloppnet --step prep

# 4. Run step (slow, ~5 days)
./submit_alloppnet_pipeline.sh conf_ils_low_10M --run-only

# 5. Monitor
squeue -u $USER | grep alloppnet
tail -f ../logs/alloppnet_run_*.out

# 6. Post-process results (after BEAST completes or times out)
cd ../jobs
./submit_alloppnet_postprocess.sh conf_ils_low_10M

# 7. Wait for post-processing jobs to complete, then validate
cd ../scripts
python check_pipeline_status.py conf_ils_low_10M --method alloppnet --step run --verbose
```

**Expected timeline (three-step workflow):**
- Day 0: Submit prep jobs (~minutes to complete)
- Day 0: Validate prep → Submit run jobs
- Day 5: First BEAST runs complete (or timeout)
- Day 5+: Post-process results (can run even if BEAST timed out)
- Day 20-25: All BEAST runs complete (if running 10 parallel)

**Benefits of separated workflow:**
- Validate prep outputs before expensive BEAST run
- Catch errors early (ploidy detection, XML generation)
- Post-processing can run independently (even if BEAST times out)
- Consistent with other methods in the pipeline
- Can re-run post-processing without re-running BEAST

**Final outputs:**
- 40 consensus trees (8 networks × 5 replicates)
- Location: `results/<Network>/<config>/alloppnet/replicate_N/alloppnet_result.tre`
