# GRANDMA_SPLIT Setup Guide

## Overview

GRANDMA_SPLIT has been integrated into your simulation pipeline. It uses the **same inputs as GRAMPA** (gene trees and species tree), so it leverages the existing GRAMPA prep and ASTRAL steps.

## Files Created

1. **`run_grandma_split.sh`** - Main SLURM job script
   - Runs 105 array jobs (21 networks × 5 replicates)
   - Uses inputs from `grampa_input/` directory
   - Outputs to `grandma_split/` directory
   - Calls: `/groups/itay_mayrose/ronenshtein/Grampa_revamp/run.py` with `-m split`

2. **`submit_grandma_split_pipeline.sh`** - Submission wrapper
   - Handles prep, ASTRAL, and GRANDMA_SPLIT submission with dependencies
   - Automatically integrates with `submit_all_methods.sh`

3. **Updated `submit_all_methods.sh`**
   - Added `grandma_split` to available methods list

## Usage

### Option 1: Run via submit_all_methods.sh (Recommended)

```bash
# Run all methods including grandma_split
cd simulations/jobs
./submit_all_methods.sh conf_ils_low_10M

# Run only grandma_split
./submit_all_methods.sh conf_ils_low_10M --methods grandma_split

# Run multiple specific methods
./submit_all_methods.sh conf_ils_low_10M --methods grampa,grandma_split

# Dry run to preview
./submit_all_methods.sh conf_ils_low_10M --methods grandma_split --dry-run
```

### Option 2: Run via submit_grandma_split_pipeline.sh

```bash
# Full pipeline (prep + ASTRAL + grandma_split)
./submit_grandma_split_pipeline.sh conf_ils_low_10M

# Only prep and ASTRAL
./submit_grandma_split_pipeline.sh conf_ils_low_10M --prep-only

# Only grandma_split (assumes prep and ASTRAL already done)
./submit_grandma_split_pipeline.sh conf_ils_low_10M --run-only

# Custom number of replicates
./submit_grandma_split_pipeline.sh conf_ils_low_10M --replicates 3

# Dry run
./submit_grandma_split_pipeline.sh conf_ils_low_10M --dry-run
```

### Option 3: Direct SLURM submission

```bash
# Submit directly (assumes prep and ASTRAL done)
sbatch --export=CONFIG=conf_ils_low_10M simulations/jobs/run_grandma_split.sh

# Custom replicates (63 jobs for 3 replicates)
sbatch --array=1-63 --export=CONFIG=conf_ils_low_10M,NUM_REPLICATES=3 run_grandma_split.sh
```

## Input/Output Structure

### Input (same as GRAMPA)
```
simulations/
└── {network}/
    └── processed/
        └── {config}/
            └── grampa_input/
                └── replicate_{N}/
                    ├── grampa_trees.tre  # Gene trees
                    └── species.tre       # Species tree
```

### Output (grandma_split specific)
```
simulations/
└── {network}/
    └── results/
        └── {config}/
            └── grandma_split/
                └── replicate_{N}/
                    └── [GRANDMA_SPLIT output files]
```

## Dependencies

GRANDMA_SPLIT requires:
1. **GRAMPA prep** - Process gene trees (clean, fix substrings, add copy numbers)
2. **ASTRAL** - Create species tree from gene trees

These are automatically handled by the submission scripts with job dependencies.

## Monitoring

```bash
# Check all your jobs
squeue -u $USER

# Check grandma_split specific logs
tail -f /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/run_grandma_split_*.out

# Check pipeline status
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run
```

## Next Steps After Job Completion

1. **Validate results:**
   ```bash
   python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run --verbose
   ```

2. **Post-process results:**
   ```bash
   python simulations/scripts/postprocess_results.py conf_ils_low_10M
   ```

3. **Run summary analysis:**
   ```bash
   python simulations/scripts/run_full_summary.py conf_ils_low_10M
   ```

## Before Pushing to Cluster

When you push to GitHub and pull on the cluster, make sure the scripts are executable:

```bash
chmod +x simulations/jobs/run_grandma_split.sh
chmod +x simulations/jobs/submit_grandma_split_pipeline.sh
```

Or do it all at once:
```bash
chmod +x simulations/jobs/*.sh
```

## Troubleshooting

### Issue: Input files not found
**Solution:** Make sure GRAMPA prep and ASTRAL have completed first:
```bash
./submit_grandma_split_pipeline.sh conf_ils_low_10M --prep-only
# Wait for completion, then:
./submit_grandma_split_pipeline.sh conf_ils_low_10M --run-only
```

### Issue: GRANDMA_SPLIT script not found
**Solution:** Verify the script path exists:
```bash
ls -l /groups/itay_mayrose/ronenshtein/Grampa_revamp/run.py
```

### Issue: Python/conda errors
**Solution:** Ensure gene2net conda environment has all required dependencies.
