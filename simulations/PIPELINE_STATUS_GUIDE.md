# Pipeline Status Checker - Complete Guide

## Overview

`check_pipeline_status.py` validates all stages of your phylogenetic network inference pipeline:
- ✅ SimPhy simulations (1000 trees per replicate)
- ✅ Input file preparation for GRAMPA, Polyphest, MPSUGAR, PADRE
- ✅ Output file generation after running methods
- ✅ Support for multiple configurations and parameter variations

## Your Setup

- **Configurations:** `conf_ils_low_10M`, `conf_ils_medium_10M`, `conf_ils_high_10M`
- **Networks:** 21 networks
- **Replicates:** 5 per network
- **Total checks per step:** 21 × 5 = 105 combinations
- **Polyphest runs:** 3 percentiles (50, 70, 90)

---

## Quick Reference Commands

### Check Everything for a Configuration
```bash
python scripts/check_pipeline_status.py conf_ils_low_10M
```

### Check SimPhy Simulations (1000 trees per replicate)
```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --step simphy
```

### Check Input Preparation (Before Running Methods)
```bash
# All methods
python scripts/check_pipeline_status.py conf_ils_low_10M --step prep

# Specific method
python scripts/check_pipeline_status.py conf_ils_low_10M --method grampa --step prep
```

### Check Output Generation (After Running Methods)
```bash
# GRAMPA
python scripts/check_pipeline_status.py conf_ils_low_10M --method grampa --step run

# Polyphest (MUST specify percentile!)
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 50

# MPSUGAR
python scripts/check_pipeline_status.py conf_ils_low_10M --method mpsugar --step run

# PADRE
python scripts/check_pipeline_status.py conf_ils_low_10M --method padre --step run
```

### Verbose Mode (Show All Details)
```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --verbose
```

### Export to CSV
```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --export results.csv
```

---

## Step-by-Step Usage

### 1. After SimPhy Completes

Verify all gene trees were generated:

```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --step simphy
```

**What it checks:**
- Each of 21 networks
- Each of 5 replicates per network
- Exactly 1000 gene tree files (`g_*`) per replicate
- Total: 105 combinations

**Expected output:**
```
Checking: 21 networks × 5 replicates = 105 combinations

Success:  105 / 105 (100.0%)
```

### 2. Before Running Methods

Check that all input files were prepared correctly:

```bash
# Check all methods at once
python scripts/check_pipeline_status.py conf_ils_low_10M --step prep

# Or check each method individually
python scripts/check_pipeline_status.py conf_ils_low_10M --method grampa --step prep
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step prep
python scripts/check_pipeline_status.py conf_ils_low_10M --method mpsugar --step prep
python scripts/check_pipeline_status.py conf_ils_low_10M --method padre --step prep
```

**Why check before running?** Avoid wasting compute time if inputs are missing or incomplete.

### 3. After Methods Complete

Verify outputs were generated successfully:

```bash
# GRAMPA
python scripts/check_pipeline_status.py conf_ils_low_10M --method grampa --step run

# Polyphest - check all three percentile runs
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 50
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 70
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 90

# MPSUGAR
python scripts/check_pipeline_status.py conf_ils_low_10M --method mpsugar --step run

# PADRE
python scripts/check_pipeline_status.py conf_ils_low_10M --method padre --step run
```

---

## Expected Files

### SimPhy Simulation Output
- **Location:** `{network}/data/{CONFIG}/replicate_{1-5}/1/`
- **Files:** 1000 gene tree files named `g_*`

### Method Input Files

**GRAMPA** (`{network}/processed/{CONFIG}/grampa_input/replicate_{X}/`):
- ✅ Required: `grampa_trees.tre`, `clean_trees.tre`, `species.tre`
- ℹ️ Optional: `taxa_map.txt` (only if substring fixes were needed)

**Polyphest** (`{network}/processed/{CONFIG}/polyphest_input/replicate_{X}/`):
- ✅ Required: `polyphest_trees.tre`, `multi_set.txt`

**MPSUGAR** (`{network}/processed/{CONFIG}/mpsugar_input/replicate_{X}/`):
- ✅ Required: `mpsugar_trees.nex` (NEXUS format), `taxon_map.json` (JSON format)

**PADRE** (`{network}/processed/{CONFIG}/padre_input/replicate_{X}/`):
- ✅ Required: `padre_trees.tre`

### Method Output Files

**GRAMPA** (`{network}/results/{CONFIG}/grampa/replicate_{X}/`):
- Success indicator: `analysis/grampa-scores.txt` exists

**Polyphest** (`{network}/results/{CONFIG}/polyphest_p{PERCENTILE}/replicate_{X}/`):
- Success indicator: `polyphest_trees-polyphest.txt` exists

**MPSUGAR** (`{network}/results/{CONFIG}/mpsugar/replicate_{X}/`):
- Success indicator: `mpsugar_results.txt` exists

**PADRE** (`{network}/results/{CONFIG}/padre/replicate_{X}/`):
- Success indicator: `padre_tree-result.tre` exists

---

## Complete Validation Script

Use this bash script to check everything at once:

```bash
#!/bin/bash
# validate_pipeline.sh - Check all stages of the pipeline

CONFIG="conf_ils_low_10M"

echo "========================================"
echo "Pipeline Validation: $CONFIG"
echo "========================================"

echo ""
echo "=== 1. SimPhy Simulations ==="
python scripts/check_pipeline_status.py $CONFIG --step simphy

echo ""
echo "=== 2. GRAMPA Pipeline ==="
python scripts/check_pipeline_status.py $CONFIG --method grampa --step prep
python scripts/check_pipeline_status.py $CONFIG --method grampa --step run

echo ""
echo "=== 3. Polyphest Pipeline (3 percentiles) ==="
python scripts/check_pipeline_status.py $CONFIG --method polyphest --step prep
python scripts/check_pipeline_status.py $CONFIG --method polyphest --step run --percentile 50
python scripts/check_pipeline_status.py $CONFIG --method polyphest --step run --percentile 70
python scripts/check_pipeline_status.py $CONFIG --method polyphest --step run --percentile 90

echo ""
echo "=== 4. MPSUGAR Pipeline ==="
python scripts/check_pipeline_status.py $CONFIG --method mpsugar --step prep
python scripts/check_pipeline_status.py $CONFIG --method mpsugar --step run

echo ""
echo "=== 5. PADRE Pipeline ==="
python scripts/check_pipeline_status.py $CONFIG --method padre --step prep
python scripts/check_pipeline_status.py $CONFIG --method padre --step run

echo ""
echo "========================================"
echo "Validation Complete!"
echo "========================================"
```

---

## Advanced Usage

### Export Results for Multiple Configurations

```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --export low_results.csv
python scripts/check_pipeline_status.py conf_ils_medium_10M --export medium_results.csv
python scripts/check_pipeline_status.py conf_ils_high_10M --export high_results.csv
```

CSV includes: Network, Replicate, Status, Details, Missing Files, Errors

### Verbose Mode for Troubleshooting

```bash
# Show all networks including successful ones
python scripts/check_pipeline_status.py conf_ils_low_10M --method grampa --step run --verbose
```

### Check Specific Method with Custom Parameters

```bash
# Polyphest with specific percentile
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 98

# MPSUGAR with specific parameters (if you changed defaults)
python scripts/check_pipeline_status.py conf_ils_low_10M --method mpsugar --step run --iterations 1000 --chains 2
```

---

## Understanding Output

### Summary Output
```
================================================================================
GRAMPA Output Files - conf_ils_low_10M
================================================================================
Checking: 21 networks × 5 replicates = 105 combinations

Success:  105 / 105 (100.0%)
================================================================================
```

**Meaning:**
- Checking all 21 networks
- Each network has 5 replicates
- Total: 105 combinations to check
- All 105 passed successfully

### Status Codes

- ✅ **SUCCESS** - All expected files present and valid
- ~ **PARTIAL** - Some files present but incomplete (e.g., 500/1000 trees)
- ✗ **FAILED** - Files/directories exist but are empty or have errors
- ? **MISSING** - Expected files/directories do not exist

### Detailed Output (with --verbose)

```
Bendiksby_2011:
--------------------------------------------------------------------------------
  Replicate 1: [✓] SUCCESS  - All 3 required files present (+1 optional)
  Replicate 2: [✓] SUCCESS  - All 3 required files present
  Replicate 3: [✗] FAILED   - Missing 1/3 required files
    Missing: species.tre
  Replicate 4: [✓] SUCCESS  - All 3 required files present
  Replicate 5: [✓] SUCCESS  - All 3 required files present
```

---

## Common Issues & Solutions

### Issue: "Output directory does not exist" for Polyphest

**Cause:** You didn't specify the percentile, and the script is looking for `polyphest_p60/` but you have `polyphest_p50/`.

**Solution:** Always specify the percentile you used:
```bash
python scripts/check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 50
```

### Issue: taxa_map.txt marked as missing for GRAMPA

**Status:** Fixed! taxa_map.txt is now optional.

**Explanation:** This file is only created if substring fixes were needed. The script now shows:
- If present: "All 3 required files present (+1 optional)"
- If absent: "All 3 required files present"

### Issue: Unclear what "Total: 105" means

**Status:** Fixed! Output now shows:
```
Checking: 21 networks × 5 replicates = 105 combinations
```

---

## Integration with Workflow

### Typical Research Workflow

```bash
# 1. Run SimPhy simulations
sbatch simulations/jobs/submit_simphy.sh conf_ils_low_10M

# 2. Verify simulations completed
python scripts/check_pipeline_status.py conf_ils_low_10M --step simphy

# 3. Prepare inputs for methods
sbatch simulations/jobs/submit_grampa_pipeline.sh conf_ils_low_10M --prep-only
sbatch simulations/jobs/submit_polyphest_pipeline.sh conf_ils_low_10M --prep-only
# ... (other methods)

# 4. Verify inputs are ready
python scripts/check_pipeline_status.py conf_ils_low_10M --step prep

# 5. Run methods
sbatch simulations/jobs/submit_grampa_pipeline.sh conf_ils_low_10M --run-only
sbatch simulations/jobs/submit_polyphest_pipeline.sh conf_ils_low_10M --run-only
# ... (other methods)

# 6. Verify outputs were generated
python scripts/check_pipeline_status.py conf_ils_low_10M --step run --verbose

# 7. Export results for analysis
python scripts/check_pipeline_status.py conf_ils_low_10M --export results.csv
```

---

## Technical Details

### Directory Structure
```
simulations/simulations/
├── {network}/
│   ├── data/{CONFIG}/
│   │   └── replicate_{1-5}/
│   │       └── 1/
│   │           └── g_* (1000 files)
│   ├── processed/{CONFIG}/
│   │   ├── grampa_input/replicate_{1-5}/
│   │   ├── polyphest_input/replicate_{1-5}/
│   │   ├── mpsugar_input/replicate_{1-5}/
│   │   └── padre_input/replicate_{1-5}/
│   └── results/{CONFIG}/
│       ├── grampa/replicate_{1-5}/
│       ├── polyphest_p{PERCENTILE}/replicate_{1-5}/
│       ├── mpsugar/replicate_{1-5}/
│       └── padre/replicate_{1-5}/
```

### Method-Specific Notes

**Polyphest:**
- Output directory includes percentile: `polyphest_p{PERCENTILE}/`
- You MUST specify `--percentile XX` when checking
- Default is 60, but you use 50, 70, and 90

**MPSUGAR:**
- Output directory does NOT include parameters: just `mpsugar/`
- If you run with different iterations/chains, outputs overwrite
- Input files are NEXUS (.nex) and JSON (.json), not text/tre

**GRAMPA:**
- Success indicator is in subdirectory: `analysis/grampa-scores.txt`
- taxa_map.txt is optional (only for substring fixes)

**PADRE:**
- Uses standard .tre format for both input and output

### SimPhy Output Structures Supported

The script handles both:
1. **Single batch:** `replicate_X/1/g_*`
2. **Multiple batches:** `replicate_X/batch_Y/1/g_*`

Your setup uses single batch structure.

---

## Command-Line Options

```bash
usage: check_pipeline_status.py [-h] [--step {simphy,prep,run,all}]
                                [--method {grampa,polyphest,mpsugar,padre,all}]
                                [--verbose] [--export FILE]
                                [--percentile PERCENTILE]
                                [--iterations ITERATIONS] [--chains CHAINS]
                                config

Arguments:
  config                Configuration name (e.g., conf_ils_low_10M)

Options:
  --step {simphy,prep,run,all}
                        Which step to check (default: all)
  --method {grampa,polyphest,mpsugar,padre,all}
                        Which method to check (default: all)
  --verbose, -v         Show detailed results including successful runs
  --export FILE         Export results to CSV file
  --percentile PERCENTILE
                        Polyphest percentile parameter (default: 60)
  --iterations ITERATIONS
                        MPSUGAR iterations parameter (default: 500)
  --chains CHAINS       MPSUGAR chains parameter (default: 1)
```

---

## Tips & Best Practices

1. **Always check SimPhy first** - No point checking methods if simulations failed
2. **Check inputs before running** - Avoid wasting compute time on missing files
3. **Use --verbose for troubleshooting** - See exactly which networks/replicates failed
4. **Export to CSV for tracking** - Easier to track progress across configurations
5. **Specify percentile for Polyphest** - Don't rely on the default (60)
6. **Run checks after each stage** - Catch issues early in the pipeline

---

## Getting Help

For issues or questions about the script:
1. Check this guide first
2. Run with `--verbose` to see detailed output
3. Check SLURM logs in `simulations/logs/` for method-specific errors
4. Use `python scripts/check_pipeline_status.py --help` for command options
