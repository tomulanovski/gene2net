# Pipeline Monitoring Quick Reference

## Overview

Use `check_pipeline_status.py` to validate your simulation and inference pipeline at every stage.

## Basic Usage

```bash
# Check everything for a configuration
python check_pipeline_status.py conf_ils_low_10M
```

## Step-by-Step Validation

### 1. After SimPhy Simulations

Verify that all SimPhy simulations completed with 1000 trees per replicate:

```bash
python check_pipeline_status.py conf_ils_low_10M --step simphy
```

**What it checks:**
- 21 networks × 5 replicates = 105 total combinations
- Each replicate should have exactly 1000 gene trees
- Handles both single-batch and multi-batch SimPhy output structures

### 2. Before Running Methods (Input Validation)

Check that all input files are prepared:

```bash
# Check all methods
python check_pipeline_status.py conf_ils_low_10M --step prep

# Check specific method
python check_pipeline_status.py conf_ils_low_10M --method grampa --step prep
python check_pipeline_status.py conf_ils_low_10M --method polyphest --step prep
python check_pipeline_status.py conf_ils_low_10M --method mpsugar --step prep
python check_pipeline_status.py conf_ils_low_10M --method padre --step prep
```

**Expected input files per method:**

**GRAMPA:**
- `grampa_trees.tre` - Gene trees formatted for GRAMPA
- `clean_trees.tre` - Clean gene trees (for ASTRAL)
- `taxa_map.txt` - Taxa mapping file
- `species.tre` - Species tree from ASTRAL

**Polyphest:**
- `polyphest_trees.tre` - Gene trees formatted for Polyphest
- `multi_set.txt` - Consensus multiset file

**MPSUGAR:**
- `mpsugar_trees.nex` - Gene trees in NEXUS format
- `taxon_map.json` - Taxa mapping in JSON format

**PADRE:**
- `padre_trees.tre` - Gene trees formatted for PADRE

**Expected output files per method:**

**GRAMPA:**
- `analysis/grampa-scores.txt` - GRAMPA completion indicator

**Polyphest:**
- `polyphest_trees-polyphest.txt` - Polyphest network output

**MPSUGAR:**
- `mpsugar_results.txt` - MPSUGAR results file

**PADRE:**
- `padre_tree-result.tre` - PADRE species tree

### 3. After Running Methods (Output Validation)

Check that methods produced output files:

```bash
# Check all methods
python check_pipeline_status.py conf_ils_low_10M --step run

# Check specific method
python check_pipeline_status.py conf_ils_low_10M --method grampa --step run
python check_pipeline_status.py conf_ils_low_10M --method polyphest --step run
python check_pipeline_status.py conf_ils_low_10M --method mpsugar --step run
python check_pipeline_status.py conf_ils_low_10M --method padre --step run
```

**What it checks:**
- Output directories exist
- Result files (.nwk, .tre, network files) are present
- Reports which network-replicate combinations succeeded/failed

### 4. With Custom Parameters

For methods with variable parameters:

```bash
# Polyphest with different percentile
python check_pipeline_status.py conf_ils_low_10M --method polyphest --step run --percentile 98

# MPSUGAR with different iterations/chains
python check_pipeline_status.py conf_ils_low_10M --method mpsugar --step run --iterations 1000 --chains 2
```

## Output Modes

### Summary Mode (Default)

Shows overview statistics and only failures/missing data:

```bash
python check_pipeline_status.py conf_ils_low_10M
```

### Verbose Mode

Shows detailed information including successful runs:

```bash
python check_pipeline_status.py conf_ils_low_10M --verbose
```

### Export to CSV

Export results for analysis in Excel/R/Python:

```bash
python check_pipeline_status.py conf_ils_low_10M --export results.csv
```

CSV includes:
- Step (e.g., "SimPhy", "GRAMPA_prep", "GRAMPA_run")
- Network name
- Replicate number
- Status (SUCCESS/FAILED/MISSING/PARTIAL)
- Details
- Missing files
- Error messages

## Common Workflows

### After running SimPhy for all 3 configurations:

```bash
python check_pipeline_status.py conf_ils_low_10M --step simphy
python check_pipeline_status.py conf_ils_medium_10M --step simphy
python check_pipeline_status.py conf_ils_high_10M --step simphy
```

### Before submitting pipeline jobs:

```bash
# Make sure inputs are ready to avoid wasting compute time
python check_pipeline_status.py conf_ils_low_10M --step prep
```

### After pipeline completes:

```bash
# Check which jobs succeeded and which failed
python check_pipeline_status.py conf_ils_low_10M --step run --verbose

# Export for further analysis
python check_pipeline_status.py conf_ils_low_10M --export low_results.csv
python check_pipeline_status.py conf_ils_medium_10M --export medium_results.csv
python check_pipeline_status.py conf_ils_high_10M --export high_results.csv
```

### Troubleshooting specific methods:

```bash
# Check if GRAMPA inputs are ready
python check_pipeline_status.py conf_ils_low_10M --method grampa --step prep --verbose

# Check which GRAMPA jobs failed and why
python check_pipeline_status.py conf_ils_low_10M --method grampa --step run --verbose
```

## Understanding Status Codes

- **SUCCESS**: All expected files present and valid
- **PARTIAL**: Some files present but incomplete (e.g., 500/1000 trees for SimPhy)
- **FAILED**: Files or directories exist but are empty or corrupted
- **MISSING**: Expected files/directories do not exist

## Tips

1. **Run checks before submitting jobs** to avoid wasting compute time
2. **Use `--step prep`** after prep jobs to ensure inputs are complete
3. **Use `--step run`** after method runs to identify failures quickly
4. **Export to CSV** for tracking progress across multiple configurations
5. **Use `--verbose`** when troubleshooting specific network-replicate combinations
6. **Check SimPhy first** - no point in running methods if simulations failed

## Integration with Existing Scripts

This replaces `check_simulation_status.sh` with extended functionality:

**Old way (only SimPhy):**
```bash
export CHECK_CONFIG="conf_ils_low_10M"
./check_simulation_status.sh
```

**New way (SimPhy + all methods):**
```bash
python check_pipeline_status.py conf_ils_low_10M
```
