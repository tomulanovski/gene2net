# Post-Processing Guide

This guide explains how to post-process raw program outputs into clean MUL-tree files required by the summary pipeline.

## Purpose

After running phylogenetic network inference methods (GRAMPA, Polyphest, MPSUGAR, PADRE), their raw output files need to be post-processed to extract clean MUL-tree strings. The `postprocess_results.py` script automates this extraction.

---

## Why Post-Processing?

Each program outputs results in different formats:

| Program | Raw Output File | Location | Format | What We Extract |
|---------|----------------|----------|--------|-----------------|
| **Polyphest** | `polyphest_trees-polyphest.txt` | `results/{method}/replicate_N/` | Text with "multree:" and "network:" lines | MUL-tree from "multree:" line |
| **GRAMPA** | `grampa-scores.txt` | `results/grampa/replicate_N/` | TSV with multiple ranked trees | Best tree (first row, last column) |
| **MPSUGAR** | `mpsugar_results.txt` | `results/mpsugar/replicate_N/` | Detailed report with "Newick:" line | Newick tree from "Newick:" line |
| **PADRE** | `padre_trees-result.tre` | `processed/padre_input/replicate_N/` ⚠️ | Already clean MUL-tree | Copy to results directory |
| **AlloppNET** | `sampledmultrees.txt` | `results/alloppnet/replicate_N/` | BEAST MCMC output (NEXUS trees) | TreeAnnotator consensus + copy number removal |

All programs output standardized `{program}_result.tre` files after post-processing.

---

## Expected Input Files

### Directory Structure

```
/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/
└── {NetworkName}/
    ├── processed/                                      # Input files
    │   └── {CONFIG}/
    │       └── padre_input/
    │           ├── replicate_1/
    │           │   ├── padre_trees.tre                 # Input
    │           │   └── padre_trees-result.tre          # PADRE writes here!
    │           └── ... (replicates 2-5)
    │
    └── results/                                        # Output files
        └── {CONFIG}/
            ├── grampa/
            │   ├── replicate_1/
            │   │   └── grampa-scores.txt              # Raw output (TSV)
            │   ├── replicate_2/
            │   └── ... (replicates 3-5)
            │
            ├── polyphest_p50/
            │   ├── replicate_1/
            │   │   └── polyphest_trees-polyphest.txt  # Raw output (text)
            │   ├── replicate_2/
            │   └── ... (replicates 3-5)
            │
            ├── mpsugar/
            │   ├── replicate_1/
            │   │   └── mpsugar_results.txt            # Raw output (report)
            │   ├── replicate_2/
            │   └── ... (replicates 3-5)
            │
            ├── padre/
            │   ├── replicate_1/
            │   │   └── padre_result.tre               # Created by post-processing
            │   ├── replicate_2/
            │   └── ... (replicates 3-5)
            │
            └── alloppnet/
                ├── replicate_1/
                │   ├── sampledmultrees.txt            # Raw BEAST output
                │   └── alloppnet_result.tre           # Created by post-processing
                ├── replicate_2/
                └── ... (replicates 3-5, only 8 compatible networks)
```

**Note**: PADRE uniquely writes output to its input directory (`processed/padre_input/`), not the results directory. Post-processing copies it to `results/padre/` where the summary pipeline expects it.

### Example Raw Outputs

#### Polyphest Output (`polyphest_trees-polyphest.txt`)

```
multree: (((A,B),C),D);
network: (((A,B)#H1,C),#H1);
Time taken: 9.92 seconds
```

**Extracted**: `(((A,B),C),D);`

#### GRAMPA Output (`grampa-scores.txt`)

```
mul.tree    h1.node    h2.node    score    labeled.tree
5137        <29>       <39>       3309     (((A,B),C),D);
5116        <28>       <39>       3326     (((A,C),B),D);
1506        ...        ...        3375     ...
```

**Extracted**: `(((A,B),C),D);` (first data row, last column)

#### MPSUGAR Output (`mpsugar_results.txt`)

```
MP-SUGAR Results
================================================================================
Input file: /path/to/input.nex
Date: 2025-12-08 20:23:05

Network 1:
Newick: (((A,B)UID_1,C)UID_2,D)UID_3;
Score: -80228
```

**Extracted**: `(((A,B)UID_1,C)UID_2,D)UID_3;`

#### PADRE Output (`padre_trees-result.tre`)

**Location**: `processed/conf_ils_low_10M/padre_input/replicate_1/padre_trees-result.tre`

**Note**: PADRE writes output to input directory, not results directory!

```
(((A,B),C),D);
```

**Extracted**: Same (already clean), but **copied to results directory**

#### AlloppNET Output (`sampledmultrees.txt`)

**Location**: `results/conf_ils_low_10M/alloppnet/replicate_1/sampledmultrees.txt`

**Note**: AlloppNET post-processing runs TreeAnnotator and removes copy number suffixes.

```
#NEXUS

BEGIN TREES;
tree STATE_0 = [initial tree];
tree STATE_1000 = [sampled tree];
tree STATE_2000 = [sampled tree];
...
tree STATE_100000000 = [final sampled tree];
END;
```

**Extraction Process**:
1. Count trees in `sampledmultrees.txt`
2. Calculate 10% burnin (exclude STATE_0)
3. Run TreeAnnotator to create consensus tree
4. Remove `_0` and `_1` suffixes from tip labels
5. Write to `alloppnet_result.tre`

**Extracted**: Clean MUL-tree without copy number suffixes

---

## Usage

### Basic Usage

```bash
# Navigate to gene2net root
cd /groups/itay_mayrose/tomulanovski/gene2net

# Post-process all methods for one configuration
python simulations/scripts/postprocess_results.py conf_ils_low_10M
```

**Expected Output**:
- Processes files for all methods (varies by method - AlloppNET only runs on 8 networks)
- Creates `{program}_result.tre` in each replicate directory
- Shows success/failure statistics

### Process Specific Methods

```bash
# Only GRAMPA
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods grampa

# GRAMPA and Polyphest p50
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods grampa polyphest_p50

# AlloppNET (only runs on 8 compatible networks)
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods alloppnet

# All Polyphest variants
python simulations/scripts/postprocess_results.py conf_ils_low_10M \
    --methods polyphest_p50 polyphest_p70 polyphest_p90
```

### Dry Run

```bash
# Preview what would be done without writing files
python simulations/scripts/postprocess_results.py conf_ils_low_10M --dry-run
```

### Process Multiple Configurations

```bash
# Process all configurations
for config in conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M; do
    echo "Post-processing $config..."
    python simulations/scripts/postprocess_results.py $config
done
```

---

## Output Files

After post-processing, each replicate directory will contain:

```
results/{CONFIG}/{method}/replicate_N/
├── {original_output_file}      # Raw program output (unchanged)
└── {program}_result.tre         # Extracted clean MUL-tree (NEW)
```

### Examples

```
# GRAMPA
results/conf_ils_low_10M/grampa/replicate_1/
├── grampa-scores.txt           # Original (1000+ lines)
└── grampa_result.tre           # Extracted best tree (1 line)

# Polyphest
results/conf_ils_low_10M/polyphest_p50/replicate_1/
├── polyphest_trees-polyphest.txt    # Original (3 lines)
└── polyphest_result.tre             # Extracted MUL-tree (1 line)

# MPSUGAR
results/conf_ils_low_10M/mpsugar/replicate_1/
├── mpsugar_results.txt         # Original (50+ lines)
└── mpsugar_result.tre          # Extracted Newick (1 line)

# PADRE
results/conf_ils_low_10M/padre/replicate_1/
├── padre_tree-result.tre       # Original (1 line)
└── padre_result.tre            # Copy (1 line)

# AlloppNET
results/conf_ils_low_10M/alloppnet/replicate_1/
├── sampledmultrees.txt         # Original (BEAST MCMC output, 100K+ trees)
└── alloppnet_result.tre        # Extracted consensus tree (1 line)
```

---

## Complete Workflow

### Typical Workflow: From Raw Results to Summary

```bash
# 1. Navigate to project
cd /groups/itay_mayrose/tomulanovski/gene2net

# 2. Check which methods completed
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run

# 3. Post-process raw outputs to extract MUL-trees
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# Expected output:
#   Post-processing Statistics
#   ==================================================
#   Total files:        630
#   Successfully processed: 542 (86.0%)
#   Skipped (missing or already done): 76 (12.1%)
#   Failed:             12 (1.9%)

# 4. Verify post-processed files exist
python simulations/scripts/check_pipeline_status.py conf_ils_low_10M --step run

# 5. Run summary pipeline
python simulations/scripts/run_full_summary.py conf_ils_low_10M
```

---

## Statistics and Reporting

### Output Statistics

The script reports:
- **Total files**: Expected files based on config (21 networks × N methods × 5 replicates)
- **Successfully processed**: Files successfully extracted and written
- **Skipped**: Files missing or already processed
- **Failed**: Files that exist but extraction failed

### Example Output

```
================================================================================
Post-processing Results: conf_ils_low_10M
================================================================================
Methods: grampa, polyphest_p50, polyphest_p70, polyphest_p90, mpsugar, padre
Networks: 21
Replicates: 5
Dry run: False
================================================================================

Processing grampa...
  Processed 50 files (Success: 45, Skipped: 3, Failed: 2)
  Processed 100 files (Success: 92, Skipped: 5, Failed: 3)
  Processed 105 files (Success: 98, Skipped: 5, Failed: 2)

Processing polyphest_p50...
  ...

================================================================================
Post-processing Statistics
================================================================================
Total files:        630
Successfully processed: 542 (86.0%)
Skipped (missing or already done): 76 (12.1%)
Failed:             12 (1.9%)

Errors (12):
  1. Bendiksby_2011 / grampa / rep_3
     Failed to extract tree from grampa-scores.txt
  2. Wu_2015 / mpsugar / rep_1
     Failed to extract tree from mpsugar_results.txt
  ...
================================================================================
```

---

## Extraction Details

### Polyphest Extraction

**Input Format**:
```
multree: (tree_string_here);
network: (network_string_here);
Time taken: X.XXX seconds
```

**Extraction Logic**:
1. Read file line by line
2. Find line starting with `multree:`
3. Extract everything after `multree: ` (including semicolon)
4. Write to `polyphest_result.tre`

### GRAMPA Extraction

**Input Format** (TSV):
```
mul.tree    h1.node    h2.node    score    labeled.tree
5137        <29>       <39>       3309     (best_tree);
5116        <28>       <39>       3326     (second_best);
...
```

**Extraction Logic**:
1. Skip header line (line 1)
2. Read first data line (line 2) - this is the best-scoring tree
3. Split by tab (`\t`)
4. Extract last column (`labeled.tree`)
5. Ensure ends with semicolon
6. Write to `grampa_result.tre`

### MPSUGAR Extraction

**Input Format**:
```
MP-SUGAR Results
...
Network 1:
Newick: (tree_string_here);
Score: -12345
...
```

**Extraction Logic**:
1. Read file line by line
2. Find line starting with `Newick:`
3. Extract everything after `Newick: ` (including semicolon)
4. Write to `mpsugar_result.tre`

### PADRE Extraction

**Input Location**: `processed/{config}/padre_input/replicate_{N}/padre_trees-result.tre`

**Output Location**: `results/{config}/padre/replicate_{N}/padre_result.tre`

**Input Format**: Already clean MUL-tree
```
(tree_string_here);
```

**Extraction Logic**:
1. Read from `processed/padre_input/` directory (PADRE's output location)
2. Strip whitespace
3. Ensure ends with semicolon
4. Create `results/padre/replicate_N/` directory if needed
5. Write to `padre_result.tre` in results directory

**Why copy?** PADRE writes output to input directory, but summary pipeline expects results in standard location.

### AlloppNET Extraction

**Input Location**: `results/{config}/alloppnet/replicate_{N}/sampledmultrees.txt`

**Output Location**: `results/{config}/alloppnet/replicate_{N}/alloppnet_result.tre`

**Input Format**: BEAST MCMC output (NEXUS format with multiple trees)
```
#NEXUS
BEGIN TREES;
tree STATE_0 = [tree];
tree STATE_1000 = [tree];
...
END;
```

**Extraction Logic**:
1. Count trees in `sampledmultrees.txt` (count "tree STATE_" occurrences)
2. Calculate 10% burnin: `(total_trees - 1) // 10` (exclude STATE_0)
3. Run TreeAnnotator with burnin and mean heights:
   - Input: `sampledmultrees.txt`
   - Output: `alloppnet_consensus.tre` (temporary)
4. Remove copy number suffixes (`_0`, `_1`) from tip labels using `remove_copy_numbers.py`
5. Write cleaned tree to `alloppnet_result.tre`

**Requirements**:
- Conda environment `alloppnet` (for TreeAnnotator)
- Conda environment `gene2net` (for Python/Biopython)
- `remove_copy_numbers.py` script in `scripts/alloppnet/`

**Note**: Works even if BEAST job timed out (uses whatever trees were sampled)

---

## Troubleshooting

### "Failed to extract tree"

**Possible Causes**:
1. Output file is empty or corrupted
2. Unexpected format (program version changed?)
3. Incomplete program run (crashed mid-execution)

**Solutions**:
```bash
# Check raw output file manually
cat results/conf_ils_low_10M/grampa/replicate_1/grampa-scores.txt

# Re-run the program if output is incomplete
# Then re-run post-processing
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods grampa
```

### "Skipped (missing or already done)"

**Normal Behavior**: Script skips:
1. Files that don't exist (method didn't run)
2. Files already processed (non-empty `*_result.tre` exists)

**To Force Reprocess**:
```bash
# Delete existing result files
find results/conf_ils_low_10M/ -name "*_result.tre" -delete

# Re-run post-processing
python simulations/scripts/postprocess_results.py conf_ils_low_10M
```

### High Failure Rate

**Diagnosis**:
- Check error messages in output
- Verify raw output files are complete
- Check for program-specific issues

**Example Debug**:
```bash
# Post-process with dry run to see what would happen
python simulations/scripts/postprocess_results.py conf_ils_low_10M --dry-run

# Process just one method to debug
python simulations/scripts/postprocess_results.py conf_ils_low_10M --methods grampa

# Manually check a failed file
head -20 results/conf_ils_low_10M/grampa/replicate_1/grampa-scores.txt
```

---

## Integration with Summary Pipeline

The post-processing script **must** be run before the summary pipeline:

```bash
# STEP 1: Run phylogenetic methods (SLURM jobs)
sbatch simulations/jobs/submit_grampa_pipeline.sh conf_ils_low_10M
sbatch simulations/jobs/submit_polyphest_pipeline.sh conf_ils_low_10M
# ... etc

# STEP 2: Post-process raw outputs
python simulations/scripts/postprocess_results.py conf_ils_low_10M

# STEP 3: Run summary pipeline
python simulations/scripts/run_full_summary.py conf_ils_low_10M
```

**Why?** The summary pipeline expects `{program}_result.tre` files, not raw program outputs.

---

## Advanced Usage

### Custom Configuration

```bash
# Use custom config file
python simulations/scripts/postprocess_results.py conf_ils_low_10M \
    --config /path/to/custom_config.yaml
```

### Re-process Specific Networks

For fine-grained control, edit the script or use individual extractions:

```bash
# Example: Manually extract from one file
python -c "
with open('results/conf_ils_low_10M/grampa/replicate_1/grampa-scores.txt') as f:
    lines = f.readlines()
    best_tree = lines[1].strip().split('\t')[-1]
    with open('results/conf_ils_low_10M/grampa/replicate_1/grampa_result.tre', 'w') as out:
        out.write(best_tree + '\n')
"
```

---

## Getting Help

```bash
# Built-in help
python simulations/scripts/postprocess_results.py --help
```

---

## Summary

**What it does**: Extracts clean MUL-trees from program-specific output formats

**When to use**: After methods complete, before running summary pipeline

**Input**: Raw output files (`grampa-scores.txt`, `polyphest_trees-polyphest.txt`, `sampledmultrees.txt`, etc.)

**Output**: Standardized `{program}_result.tre` files

**Typical command**:
```bash
python simulations/scripts/postprocess_results.py conf_ils_low_10M
```
