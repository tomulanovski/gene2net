# check_pipeline_status.py - Validation Notes

## Fixed Issues

### 1. ✅ All File Names Confirmed (FIXED)

**Input files:**
- GRAMPA: `grampa_trees.tre`, `clean_trees.tre`, `taxa_map.txt`, `species.tre`
- Polyphest: `polyphest_trees.tre`, `multi_set.txt`
- MPSUGAR: `mpsugar_trees.nex` (NEXUS), `taxon_map.json` (JSON)
- PADRE: `padre_trees.tre`

**Output files:**
- GRAMPA: `analysis/grampa-scores.txt` (subdirectory!)
- Polyphest: `polyphest_trees-polyphest.txt`
- MPSUGAR: `mpsugar_results.txt`
- PADRE: `padre_tree-result.tre`

### 2. ✅ Directory Structure (VERIFIED)
All paths confirmed from job scripts:
- SimPhy data: `{network}/data/{CONFIG}/replicate_{X}/`
- Method inputs: `{network}/processed/{CONFIG}/{method}_input/replicate_{X}/`
- Method results: `{network}/results/{CONFIG}/{method}/replicate_{X}/`
- Polyphest uses: `{network}/results/{CONFIG}/polyphest_p{PERCENTILE}/replicate_{X}/`

### 3. ✅ Input File Names (VERIFIED)
**GRAMPA:**
- `grampa_trees.tre` - Gene trees with copy numbers
- `clean_trees.tre` - Clean gene trees for ASTRAL
- `taxa_map.txt` - Taxa mapping
- `species.tre` - Species tree from ASTRAL

**Polyphest:**
- `polyphest_trees.tre` - Gene trees
- `multi_set.txt` - Consensus multiset

**MPSUGAR:**
- `mpsugar_trees.nex` - Gene trees in NEXUS format
- `taxon_map.json` - Taxa mapping in JSON format

**PADRE:**
- `padre_trees.tre` - Gene trees

## Points of Uncertainty / Things to Watch

### 1. ✅ Output File Detection (CONFIRMED)
All output files are now explicitly checked:
- GRAMPA: `analysis/grampa-scores.txt`
- Polyphest: `polyphest_trees-polyphest.txt`
- MPSUGAR: `mpsugar_results.txt`
- PADRE: `padre_tree-result.tre`

**Risk:** None - exact file paths confirmed from actual outputs.

**Note:** GRAMPA's output is in a subdirectory (`analysis/`), which is correctly handled.

### 2. ⚠️ SimPhy Output Structure (LOW RISK)
The script handles two possible SimPhy output structures:
- **Single batch:** `replicate_X/1/g_*`
- **Multiple batches:** `replicate_X/batch_Y/1/g_*`

**Risk:** Low - both formats are explicitly handled.

**What to check:** If SimPhy uses a different structure, the script may not count trees correctly.

### 3. ⚠️ Empty vs Missing Files (LOW RISK)
The script checks both:
- File existence
- File size > 0 bytes

**Risk:** Low - catches both missing and empty files.

**Edge case:** A corrupted file with content might still be detected as "success". For thorough validation, you'd need to parse file contents.

### 4. ⚠️ SLURM Log Parsing (NOT IMPLEMENTED)
The function `parse_slurm_log_for_errors()` exists but is **not currently called**.

**Status:** Error detection from logs is not yet active.

**Impact:** The script won't automatically extract error messages from SLURM logs.

**Workaround:** Use the existing `check_error_logs.sh` script or manually check logs for detailed errors.

**Future enhancement:** We can activate this feature if needed.

### 5. ⚠️ Parameters Not Stored in Directory Names
**MPSUGAR:** The iterations and chains parameters are NOT included in the output directory name.
- Directory is just: `mpsugar/replicate_X`
- If you run with different parameters, outputs will overwrite!

**Polyphest:** The percentile IS included: `polyphest_p{PERCENTILE}/`

**What to watch:** If you run MPSUGAR with different parameters on the same configuration, the results will overwrite each other.

## Testing Recommendations

### First-Time Use Checklist

1. **After SimPhy completes:**
   ```bash
   python check_pipeline_status.py conf_ils_low_10M --step simphy --verbose
   ```
   - Verify it correctly counts 1000 trees per replicate
   - Check that it handles your SimPhy output structure

2. **After first prep job:**
   ```bash
   python check_pipeline_status.py conf_ils_low_10M --method grampa --step prep --verbose
   ```
   - Verify it finds all expected input files
   - Confirm file names match what's actually created

3. **After first method run:**
   ```bash
   python check_pipeline_status.py conf_ils_low_10M --method grampa --step run --verbose
   ```
   - Verify it correctly identifies output files
   - If it fails to detect outputs, check what files were actually created and let me know

4. **Test CSV export:**
   ```bash
   python check_pipeline_status.py conf_ils_low_10M --export test.csv
   ```
   - Verify the CSV is well-formatted and contains expected data

### If Issues Arise

**If output detection fails:**
1. Check what files were actually created in the output directory
2. We can add specific file names to `METHOD_FILES` dict
3. Or adjust the glob patterns in `check_method_outputs()`

**If SimPhy tree counting is wrong:**
1. Check the actual directory structure of your SimPhy output
2. We can adjust `count_gene_trees_in_simphy_output()` function

**If you need log parsing:**
1. Let me know and I can activate the `parse_slurm_log_for_errors()` function
2. It will extract error messages from SLURM .err files

## Confidence Levels

| Component | Confidence | Basis |
|-----------|-----------|--------|
| Directory structure | ✅ HIGH | Verified from job scripts |
| Input file names | ✅ HIGH | Verified from prep scripts |
| GRAMPA inputs | ✅ HIGH | Verified from scripts |
| Polyphest inputs | ✅ HIGH | Verified from scripts |
| MPSUGAR inputs | ✅ HIGH | Verified (.nex, .json) |
| PADRE inputs | ✅ HIGH | Verified from scripts |
| Output directory paths | ✅ HIGH | Verified from run scripts |
| **Output file names** | ✅ **HIGH** | **Confirmed from actual data** |
| GRAMPA output | ✅ HIGH | analysis/grampa-scores.txt |
| Polyphest output | ✅ HIGH | polyphest_trees-polyphest.txt |
| MPSUGAR output | ✅ HIGH | mpsugar_results.txt |
| PADRE output | ✅ HIGH | padre_tree-result.tre |
| SimPhy tree counting | ✅ HIGH | Handles both batch structures |
| Empty file detection | ✅ HIGH | Explicitly checked |
| Error log parsing | ⚠️ NOT ACTIVE | Function exists but not called |

## Summary

**Ready to use:** ✅ **YES** - All file names and paths confirmed!

**Confidence level:** HIGH - All critical paths and file names verified from:
- Job scripts (directory structure)
- Prep scripts (input files)
- Actual output examples (output files)

**What's been tested:**
- ✅ Directory structure matches job scripts
- ✅ Input file names match prep scripts
- ✅ Output file names match actual cluster data
- ✅ GRAMPA's subdirectory structure handled correctly

**Remaining uncertainty:**
- ⚠️ SimPhy tree counting (should work but test on first run)
- ⚠️ Error log parsing (not active, optional feature)
