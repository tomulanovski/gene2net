# MUL-tree Metrics Update - Complete Guide

## What Changed?

The analysis pipeline now uses **MUL-tree based metrics** as the PRIMARY accuracy measures, replacing network-based edit distance.

---

## 🎯 Primary Metrics (MUL-tree Based)

### 1. **Edit Distance on MUL-trees** (`edit_distance_multree`)
**What it does:** Compares tree structures directly before folding to networks

**Algorithm:**
- Converts both MUL-trees to directed graphs (preserving tree structure)
- Computes graph edit distance (node/edge insertions/deletions)
- Normalizes by `max(nodes₁ + edges₁, nodes₂ + edges₂)`

**Normalization:** YES - `distance / max(tree_size₁, tree_size₂)`

**Range:** 0.0 (identical) to ~1.0 (very different)

**Advantages:**
- ✅ More stable than network edit distance
- ✅ Directly comparable structures
- ✅ Handles polyploids naturally
- ✅ Rarely exceeds 1.0

---

### 2. **Robinson-Foulds Distance on MUL-trees** (`rf_distance`)
**What it does:** Counts bipartitions (splits) that differ between MUL-trees

**Algorithm:**
```python
# For each internal node, create a bipartition
# Example: Node splits tree into {A, A, B} | {C, D}
bipartitions = set of all splits in tree

# Compare trees
unique_to_tree1 = bipartitions1 - bipartitions2
unique_to_tree2 = bipartitions2 - bipartitions1
rf_distance = len(unique_to_tree1) + len(unique_to_tree2)

# Normalize
normalized_rf = rf_distance / (len(bipartitions1) + len(bipartitions2))
```

**Key Innovation:** Uses `frozenset` to handle **duplicate leaf labels** (polyploids)

**Normalization:** YES - `RF / total_bipartitions_in_both_trees`

**Range:** 0.0 (identical topology) to 1.0 (completely different)

**Advantages:**
- ✅ Classic phylogenetic metric adapted for polyploids
- ✅ Focuses on topology (splits/clades)
- ✅ Computationally efficient
- ✅ Always normalized to [0, 1]

**Difference from Classic RF:**
- Classic RF: Requires unique leaf labels, fails on polyploids
- Our RF: Uses frozensets, works with duplicated labels

---

## 📊 Secondary/Legacy Metric

### 3. **Edit Distance on Networks** (`edit_distance`)
**What it does:** Compares folded networks (DAGs with reticulation nodes)

**Issues:**
- ⚠️ Can exceed 1.0 in edge cases (e.g., MPSUGAR often shows > 1.0)
- ⚠️ Depends on folding algorithm
- ⚠️ Less stable normalization

**When it exceeds 1.0:**
```
Example:
- Ground truth: 10 nodes + 9 edges = 19 total
- Inferred: 8 nodes + 7 edges = 15 total
- Normalization = max(19, 15) = 19
- Edit operations needed: 25
- Result: 25/19 = 1.32 (exceeds 1.0!)
```

**Status:** Kept for backward compatibility and comparison

---

## 📈 What You Get Now

### New Plots

1. **`08a_distance_metrics_comparison.pdf`** ⭐ NEW!
   - Side-by-side comparison of all three distance metrics
   - Green borders highlight primary metrics (MUL-tree based)
   - Shows mean values for each method

2. **Updated Accuracy Plots** (Category 2):
   - All now use `edit_distance_multree` instead of `edit_distance`
   - Examples:
     - `11_combined_editdist_multree_vs_num_species.pdf`
     - `12_combined_editdist_multree_vs_h_strict.pdf`
     - etc.

3. **New RF Distance Plots**:
   - `15_combined_rf_vs_num_species.pdf`
   - `16_combined_rf_vs_h_strict.pdf`
   - Plus faceted versions

4. **Updated Distribution Plots**:
   - `08_edit_distance_multree_boxplot.pdf` - Now shows MUL-tree metric
   - Method summary now includes MUL-tree edit distance

### Updated Summary Tables

**`01_method_performance_summary.csv`** now includes:

| Column | Description |
|--------|-------------|
| `Mean_Edit_Distance_MULtree` | Primary accuracy metric |
| `Median_Edit_Distance_MULtree` | Median MUL-tree edit distance |
| `Std_Edit_Distance_MULtree` | Variability |
| `Mean_RF_Distance` | Primary topology metric |
| `Median_RF_Distance` | Median RF distance |
| `Mean_Reticulation_MAE` | Absolute reticulation error |
| `Mean_Reticulation_Bias` | Signed reticulation error |

**`02_per_network_performance.csv`** now includes:
- `{method}_EditDist_MULtree` - MUL-tree edit distance per method
- `{method}_RF` - RF distance per method

### New Pivot Tables (Level 1)

In `level1_detailed_per_network/`:
- `edit_distance_multree.csv` - PRIMARY accuracy metric
- `rf_distance.csv` - PRIMARY topology metric
- `edit_distance_network.csv` - Legacy network metric (for comparison)

---

## 🔄 Re-running the Analysis

You've already re-run with `--force-recompute`, so all metrics are computed! 

The new outputs are automatically generated when you run:

```bash
python simulations/scripts/create_analysis_figures.py --config conf_ils_low_10M
```

Or re-run the full pipeline:

```bash
python simulations/scripts/run_full_summary.py conf_ils_low_10M
```

---

## 📖 Interpreting the Metrics

### When to Use Which Metric?

| Question | Use This Metric |
|----------|----------------|
| Overall accuracy of reconstruction? | `edit_distance_multree` |
| Topology/clade accuracy? | `rf_distance` |
| How many reticulations off? | `num_rets_diff` (MAE) |
| Over- or under-estimating reticulations? | `num_rets_bias` |
| Why is network edit distance > 1? | See explanation above; use MUL-tree metrics instead |

### Example Interpretation

If a method has:
- `edit_distance_multree = 0.15` → Good reconstruction (15% different)
- `rf_distance = 0.25` → Topology 75% correct
- `num_rets_bias = +2.3` → Over-estimates reticulations by ~2.3

---

## 🔬 Technical Details

### Why MUL-tree Metrics Are Better

1. **Stability**: MUL-trees have fixed structure (no folding ambiguity)
2. **Normalization**: More reliable bounds (less likely to exceed 1.0)
3. **Comparability**: Direct tree-to-tree comparison
4. **Polyploid-friendly**: Naturally handles duplicated labels
5. **Standard metrics**: RF is a classic phylogenetic distance

### Code Locations

- **Metric computation**: `simulations/scripts/reticulate_tree.py`
  - `get_edit_distance_multree()` (line 815-868)
  - `get_rf_distance()` (line 870-916)

- **Comparison**: `simulations/scripts/compare_reticulations.py`
  - `pairwise_compare()` includes all three metrics

- **Plotting**: `simulations/scripts/create_analysis_figures.py`
  - `plot_distance_metrics_comparison()` - 3-way comparison
  - All accuracy plots updated to use MUL-tree metrics

- **Aggregation**: `simulations/scripts/aggregate_and_summarize.py`
  - MUL-tree metrics prioritized in pivot tables

---

## ✅ Summary

### Before:
- Primary metric: Network edit distance (unstable, can exceed 1.0)
- No RF distance
- Limited to network-based comparisons

### After:
- Primary metrics: MUL-tree edit distance + RF distance (stable, normalized)
- 3-way comparison plot showing all metrics
- All plots updated to use MUL-tree metrics
- Summary tables include both primary metrics
- Network edit distance kept for comparison

### Impact on Your Thesis:
You can now confidently report:
1. **Robust accuracy**: MUL-tree edit distance (stable, interpretable)
2. **Topology correctness**: RF distance (classic phylogenetic metric)
3. **Method bias**: Signed reticulation error (over/under-estimation)
4. **Comparison**: Show that network edit distance is unreliable (>1.0 for some methods)

This gives you **stronger, more defensible results** for your thesis! 🎓

