# Bias Metrics Guide

## Overview

The analysis pipeline now includes **signed error (bias)** metrics in addition to absolute error metrics. This allows you to determine whether methods are **over-estimating** or **under-estimating** network properties.

## What Changed?

### 1. New Metric: `num_rets_bias`

In addition to `num_rets_diff` (absolute error), there is now **`num_rets_bias`** (signed error):

- **`num_rets_diff`**: Absolute difference `|predicted - truth|` (always positive)
- **`num_rets_bias`**: Signed difference `predicted - truth` (can be positive or negative)
  - **Positive bias** (+5): Method over-estimated by 5 reticulations
  - **Negative bias** (-5): Method under-estimated by 5 reticulations
  - **Zero bias** (0): Perfect estimate on average

### 2. Updated Files

#### `compare_reticulations.py`
- Modified `compare_num_rets()` to return both signed and absolute differences
- Added `num_rets_bias` to the metrics dictionary in `pairwise_compare()`

#### `aggregate_and_summarize.py`
- Added `num_rets_bias` to the list of metrics to aggregate
- Generates pivot tables for bias alongside absolute error
- Includes bias in correlation analyses

#### `create_analysis_figures.py`
- **Plot 06**: `reticulation_bias_histogram.pdf` - Shows bias distribution with mean bias
- **Plot 07**: `reticulation_bias_boxplot.pdf` - Boxplot with mean bias annotations
- **Plot 10**: `method_summary.pdf` - Now includes a 4th panel showing bias (color-coded: red = over-estimation, blue = under-estimation)
- **Tables**: Summary tables now include both MAE and bias columns

## How to Interpret Bias

### In Plots

1. **Boxplots** (Plot 07):
   - Y-axis crosses zero at perfect accuracy
   - Boxes above zero = over-estimation tendency
   - Boxes below zero = under-estimation tendency
   - Mean bias shown as text annotation above each box

2. **Histograms** (Plot 06):
   - Distribution centered at zero = unbiased method
   - Distribution shifted right = over-estimation bias
   - Distribution shifted left = under-estimation bias
   - Vertical red line shows mean bias

3. **Summary Bar Chart** (Plot 10):
   - 4th panel shows mean bias
   - Red bars = over-estimation
   - Blue bars = under-estimation
   - Height shows magnitude of bias

### In Summary Tables

The `01_method_performance_summary.csv` table now includes:

| Method | ... | Mean_Reticulation_MAE | Median_Reticulation_MAE | Mean_Reticulation_Bias | Median_Reticulation_Bias |
|--------|-----|-----------------------|-------------------------|------------------------|--------------------------|
| grampa | ... | 2.45                 | 2.00                   | +1.23                  | +1.00                   |
| ...    | ... | ...                  | ...                    | ...                    | ...                     |

- **MAE (Mean Absolute Error)**: Average magnitude of error (always positive)
- **Bias**: Average signed error
  - Positive = over-estimation
  - Negative = under-estimation

### In Level 1 Pivot Tables

New file: `summary/{config}/level1_detailed_per_network/num_rets_bias.csv`

Shows mean ± std bias for each (network, method) combination. Example:
```
network,grampa,polyphest,mpsugar
net1,+2.34 ± 0.45 (5),+1.12 ± 0.23 (5),-0.56 ± 0.12 (5)
net2,+3.45 ± 0.67 (5),...
```

Positive values = over-estimation, negative = under-estimation.

## Example Interpretation

If `grampa` has:
- `Mean_Reticulation_MAE = 2.5`
- `Mean_Reticulation_Bias = +2.3`

This tells you:
- On average, GRAMPA's predictions are **2.5 reticulations off** from the truth (magnitude)
- GRAMPA tends to **over-estimate** by **2.3 reticulations** (direction)
- The bias is close to the MAE, suggesting consistent over-estimation (not random errors)

If another method has:
- `Mean_Reticulation_MAE = 2.5`
- `Mean_Reticulation_Bias = +0.1`

This tells you:
- Similar magnitude of error (2.5)
- But nearly **unbiased** (+0.1 ≈ 0), meaning errors are roughly balanced between over- and under-estimation
- Errors are more random, not systematic

## Re-running the Analysis

The new bias metrics are automatically included when you re-run the pipeline:

```bash
# Step 1: Compute comparisons (this will generate num_rets_bias)
python simulations/scripts/compute_comparisons.py \
    simulations/analysis/summary/conf_ils_low_10M/inventory.csv \
    simulations/analysis/cache/

# Step 2: Aggregate and summarize (includes bias in pivot tables)
python simulations/scripts/aggregate_and_summarize.py \
    simulations/analysis/summary/conf_ils_low_10M/comparisons_raw.csv \
    simulations/analysis/summary/conf_ils_low_10M/

# Step 3: Create figures (generates bias plots)
python simulations/scripts/create_analysis_figures.py \
    --config conf_ils_low_10M
```

Or use the full pipeline script:
```bash
python simulations/scripts/run_full_summary.py --config conf_ils_low_10M
```

## Backward Compatibility

- Existing `num_rets_diff` metric is **unchanged** (still absolute error)
- Old scripts will continue to work
- New bias metrics are **additions**, not replacements
- If `num_rets_bias` is not available, plots fall back to `num_rets_diff`

## Questions?

The bias metrics help answer:
1. **Are methods systematically over-estimating or under-estimating?**
2. **Which methods have the most balanced errors?**
3. **Do errors vary by network complexity?**

Compare MAE (magnitude) with bias (direction) to understand error patterns!

