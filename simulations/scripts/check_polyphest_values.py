#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

POLYPHEST_THRESHOLDS = ['polyphest_p50', 'polyphest_p70', 'polyphest_p90']
BASE = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/summary')

configs = [
    'conf_ils_high_10M',
    'conf_dup_loss_low_10M_ne2M',
    'conf_dup_loss_medium_10M_ne2M',
    'conf_dup_loss_high_10M_ne2M',
]

# Diagnostic: show available metrics and methods for first config
diag_csv = BASE / configs[0] / 'comparisons_raw.csv'
if diag_csv.exists():
    diag = pd.read_csv(diag_csv)
    print("Columns:", list(diag.columns))
    print("Methods:", sorted(diag['method'].unique()))
    print("Metrics:", sorted(diag['metric'].unique()))
    print()

for metric in ['edit_distance_multree', 'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist']:
    print(f"\n--- {metric} ---")
    for config in configs:
        csv = BASE / config / 'comparisons_raw.csv'
        if not csv.exists():
            print(f"  {config}: file not found")
            continue
        df = pd.read_csv(csv)

        # Polyphest: per-dataset lowest threshold
        poly = df[df['method'].isin(POLYPHEST_THRESHOLDS) & (df['metric'] == metric)].copy()
        merged = []
        for key, group in poly.groupby(['network', 'replicate']):
            for method in POLYPHEST_THRESHOLDS:
                rows = group[(group['method'] == method) & (group['status'] == 'SUCCESS')]
                if len(rows) > 0:
                    merged.append(rows.iloc[0])
                    break
        if merged:
            result = pd.DataFrame(merged)
            print(f"  {config} [Polyphest]:   mean={result['value'].mean():.3f}  (n={len(result)})")
        else:
            print(f"  {config} [Polyphest]:   no data")

        # GRAMPAIter: straightforward
        gi = df[(df['method'] == 'grandma_split') & (df['metric'] == metric) & (df['status'] == 'SUCCESS')]
        if len(gi) > 0:
            print(f"  {config} [GRAMPAIter]: mean={gi['value'].mean():.3f}  (n={len(gi)})")
        else:
            print(f"  {config} [GRAMPAIter]: no data")
