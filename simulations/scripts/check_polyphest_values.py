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

for metric in ['ret_leaf', 'ret_sisters']:
    print(f"\n--- {metric} ---")
    for config in configs:
        csv = BASE / config / 'comparisons_raw.csv'
        if not csv.exists():
            print(f"  {config}: file not found")
            continue
        df = pd.read_csv(csv)
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
            print(f"  {config}: mean={result['value'].mean():.3f}  (n={len(result)})")
        else:
            print(f"  {config}: no data")
