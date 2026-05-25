#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

POLYPHEST_THRESHOLDS = ['polyphest_p50', 'polyphest_p70', 'polyphest_p90']
BASE = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/summary')

configs = [
    'conf_ils_low_10M',
    'conf_ils_medium_10M',
    'conf_ils_high_10M',
    'conf_dup_loss_low_10M',
    'conf_dup_loss_medium_10M',
    'conf_dup_loss_high_10M',
    'conf_dup_loss_low_10M_ne1M',
    'conf_dup_loss_medium_10M_ne1M',
    'conf_dup_loss_high_10M_ne1M',
    'conf_dup_loss_low_10M_ne2M',
    'conf_dup_loss_medium_10M_ne2M',
    'conf_dup_loss_high_10M_ne2M',
]

METRICS = ['edit_distance_multree', 'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist']

for config in configs:
    inv_path = BASE / config / 'inventory.csv'
    comp_path = BASE / config / 'comparisons_raw.csv'
    if not comp_path.exists():
        print(f"\n{config}: files not found")
        continue

    inv = pd.read_csv(inv_path)
    df = pd.read_csv(comp_path)
    print(f"\n{'='*60}")
    print(f"{config}")
    print(f"{'='*60}")

    # Completion rates per threshold
    print("  Completion rates:")
    for t in POLYPHEST_THRESHOLDS:
        rows = inv[inv['method'] == t]
        if len(rows) == 0:
            print(f"    {t}: not found")
            continue
        rate = rows['inferred_exists'].mean() * 100
        print(f"    {t}: {rate:.1f}%  ({rows['inferred_exists'].sum()}/{len(rows)})")

    # Per-threshold accuracy
    for metric in METRICS:
        print(f"  {metric}:")
        for t in POLYPHEST_THRESHOLDS:
            rows = df[(df['method'] == t) & (df['metric'] == metric) & (df['status'] == 'SUCCESS')]
            if len(rows) > 0:
                print(f"    {t}: {rows['value'].mean():.3f}  (n={len(rows)})")
            else:
                print(f"    {t}: no data")
