#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

POLYPHEST_THRESHOLDS = ['polyphest_p50', 'polyphest_p70', 'polyphest_p90']
BASE = Path('/groups/itay_mayrose/tomulanovski/gene2net/simulations/analysis/summary')

CONFIGS = [
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
METRIC_LABELS = {
    'edit_distance_multree': 'Edit Distance',
    'ret_leaf_jaccard.dist': 'Ret. Descendants',
    'ret_sisters_jaccard.dist': 'Ret. Sister',
}

for config in CONFIGS:
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

    # --- Polyphest: completion rates ---
    print("  [Polyphest] Completion rates:")
    for t in POLYPHEST_THRESHOLDS:
        rows = inv[inv['method'] == t]
        if len(rows) == 0:
            print(f"    {t}: not found")
            continue
        rate = rows['inferred_exists'].mean() * 100
        print(f"    {t}: {rate:.1f}%  ({rows['inferred_exists'].sum()}/{len(rows)})")

    # --- Polyphest: per-threshold accuracy ---
    print("  [Polyphest] Accuracy:")
    for metric in METRICS:
        label = METRIC_LABELS[metric]
        parts = []
        for t in POLYPHEST_THRESHOLDS:
            rows = df[(df['method'] == t) & (df['metric'] == metric) & (df['status'] == 'SUCCESS')]
            if len(rows) > 0:
                parts.append(f"{t.split('_')[-1]}={rows['value'].mean():.3f}(n={len(rows)})")
            else:
                parts.append(f"{t.split('_')[-1]}=n/a")
        print(f"    {label}: {' | '.join(parts)}")

    # --- GRAMPAIter (no prior): accuracy ---
    print("  [GRAMPAIter] Accuracy:")
    gi_rows = df[df['method'] == 'grandma_split']
    for metric in METRICS:
        label = METRIC_LABELS[metric]
        rows = gi_rows[(gi_rows['metric'] == metric) & (gi_rows['status'] == 'SUCCESS')]
        if len(rows) > 0:
            print(f"    {label}: {rows['value'].mean():.3f}  (n={len(rows)})")
        else:
            print(f"    {label}: no data")

    # --- GRAMPAIter (with prior): completion + accuracy ---
    prior_inv = inv[inv['method'] == 'grandma_split_prior']
    if len(prior_inv) > 0:
        rate = prior_inv['inferred_exists'].mean() * 100
        print(f"  [GRAMPAIter+prior] Completion: {rate:.1f}%  ({prior_inv['inferred_exists'].sum()}/{len(prior_inv)})")
    print("  [GRAMPAIter+prior] Accuracy:")
    gp_rows = df[df['method'] == 'grandma_split_prior']
    for metric in METRICS:
        label = METRIC_LABELS[metric]
        rows = gp_rows[(gp_rows['metric'] == metric) & (gp_rows['status'] == 'SUCCESS')]
        if len(rows) > 0:
            print(f"    {label}: {rows['value'].mean():.3f}  (n={len(rows)})")
        else:
            print(f"    {label}: no data")
