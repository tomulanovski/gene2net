#!/usr/bin/env python3
"""
report_section33_numbers.py - Pull the manuscript Section 3.3 numbers out of the
empirical summary tables.

Section 3.3 of the manuscript quotes three sets of values that changed when the
approximate MUL-tree edit distance was replaced by the normalized mu-distance.
Hunting them down by hand across four CSVs is error prone, so this prints them
in one block, labelled by the sentence they belong to.

Run it after jobs/resummary_empirical_mu.sh has finished:

    python scripts/report_section33_numbers.py
    python scripts/report_section33_numbers.py --summary-dir analysis/real_data_summary

Nothing here is silent. A missing file, a missing method pair or an absent
mu_distance column raises, because a quietly skipped number would end up as a
wrong number in the paper.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Manuscript names for the pipeline's method keys
DISPLAY = {
    'grampa': 'GRAMPA',
    'grandma_split': 'GRAMPAIter',
    'polyphest': 'Polyphest',
    'padre': 'PADRE',
    'mpallop': 'MPAllopp',
    'alloppnet': 'AlloppNET',
}

POPP = 'Popp_2005'
EPHEDRA = 'Wu_2015'


def require(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(
            f"Missing {path}.\n"
            "Run jobs/resummary_empirical_mu.sh on the cluster first. It must use "
            "--force-recompute, because the comparison cache is keyed on network "
            "file hashes only and will otherwise return pre-mu results."
        )
    return pd.read_csv(path)


def find_pair(df: pd.DataFrame, m1: str, m2: str) -> pd.DataFrame:
    """Rows for a method pair, in whichever order the pipeline stored it."""
    hit = df[((df['method1'] == m1) & (df['method2'] == m2)) |
             ((df['method1'] == m2) & (df['method2'] == m1))]
    if hit.empty:
        available = sorted({tuple(sorted(p)) for p in
                            zip(df['method1'], df['method2'])})
        raise SystemExit(
            f"No comparison found for {m1} vs {m2}.\n"
            f"Pairs present: {available}"
        )
    return hit


def pair_column(per_network: pd.DataFrame, m1: str, m2: str, metric: str) -> str:
    """Column name for a pair in the wide per-network table, either order."""
    for a, b in ((m1, m2), (m2, m1)):
        col = f'{a}_vs_{b}_{metric}'
        if col in per_network.columns:
            return col
    raise SystemExit(
        f"No column for {m1} vs {m2} / {metric} in per_network_comparisons.csv.\n"
        f"Columns carrying {metric}: "
        f"{[c for c in per_network.columns if c.endswith(metric)]}"
    )


def fmt(x) -> str:
    return 'NOT SCORED' if pd.isna(x) else f'{x:.2f}'


def main():
    parser = argparse.ArgumentParser(
        description='Report manuscript Section 3.3 numbers under the mu-distance')
    parser.add_argument('--summary-dir', default='analysis/real_data_summary',
                        help='Directory written by run_analysis.py '
                             '(default: analysis/real_data_summary)')
    args = parser.parse_args()

    d = Path(args.summary_dir)
    pairwise = require(d / 'pairwise_summary.csv')
    per_network = require(d / 'per_network_comparisons.csv')
    mu_report = require(d / 'mu_scoring_report.csv')

    if 'mu_distance' not in set(pairwise['metric']):
        raise SystemExit(
            f"{d / 'pairwise_summary.csv'} has no mu_distance rows, only "
            f"{sorted(set(pairwise['metric']))}.\n"
            "These summaries predate the mu-distance switch. Re-run "
            "run_analysis.py with --force-recompute."
        )

    out = []
    w = out.append

    w('=' * 78)
    w('Manuscript Section 3.3 numbers under the normalized mu-distance')
    w(f'Source: {d.resolve()}')
    w('=' * 78)

    # ---------------------------------------------------------------- 3.3.1 Popp
    w('')
    w('--- 3.3.1, Popp et al. paragraph -------------------------------------------')
    w('Sentence: "They disagreed on the origins of S. sorensenis and S. involucrata,')
    w('           reflected in a high reticulation sister measure (X) and')
    w('           mu-distance (Y)."   [Polyphest vs PADRE on Popp_2005]')
    w('')
    if POPP not in set(per_network['network']):
        w(f'  !! {POPP} absent from per_network_comparisons.csv. Networks present:')
        w(f'     {sorted(per_network["network"].unique())}')
    else:
        row = per_network[per_network['network'] == POPP].iloc[0]
        mu_col = pair_column(per_network, 'polyphest', 'padre', 'mu_distance')
        sis_col = pair_column(per_network, 'polyphest', 'padre',
                              'ret_sisters_jaccard.dist')
        scored_col = pair_column(per_network, 'polyphest', 'padre', 'mu_scored')
        w(f'  reticulation sister measure : {fmt(row[sis_col])}   (was 0.62)')
        w(f'  mu-distance                 : {fmt(row[mu_col])}   (edit distance was 0.85)')
        if row[scored_col] != 1.0:
            w('  !! mu_scored is not 1.0 for this pair, so the mu-distance could not')
            w('     be computed here and the sentence needs rewording, not a new number.')

    # ------------------------------------------------------------- 3.3.1 Ephedra
    w('')
    w('--- 3.3.1, Ephedra paragraph -----------------------------------------------')
    w('Sentence: "mu-distances were high across all pairs (LO to HI)."')
    w(f'          [range over every method pair on {EPHEDRA}]')
    w('')
    eph_cols = [c for c in per_network.columns if c.endswith('_mu_distance')]
    if EPHEDRA not in set(per_network['network']):
        w(f'  !! {EPHEDRA} absent from per_network_comparisons.csv.')
    else:
        row = per_network[per_network['network'] == EPHEDRA].iloc[0]
        vals = {c[:-len('_mu_distance')]: row[c] for c in eph_cols}
        present = {k: v for k, v in vals.items() if not pd.isna(v)}
        missing = sorted(k for k, v in vals.items() if pd.isna(v))
        for pair, v in sorted(present.items(), key=lambda kv: kv[1]):
            w(f'    {pair:<34} {v:.4f}')
        if present:
            w('')
            w(f'  range : {min(present.values()):.2f} to {max(present.values()):.2f}'
              f'   (edit distance was 0.76 to 0.92)')
            w(f'  pairs scored: {len(present)} of {len(vals)}')
        if missing:
            w(f'  !! not scored: {", ".join(missing)}')

    # ---------------------------------------------------------------- 3.3.2
    w('')
    w('--- 3.3.2, across all datasets ---------------------------------------------')
    w('Sentence: "mu-distances were uniformly high, from LO to HI across the N')
    w('           datasets for which they could be computed."')
    w('          [Polyphest vs GRAMPAIter across datasets]')
    w('')
    mu_rows = pairwise[pairwise['metric'] == 'mu_distance']
    pair = find_pair(mu_rows, 'polyphest', 'grandma_split').iloc[0]
    w(f'  mu-distance   min {pair["min"]:.2f}   max {pair["max"]:.2f}   '
      f'mean {pair["mean"]:.2f}   sd {pair["std"]:.2f}')
    w(f'  datasets scored (n_networks): {int(pair["n_networks"])}'
      '   (edit distance covered 10)')

    scoring = find_pair(mu_report, 'polyphest', 'grandma_split').iloc[0]
    w(f'  from mu_scoring_report: {int(scoring["n_scored"])} of '
      f'{int(scoring["n_pairs"])} scored '
      f'({100 * scoring["frac_scored"]:.0f}%)')
    unscored = scoring['networks_unscored']
    if isinstance(unscored, str) and unscored:
        w(f'  !! not scored on: {unscored.replace(";", ", ")}')
        w('     Section 3.3.2 must say which datasets these are, or why they drop out.')
    else:
        w('  every dataset with both results was scored')

    # unchanged reticulation measures, reprinted so the rest of 3.3.2 can be checked
    w('')
    w('  unchanged by the metric switch, reprinted to confirm:')
    for metric, was in (('ret_leaf_jaccard.dist', '0.33 to 0.95'),
                        ('ret_sisters_jaccard.dist', '0.72 to 0.97')):
        sub = pairwise[pairwise['metric'] == metric]
        if sub.empty:
            w(f'    !! {metric} absent from pairwise_summary.csv')
            continue
        r = find_pair(sub, 'polyphest', 'grandma_split').iloc[0]
        w(f'    {metric:<26} {r["min"]:.2f} to {r["max"]:.2f}   (was {was})')

    # ------------------------------------------------------- overall coverage
    w('')
    w('--- mu-distance coverage, all pairs ----------------------------------------')
    w('')
    for _, r in mu_report.iterrows():
        n1 = DISPLAY.get(r['method1'], r['method1'])
        n2 = DISPLAY.get(r['method2'], r['method2'])
        note = ''
        if isinstance(r['networks_unscored'], str) and r['networks_unscored']:
            note = f'   unscored: {r["networks_unscored"].replace(";", ", ")}'
        w(f'  {n1:<12} vs {n2:<12} {int(r["n_scored"])}/{int(r["n_pairs"])}{note}')

    w('')
    w('=' * 78)

    text = '\n'.join(out)
    print(text)
    return text


if __name__ == '__main__':
    main()
