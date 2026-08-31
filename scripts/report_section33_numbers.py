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


# The three measures the manuscript reports side by side. All are distances, so
# lower is closer, which is easy to invert when reading "highest agreement".
TRIO = ['mu_distance', 'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist']


def pairs_present(per_network: pd.DataFrame, metrics) -> list:
    """Method-pair prefixes appearing in the wide table for any of `metrics`."""
    names = set()
    for col in per_network.columns:
        for metric in metrics:
            suffix = f'_{metric}'
            if col.endswith(suffix):
                names.add(col[:-len(suffix)])
    return sorted(names)


def pretty_pair(pair: str) -> str:
    """'polyphest_vs_padre' -> 'Polyphest vs PADRE', for manuscript names."""
    if '_vs_' not in pair:
        return pair
    a, b = pair.split('_vs_', 1)
    return f'{DISPLAY.get(a, a)} vs {DISPLAY.get(b, b)}'


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
    median = f'{pair["median"]:.2f}' if 'median' in pair else 'RERUN FOR MEDIAN'
    w(f'  mu-distance   min {pair["min"]:.2f}   max {pair["max"]:.2f}   '
      f'mean {pair["mean"]:.2f}   median {median}   sd {pair["std"]:.2f}')
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
    for metric, was in (('ret_leaf_jaccard.dist', '0.33 to 0.95, median 0.73'),
                        ('ret_sisters_jaccard.dist', '0.72 to 0.97, median 0.87')):
        sub = pairwise[pairwise['metric'] == metric]
        if sub.empty:
            w(f'    !! {metric} absent from pairwise_summary.csv')
            continue
        r = find_pair(sub, 'polyphest', 'grandma_split').iloc[0]
        med = f', median {r["median"]:.2f}' if 'median' in r else ''
        w(f'    {metric:<26} {r["min"]:.2f} to {r["max"]:.2f}{med}   (was {was})')

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

    # ------------------------------------------- full Ephedra pair matrix
    # The paragraph claims Polyphest and PADRE had the HIGHEST agreement on
    # reticulation descendants. That is a claim over every pair, not only the
    # pairs involving Polyphest, and it is checkable by any referee in five
    # minutes, so print every pair and state outright whether it holds.
    w('')
    w(f'--- {EPHEDRA}, every method pair, all three measures ------------------------')
    w('    (all three are DISTANCES: lower means better agreement)')
    w('')
    pairs = pairs_present(per_network, TRIO)
    if EPHEDRA not in set(per_network['network']):
        w(f'  !! {EPHEDRA} absent from per_network_comparisons.csv')
    else:
        row = per_network[per_network['network'] == EPHEDRA].iloc[0]
        table = {}
        for p in pairs:
            vals = {}
            for metric in TRIO:
                col = f'{p}_{metric}'
                vals[metric] = row[col] if col in row.index else float('nan')
            if not all(pd.isna(v) for v in vals.values()):
                table[p] = vals

        best = {m: min((v[m] for v in table.values() if not pd.isna(v[m])),
                       default=float('nan')) for m in TRIO}
        w(f'  {"pair":<32}{"mu":>10}{"ret.desc":>10}{"ret.sis":>10}')
        w('  ' + '-' * 62)
        for p, vals in sorted(table.items(),
                              key=lambda kv: (kv[1][TRIO[1]]
                                              if not pd.isna(kv[1][TRIO[1]]) else 9)):
            line = f'  {pretty_pair(p):<32}'
            for metric in TRIO:
                v = vals[metric]
                mark = '*' if (not pd.isna(v) and v == best[metric]) else ' '
                line += (f'{"n/a":>9} ' if pd.isna(v) else f'{v:>9.3f}{mark}')
            w(line)
        w('  (* marks the best value in that column, i.e. the closest pair)')
        w('')
        for metric, label in zip(TRIO, ['mu-distance', 'reticulation descendants',
                                        'reticulation sister']):
            avail = {p: v[metric] for p, v in table.items() if not pd.isna(v[metric])}
            if not avail:
                continue
            winner = min(avail, key=avail.get)
            w(f'  closest pair on {label:<26}: {pretty_pair(winner)} '
              f'({avail[winner]:.3f})')
        desc_avail = {p: v[TRIO[1]] for p, v in table.items()
                      if not pd.isna(v[TRIO[1]])}
        if desc_avail:
            winner = min(desc_avail, key=desc_avail.get)
            claim = 'HOLDS' if set(winner.split('_vs_')) == {'polyphest', 'padre'} \
                    else 'DOES NOT HOLD'
            w('')
            w(f'  >> "highest agreement on reticulation descendants was between')
            w(f'      Polyphest and PADRE" across ALL pairs: {claim}')

    # ------------------------------- per-dataset Polyphest vs GRAMPAIter
    w('')
    w('--- Polyphest vs GRAMPAIter, per dataset -----------------------------------')
    w('')
    cols = {}
    for metric in TRIO:
        try:
            cols[metric] = pair_column(per_network, 'polyphest', 'grandma_split', metric)
        except SystemExit:
            cols[metric] = None
    sub = per_network[['network'] + [c for c in cols.values() if c]].copy()
    sub = sub.dropna(how='all', subset=[c for c in cols.values() if c])
    w(f'  {"dataset":<24}{"mu":>10}{"ret.desc":>10}{"ret.sis":>10}')
    w('  ' + '-' * 54)
    for _, r in sub.sort_values(cols[TRIO[0]] or 'network').iterrows():
        line = f'  {r["network"]:<24}'
        for metric in TRIO:
            c = cols[metric]
            v = r[c] if c else float('nan')
            line += f'{"n/a":>10}' if pd.isna(v) else f'{v:>10.3f}'
        w(line)
    w('')
    for metric, label in zip(TRIO, ['mu', 'ret.desc', 'ret.sis']):
        c = cols[metric]
        if not c:
            continue
        s = sub[['network', c]].dropna()
        lo, hi = s.loc[s[c].idxmin()], s.loc[s[c].idxmax()]
        w(f'  {label:<9} min {lo[c]:.3f} ({lo["network"]})   '
          f'max {hi[c]:.3f} ({hi["network"]})')
    # Do the global measure and the reticulation measures pick the same dataset
    # as the one the two methods agreed on most? If not, that is a real-data
    # instance of the framework argument and is worth a sentence.
    picks = {}
    for metric in TRIO:
        c = cols[metric]
        if c:
            s = sub[['network', c]].dropna()
            picks[metric] = s.loc[s[c].idxmin(), 'network']
    if len(set(picks.values())) > 1:
        w('')
        w('  >> the measures DISAGREE on which dataset the two methods agreed on most:')
        for metric, net in picks.items():
            w(f'       {metric:<26} -> {net}')
    elif picks:
        w(f'  >> all three measures agree the closest dataset is '
          f'{next(iter(picks.values()))}')

    # --------------------------------- reticulation counts per method
    stats_path = d / 'method_network_stats.csv'
    w('')
    w('--- inferred reticulation count, per method per dataset --------------------')
    w('')
    if not stats_path.exists():
        w(f'  !! {stats_path} not found, so reticulation counts cannot be reported.')
        w('     It is written by run_analysis.py phase 1b.')
    else:
        stats = pd.read_csv(stats_path)
        piv = stats.pivot_table(index='network', columns='method',
                                values='reticulation_count', aggfunc='first')
        piv = piv.rename(columns=DISPLAY)
        w(piv.to_string(na_rep='.'))
        if {'Polyphest', 'GRAMPAIter'} <= set(piv.columns):
            both = piv[['Polyphest', 'GRAMPAIter']].dropna()
            more = (both['GRAMPAIter'] > both['Polyphest']).sum()
            fewer = (both['GRAMPAIter'] < both['Polyphest']).sum()
            same = (both['GRAMPAIter'] == both['Polyphest']).sum()
            w('')
            w(f'  GRAMPAIter inferred MORE than Polyphest on {more}/{len(both)} '
              f'datasets, FEWER on {fewer}, equal on {same}.')
            w('  This is the sentence "neither method consistently inferred more')
            w('  reticulations than the other" - replace it with these counts.')

    w('')
    w('=' * 78)

    text = '\n'.join(out)
    print(text)
    return text


if __name__ == '__main__':
    main()
