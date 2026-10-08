#!/usr/bin/env python3
"""
empirical_stats_section33.py - paired statistics for the Polyphest vs GRAMPAIter
empirical comparison (manuscript section "Comparison of Polyphest and GRAMPAIter
across all empirical datasets").

Answers the PI comments on that section:
  1. polyploid counts per method: means and a paired test
  2. reticulation counts per method: means and a paired test
  3. "reticulations at deeper nodes": mean number of descendant species per
     reticulation, per method, with a paired test
  4. whether the mu-distance tracks the reticulation measures across datasets:
     Spearman correlations

All paired tests are two-sided Wilcoxon signed-rank tests over the datasets on
which both methods produced a network. Sign counts are printed alongside.

Usage (on the cluster, after run_analysis.py has written the summary):
    python scripts/empirical_stats_section33.py
    python scripts/empirical_stats_section33.py --summary-dir analysis/real_data_summary
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr, wilcoxon

scripts_dir = Path(__file__).parent.resolve()
sys.path.insert(0, str(scripts_dir))

M1, M2 = 'polyphest', 'grandma_split'
NAME = {M1: 'Polyphest', M2: 'GRAMPAIter'}
TRIO = ['mu_distance', 'ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist']
LABEL = {'mu_distance': 'mu-distance',
         'ret_leaf_jaccard.dist': 'ret. descendants',
         'ret_sisters_jaccard.dist': 'ret. sister'}


def import_from_file(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Missing {path}. Run run_analysis.py first.")
    return pd.read_csv(path)


def pair_column(per_network: pd.DataFrame, metric: str) -> str:
    for a, b in ((M1, M2), (M2, M1)):
        col = f'{a}_vs_{b}_{metric}'
        if col in per_network.columns:
            return col
    raise SystemExit(f"No {M1}/{M2} column for {metric} in per_network_comparisons.csv")


def mean_descendants_per_reticulation(network_path: str, engine) -> float:
    rt = engine.load_network(network_path)
    if rt is None:
        raise SystemExit(f"Could not load network {network_path}")
    leaves = rt.get_reticulation_leaves()
    if not leaves:
        return float('nan')
    sizes = [len(set(v)) for v in leaves.values()]
    return sum(sizes) / len(sizes)


def paired_report(df: pd.DataFrame, col1: str, col2: str, label: str, w):
    both = df[[col1, col2]].dropna()
    n = len(both)
    diff = both[col2] - both[col1]
    w(f'--- {label} (n = {n} datasets) ---')
    w(f'  {NAME[M1]:<11} mean {both[col1].mean():7.2f}   median {both[col1].median():7.2f}')
    w(f'  {NAME[M2]:<11} mean {both[col2].mean():7.2f}   median {both[col2].median():7.2f}')
    w(f'  {NAME[M2]} higher on {(diff > 0).sum()}, lower on {(diff < 0).sum()}, '
      f'equal on {(diff == 0).sum()}')
    if (diff != 0).sum() < 2:
        w('  Wilcoxon: not computable (fewer than two non-zero differences)')
    else:
        stat, p = wilcoxon(both[col1], both[col2], alternative='two-sided')
        w(f'  Wilcoxon signed-rank, two-sided: W = {stat:.1f}, p = {p:.4f}')
    w('')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--summary-dir', default='analysis/real_data_summary')
    parser.add_argument('--cache-dir', default=None,
                        help='network cache dir for ComparisonEngine '
                             '(default: <summary-dir>/cache)')
    args = parser.parse_args()

    d = Path(args.summary_dir)
    stats = require(d / 'method_network_stats.csv')
    per_network = require(d / 'per_network_comparisons.csv')
    inventory = require(d / 'inventory.csv')

    out = []
    w = out.append

    # ---- per-method counts, paired over datasets where both methods exist
    piv_poly = stats.pivot_table(index='network', columns='method',
                                 values='num_polyploids', aggfunc='first')
    piv_ret = stats.pivot_table(index='network', columns='method',
                                values='reticulation_count', aggfunc='first')
    for piv, what in ((piv_poly, 'num_polyploids'), (piv_ret, 'reticulation_count')):
        missing = {M1, M2} - set(piv.columns)
        if missing:
            raise SystemExit(f"{what}: no rows for {missing} in method_network_stats.csv")

    paired_report(piv_poly, M1, M2, 'Number of polyploid species', w)
    paired_report(piv_ret, M1, M2, 'Number of reticulation events', w)

    # ---- depth: mean descendant species per reticulation, from the network files
    compute_comparisons = import_from_file('compute_comparisons',
                                           scripts_dir / 'compute_comparisons.py')
    cache_dir = args.cache_dir or str(d / 'cache')
    engine = compute_comparisons.ComparisonEngine(cache_dir)

    avail = inventory[(inventory['exists']) & (inventory['method'].isin([M1, M2]))]
    depth_rows = []
    for _, r in avail.iterrows():
        depth_rows.append({'network': r['network'], 'method': r['method'],
                           'mean_desc': mean_descendants_per_reticulation(
                               r['network_path'], engine)})
    depth = pd.DataFrame(depth_rows).pivot_table(index='network', columns='method',
                                                 values='mean_desc', aggfunc='first')
    paired_report(depth, M1, M2,
                  'Mean descendant species per reticulation (larger = deeper)', w)

    # ---- per-dataset table and correlations between the three measures
    cols = {m: pair_column(per_network, m) for m in TRIO}
    sub = per_network[['network'] + list(cols.values())].dropna(
        how='all', subset=list(cols.values()))
    w(f'--- Measures per dataset, {NAME[M1]} vs {NAME[M2]} ---')
    w(f'  {"dataset":<24}' + ''.join(f'{LABEL[m]:>18}' for m in TRIO))
    for _, r in sub.sort_values(cols['mu_distance']).iterrows():
        w(f'  {r["network"]:<24}' + ''.join(
            f'{"n/a":>18}' if pd.isna(r[cols[m]]) else f'{r[cols[m]]:>18.3f}'
            for m in TRIO))
    w('')
    for m in TRIO:
        s = sub[cols[m]].dropna()
        w(f'  {LABEL[m]:<18} mean {s.mean():.2f}  median {s.median():.2f}  '
          f'min {s.min():.2f}  max {s.max():.2f}  (n = {len(s)})')
    w('')
    w('--- Spearman correlations across datasets ---')
    for a, b in (('mu_distance', 'ret_leaf_jaccard.dist'),
                 ('mu_distance', 'ret_sisters_jaccard.dist'),
                 ('ret_leaf_jaccard.dist', 'ret_sisters_jaccard.dist')):
        pair = sub[[cols[a], cols[b]]].dropna()
        rho, p = spearmanr(pair[cols[a]], pair[cols[b]])
        w(f'  {LABEL[a]:<18} vs {LABEL[b]:<18} rho = {rho:6.3f}, p = {p:.4f}, '
          f'n = {len(pair)}')

    text = '\n'.join(out)
    print(text)
    (d / 'empirical_stats_section33.txt').write_text(text + '\n')
    print(f'\nSaved to {d / "empirical_stats_section33.txt"}')


if __name__ == '__main__':
    main()
