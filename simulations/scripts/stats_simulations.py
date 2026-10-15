#!/usr/bin/env python3
"""
stats_simulations.py - significance tests for the simulation results.

Every test is a two-sided Wilcoxon signed-rank test paired by network. Replicate
values are first averaged per network, as in the figures, because replicates of
one network are not independent. Holm correction is applied within each block.

Blocks (one per manuscript claim):
  A  Polyphest vs GRAMPAIter in each of the 12 baseline configurations
     (section 3.3.2).
  B  Polyphest vs GRAMPAIter at each retention level of the fractionation
     series (fractionation paragraph, "at 25% retention ...").
  C  Sensitivity: per network, the slope of each measure across a series of
     configurations, compared between Polyphest and GRAMPAIter. Series are
     duplication/loss None -> High at each ILS level (x = 0, 1, 2, 3), and
     retention 1.0 -> 0.25 (x = 1 - retention). Only networks with a value at
     every level of the series for both methods are used.
  D  Tetraploid subset: MPAllopp vs GRAMPAIter on the 8 AlloppNET-compatible
     networks, values averaged per network over replicates and the three
     ILS-only configurations, as in Table 4.

Usage (on the cluster, after the summaries carry mu_distance):
    cd simulations/scripts
    python stats_simulations.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from create_analysis_figures import merge_polyphest_comparisons

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"
OUT_DIR = SUMMARY_BASE / "significance_tests"

NAME = {'polyphest': 'Polyphest', 'grandma_split': 'GRAMPAIter',
        'mpsugar': 'MPAllopp'}
TRIO = [('mu_distance', 'mu-distance'),
        ('ret_leaf_jaccard.dist', 'ret. descendants'),
        ('ret_sisters_jaccard.dist', 'ret. sister')]

BASELINE = [
    'conf_ils_low_10M', 'conf_ils_medium_10M', 'conf_ils_high_10M',
    'conf_dup_loss_low_10M', 'conf_dup_loss_medium_10M', 'conf_dup_loss_high_10M',
    'conf_dup_loss_low_10M_ne1M', 'conf_dup_loss_medium_10M_ne1M', 'conf_dup_loss_high_10M_ne1M',
    'conf_dup_loss_low_10M_ne2M', 'conf_dup_loss_medium_10M_ne2M', 'conf_dup_loss_high_10M_ne2M',
]
# (x value, config). Retention series uses x = 1 - retention, so a positive
# slope means the measure rises (accuracy drops) as fractionation increases.
FRACTIONATION = [(0.00, 'conf_dup_loss_medium_10M_ne1M'),
                 (0.25, 'conf_dup_loss_medium_10M_ne1M_fix075'),
                 (0.50, 'conf_dup_loss_medium_10M_ne1M_fix050'),
                 (0.75, 'conf_dup_loss_medium_10M_ne1M_fix025')]
DUPLOSS_SERIES = {
    'dup/loss, low ILS': [(0, 'conf_ils_low_10M'), (1, 'conf_dup_loss_low_10M'),
                          (2, 'conf_dup_loss_medium_10M'), (3, 'conf_dup_loss_high_10M')],
    'dup/loss, medium ILS': [(0, 'conf_ils_medium_10M'), (1, 'conf_dup_loss_low_10M_ne1M'),
                             (2, 'conf_dup_loss_medium_10M_ne1M'), (3, 'conf_dup_loss_high_10M_ne1M')],
    'dup/loss, high ILS': [(0, 'conf_ils_high_10M'), (1, 'conf_dup_loss_low_10M_ne2M'),
                           (2, 'conf_dup_loss_medium_10M_ne2M'), (3, 'conf_dup_loss_high_10M_ne2M')],
}
# Never simulated under fractionation; dropped from blocks B and C2, including
# the 1.0 baseline, so every retention level covers the same networks.
FRACTIONATION_EXCLUDE_NETWORKS = ['Morales-Briones_2021']
TETRAPLOID_NETWORKS = ['Bendiksby_2011', 'Ding_2023', 'Koenen_2020', 'Liu_2023',
                       'Shahrestani_2015', 'Wisecaver_2023', 'Wu_2015', 'Zhao_2021']
ILS_ONLY = ['conf_ils_low_10M', 'conf_ils_medium_10M', 'conf_ils_high_10M']


def load(configs):
    frames = []
    for cfg in configs:
        f = SUMMARY_BASE / cfg / 'comparisons_raw.csv'
        if not f.exists():
            raise SystemExit(f"Missing {f}")
        df = pd.read_csv(f)
        df['config'] = cfg
        frames.append(df)
    comp = merge_polyphest_comparisons(pd.concat(frames, ignore_index=True))
    comp = comp[comp['status'] == 'SUCCESS']
    absent = [k for k, _ in TRIO if k not in set(comp['metric'])]
    if absent:
        raise SystemExit(f"No rows for {absent}. Re-run the mu re-summary first.")
    return comp


def per_network(comp, configs, metric):
    """network x method table of means over replicates (and over the configs given)."""
    d = comp[comp['config'].isin(configs) & (comp['metric'] == metric)]
    return d.groupby(['network', 'method'])['value'].mean().unstack('method')


def holm(pvals):
    p = np.asarray(pvals, dtype=float)
    adj = np.empty(len(p))
    running = 0.0
    for rank, i in enumerate(np.argsort(p)):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(running, 1.0)
    return adj


def paired(a, b, label):
    """Wilcoxon on two aligned Series (same networks). Fails loud if untestable."""
    both = pd.concat([a, b], axis=1, keys=['a', 'b']).dropna()
    diff = both['b'] - both['a']
    if (diff != 0).sum() < 2:
        raise SystemExit(f"{label}: fewer than two non-zero differences")
    stat, p = wilcoxon(both['a'], both['b'], alternative='two-sided')
    return {'n': len(both), 'mean_1': both['a'].mean(), 'mean_2': both['b'].mean(),
            'first_lower_on': int((diff > 0).sum()), 'W': stat, 'p': p}


def report(res, title, name1, name2, w):
    res['p_holm'] = holm(res['p'])
    w(f'=== {title} ({len(res)} tests, Holm within block) ===')
    w(f'  {"test":<46}{"measure":<18}{"n":>3}{name1:>12}{name2:>12}'
      f'{"1 lower on":>11}{"W":>7}{"p":>10}{"p_Holm":>10}')
    for _, r in res.iterrows():
        w(f'  {r["test"]:<46}{r["measure"]:<18}{r["n"]:>3}{r["mean_1"]:>12.3f}'
          f'{r["mean_2"]:>12.3f}{f"{r.first_lower_on}/{r.n}":>11}{r["W"]:>7.1f}'
          f'{r["p"]:>10.2e}{r["p_holm"]:>10.2e}')
    w(f'  significant after Holm: {(res["p_holm"] < 0.05).sum()}/{len(res)}, '
      f'max p_Holm {res["p_holm"].max():.2e}')
    w('')
    return res.assign(block=title)


def block_per_config(comp, configs, m1, m2, title, w):
    rows = []
    for cfg in configs:
        for key, label in TRIO:
            piv = per_network(comp, [cfg], key)
            for m in (m1, m2):
                if m not in piv.columns:
                    raise SystemExit(f"{cfg} {key}: no rows for {m}")
            rows.append({'test': cfg, 'measure': label,
                         **paired(piv[m1], piv[m2], f'{cfg} {key}')})
    return report(pd.DataFrame(rows), title, NAME[m1], NAME[m2], w)


def slopes(comp, series, metric, method):
    """Per-network OLS slope of the measure against x, networks complete on the series."""
    cols = {}
    for x, cfg in series:
        piv = per_network(comp, [cfg], metric)
        if method not in piv.columns:
            raise SystemExit(f"{cfg} {metric}: no rows for {method}")
        cols[x] = piv[method]
    tab = pd.DataFrame(cols).dropna()
    xs = np.array(list(cols.keys()), dtype=float)
    return tab.apply(lambda row: np.polyfit(xs, row.values, 1)[0], axis=1)


def block_slopes(comp, series_map, m1, m2, title, w):
    rows = []
    for sname, series in series_map.items():
        for key, label in TRIO:
            rows.append({'test': f'slope, {sname}', 'measure': label,
                         **paired(slopes(comp, series, key, m1),
                                  slopes(comp, series, key, m2), f'{sname} {key}')})
    return report(pd.DataFrame(rows), title, NAME[m1], NAME[m2], w)


def block_tetraploid(comp, m1, m2, title, w):
    rows = []
    for key, label in TRIO:
        piv = per_network(comp, ILS_ONLY, key)
        piv = piv[piv.index.isin(TETRAPLOID_NETWORKS)]
        for m in (m1, m2):
            if m not in piv.columns:
                raise SystemExit(f"tetraploid {key}: no rows for {m}")
        rows.append({'test': 'tetraploid subset, ILS-only configs', 'measure': label,
                     **paired(piv[m1], piv[m2], f'tetraploid {key}')})
    return report(pd.DataFrame(rows), title, NAME[m1], NAME[m2], w)


def main():
    out = []
    w = out.append
    w('Means are over networks. "1 lower on" counts networks where the first method')
    w('has the lower value. For slope tests the means are mean slopes per step of x.')
    w('')

    comp = load(BASELINE + [c for _, c in FRACTIONATION[1:]])
    P, G, M = 'polyphest', 'grandma_split', 'mpsugar'
    excluded = comp['network'].isin(FRACTIONATION_EXCLUDE_NETWORKS)
    if (excluded & comp['config'].isin([c for _, c in FRACTIONATION[1:]])).any():
        raise SystemExit(f"{FRACTIONATION_EXCLUDE_NETWORKS} have fractionation results, "
                         f"but are listed as never run")
    comp_frac = comp[~excluded]
    blocks = [
        block_per_config(comp, BASELINE, P, G, 'A. Baseline configurations', w),
        block_per_config(comp_frac, [c for _, c in FRACTIONATION], P, G,
                         'B. Fractionation levels (Morales-Briones excluded)', w),
        block_slopes(comp, DUPLOSS_SERIES, P, G,
                     'C1. Sensitivity to duplication and loss (x = level 0..3)', w),
        block_slopes(comp_frac, {'retention 1.0 -> 0.25': FRACTIONATION}, P, G,
                     'C2. Sensitivity to fractionation (x = 1 - retention)', w),
        block_tetraploid(comp, M, G, 'D. Tetraploid subset, MPAllopp vs GRAMPAIter', w),
    ]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(blocks, ignore_index=True).to_csv(OUT_DIR / 'simulation_tests.csv', index=False)
    text = '\n'.join(out)
    print(text)
    (OUT_DIR / 'simulation_tests.txt').write_text(text + '\n')
    print(f'Saved to {OUT_DIR}')


if __name__ == '__main__':
    main()
