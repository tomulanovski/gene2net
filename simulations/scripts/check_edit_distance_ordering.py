#!/usr/bin/env python3
"""
check_edit_distance_ordering.py - Does Newick child order affect our edit distances?

READ-ONLY DIAGNOSTIC. Reads the MUL-trees already saved by postprocess_results.py
and recomputes `edit_distance_multree` two ways:

  ed_current    exactly what ReticulateTree.get_edit_distance_multree does today
                (nodes indexed in tree.traverse() order = the child order the
                 program happened to write into the file)

  ed_canonical  identical in every other respect, but children are visited in
                canonical order (sorted by the recursive canonical form of their
                subtree) so the result cannot depend on file ordering

Why this matters: `next(nx.optimize_graph_edit_distance(...))` returns the FIRST
complete edit path a greedy search finds, not the minimum. That search enumerates
candidates in node-insertion order, so the approximation error depends on how well
the two trees' child orderings happen to line up. Child order is phylogenetically
meaningless, and each program writes its own. The exact GED would be unaffected;
only the approximation is.

Nothing is modified: no pipeline file is touched, no .tre file is rewritten.

Usage:
    # from an inventory produced by collect_results.py (preferred)
    python check_edit_distance_ordering.py --inventory inventory.csv

    # or straight from summary_config.yaml
    python check_edit_distance_ordering.py conf_ils_low_10M

    # limit the work while you sanity-check it
    python check_edit_distance_ordering.py conf_ils_low_10M --methods grampa polyphest_p50 --replicates 1
"""

import argparse
import math
import os
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml
from ete3 import Tree

# Matches ReticulateTree.get_edit_distance_multree
MAX_NODES_FOR_GED = 500
DEFAULT_TIMEOUT = 300

NODE_MATCH = lambda u, v: u.get('label') == v.get('label')


# ───── graph construction ─────

def build_graph(tree, canonical):
    """ete3 tree -> DiGraph, labels on leaves only (internal nodes get None).

    canonical=False reproduces the current tree_to_simple_graph exactly: nodes are
    inserted in tree.traverse() order, which follows the child order of the input
    Newick.

    canonical=True visits children sorted by the recursive canonical form of their
    subtree, so insertion order depends only on topology. Does not modify the tree.
    """
    G = nx.DiGraph()

    if not canonical:
        for node in tree.traverse():
            G.add_node(id(node), label=node.name if node.is_leaf() else None)
            if not node.is_root():
                G.add_edge(id(node.up), id(node))
        return G

    canon = {}

    def canonical_form(node):
        if node.is_leaf():
            form = node.name
        else:
            form = '(' + ','.join(sorted(canonical_form(c) for c in node.children)) + ')'
        canon[node] = form
        return form

    canonical_form(tree)

    def add_subtree(node):
        G.add_node(id(node), label=node.name if node.is_leaf() else None)
        if not node.is_root():
            G.add_edge(id(node.up), id(node))
        # ties (identical sibling subtrees, i.e. autopolyploidy) keep input order,
        # which is harmless because the subtrees are identical
        for child in sorted(node.children, key=lambda c: canon[c]):
            add_subtree(child)

    add_subtree(tree)
    return G


# ───── the metric ─────

def _alarm_supported():
    return hasattr(__import__('signal'), 'SIGALRM')


def edit_distance(t1, t2, canonical, timeout=DEFAULT_TIMEOUT):
    """Same metric, same normalization, same node cap as the pipeline.

    Returns (value, status). value is NaN when skipped or timed out.
    """
    g1 = build_graph(t1, canonical)
    g2 = build_graph(t2, canonical)

    n1, n2 = len(g1.nodes), len(g2.nodes)
    if max(n1, n2) > MAX_NODES_FOR_GED:
        return float('nan'), f'skipped_too_large({n1},{n2})'

    use_alarm = timeout and _alarm_supported()
    if use_alarm:
        import signal

        def _handler(signum, frame):
            raise TimeoutError()

        old = signal.signal(signal.SIGALRM, _handler)
        signal.alarm(timeout)

    try:
        distance = next(nx.optimize_graph_edit_distance(g1, g2, node_match=NODE_MATCH))
    except TimeoutError:
        return float('nan'), 'timeout'
    except Exception as e:
        return float('nan'), f'error({type(e).__name__})'
    finally:
        if use_alarm:
            import signal
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)

    normalization = max(len(g1.nodes) + len(g1.edges), len(g2.nodes) + len(g2.edges))
    if normalization > 0:
        distance /= normalization
    return distance, 'ok'


def load_tree(path):
    try:
        with open(path) as f:
            newick = f.read().strip()
        if not newick:
            return None
        return Tree(newick, format=1)
    except Exception:
        return None


# ───── inventory ─────

def inventory_from_config(config, config_path, methods_filter, replicates):
    cfg = yaml.safe_load(open(config_path))
    base_dir = Path(cfg['base_dir'])
    networks_dir = Path(cfg['networks_dir'])
    n_reps = replicates or cfg['num_replicates']

    methods = cfg['methods']
    if methods_filter:
        methods = {k: v for k, v in methods.items() if k in methods_filter}

    rows = []
    for method, mcfg in methods.items():
        for network in cfg['networks']:
            for rep in range(1, n_reps + 1):
                if method == 'padre':
                    inf = (base_dir / network / 'results' / config /
                           mcfg['directory'] / f'replicate_{rep}' / mcfg['output_file'])
                else:
                    inf = (base_dir / network / 'results' / config /
                           mcfg['directory'] / f'replicate_{rep}' / mcfg['output_file'])
                gt = networks_dir / f'{network}.tre'
                rows.append({
                    'network': network, 'config': config, 'method': method,
                    'replicate': rep, 'gt_path': str(gt), 'inferred_path': str(inf),
                    'gt_exists': gt.exists(), 'inferred_exists': inf.exists(),
                })
    return pd.DataFrame(rows)


# ───── main ─────

def main():
    p = argparse.ArgumentParser(
        description='Check whether Newick child ordering affects edit_distance_multree',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('configuration', nargs='?', help='e.g. conf_ils_low_10M')
    p.add_argument('--inventory', help='inventory CSV from collect_results.py (preferred)')
    p.add_argument('--config', default='simulations/summary_config.yaml')
    p.add_argument('--methods', nargs='+', help='subset of methods')
    p.add_argument('--replicates', type=int, help='use only replicates 1..N')
    p.add_argument('--timeout', type=int, default=DEFAULT_TIMEOUT, help='seconds per GED (POSIX only)')
    p.add_argument('--out', default=None, help='output CSV (default: edit_ordering_{config}.csv)')
    args = p.parse_args()

    if args.inventory:
        inv = pd.read_csv(args.inventory)
        if args.methods:
            inv = inv[inv['method'].isin(args.methods)]
        if args.replicates:
            inv = inv[inv['replicate'] <= args.replicates]
        config = str(inv['config'].iloc[0]) if 'config' in inv.columns and len(inv) else 'inventory'
    elif args.configuration:
        config = args.configuration
        inv = inventory_from_config(config, args.config, args.methods, args.replicates)
    else:
        p.error('provide a configuration name or --inventory')

    valid = inv[inv['gt_exists'] & inv['inferred_exists']].copy()
    print(f'{len(valid)} of {len(inv)} combinations have both trees present\n')
    if valid.empty:
        print('Nothing to compare. Check base_dir in the config, or pass --inventory.')
        return

    out_path = Path(args.out or f'edit_ordering_{config}.csv')
    results = []
    gt_cache = {}

    for i, (_, row) in enumerate(valid.iterrows(), 1):
        gt_path = row['gt_path']
        if gt_path not in gt_cache:
            gt_cache[gt_path] = load_tree(gt_path)
        t_true = gt_cache[gt_path]
        t_inf = load_tree(row['inferred_path'])

        if t_true is None or t_inf is None:
            results.append({**{k: row[k] for k in ('network', 'method', 'replicate')},
                            'ed_current': float('nan'), 'ed_canonical': float('nan'),
                            'delta': float('nan'), 'status': 'parse_failed'})
            continue

        cur, s1 = edit_distance(t_true, t_inf, canonical=False, timeout=args.timeout)
        can, s2 = edit_distance(t_true, t_inf, canonical=True, timeout=args.timeout)

        results.append({
            'network': row['network'], 'method': row['method'], 'replicate': row['replicate'],
            'n_nodes_gt': len(list(t_true.traverse())), 'n_nodes_inf': len(list(t_inf.traverse())),
            'ed_current': cur, 'ed_canonical': can,
            'delta': (can - cur) if (not math.isnan(cur) and not math.isnan(can)) else float('nan'),
            'status': s1 if s1 != 'ok' else s2,
        })

        if i % 25 == 0 or i == len(valid):
            print(f'  {i}/{len(valid)}', flush=True)
            pd.DataFrame(results).to_csv(out_path, index=False)

    df = pd.DataFrame(results)
    df.to_csv(out_path, index=False)
    print(f'\nPer-pair results written to {out_path}\n')

    ok = df.dropna(subset=['ed_current', 'ed_canonical'])
    if ok.empty:
        print('No pair produced both values.')
        return

    summary = (ok.groupby('method')
                 .agg(n=('ed_current', 'size'),
                      ed_current=('ed_current', 'mean'),
                      ed_canonical=('ed_canonical', 'mean'))
                 .assign(delta=lambda d: d['ed_canonical'] - d['ed_current'])
                 .sort_values('ed_current'))

    summary['rank_current'] = summary['ed_current'].rank().astype(int)
    summary['rank_canonical'] = summary['ed_canonical'].rank().astype(int)

    print('=' * 78)
    print(f'Per-method mean edit distance  ({config})')
    print('=' * 78)
    print(summary.to_string(float_format=lambda v: f'{v:.4f}'))
    print('=' * 78)

    flipped = summary[summary['rank_current'] != summary['rank_canonical']]
    best_cur = summary['ed_current'].idxmin()
    best_can = summary['ed_canonical'].idxmin()

    print(f'\nBest method under current metric:   {best_cur}')
    print(f'Best method under canonical order:  {best_can}')
    if best_cur != best_can:
        print('\n  *** THE WINNER CHANGES. The current ranking is affected by Newick')
        print('      child ordering, which carries no phylogenetic meaning.')
    elif len(flipped):
        print(f'\n  Winner unchanged, but {len(flipped)} method(s) move in the ranking:')
        print('  ' + ', '.join(flipped.index))
    else:
        print('\n  Ranking is identical. Your existing ordering was not a confound;')
        print('  the fix is then a one-sentence methods addition, not a revision.')

    shift = ok['delta'].abs().mean()
    print(f'\nMean absolute change per pair: {shift:.4f}')
    print(f'Pairs where canonical is lower (tighter bound): '
          f'{(ok["delta"] < 0).sum()}/{len(ok)}')


if __name__ == '__main__':
    main()
