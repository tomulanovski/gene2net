#!/usr/bin/env python3
"""
Tests for the mu-distance entries that pairwise_compare reports.

A pair whose networks fall outside the class of Theorem 1 must not take the
whole comparison down with it, because the other metrics for that pair are
still valid. Instead the mu columns carry NaN and mu_status records the reason,
so failures are counted and visible in the results rather than skipped.
"""

import math
import importlib.util
import sys
from pathlib import Path

import networkx as nx
import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from compare_reticulations import pairwise_compare
from reticulate_tree import ReticulateTree

requires_reference = pytest.mark.skipif(
    importlib.util.find_spec('phylonetwork') is None,
    reason='mu-vectors come from phylonetwork, which is not installed here'
)


HEXAPLOID = "(((A,X),(B,X)),(C,X));"
TETRAPLOID = "(((A,X),(B,X)),C);"
RENAMED = "(((A,X),(B,X)),(Z,X));"


def tree(newick):
    return ReticulateTree(newick, is_multree=True)


@requires_reference
def test_pairwise_compare_reports_mu_columns():
    metrics = pairwise_compare(tree(HEXAPLOID), tree(HEXAPLOID))
    assert {'mu_distance', 'mu_distance_raw', 'mu_scored'} <= set(metrics)


@requires_reference
def test_edit_distance_is_no_longer_reported():
    """
    Commented out in pairwise_compare, superseded by mu_distance. The method
    itself is untouched on ReticulateTree and still callable.
    """
    metrics = pairwise_compare(tree(HEXAPLOID), tree(HEXAPLOID))
    assert 'edit_distance_multree' not in metrics
    assert tree(HEXAPLOID).get_edit_distance_multree(tree(HEXAPLOID)) == 0.0


@requires_reference
def test_identical_networks_score_zero_and_ok():
    metrics = pairwise_compare(tree(HEXAPLOID), tree(HEXAPLOID))
    assert metrics['mu_distance'] == 0.0
    assert metrics['mu_distance_raw'] == 0.0
    assert metrics['mu_scored'] == 1.0


@requires_reference
def test_different_networks_score_above_zero():
    metrics = pairwise_compare(tree(HEXAPLOID), tree(TETRAPLOID))
    assert metrics['mu_distance'] > 0.0
    assert metrics['mu_scored'] == 1.0


@requires_reference
def test_raw_distance_is_the_published_unnormalized_count():
    metrics = pairwise_compare(tree(HEXAPLOID), tree(TETRAPLOID))
    raw = metrics['mu_distance_raw']
    assert raw == int(raw) and raw > metrics['mu_distance']


@requires_reference
def test_mismatched_taxa_report_a_reason_without_losing_other_metrics():
    metrics = pairwise_compare(tree(HEXAPLOID), tree(RENAMED))

    assert math.isnan(metrics['mu_distance'])
    assert metrics['mu_scored'] == 0.0
    # the rest of the comparison is unaffected and still reported
    assert metrics['rf_distance'] is not None
    assert 'num_rets_diff' in metrics


@requires_reference
def test_network_outside_the_class_reports_the_reason():
    G = nx.DiGraph()
    G.add_edges_from([('root', 'p1'), ('root', 'p2'), ('root', 'p3'),
                      ('p1', 'r1'), ('p2', 'r1'), ('p3', 'r2'), ('r1', 'r2'),
                      ('r2', 'X'), ('p1', 'A'), ('p2', 'B'), ('p3', 'C')])
    for leaf in ['A', 'B', 'C', 'X']:
        G.nodes[leaf]['label'] = leaf
    stacked = ReticulateTree(G)

    metrics = pairwise_compare(stacked, stacked)

    assert math.isnan(metrics['mu_distance'])
    assert metrics['mu_scored'] == 0.0


@requires_reference
def test_every_reported_metric_is_numeric():
    """
    Metrics land in a long-format 'value' column that aggregate_replicates runs
    mean/std/min/max over. A string there raises TypeError and takes the whole
    summary down, so nothing reported may be text.
    """
    for pair in [(tree(HEXAPLOID), tree(HEXAPLOID)), (tree(HEXAPLOID), tree(RENAMED))]:
        for name, value in pairwise_compare(*pair).items():
            flat = value.values() if isinstance(value, dict) else [value]
            for v in flat:
                assert isinstance(v, (int, float)), f'{name} is not numeric: {v!r}'


@requires_reference
def test_unscorable_pair_warns_rather_than_passing_silently():
    """A pair we could not score must say so, not just leave a NaN behind."""
    import warnings
    import compare_reticulations
    compare_reticulations._MU_REPORTED.clear()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        pairwise_compare(tree(HEXAPLOID), tree(RENAMED))
    assert any('taxa' in str(w.message) for w in caught)


@requires_reference
def test_scored_flag_aggregates_to_a_success_fraction():
    """mu_scored must survive numeric aggregation so losses are countable."""
    import pandas as pd
    rows = [pairwise_compare(tree(HEXAPLOID), tree(HEXAPLOID))['mu_scored'],
            pairwise_compare(tree(HEXAPLOID), tree(RENAMED))['mu_scored']]
    assert pd.Series(rows).mean() == 0.5
