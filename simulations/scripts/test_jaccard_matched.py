#!/usr/bin/env python3
"""
Tests for the unpenalised companions to the reticulation Jaccard metrics.

ret_leaf_jaccard and ret_sisters_jaccard normalise by max(len1, len2) and charge
for every reticulation that found no partner, so a method that recovers three of
ten reticulations perfectly still scores badly. That answers "did it find the
reticulations". It does not answer "were the ones it found right".

The *_matched columns answer the second question by normalising over matched
pairs only. They use partial_match=True for EVERY method, so unlike
ret_leaf_jaccard - which is partial only for the methods in
SINGLE_RETICULATION_METHODS - they mean the same thing in every row and can be
compared across methods.

The penalised columns are untouched, so no existing number moves.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from compare_reticulations import pairwise_compare
from reticulate_tree import ReticulateTree

# pairwise_compare computes mu, which needs the reference implementation.
pytestmark = pytest.mark.skipif(
    importlib.util.find_spec('phylonetwork') is None,
    reason='pairwise_compare needs phylonetwork for the mu columns'
)

# X and Y are both tetraploid -> two reticulations
TRUTH = "(((A,X),(B,X)),((C,Y),(D,Y)));"
# X tetraploid, Y not recovered -> one reticulation, one unmatched in the truth
PARTIAL = "(((A,X),(B,X)),(C,(D,Y)));"

MATCHED = ['ret_leaf_jaccard_matched', 'ret_sisters_jaccard_matched']


def tree(newick):
    return ReticulateTree(newick, is_multree=True)


def test_matched_columns_are_reported():
    metrics = pairwise_compare(tree(TRUTH), tree(PARTIAL))
    assert all(name in metrics for name in MATCHED)


def test_fixture_really_has_an_unmatched_reticulation():
    assert tree(TRUTH).get_reticulation_count() == 2
    assert tree(PARTIAL).get_reticulation_count() == 1


def test_matched_scores_better_than_penalised_when_one_is_missed():
    """The whole point: no charge for the reticulation that was never found."""
    m = pairwise_compare(tree(TRUTH), tree(PARTIAL))
    assert m['ret_leaf_jaccard_matched']['dist'] < m['ret_leaf_jaccard']['dist']


def test_penalised_columns_are_unchanged_for_a_single_reticulation_method():
    """partial_match still drives the original column exactly as before."""
    plain = pairwise_compare(tree(TRUTH), tree(PARTIAL), partial_match=False)
    carved = pairwise_compare(tree(TRUTH), tree(PARTIAL), partial_match=True)
    assert carved['ret_leaf_jaccard']['dist'] < plain['ret_leaf_jaccard']['dist']


def test_matched_column_ignores_the_partial_match_flag():
    """
    It must mean the same thing for every method, otherwise it is no more
    comparable across rows than the penalised column it exists to supplement.
    """
    plain = pairwise_compare(tree(TRUTH), tree(PARTIAL), partial_match=False)
    carved = pairwise_compare(tree(TRUTH), tree(PARTIAL), partial_match=True)
    for name in MATCHED:
        assert plain[name] == carved[name]


def test_matched_columns_are_numeric():
    metrics = pairwise_compare(tree(TRUTH), tree(PARTIAL))
    for name in MATCHED:
        for v in metrics[name].values():
            assert isinstance(v, (int, float))
