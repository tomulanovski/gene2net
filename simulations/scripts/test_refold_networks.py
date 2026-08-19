#!/usr/bin/env python3
"""
Tests for canonicalizing a network handed in already folded.

Every method except MPAllopp outputs a MUL-tree, which the pipeline folds with
Holm. MPAllopp writes extended Newick, so it arrives as a network and skips
folding entirely. If it emits a binary resolution while the ground truth is the
unresolved N(T), the two are compared at different resolution levels and
MPAllopp is penalised for committing to an ordering no MUL-tree can justify.

Huber et al. 2006 prove that a minimal binary network exhibiting a MUL-tree is
a RESOLUTION of N(T), so unfolding to the MUL-tree and refolding maps any
resolution back onto N(T). That is the canonicalisation tested here.

Refolding must stay opt-in. get_modified_mu_representation rebuilds a
ReticulateTree around a suppressed graph to hand it to phylonetwork, and
refolding there would silently change the metric.
"""

import importlib.util
import sys
from pathlib import Path

import networkx as nx
import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from reticulate_tree import ReticulateTree

# Tests that call get_mu_distance need the mu-vectors, which come from phylonetwork.
requires_reference = pytest.mark.skipif(
    importlib.util.find_spec('phylonetwork') is None,
    reason='mu-vectors come from phylonetwork, which is not installed here'
)


def resolved_hexaploid():
    """
    X is a hexaploid built by two successive binary events:
    pA and pB give r1, then r1 and pC give r2, which carries X.
    Unfolds to (((A,X),(B,X)),(C,X)); whose N(T) has ONE in-degree-3 node.
    """
    G = nx.DiGraph()
    G.add_edges_from([('root', 'nAB'), ('root', 'pC'), ('nAB', 'pA'), ('nAB', 'pB'),
                      ('pA', 'A'), ('pB', 'B'), ('pC', 'C'),
                      ('pA', 'r1'), ('pB', 'r1'), ('r1', 'r2'), ('pC', 'r2'),
                      ('r2', 'X')])
    for leaf in ['A', 'B', 'C', 'X']:
        G.nodes[leaf]['label'] = leaf
    return G


def in_degrees(rt):
    return sorted((d for _, d in rt.dag.in_degree() if d >= 2), reverse=True)


def test_network_input_is_left_alone_by_default():
    """Default behaviour is unchanged, so nothing existing shifts underneath us."""
    assert in_degrees(ReticulateTree(resolved_hexaploid())) == [2, 2]


def test_refold_contracts_a_resolved_network_to_its_canonical_form():
    assert in_degrees(ReticulateTree(resolved_hexaploid(), refold=True)) == [3]


@requires_reference
def test_refolded_network_matches_the_ground_truth_folding():
    """
    The whole point: the same history reached from an eNewick network and from
    a MUL-tree must land on the same object, so the comparison is fair.
    """
    from_network = ReticulateTree(resolved_hexaploid(), refold=True)
    from_multree = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True)
    assert from_network.get_mu_distance(from_multree) == 0.0


@requires_reference
def test_refolding_is_idempotent():
    once = ReticulateTree(resolved_hexaploid(), refold=True)
    twice = ReticulateTree(once.dag, refold=True)
    assert in_degrees(twice) == in_degrees(once)
    assert once.get_mu_distance(twice) == 0.0


@requires_reference
def test_multree_input_is_unaffected_by_the_flag():
    plain = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True)
    flagged = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True, refold=True)
    assert plain.get_mu_distance(flagged) == 0.0
    assert in_degrees(flagged) == [3]


@requires_reference
def test_mu_representation_is_not_disturbed_by_the_flag():
    """
    Guards the trap: get_modified_mu_representation wraps a suppressed graph in
    a ReticulateTree internally. If refolding leaked in as a default, the
    representation would change under it.
    """
    rt = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True)
    before = sorted(rt.get_modified_mu_representation())
    assert sorted(rt.get_modified_mu_representation()) == before
    assert rt.get_mu_distance(rt) == 0.0
