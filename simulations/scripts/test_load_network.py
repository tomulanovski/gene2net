#!/usr/bin/env python3
"""
Tests for the loader the comparison engine actually uses.

ComparisonEngine has two paths. _compare_in_subprocess runs only when a
positive --timeout is given; the default is 0, so compare_pair -> load_network
is what executes in practice. Both must agree, or a fix applied to one silently
does nothing.
"""

import sys
import tempfile
from pathlib import Path

import networkx as nx
import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from compute_comparisons import ComparisonEngine
from reticulate_tree import ReticulateTree


def engine(tmp):
    return ComparisonEngine(str(tmp))


def write(tmp, name, text):
    p = Path(tmp) / name
    p.write_text(text, encoding='utf-8')
    return str(p)


def test_default_timeout_uses_the_load_network_path():
    """Guards the assumption the other tests rest on."""
    with tempfile.TemporaryDirectory() as tmp:
        assert engine(tmp).timeout == 0


def test_load_network_reads_utf8_regardless_of_locale():
    """A ground truth with an accented taxon must load, not raise."""
    with tempfile.TemporaryDirectory() as tmp:
        p = write(tmp, 'gt.tre', '((A,jalapa\u00ebnsis),(B,jalapa\u00ebnsis));')
        rt = engine(tmp).load_network(p)
        assert rt is not None
        assert 'jalapa\u00ebnsis' in rt._mu_taxa()


def test_load_network_refolds_a_network_given_as_enewick():
    """
    MPAllopp writes extended Newick. It must be refolded here too, or it is
    compared at a different resolution than every MUL-tree method.
    """
    resolved = ReticulateTree(_resolved_hexaploid(), refold=False)
    with tempfile.TemporaryDirectory() as tmp:
        p = write(tmp, 'inf.tre', resolved.to_enewick())
        loaded = engine(tmp).load_network(p)
        assert sorted((d for _, d in loaded.dag.in_degree() if d >= 2), reverse=True) == [3]


def test_load_network_reports_why_it_failed():
    """A load failure must carry its reason, not vanish into None."""
    with tempfile.TemporaryDirectory() as tmp:
        result = engine(tmp).compare_pair(
            write(tmp, 'gt.tre', 'this is not newick'),
            write(tmp, 'inf.tre', '((A,B),C);'), 'N', 'm')
        assert result['status'] == 'ERROR'
        assert len(result['error']) > len('Failed to load ground truth: ') + len(tmp)


def _resolved_hexaploid():
    G = nx.DiGraph()
    G.add_edges_from([('root', 'nAB'), ('root', 'pC'), ('nAB', 'pA'), ('nAB', 'pB'),
                      ('pA', 'A'), ('pB', 'B'), ('pC', 'C'),
                      ('pA', 'r1'), ('pB', 'r1'), ('r1', 'r2'), ('pC', 'r2'),
                      ('r2', 'X')])
    for leaf in ['A', 'B', 'C', 'X']:
        G.nodes[leaf]['label'] = leaf
    return G
