#!/usr/bin/env python3
"""
Tests for the extended mu-distance on folded networks.

Reference: Cardona, Pons, Ribas, Martinez Coronado (2024), "Comparison of
orchard networks using their extended mu-representation", IEEE/ACM TCBB.
Semi-binary extension: arXiv:2412.05107 (Advances in Applied Mathematics).

The ground truth MUL-trees fold to SEMI-BINARY networks, with reticulations of
in-degree up to 6. The extended mu-representation of Cardona et al. encodes only
BINARY orchard networks - arXiv:2412.05107 Figure 1 exhibits two non-isomorphic
semi-binary stack-free networks sharing one. So the metric here is the MODIFIED
mu-representation of that paper:

    mu-bar(v) = ( in_degree(v), mu_1(v), ..., mu_n(v) )       over ALL nodes
    d(N1, N2) = |mu-bar(N1) symmetric-difference mu-bar(N2)|

Theorem 1: for N1 semi-binary stack-free orchard and N2 semi-binary stack-free,
mu-bar(N1) = mu-bar(N2) if and only if N1 is isomorphic to N2.
"""

import glob
import os
import importlib.util
import sys
from pathlib import Path

import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from reticulate_tree import ReticulateTree

requires_reference = pytest.mark.skipif(
    importlib.util.find_spec('phylonetwork') is None,
    reason='mu-vectors come from phylonetwork, which is not installed here'
)


NETWORKS_DIR = scripts_dir.parent / 'networks'

# X is a hexaploid with subgenomes in the A, B and C lineages.
# Folding contracts the two successive events into one in-degree-3 reticulation.
HEXAPLOID = "(((A,X),(B,X)),(C,X));"
# Same taxa {A,B,C,X}, but X is only a tetraploid from A and B.
TETRAPLOID = "(((A,X),(B,X)),C);"
# Same shape as HEXAPLOID, but taxon C is renamed to Z.
HEXAPLOID_RENAMED = "(((A,X),(B,X)),(Z,X));"
# Same shape as HEXAPLOID plus an extra taxon.
HEXAPLOID_EXTRA_TAXON = "(((A,X),(B,X)),((C,D),X));"


def tree(newick):
    return ReticulateTree(newick, is_multree=True)


@requires_reference
def test_identical_network_has_zero_distance():
    t = tree(TETRAPLOID)
    assert t.get_mu_distance(tree(TETRAPLOID)) == 0.0


@requires_reference
def test_identical_semibinary_network_has_zero_distance():
    """The in-degree-3 reticulation must not break the representation."""
    t = tree(HEXAPLOID)
    reticulation_in_degrees = [d for _, d in t.dag.in_degree() if d >= 2]
    assert reticulation_in_degrees == [3], "fixture no longer exercises semi-binary"
    assert t.get_mu_distance(tree(HEXAPLOID)) == 0.0


@requires_reference
def test_different_networks_have_positive_distance():
    assert tree(HEXAPLOID).get_mu_distance(tree(TETRAPLOID)) > 0.0


@requires_reference
def test_distance_is_symmetric():
    a, b = tree(HEXAPLOID), tree(TETRAPLOID)
    assert a.get_mu_distance(b) == b.get_mu_distance(a)


@requires_reference
def test_renamed_taxon_raises_instead_of_scoring_zero():
    """
    The reference implementation indexes mu-vectors by each network's own
    sorted taxa, so two networks of identical shape on {A,B,C,X} and {A,B,Z,X}
    score 0 - a perfect match for different species. That must be rejected.
    """
    with pytest.raises(ValueError, match='taxa'):
        tree(HEXAPLOID).get_mu_distance(tree(HEXAPLOID_RENAMED))


@requires_reference
def test_different_taxon_counts_raise():
    with pytest.raises(ValueError, match='taxa'):
        tree(HEXAPLOID).get_mu_distance(tree(HEXAPLOID_EXTRA_TAXON))


@requires_reference
def test_normalized_distance_is_in_unit_interval():
    d = tree(HEXAPLOID).get_mu_distance(tree(TETRAPLOID), normalize=True)
    assert 0.0 < d <= 1.0


@requires_reference
def test_unnormalized_distance_is_a_whole_count():
    d = tree(HEXAPLOID).get_mu_distance(tree(TETRAPLOID), normalize=False)
    assert d == int(d) and d > 0


@requires_reference
def test_normalization_divides_by_internal_nodes_only():
    """
    The denominator is the internal nodes of both networks, not every node.

    Every leaf has in-degree 1 by the semi-binary precondition, so its modified
    mu-vector is (1, e_x) in both networks and can never enter the symmetric
    difference. Including the 2n leaf vectors in the denominator would cap the
    distance at 1 - 2n/(|V1|+|V2|), a bound that varies with each pair's
    leaf-to-node ratio, so distances on differently shaped networks would sit on
    different scales and could not be averaged. This pins the fix.
    """
    a, b = tree(HEXAPLOID), tree(TETRAPLOID)
    mu_a, mu_b = a.get_modified_mu_representation(), b.get_modified_mu_representation()
    n_taxa = len(a._mu_taxa(a._mu_network()))

    raw = a.get_mu_distance(b, normalize=False)
    normalized = a.get_mu_distance(b, normalize=True)

    internal = (len(mu_a) - n_taxa) + (len(mu_b) - n_taxa)
    assert normalized == pytest.approx(raw / internal)

    # The old all-node denominator was strictly larger, so it reported a
    # strictly smaller distance and could never reach 1.
    assert normalized > raw / (len(mu_a) + len(mu_b))


@requires_reference
@pytest.mark.parametrize('path', sorted(glob.glob(str(NETWORKS_DIR / '*.tre'))),
                         ids=lambda p: os.path.basename(p).replace('.tre', ''))
def test_ground_truth_has_zero_self_distance(path):
    """Every real ground truth must be at distance 0 from itself."""
    newick = open(path, encoding='utf-8').read().strip()
    assert tree(newick).get_mu_distance(tree(newick)) == 0.0


@requires_reference
def test_modified_mu_representation_covers_every_node():
    """Unlike the extended representation, reticulations are NOT excluded."""
    t = tree(HEXAPLOID)
    assert len(t.get_modified_mu_representation()) == t.dag.number_of_nodes()


@requires_reference
def test_modified_mu_vector_leads_with_in_degree():
    """
    The first coordinate is the in-degree. The in-degree-3 reticulation of
    HEXAPLOID carries a single copy of X below it, so its vector is (3,0,0,0,1)
    against taxa [A,B,C,X]. A binary resolution would put 2 there instead, which
    is what makes mu-bar separate semi-binary networks that mu conflates.
    """
    t = tree(HEXAPLOID)
    assert (3, 0, 0, 0, 1) in t.get_modified_mu_representation()


@requires_reference
def test_stacked_reticulations_raise():
    """Theorem 1 assumes stack-free; a reticulation under a reticulation is not."""
    import networkx as nx
    G = nx.DiGraph()
    G.add_edges_from([('root', 'p1'), ('root', 'p2'), ('root', 'p3'),
                      ('p1', 'r1'), ('p2', 'r1'), ('p3', 'r2'), ('r1', 'r2'),
                      ('r2', 'X'), ('p1', 'A'), ('p2', 'B'), ('p3', 'C')])
    for leaf in ['A', 'B', 'C', 'X']:
        G.nodes[leaf]['label'] = leaf
    stacked = ReticulateTree(G)
    with pytest.raises(ValueError, match='stack-free'):
        stacked.get_mu_distance(stacked)


@requires_reference
def test_mu_vector_of_root_is_the_copy_number_vector():
    """
    mu_i(root) counts root-to-leaf-i paths, which is the number of copies of
    taxon i in the MUL-tree. For HEXAPLOID that is A:1 B:1 C:1 X:3.
    """
    t = tree(HEXAPLOID)
    taxa, root_vector = t.get_mu_root_vector()
    assert dict(zip(taxa, root_vector)) == {'A': 1, 'B': 1, 'C': 1, 'X': 3}


@requires_reference
@pytest.mark.parametrize('path', sorted(glob.glob(str(NETWORKS_DIR / '*.tre'))),
                         ids=lambda p: os.path.basename(p).replace('.tre', ''))
def test_distance_is_independent_of_newick_child_order(path):
    """
    Characterisation test for a property the metric has by construction, and
    the reason no canonical ordering step is needed here.

    mu-bar is built from path counts and in-degrees, which are properties of
    the graph rather than of any traversal, and the taxa are sorted, so the
    coordinates are canonical already. The comparison is between multisets, so
    it cannot depend on the order children happen to appear in the input file.

    This is exactly what the graph edit distance did not have. Shuffling
    children moved raw-order GED to as much as 0.67 on these same networks.
    """
    import random
    from ete3 import Tree

    newick = open(path, encoding='utf-8').read().strip()
    shuffled = Tree(newick, format=1)
    rnd = random.Random(1)
    for node in shuffled.traverse():
        if not node.is_leaf():
            rnd.shuffle(node.children)

    original = ReticulateTree(newick, is_multree=True)
    reordered = ReticulateTree(shuffled.write(format=1), is_multree=True)

    assert original.get_mu_distance(reordered) == 0.0


def test_mu_vectors_require_the_reference_implementation(monkeypatch):
    """
    The mu-vectors come from phylonetwork. If it is missing we say so, never
    silently substitute our own path counting.
    """
    monkeypatch.setitem(sys.modules, 'phylonetwork', None)
    with pytest.raises(ImportError, match='phylonetwork'):
        tree(HEXAPLOID).get_modified_mu_representation()


@requires_reference
def test_internal_node_names_do_not_affect_the_distance():
    """
    Regression. phylonetwork indexes mu-vectors by every node LABEL, not only
    leaves, so a network whose internal nodes carry names produced longer
    vectors. No tuple could then compare equal and the distance pinned at
    exactly 1.0 while the taxon guard, which looks at leaves only, saw nothing
    wrong. That is how mpsugar scored 1.0000 on every network.
    """
    plain = ReticulateTree("((A,X),(B,X));", is_multree=True)
    named = ReticulateTree("((A,X)ancestor1,(B,X)ancestor2);", is_multree=True)

    assert plain._mu_taxa() == named._mu_taxa()
    assert plain.get_mu_distance(named) == 0.0


@requires_reference
def test_vector_width_is_one_plus_the_leaf_taxa():
    """The only coordinates are the in-degree and one per leaf taxon."""
    rt = ReticulateTree("((A,X)ancestor1,(B,X)ancestor2);", is_multree=True)
    width = len(rt._mu_taxa()) + 1
    assert all(len(v) == width for v in rt.get_modified_mu_representation())
