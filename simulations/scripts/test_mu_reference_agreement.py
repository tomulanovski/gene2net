#!/usr/bin/env python3
"""
Cross-validates our mu-vectors against the reference implementation.

We do NOT call `phylonetwork` for the distance. Its `distance()` computes the
EXTENDED mu-representation of Cardona et al. 2024, which encodes only BINARY
orchard networks. Figure 1 of arXiv:2412.05107 exhibits two non-isomorphic
semi-binary stack-free networks that share an extended mu-representation, and
MUL-tree folding produces semi-binary networks with reticulations of in-degree
up to 6, so that metric can report 0 for genuinely different networks here.

The metric we report is the MODIFIED mu-representation of arXiv:2412.05107,
which prepends each node's in-degree to its standard mu-vector and ranges over
all nodes. The package does not implement it.

What the package can still do is validate the part we share: the standard
mu-vectors, mu_i(v) = number of paths from v to leaf i. These tests check ours
against theirs on every ground truth. They skip when phylonetwork is not
installed, so they run on the cluster conda environments and are inert
elsewhere.

Two further reference behaviours are documented here because they are the
reasons our guards exist.
"""

import glob
import os
import sys
from pathlib import Path

import pytest

scripts_dir = Path(__file__).parent
sys.path.insert(0, str(scripts_dir))

from reticulate_tree import ReticulateTree

NETWORKS_DIR = scripts_dir.parent / 'networks'

pn = pytest.importorskip('phylonetwork',
                         reason='reference implementation not installed')


def as_reference(rt):
    """Our folded network, handed to the reference implementation."""
    return pn.PhylogeneticNetwork(eNewick=rt.to_enewick())


@pytest.mark.parametrize('path', sorted(glob.glob(str(NETWORKS_DIR / '*.tre'))),
                         ids=lambda p: os.path.basename(p).replace('.tre', ''))
def test_standard_mu_vectors_match_reference(path):
    """
    Our path counts must equal the reference's. Compared as a sorted multiset
    of vectors, since node identities differ between the two representations.
    """
    rt = ReticulateTree(open(path, encoding='utf-8').read().strip(), is_multree=True)
    reference = as_reference(rt)

    ours = sorted(tuple(v) for v in rt._mu_path_counts(rt.dag)[0].values())
    theirs = sorted(tuple(int(x) for x in vector)
                    for vector in reference.mu_dict.values())

    assert ours == theirs


def test_our_taxon_order_matches_reference():
    """mu-vector coordinates are only comparable under a shared taxon order."""
    rt = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True)
    assert rt._mu_taxa() == list(as_reference(rt).ordered_taxa)


def test_reference_extended_representation_cannot_see_in_degree():
    """
    Why we do not use the package's distance.

    The extended mu-vector of a reticulation is dropped from its representation
    entirely, so nothing in it records that a reticulation has three parents
    rather than two. Our modified representation keeps every node and leads
    with the in-degree.
    """
    rt = ReticulateTree("(((A,X),(B,X)),(C,X));", is_multree=True)
    reference = as_reference(rt)

    n_reticulations = sum(1 for _, d in rt.dag.in_degree() if d >= 2)
    assert n_reticulations == 1
    assert len(reference.emu_representation) == rt.dag.number_of_nodes() - 1

    in_degrees_we_record = {vector[0] for vector in rt.get_modified_mu_representation()}
    assert 3 in in_degrees_we_record


def test_reference_scores_renamed_taxon_as_identical():
    """
    Why we reject differing taxon sets.

    Renaming B to C keeps the renamed taxon in the same position of the sorted
    taxa list, so every vector lines up and two networks on {A,B,X} and {A,C,X}
    are reported as distance 0 despite describing different species. A rename
    that shifts the sort order instead yields an arbitrary nonzero count.
    """
    original = ReticulateTree("((A,X),(B,X));", is_multree=True)
    renamed = ReticulateTree("((A,X),(C,X));", is_multree=True)

    assert as_reference(original).distance(as_reference(renamed)) == 0

    with pytest.raises(ValueError, match='taxa'):
        original.get_mu_distance(renamed)
