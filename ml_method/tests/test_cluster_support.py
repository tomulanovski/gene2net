"""Tests for the copy-aware cluster-support partner feature.

The signal: for a duplicated species, one copy stays "home" (near its own
species-tree clade) and the other lands "away" near its allopolyploid partner.
We measure, per away copy, what fraction of its local gene-tree neighborhood is
made of partner-clade species, accumulated into an [E, E, 2] (sum, max) tensor.

The ``sample_gene_trees`` fixture is exactly this case: species D is an
allotetraploid whose away copy sits next to B (the A/B clade) in 7 of 10 trees.
"""
import pytest
import torch

from gene2net_gnn.data.tree_io import tree_to_edge_index, reorder_edge_index_preorder
from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.features import (
    gene_tree_copy_neighborhoods,
    copy_aware_cluster_support,
    species_coclustering_matrix,
    edge_clades_species,
    pairwise_partner_features,
)
from gene2net_gnn.training.trainer_reconstruct import build_pairwise_feat

# Fixed species indexing used throughout these tests.
SP = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
N_SPECIES = 5


def _to_tensors(tree):
    """Convert one ETE3 gene tree to (edge_index, species_ids, leaf_mask),
    matching the on-disk format produced by Gene2NetSample.from_trees."""
    edge_index, names = tree_to_edge_index(tree)
    species_ids, leaf_mask = [], []
    for name in names:
        if name in SP:
            species_ids.append(SP[name])
            leaf_mask.append(True)
        else:
            species_ids.append(-1)
            leaf_mask.append(False)
    return (
        edge_index,
        torch.tensor(species_ids, dtype=torch.long),
        torch.tensor(leaf_mask, dtype=torch.bool),
    )


# Species-tree ((A,B),(C,(D,E))) edge clades, in a fixed order for the tests.
# 8 edges for a 5-leaf rooted binary tree (2n - 2).
EDGE_CLADES = [
    {SP["A"], SP["B"]},          # 0: (A,B)
    {SP["A"]},                   # 1: A
    {SP["B"]},                   # 2: B
    {SP["C"], SP["D"], SP["E"]}, # 3: (C,(D,E))
    {SP["C"]},                   # 4: C
    {SP["D"], SP["E"]},          # 5: (D,E)
    {SP["D"]},                   # 6: D
    {SP["E"]},                   # 7: E
]


def test_copy_neighborhoods_splits_duplicated_species(sample_gene_trees):
    """A two-copy gene tree yields one neighborhood per copy of the duplicated
    species; with k=2 each neighborhood is just the copy's sister leaf."""
    # ((A,(D,B)),(C,(D,E))): D has an away copy (sister B) and a home copy (sister E).
    ei, sp, lm = _to_tensors(sample_gene_trees[0])
    result = gene_tree_copy_neighborhoods(ei, sp, lm, k=2)

    # Only D is duplicated -> exactly two entries, both for species D.
    assert len(result) == 2
    assert all(a == SP["D"] for a, _ in result)

    neighborhoods = sorted(tuple(sorted(n.items())) for _, n in result)
    assert neighborhoods == sorted([((SP["B"], 1),), ((SP["E"], 1),)])


def test_copy_neighborhoods_ignores_non_duplicated(sample_gene_trees):
    """A gene tree with no duplicated species yields no neighborhoods."""
    # sample_gene_trees[7] is ((A,B),(C,(D,E))) -> every species appears once.
    ei, sp, lm = _to_tensors(sample_gene_trees[7])
    assert gene_tree_copy_neighborhoods(ei, sp, lm, k=2) == []


def test_cluster_support_shape(sample_gene_trees):
    tensors = [_to_tensors(t) for t in sample_gene_trees]
    feat = copy_aware_cluster_support(
        [t[0] for t in tensors], [t[1] for t in tensors], [t[2] for t in tensors],
        EDGE_CLADES, N_SPECIES, k=2, away_threshold=0.5,
    )
    assert feat.shape == (len(EDGE_CLADES), len(EDGE_CLADES), 2)
    assert torch.isfinite(feat).all()


def test_cluster_support_captures_allo_partner(sample_gene_trees):
    """For the (D,E) clade edge, the away copy of D routes its support to the
    (A,B) clade edge (its allo partner), and never to its own clade."""
    tensors = [_to_tensors(t) for t in sample_gene_trees]
    feat = copy_aware_cluster_support(
        [t[0] for t in tensors], [t[1] for t in tensors], [t[2] for t in tensors],
        EDGE_CLADES, N_SPECIES, k=2, away_threshold=0.5,
    )
    de, ab = 5, 0  # edge indices for (D,E) and (A,B)

    # 7 of 10 trees have an away D copy whose neighborhood is 100% the (A,B)
    # clade -> sum = 7/10, max = 1.0.
    assert feat[de, ab, 0].item() == pytest.approx(0.7, abs=1e-6)
    assert feat[de, ab, 1].item() == pytest.approx(1.0, abs=1e-6)

    # The home copy (sister E) is filtered out, so the (D,E) edge gets no
    # self-support from it.
    assert feat[de, de, 0].item() == pytest.approx(0.0, abs=1e-6)


def test_cluster_support_zero_for_non_duplicated_clade(sample_gene_trees):
    """An edge whose clade species are never duplicated contributes nothing."""
    tensors = [_to_tensors(t) for t in sample_gene_trees]
    feat = copy_aware_cluster_support(
        [t[0] for t in tensors], [t[1] for t in tensors], [t[2] for t in tensors],
        EDGE_CLADES, N_SPECIES, k=2, away_threshold=0.5,
    )
    a_edge = 1  # clade {A}; A is never duplicated
    assert feat[a_edge].sum().item() == pytest.approx(0.0, abs=1e-6)


def test_build_pairwise_feat_concatenates_cluster_support(
    simple_species_tree, sample_gene_trees
):
    """build_pairwise_feat augments the 2-dim co-clustering feature with the new
    2-dim cluster-support feature, giving [E, E, 4]; the two halves match the
    standalone functions."""
    species_list = ["A", "B", "C", "D", "E"]
    sample = Gene2NetSample.from_trees(simple_species_tree, sample_gene_trees, species_list)
    sample._edge_index_pre = reorder_edge_index_preorder(sample.species_tree_edge_index)

    feat = build_pairwise_feat(sample)
    E = feat.shape[0]
    assert feat.shape == (E, E, 4)

    # Recompute the two halves independently and compare.
    n_species = len(species_list)
    C = species_coclustering_matrix(
        sample.gene_tree_edge_indices, sample.gene_tree_species_ids,
        sample.gene_tree_leaf_masks, n_species,
    )
    sp_to_idx = {sp: i for i, sp in enumerate(species_list)}
    node_species = [sp_to_idx.get(name, -1) for name in sample.species_tree_node_names]
    clades = edge_clades_species(sample._edge_index_pre, sample.species_tree_is_leaf, node_species)

    expected_old = pairwise_partner_features(C, clades)
    expected_new = copy_aware_cluster_support(
        sample.gene_tree_edge_indices, sample.gene_tree_species_ids,
        sample.gene_tree_leaf_masks, clades, n_species,
    )
    assert torch.allclose(feat[..., :2], expected_old)
    assert torch.allclose(feat[..., 2:], expected_new)
