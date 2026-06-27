"""Equivalence guards for the vectorized co-clustering / pairwise functions.

These functions run at load time on every sample, so they were vectorized for
speed. Each test pins the optimized function against a deliberately naive
reference so the rewrite is provably behavior-preserving.
"""
import random

import torch
from ete3 import Tree

from gene2net_gnn.data.tree_io import tree_to_edge_index
from gene2net_gnn.data.features import (
    species_coclustering_matrix,
    pairwise_partner_features,
)


def _random_numeric_gene_tree(n_species, n_leaves, rng):
    t = Tree()
    t.populate(n_leaves)
    for leaf in t.get_leaves():
        leaf.name = str(rng.randint(0, n_species - 1))
    return t


def _numeric_to_tensors(tree):
    ei, _ = tree_to_edge_index(tree)
    sp, lm = [], []
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            sp.append(int(node.name)); lm.append(True)
        else:
            sp.append(-1); lm.append(False)
    return ei, torch.tensor(sp, dtype=torch.long), torch.tensor(lm, dtype=torch.bool)


def _reference_coclust(eis, sps, lms, n_species):
    C = torch.zeros(n_species, n_species)
    nt = len(eis)
    if nt == 0:
        return C
    for ei, sp, lm in zip(eis, sps, lms):
        parent_leaf_species = {}
        for k in range(0, ei.shape[1], 2):
            p = int(ei[0, k]); c = int(ei[1, k])
            if bool(lm[c]) and int(sp[c]) >= 0:
                parent_leaf_species.setdefault(p, set()).add(int(sp[c]))
        seen = set()
        for specs in parent_leaf_species.values():
            specs = sorted(specs)
            for a_i in range(len(specs)):
                for b_i in range(a_i + 1, len(specs)):
                    seen.add((specs[a_i], specs[b_i]))
        for a, b in seen:
            C[a, b] += 1.0
            C[b, a] += 1.0
    return C / nt


def _reference_pairwise(coclust, edge_clades):
    E = len(edge_clades)
    feat = torch.zeros(E, E, 2)
    idx_lists = [torch.tensor(sorted(c), dtype=torch.long) for c in edge_clades]
    for i in range(E):
        ci = idx_lists[i]
        if ci.numel() == 0:
            continue
        for j in range(E):
            cj = idx_lists[j]
            if cj.numel() == 0:
                continue
            sub = coclust[ci][:, cj]
            feat[i, j, 0] = sub.mean()
            feat[i, j, 1] = sub.max()
    return feat


def test_coclustering_matches_reference_randomized():
    rng = random.Random(7)
    for _ in range(8):
        n_species = rng.randint(3, 9)
        eis, sps, lms = [], [], []
        for _ in range(rng.randint(2, 10)):
            ei, sp, lm = _numeric_to_tensors(
                _random_numeric_gene_tree(n_species, rng.randint(n_species, 2 * n_species), rng)
            )
            eis.append(ei); sps.append(sp); lms.append(lm)
        fast = species_coclustering_matrix(eis, sps, lms, n_species)
        ref = _reference_coclust(eis, sps, lms, n_species)
        assert torch.allclose(fast, ref, atol=1e-6)


def test_pairwise_partner_matches_reference_randomized():
    rng = random.Random(11)
    for _ in range(8):
        n_species = rng.randint(3, 10)
        m = torch.rand(n_species, n_species)
        coclust = (m + m.t()) / 2  # symmetric, like a real co-clustering matrix
        E = rng.randint(2, 7)
        edge_clades = [set(rng.sample(range(n_species), rng.randint(1, n_species))) for _ in range(E)]
        fast = pairwise_partner_features(coclust, edge_clades)
        ref = _reference_pairwise(coclust, edge_clades)
        assert torch.allclose(fast, ref, atol=1e-6)
