"""Smoke test for ablations B (duplication-conditioned co-clustering) and
D (n_eff node feature). Fully synthetic -- no packaged data needed.

Builds a tiny case where polyploid X is duplicated with one copy sistering
partner P and the other sistering C, plus single-copy trees that should DILUTE
the unconditioned matrix but NOT the duplication-conditioned one. Verifies:
  - unconditioned C is symmetric; conditioned C is directional and sharper for X;
  - single-copy dilution is removed by conditioning;
  - n_eff is finite and ~2 for a two-way-split polyploid;
  - build_pairwise_feat runs in both modes and augment_node_features_neff widens
    node features 13 -> 14.

Run in the gene2net env:  python scripts/smoke_new_features.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import torch
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.data.features import species_coclustering_matrix, n_eff_per_species
from gene2net_gnn.training.trainer_reconstruct import (
    build_pairwise_feat, augment_node_features_neff, set_feature_opts,
)

species_list = ["A", "B", "C", "X", "P"]
idx = {s: i for i, s in enumerate(species_list)}
sp_tree = Tree("(((A,B),C),(X,P));", format=1)

# Gene trees:
#  gt1: X duplicated -> one copy sisters P, one copy sisters C (the two subgenomes)
#  gt2,gt3: X SINGLE copy sistering P -> these dilute the UNconditioned row of X
gene_trees = [
    Tree("((A,B),(C,((X,P),(X,C))));", format=1),
    Tree("(((A,B),C),(X,P));", format=1),
    Tree("(((A,B),C),(X,P));", format=1),
]

sample = Gene2NetSample.from_trees(sp_tree, gene_trees, species_list)
sample._edge_index_pre = reorder_edge_index_preorder(sample.species_tree_edge_index)
n = len(species_list)

C0 = species_coclustering_matrix(sample.gene_tree_edge_indices, sample.gene_tree_species_ids,
                                 sample.gene_tree_leaf_masks, n, condition_on_dup=False)
C1 = species_coclustering_matrix(sample.gene_tree_edge_indices, sample.gene_tree_species_ids,
                                 sample.gene_tree_leaf_masks, n, condition_on_dup=True)

xi, pi, ci = idx["X"], idx["P"], idx["C"]
print("Unconditioned C (should be symmetric):")
print(C0.round(decimals=3))
print(f"  symmetric? {torch.allclose(C0, C0.t(), atol=1e-6)}")
print(f"  C0[X,P]={C0[xi,pi]:.3f}  C0[X,C]={C0[xi,ci]:.3f}  (P diluted up by single-copy trees)")

print("\nDuplication-conditioned C (X row credited only in the 1 tree where X is duplicated):")
print(C1.round(decimals=3))
print(f"  C1[X,P]={C1[xi,pi]:.3f}  C1[X,C]={C1[xi,ci]:.3f}  (both ~1: the two subgenome peaks)")
print(f"  C1[P,X]={C1[pi,xi]:.3f}  (0: P never duplicated -> empty row, still a valid column)")

neff = n_eff_per_species(C0)
print(f"\nn_eff per species: {neff.round(decimals=3).tolist()}")
print(f"  n_eff[X]={neff[xi]:.3f}  (>1 -> X splits over more than one partner)")
print(f"  all finite? {torch.isfinite(neff).all().item()}")

# Integration: pairwise feature in both modes, and node augmentation.
set_feature_opts(coclust_condition_on_dup=False, use_n_eff=False)
sample._pairwise_feat = None
pw0 = build_pairwise_feat(sample)
set_feature_opts(coclust_condition_on_dup=True, use_n_eff=True)
pw1 = build_pairwise_feat(sample)
print(f"\npairwise feat shape {tuple(pw0.shape)} | finite base={torch.isfinite(pw0).all().item()} "
      f"cond={torch.isfinite(pw1).all().item()} | differ={not torch.allclose(pw0, pw1)}")

before = tuple(sample.species_tree_node_features.shape)
augment_node_features_neff(sample)
after = tuple(sample.species_tree_node_features.shape)
augment_node_features_neff(sample)  # idempotent
after2 = tuple(sample.species_tree_node_features.shape)
print(f"node features {before} -> {after} (idempotent: {after == after2})")
print(f"  widened by exactly 1 col? {after[1] == before[1] + 1}")

ok = (torch.allclose(C0, C0.t(), atol=1e-6)              # unconditioned is symmetric
      and float(C1[xi, ci]) > float(C0[xi, ci]) + 0.1    # conditioning lifts the diluted 2nd peak
      and float(C1[pi, xi]) == 0.0                       # never-duplicated species: empty row
      and torch.isfinite(neff).all().item() and float(neff[xi]) > 1.0
      and after[1] == before[1] + 1 and after == after2  # n_eff widens by exactly 1, idempotent
      and not torch.allclose(pw0, pw1))                  # the B toggle actually changes the feature
print("\nSMOKE TEST:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
