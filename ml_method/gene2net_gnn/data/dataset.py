"""PyG-compatible dataset for Gene2Net training examples."""
import os
import pickle
from dataclasses import dataclass, field
from typing import List, Optional

import torch
from ete3 import Tree

from gene2net_gnn.data.tree_io import tree_to_edge_index, get_species_set
from gene2net_gnn.data.features import (
    compute_copy_count_features,
    compute_clustering_summary,
    compute_species_tree_edge_features,
)
from gene2net_gnn.data.label_extractor import (
    TrainingLabels,
    decompose_mul_tree,
    map_events_to_astral,
)


@dataclass
class Gene2NetSample:
    """One training example: gene trees + species tree + optional labels."""
    # Species tree (PyG format)
    species_tree_edge_index: torch.Tensor = None        # [2, E]
    species_tree_node_features: torch.Tensor = None     # [N, F_node]
    species_tree_edge_features: torch.Tensor = None     # [E//2, F_edge]
    species_tree_node_names: List[str] = field(default_factory=list)
    species_tree_is_leaf: torch.Tensor = None           # [N] bool
    n_species: int = 0
    species_list: List[str] = field(default_factory=list)

    # Gene trees (lists of tensors - variable count)
    gene_tree_edge_indices: List[torch.Tensor] = field(default_factory=list)
    gene_tree_species_ids: List[torch.Tensor] = field(default_factory=list)
    gene_tree_branch_lengths: List[torch.Tensor] = field(default_factory=list)
    gene_tree_leaf_masks: List[torch.Tensor] = field(default_factory=list)

    # Labels (optional - only for training)
    labels: Optional[TrainingLabels] = None

    @classmethod
    def from_trees(cls, species_tree, gene_trees, species_list, mul_tree=None):
        """Build sample from ETE3 trees."""
        sample = cls()
        sample.species_list = species_list
        sample.n_species = len(species_list)
        sp_to_idx = {sp: i for i, sp in enumerate(species_list)}

        # Convert species tree
        edge_index, node_names = tree_to_edge_index(species_tree)
        sample.species_tree_edge_index = edge_index
        sample.species_tree_node_names = node_names

        # Leaf mask
        is_leaf = torch.zeros(len(node_names), dtype=torch.bool)
        for i, name in enumerate(node_names):
            if name in species_list:
                is_leaf[i] = True
        sample.species_tree_is_leaf = is_leaf

        # Node features: 8 copy count features + 5 clustering summary = 13 dims (fixed)
        node_features = []
        for name in node_names:
            if name in species_list:
                cc = compute_copy_count_features(gene_trees, name)
                cs = compute_clustering_summary(gene_trees, name, set(species_list))
                feat = [cc["mean_copies"], cc["var_copies"], cc["mode_copies"],
                        cc["p_absent"], cc["p_1_copy"], cc["p_2_copies"], cc["p_3plus_copies"],
                        cc["max_copies"]]
                feat.extend(cs)  # 5 summary stats
                node_features.append(feat)
            else:
                # Internal node: zeros (will be filled by GNN)
                node_features.append([0.0] * 13)
        sample.species_tree_node_features = torch.tensor(node_features, dtype=torch.float)

        # Edge features
        edge_feats_dict = compute_species_tree_edge_features(species_tree, gene_trees)
        edge_feat_list = []
        # Process edges in pairs (we stored undirected, take every other)
        n_undirected = edge_index.shape[1] // 2
        for i in range(n_undirected):
            ef = edge_feats_dict.get(i, {})
            edge_feat_list.append([
                ef.get("concordance_factor", 0.0),
                ef.get("branch_length", 0.0),
                ef.get("clade_size", 0.0),
                ef.get("depth", 0.0),
            ])
        sample.species_tree_edge_features = torch.tensor(edge_feat_list, dtype=torch.float)

        # Convert gene trees
        for gt in gene_trees:
            ei, names = tree_to_edge_index(gt)
            sample.gene_tree_edge_indices.append(ei)

            species_ids = []
            leaf_mask = []
            for name in names:
                if name in sp_to_idx:
                    species_ids.append(sp_to_idx[name])
                    leaf_mask.append(True)
                else:
                    species_ids.append(-1)
                    leaf_mask.append(False)
            sample.gene_tree_species_ids.append(torch.tensor(species_ids, dtype=torch.long))
            sample.gene_tree_leaf_masks.append(torch.tensor(leaf_mask, dtype=torch.bool))

            # Branch lengths: use actual lengths from tree, default 1.0 if missing
            bl = []
            node_idx_map = {}
            for idx_n, node in enumerate(gt.traverse("preorder")):
                node_idx_map[id(node)] = idx_n
            for node in gt.traverse("preorder"):
                for child in node.children:
                    dist = child.dist if child.dist and child.dist > 0 else 1.0
                    bl.append(dist)  # parent→child
                    bl.append(dist)  # child→parent (same length)
            sample.gene_tree_branch_lengths.append(torch.tensor(bl, dtype=torch.float))

        # Labels from MUL-tree if provided
        if mul_tree is not None:
            events = decompose_mul_tree(mul_tree)
            labels = map_events_to_astral(events, species_tree)
            sample.labels = labels

        return sample

    def save(self, directory):
        """Save sample to directory."""
        os.makedirs(directory, exist_ok=True)
        # Save everything as pickle
        data = {
            "species_tree_edge_index": self.species_tree_edge_index,
            "species_tree_node_features": self.species_tree_node_features,
            "species_tree_edge_features": self.species_tree_edge_features,
            "species_tree_node_names": self.species_tree_node_names,
            "species_tree_is_leaf": self.species_tree_is_leaf,
            "n_species": self.n_species,
            "species_list": self.species_list,
            "gene_tree_edge_indices": self.gene_tree_edge_indices,
            "gene_tree_species_ids": self.gene_tree_species_ids,
            "gene_tree_branch_lengths": self.gene_tree_branch_lengths,
            "gene_tree_leaf_masks": self.gene_tree_leaf_masks,
            "labels": self.labels,
        }
        with open(os.path.join(directory, "sample.pkl"), "wb") as f:
            pickle.dump(data, f)

    @classmethod
    def load(cls, directory):
        """Load sample from directory."""
        with open(os.path.join(directory, "sample.pkl"), "rb") as f:
            data = pickle.load(f)
        sample = cls()
        for key, val in data.items():
            setattr(sample, key, val)
        return sample


class Gene2NetDataset:
    """Dataset of Gene2Net training examples stored on disk."""

    def __init__(self, root_dir):
        self.root_dir = root_dir
        self.example_dirs = sorted([
            os.path.join(root_dir, d) for d in os.listdir(root_dir)
            if os.path.isdir(os.path.join(root_dir, d))
        ])

    def __len__(self):
        return len(self.example_dirs)

    def __getitem__(self, idx):
        return Gene2NetSample.load(self.example_dirs[idx])
