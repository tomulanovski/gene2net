import pytest
import torch
import os
import pickle
from ete3 import Tree
from gene2net_gnn.data.dataset import Gene2NetSample, Gene2NetDataset
from gene2net_gnn.data.collate import gene2net_collate

def test_sample_from_trees():
    """Create a sample from raw trees."""
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gene_trees = [
        Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
        for _ in range(5)
    ]
    species_list = sorted(["A", "B", "C", "D", "E"])

    sample = Gene2NetSample.from_trees(species_tree, gene_trees, species_list)
    assert sample.species_tree_edge_index is not None
    assert sample.species_tree_edge_index.shape[0] == 2
    assert len(sample.gene_tree_edge_indices) == 5
    assert sample.species_tree_node_features is not None

def test_sample_has_labels():
    """Sample with labels from a MUL-tree."""
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    mul_tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gene_trees = [
        Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
        for _ in range(5)
    ]
    species_list = sorted(["A", "B", "C", "D", "E"])

    sample = Gene2NetSample.from_trees(
        species_tree, gene_trees, species_list, mul_tree=mul_tree
    )
    assert sample.labels is not None
    assert sample.labels.n_edges > 0

def test_sample_save_load(tmp_path):
    """Save and load a sample."""
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gene_trees = [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)] * 3
    species_list = sorted(["A", "B", "C", "D", "E"])

    sample = Gene2NetSample.from_trees(species_tree, gene_trees, species_list)
    sample.save(str(tmp_path / "example_0001"))

    loaded = Gene2NetSample.load(str(tmp_path / "example_0001"))
    assert loaded.species_tree_edge_index.shape == sample.species_tree_edge_index.shape

def test_collate_batch():
    """Collate multiple samples into a batch."""
    species_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    gene_trees = [Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)] * 3
    species_list = sorted(["A", "B", "C", "D", "E"])

    samples = [
        Gene2NetSample.from_trees(species_tree, gene_trees, species_list)
        for _ in range(2)
    ]
    batch = gene2net_collate(samples)
    assert "species_trees" in batch
    assert "gene_trees" in batch
    assert len(batch["species_trees"]) == 2
