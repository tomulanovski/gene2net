import pytest
from ete3 import Tree
from gene2net_gnn.data.features import (
    compute_copy_count_features,
    compute_clustering_profile,
    compute_species_tree_edge_features,
)

def test_copy_count_features(sample_gene_trees):
    features = compute_copy_count_features(sample_gene_trees, species="D")
    assert "mean_copies" in features
    assert "var_copies" in features
    # 7 trees with 2 copies + 3 trees with 1 copy = mean 1.7
    assert features["mean_copies"] == pytest.approx(1.7, abs=0.1)

def test_copy_count_features_diploid(sample_gene_trees):
    features = compute_copy_count_features(sample_gene_trees, species="A")
    assert features["mean_copies"] == pytest.approx(1.0, abs=0.01)
    assert features["var_copies"] == pytest.approx(0.0, abs=0.01)

def test_clustering_profile(sample_gene_trees):
    profile = compute_clustering_profile(
        sample_gene_trees, species="D", all_species={"A","B","C","D","E"}
    )
    # D clusters with other species
    assert "E" in profile
    assert len(profile) == 4  # all species except D itself

def test_species_tree_edge_features(simple_species_tree, sample_gene_trees):
    edge_features = compute_species_tree_edge_features(simple_species_tree, sample_gene_trees)
    assert len(edge_features) > 0

def test_copy_count_features_absent_species(sample_gene_trees):
    """Species not in any tree should have 0 copies."""
    features = compute_copy_count_features(sample_gene_trees, species="Z")
    assert features["mean_copies"] == 0.0
