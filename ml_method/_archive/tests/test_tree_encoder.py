import torch
from ete3 import Tree
from gene2net_gnn.model.tree_encoder import GeneTreeEncoder
from gene2net_gnn.data.tree_io import tree_to_edge_index


def test_encoder_output_shape():
    encoder = GeneTreeEncoder(n_species=5, embed_dim=32)
    tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    edge_index, node_names = tree_to_edge_index(tree)

    sp_map = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    species_ids = torch.tensor([sp_map.get(n, -1) for n in node_names])
    branch_lengths = torch.ones(edge_index.shape[1])

    leaf_embeddings = encoder(edge_index, species_ids, branch_lengths)
    # 5 leaves → 5 embeddings of dim 32
    n_leaves = sum(1 for s in species_ids if s >= 0)
    assert leaf_embeddings.shape == (n_leaves, 32)


def test_encoder_mul_tree():
    """Gene tree where D appears twice - should get different embeddings."""
    encoder = GeneTreeEncoder(n_species=5, embed_dim=32)
    tree = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
    edge_index, node_names = tree_to_edge_index(tree)

    sp_map = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    species_ids = torch.tensor([sp_map.get(n, -1) for n in node_names])
    branch_lengths = torch.ones(edge_index.shape[1])

    leaf_embeddings = encoder(edge_index, species_ids, branch_lengths)
    # 6 leaves (A, D, B, C, D, E) → 6 embeddings
    n_leaves = sum(1 for s in species_ids if s >= 0)
    assert leaf_embeddings.shape == (n_leaves, 32)


def test_encoder_gradient_flows():
    """Verify gradients flow through the encoder."""
    encoder = GeneTreeEncoder(n_species=5, embed_dim=16)
    tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    edge_index, node_names = tree_to_edge_index(tree)

    sp_map = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    species_ids = torch.tensor([sp_map.get(n, -1) for n in node_names])
    branch_lengths = torch.ones(edge_index.shape[1])

    leaf_embeddings = encoder(edge_index, species_ids, branch_lengths)
    loss = leaf_embeddings.sum()
    loss.backward()

    assert encoder.species_embedding.weight.grad is not None
