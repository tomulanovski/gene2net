import torch
from ete3 import Tree
from gene2net_gnn.model.gene2net_model import Gene2NetModel
from gene2net_gnn.data.tree_io import tree_to_edge_index

def test_full_forward_pass():
    model = Gene2NetModel(
        n_species=5,
        embed_dim=32,
        n_attention_heads=4,
        n_gat_layers=2,
        handcrafted_dim=12,
    )

    # Build species tree data
    sp_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    sp_edge_index, sp_node_names = tree_to_edge_index(sp_tree)
    sp_map = {"A":0, "B":1, "C":2, "D":3, "E":4}

    n_sp_nodes = len(sp_node_names)
    species_tree_data = {
        "edge_index": sp_edge_index,
        "node_features": torch.randn(n_sp_nodes, 12),  # hand-crafted
        "is_leaf": torch.tensor([sp_map.get(n, -1) >= 0 for n in sp_node_names], dtype=torch.bool),
        "species_ids": torch.tensor([sp_map.get(n, -1) for n in sp_node_names]),
    }

    # Build gene tree data (3 trees, D appears twice in each)
    gene_tree_data = []
    for _ in range(3):
        gt = Tree("((A:1,(D:1,B:1):1):1,(C:1,(D:1,E:1):1):1);", format=1)
        ei, names = tree_to_edge_index(gt)
        sids = torch.tensor([sp_map.get(n, -1) for n in names])
        gene_tree_data.append({
            "edge_index": ei,
            "species_ids": sids,
            "branch_lengths": torch.ones(ei.shape[1]),
            "leaf_mask": sids >= 0,
        })

    wgd_logits, partner_logits, edge_embeddings = model(gene_tree_data, species_tree_data)

    n_edges = sp_edge_index.shape[1] // 2
    assert wgd_logits.shape[0] == n_edges
    assert wgd_logits.shape[1] == 4  # 0,1,2,3 events
    assert edge_embeddings.shape == (n_edges, 32)

def test_forward_no_gene_trees():
    """Should handle zero gene trees (degenerate case)."""
    model = Gene2NetModel(n_species=5, embed_dim=16, n_attention_heads=2, n_gat_layers=1, handcrafted_dim=8)

    sp_tree = Tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);", format=1)
    sp_edge_index, sp_node_names = tree_to_edge_index(sp_tree)
    sp_map = {"A":0, "B":1, "C":2, "D":3, "E":4}

    species_tree_data = {
        "edge_index": sp_edge_index,
        "node_features": torch.randn(len(sp_node_names), 8),
        "is_leaf": torch.tensor([sp_map.get(n, -1) >= 0 for n in sp_node_names], dtype=torch.bool),
        "species_ids": torch.tensor([sp_map.get(n, -1) for n in sp_node_names]),
    }

    wgd_logits, partner_logits, edge_embeddings = model([], species_tree_data)
    assert wgd_logits.shape[0] == sp_edge_index.shape[1] // 2

def test_forward_gradient():
    model = Gene2NetModel(n_species=3, embed_dim=16, n_attention_heads=2, n_gat_layers=1, handcrafted_dim=4)

    sp_tree = Tree("(A:1,(B:1,C:1):1);", format=1)
    sp_edge_index, sp_node_names = tree_to_edge_index(sp_tree)
    sp_map = {"A":0, "B":1, "C":2}

    species_tree_data = {
        "edge_index": sp_edge_index,
        "node_features": torch.randn(len(sp_node_names), 4),
        "is_leaf": torch.tensor([sp_map.get(n, -1) >= 0 for n in sp_node_names], dtype=torch.bool),
        "species_ids": torch.tensor([sp_map.get(n, -1) for n in sp_node_names]),
    }

    gt = Tree("(A:1,(B:1,C:1):1);", format=1)
    ei, names = tree_to_edge_index(gt)
    sids = torch.tensor([sp_map.get(n, -1) for n in names])
    gene_tree_data = [{"edge_index": ei, "species_ids": sids, "branch_lengths": torch.ones(ei.shape[1]), "leaf_mask": sids >= 0}]

    wgd_logits, partner_logits, edge_embeddings = model(gene_tree_data, species_tree_data)
    loss = wgd_logits.sum() + partner_logits.sum()
    loss.backward()

    # Check gradients flow to all modules
    assert model.tree_encoder.species_embedding.weight.grad is not None
