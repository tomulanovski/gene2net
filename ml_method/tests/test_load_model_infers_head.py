import torch
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from scripts.reconstruct_mul_tree import load_model


def test_load_model_infers_head_shape_from_checkpoint(tmp_path):
    # Save a two-parent, pair_dim=4 model...
    m = SpeciesTreeGNNv2(node_feat_dim=13, edge_feat_dim=9, hidden_dim=16,
                         n_gat_layers=2, n_gat_heads=2, partner_pair_feat_dim=4, n_parents=2)
    d = tmp_path / "model"
    d.mkdir()
    torch.save(m.state_dict(), d / "best_model.pt")

    # ...and load it with a config that LIES (one-partner, pair_dim=2). load_model
    # must infer the true shapes from the checkpoint and load without error.
    cfg = {"node_feat_dim": 13, "edge_feat_dim": 9, "hidden_dim": 16,
           "n_gat_layers": 2, "n_gat_heads": 2, "n_parents": 1, "partner_pair_feat_dim": 2}
    loaded = load_model(str(d), cfg, "cpu")
    assert loaded.n_parents == 2
    assert loaded.partner_pair_feat_dim == 4
