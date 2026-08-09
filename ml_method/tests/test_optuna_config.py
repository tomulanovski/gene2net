import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from optuna.trial import FixedTrial
from optuna_reconstruct import build_config_from_trial


def test_build_config_maps_all_params_and_preserves_base():
    base_mc = {"node_feat_dim": 13, "edge_feat_dim": 9, "partner_pair_feat_dim": 4,
               "n_parents": 1, "hidden_dim": 64, "n_gat_layers": 3, "dropout": 0.2,
               "conv_type": "gat"}
    base_tc = {"lr": 0.001, "weight_decay": 0.0001, "max_epochs": 200, "val_split": 0.2}
    t = FixedTrial({"conv_type": "gin", "hidden_dim": 128, "n_gat_layers": 2,
                    "dropout": 0.3, "lr": 0.0005, "weight_decay": 1e-5})
    mc, tc = build_config_from_trial(t, base_mc, base_tc)
    assert mc["conv_type"] == "gin" and mc["hidden_dim"] == 128
    assert mc["n_gat_layers"] == 2 and mc["dropout"] == 0.3
    assert tc["lr"] == 0.0005 and tc["weight_decay"] == 1e-5
    # method invariants preserved from the base config
    assert mc["n_parents"] == 1 and mc["partner_pair_feat_dim"] == 4
    assert tc["max_epochs"] == 200
