"""Generate the ablation configs from the final one-parter config by varying one
hyperparameter at a time. Base (n_gat_layers=3, hidden_dim=64) is the shipped model
reconstruct_op_away_clade, so it is not regenerated. Run once:
    python scripts/gen_ablation_configs.py
"""
import os
import yaml

BASE = os.path.join(os.path.dirname(__file__), "..", "configs", "reconstruct_oneparter.yaml")
CFG_DIR = os.path.join(os.path.dirname(__file__), "..", "configs")

with open(BASE) as f:
    base = yaml.safe_load(f)

variants = {
    # depth ablation (base is 3)
    "abl_depth1": ("n_gat_layers", 1),
    "abl_depth2": ("n_gat_layers", 2),
    "abl_depth4": ("n_gat_layers", 4),
    # width ablation (base is 64)
    "abl_hidden32": ("hidden_dim", 32),
    "abl_hidden128": ("hidden_dim", 128),
    # learning-rate ablation (base is 0.001)
    "abl_lr3e-4": ("lr", 0.0003),
}

for name, (key, val) in variants.items():
    cfg = {"model": dict(base["model"]), "training": dict(base["training"])}
    if key in cfg["model"]:
        cfg["model"][key] = val
    else:
        cfg["training"][key] = val
    out = os.path.join(CFG_DIR, f"{name}.yaml")
    with open(out, "w") as f:
        f.write(f"# Ablation: {key} = {val} (else identical to reconstruct_oneparter.yaml)\n")
        yaml.safe_dump(cfg, f, sort_keys=False)
    print("wrote", out)
