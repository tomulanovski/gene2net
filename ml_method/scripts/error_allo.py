"""Explore allo partner errors.

When the model predicts the WRONG second parent for an allopolyploid, is it predicting
the OTHER true parent (the 'home' = the polyploid's sibling clade in the backbone), or a
random lineage? If wrong predictions are mostly the home, then 0.44 is the symmetric-
parents problem (the model finds both parents but can't tell which is the labelled
"partner") -> a two-parent phaser fixes it. If they're random, it's a different problem.

Sanity: the reported 'correct' rate should land near the known allo accuracy (~0.44); if
it does, the edge<->clade alignment is right.

Run in final_project.
Usage:
  python scripts/error_allo.py --model-dir output/reconstruct_cladelabels_rooted \
      --model-config configs/reconstruct_base.yaml --data-root data/mul_trees_2k \
      --configs ils_low dup_loss_high_ne1M --max-samples 200
"""
import argparse
import os
import sys

import torch
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2
from gene2net_gnn.training.trainer_phase1 import prepare_sample
from gene2net_gnn.training.trainer_reconstruct import build_pairwise_feat, build_partner_targets
from gene2net_gnn.data.metadata_labels import sample_edge_bipartitions


def load_model(model_dir, model_config, device):
    with open(model_config) as f:
        cfg = yaml.safe_load(f)
    mc = cfg.get("model", {})
    model = SpeciesTreeGNNv2(
        node_feat_dim=int(mc.get("node_feat_dim", 13)),
        edge_feat_dim=int(mc.get("edge_feat_dim", 9)),
        hidden_dim=int(mc.get("hidden_dim", 64)),
        n_gat_layers=int(mc.get("n_gat_layers", 3)),
        n_gat_heads=int(mc.get("n_gat_heads", 4)),
        dropout=float(mc.get("dropout", 0.2)),
        partner_pair_feat_dim=int(mc.get("partner_pair_feat_dim", 2)),
    )
    ckpt = next((os.path.join(model_dir, n) for n in ("best_model.pt", "model.pt")
                 if os.path.exists(os.path.join(model_dir, n))), None)
    if ckpt is None:
        raise FileNotFoundError(f"No checkpoint in {model_dir}")
    model.load_state_dict(torch.load(ckpt, map_location="cpu", weights_only=True))
    model.to(device).eval()
    return model


def edge_clade_map(sample):
    sd = {
        "species_tree_edge_index": sample.species_tree_edge_index,
        "species_list": sample.species_list,
        "species_tree_node_names": sample.species_tree_node_names,
        "species_tree_is_leaf": sample.species_tree_is_leaf,
    }
    return {i: c for i, c in sample_edge_bipartitions(sd)}


def home_clade(x_clade, clades):
    """The polyploid's sibling clade = smallest clade strictly containing x, minus x."""
    supersets = [c for c in clades if x_clade < c]
    if not supersets:
        return None
    parent = min(supersets, key=len)
    return parent - x_clade


def analyze_config(model, ds, device, max_samples):
    n_allo = n_correct = n_home = n_other = n_home_but_correct = 0
    n_correct_masked = 0   # partner accuracy when the home edge is excluded from candidates
    n_partner_is_home = 0  # how often the LABELLED partner == the home (the labelling bug)
    for i in range(min(len(ds), max_samples)):
        try:
            sample = ds[i]
        except Exception:
            continue
        if sample.labels is None:
            continue
        with torch.no_grad():
            prepared = prepare_sample(sample, device)
            if prepared is None:
                continue
            _, wgd_targets, _ = prepared
            n_edges = wgd_targets.shape[0]
            inputs, _, _ = prepared
            wgd_logits, edge_emb = model(**inputs)
            partner_targets = build_partner_targets(sample, n_edges, device)
            query_idx = (partner_targets >= 0).nonzero(as_tuple=True)[0]
            if query_idx.numel() == 0:
                continue
            pairwise_feat = build_pairwise_feat(sample).to(device)
            scores = model.compute_partner_scores_rows(edge_emb, query_idx, pairwise_feat)
            pred = scores.argmax(dim=-1)
            tgt = partner_targets[query_idx]

        cmap = edge_clade_map(sample)
        clade_values = list(cmap.values())
        for k in range(query_idx.numel()):
            q, t, p = int(query_idx[k]), int(tgt[k]), int(pred[k])
            if t == q:
                continue  # auto
            n_allo += 1
            xc = cmap.get(q)
            h = home_clade(xc, clade_values) if xc is not None else None
            tc = cmap.get(t)
            if h is not None and tc is not None and bool(tc & h):
                n_partner_is_home += 1   # ASTRAL placed X next to its labelled partner
            pc = cmap.get(p)
            pred_is_home = h is not None and pc is not None and bool(pc & h)
            if p == t:
                n_correct += 1
                if pred_is_home:
                    n_home_but_correct += 1   # home actually was the partner
            elif pred_is_home:
                n_home += 1
            else:
                n_other += 1

            # Smart tweak: re-predict the partner with only the HOME edge excluded.
            # Keep self (edge q) as a candidate — self-partner is how AUTOpolyploidy is
            # predicted, so masking it would break auto. We only remove the home (sibling).
            home_idx = next((i for i, c in cmap.items() if h is not None and c == h), None)
            row = scores[k].clone()
            if home_idx is not None and home_idx != q:
                row[home_idx] = float("-inf")
            if int(row.argmax()) == t:
                n_correct_masked += 1
    return (n_allo, n_correct, n_home, n_other, n_home_but_correct,
            n_correct_masked, n_partner_is_home)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--model-config", required=True)
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--subdir", default="training_rooted")
    ap.add_argument("--configs", nargs="+", required=True)
    ap.add_argument("--max-samples", type=int, default=200)
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = load_model(args.model_dir, args.model_config, device)

    for cfg in args.configs:
        d = os.path.join(args.data_root, args.subdir, cfg)
        if not os.path.isdir(d):
            print(f"{cfg}: missing {d}")
            continue
        ds = Gene2NetDataset(d, clade_labels=True)
        (n_allo, n_ok, n_home, n_other, n_hbc,
         n_ok_masked, n_partner_home) = analyze_config(model, ds, device, args.max_samples)
        print(f"\n=== {cfg} — {n_allo} allo events ===")
        if not n_allo:
            continue
        wrong = n_home + n_other
        print(f"  correct (partner right):        {n_ok}/{n_allo} ({100*n_ok/n_allo:.1f}%)   <- should be ~44% (sanity)")
        print(f"  WRONG -> predicted the HOME:     {n_home}/{n_allo} ({100*n_home/n_allo:.1f}%)")
        print(f"  WRONG -> predicted OTHER:        {n_other}/{n_allo} ({100*n_other/n_allo:.1f}%)")
        if wrong:
            print(f"  of the wrong ones, home-collision: {100*n_home/wrong:.1f}%")
        print(f"\n  *** LABELLING BUG CHECK ***")
        print(f"  labelled partner == the home:    {n_partner_home}/{n_allo} ({100*n_partner_home/n_allo:.1f}%)")
        print(f"    (if ~30-40%, the target is the home a third of the time -> relabel partner = away parent)")
        print(f"\n  SMART TWEAK (exclude the home; hurts because partner==home so often):")
        print(f"  correct with home masked:        {n_ok_masked}/{n_allo} ({100*n_ok_masked/n_allo:.1f}%)"
              f"   (was {100*n_ok/n_allo:.1f}%)")


if __name__ == "__main__":
    main()
