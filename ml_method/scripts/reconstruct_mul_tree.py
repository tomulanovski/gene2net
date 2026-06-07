"""Build inferred MUL-trees from the reconstruction model's predictions.

For each packaged sample: run the model, take WGD edges (prob >= threshold) and
each one's partner edge (argmax of the partner head; self => auto, other => allo),
map edges to clades, assemble a MUL-tree with build_mul_tree, and write it as
Newick (output.tre) next to a copy of the ground-truth MUL-tree. The output
directory is laid out so the existing compare_nets.py / one_stop_compare can
score it (edit distance / Jaccard) against the ground truth.

Edge i of the model output corresponds to the i-th non-root node in preorder
(the alignment fixed in prepare_sample / reorder_edge_index_preorder), so clades
are enumerated the same way here.

Usage:
    python scripts/reconstruct_mul_tree.py \
        --model-dir output/reconstruct_aligned \
        --mul-trees-dir data/mul_trees_2k --config ils_low \
        --start 0 --n 50 --threshold 0.85 --out-dir output/reconstruct_aligned/mul_trees/ils_low
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2, propagate_to_internal
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent


def load_astral_tree(path):
    with open(path) as f:
        return Tree(f.read().strip(), format=1)


def preorder_edge_clades(tree):
    """edge i -> frozenset of leaf names below the i-th non-root preorder node."""
    clades = []
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        clades.append(frozenset(node.get_leaf_names()))
    return clades


def model_inputs_for(sample, device):
    """Build model inputs with the corrected (preorder) edge ordering."""
    ei = reorder_edge_index_preorder(sample.species_tree_edge_index)
    prop = propagate_to_internal(
        sample.species_tree_node_features, ei,
        sample.species_tree_is_leaf, sample.species_tree_node_features.shape[1],
    )
    return {
        "node_features": sample.species_tree_node_features.to(device),
        "edge_index": ei.to(device),
        "edge_features": sample.species_tree_edge_features.to(device),
        "is_leaf": sample.species_tree_is_leaf.to(device),
        "node_features_propagated": prop.to(device),
    }


def load_model(model_dir, model_config, device):
    model = SpeciesTreeGNNv2(
        node_feat_dim=int(model_config.get("node_feat_dim", 13)),
        edge_feat_dim=int(model_config.get("edge_feat_dim", 9)),
        hidden_dim=int(model_config.get("hidden_dim", 64)),
        n_gat_layers=int(model_config.get("n_gat_layers", 3)),
        n_gat_heads=int(model_config.get("n_gat_heads", 4)),
        dropout=float(model_config.get("dropout", 0.2)),
    )
    for name in ["best_model.pt", "best_partner_model.pt"]:
        p = os.path.join(model_dir, name)
        if os.path.exists(p):
            model.load_state_dict(torch.load(p, map_location=device, weights_only=True))
            print(f"Loaded {name}")
            break
    return model.to(device).eval()


def reconstruct_one(model, sample, astral_tree, threshold, device):
    """Return (mul_tree, n_auto, n_allo) for one sample."""
    clades = preorder_edge_clades(astral_tree)
    inputs = model_inputs_for(sample, device)

    with torch.no_grad():
        wgd_logits, edge_emb = model(**inputs)
        wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
        partner_scores = model.compute_partner_scores(edge_emb)

    n_edges = min(len(clades), wgd_prob.shape[0])
    events = []
    n_auto = n_allo = 0
    for i in range(n_edges):
        if wgd_prob[i].item() < threshold:
            continue
        j = int(partner_scores[i, :n_edges].argmax())
        if j == i:
            n_auto += 1
        else:
            n_allo += 1
        events.append(WGDEvent(
            wgd_edge_clade=clades[i],
            partner_edge_clade=clades[j],
            confidence=wgd_prob[i].item(),
        ))

    mul_tree = build_mul_tree(astral_tree, events)
    return mul_tree, n_auto, n_allo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True, help="config name, e.g. ils_low")
    parser.add_argument("--model-config", default=None, help="YAML (default: configs/reconstruct.yaml)")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=0.85)
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    with open(cfg_path) as f:
        model_config = yaml.safe_load(f).get("model", {})
    device = "cuda" if torch.cuda.is_available() else "cpu"

    out_dir = args.out_dir or os.path.join(args.model_dir, "mul_trees", args.config)
    os.makedirs(out_dir, exist_ok=True)

    model = load_model(args.model_dir, model_config, device)

    done = skipped = 0
    for idx in range(args.start, args.start + args.n):
        idx_str = f"{idx:04d}"
        sample_dir = os.path.join(args.mul_trees_dir, "training", args.config, f"sample_{idx_str}")
        gt_path = os.path.join(args.mul_trees_dir, f"mul_tree_{idx_str}.nex")
        astral_path = os.path.join(
            args.mul_trees_dir, "simphy", args.config, idx_str,
            f"replicate_{args.replicate}", "astral_species.tre",
        )
        if not (os.path.isdir(sample_dir) and os.path.exists(astral_path)):
            skipped += 1
            continue

        sample = Gene2NetSample.load(sample_dir)
        if sample.species_tree_edge_features is None or sample.species_tree_edge_features.shape[1] != int(model_config.get("edge_feat_dim", 9)):
            skipped += 1
            continue

        astral_tree = load_astral_tree(astral_path)
        mul_tree, n_auto, n_allo = reconstruct_one(model, sample, astral_tree, args.threshold, device)

        case_dir = os.path.join(out_dir, f"sample_{idx_str}")
        os.makedirs(case_dir, exist_ok=True)
        mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
        # Copy ground truth alongside for the comparison tool.
        if os.path.exists(gt_path):
            with open(gt_path) as f:
                gt = f.read()
            with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                f.write(gt)

        print(f"[{idx_str}] events={n_auto+n_allo} (auto={n_auto}, allo={n_allo}) -> {case_dir}/output.tre")
        done += 1

    print(f"\nDone: {done} reconstructed, {skipped} skipped. Output under {out_dir}")


if __name__ == "__main__":
    main()
