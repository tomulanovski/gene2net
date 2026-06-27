"""Build inferred MUL-trees from the reconstruction model's predictions.

For each packaged sample: run the model once (the expensive part: forward +
co-clustering), then decode at one or more thresholds (cheap), assembling a
MUL-tree with build_mul_tree per threshold and writing it as Newick (output.tre)
next to the ground-truth MUL-tree, laid out for compare_nets / one_stop_compare.

Speed:
  - --workers parallelises across samples (CPU; the model is tiny and the
    bottleneck is the per-sample co-clustering, so CPU workers scale well).
  - --thresholds runs the model once and emits every threshold in one pass,
    so a sweep doesn't re-run the model per threshold.

Edge i of the model output is the i-th non-root node in preorder (alignment
fixed in reorder_edge_index_preorder), so clades are enumerated the same way.

Usage:
    python scripts/reconstruct_mul_tree.py \
        --model-dir output/reconstruct_allo \
        --mul-trees-dir data/mul_trees_2k --config ils_low \
        --start 0 --n 50 --thresholds 0.5,0.7,0.9,0.95 --workers 8
Output: <out-base>/t{T}/sample_NNNN/{output.tre,ground_truth.nex}
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.model.species_gnn_v2 import SpeciesTreeGNNv2, propagate_to_internal
from gene2net_gnn.training.trainer_reconstruct import build_pairwise_feat
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent

# Per-worker globals (set by init_worker).
_M = {}


def load_astral_tree(path):
    with open(path) as f:
        return Tree(f.read().strip(), format=1)


def load_nexus_tree(path):
    """Load the true backbone (species_tree_NNNN.nex) from a NEXUS file."""
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree found in {path}")


def preorder_edge_clades(tree):
    clades = []
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        clades.append(frozenset(node.get_leaf_names()))
    return clades


def reroot_to_match(astral_tree, true_tree):
    """Re-root the ASTRAL tree so its root split matches the true backbone's root.

    ASTRAL output is unrooted, so ete3 roots it arbitrarily. The MUL-tree edit
    distance is rooting-sensitive, so a correct unrooted topology with the wrong
    root still scores badly. We pick the smaller side of the true root as the
    outgroup and set it on the ASTRAL tree. Falls back to the original tree if the
    outgroup is not monophyletic in ASTRAL (i.e. unrooted topologies differ)."""
    kids = true_tree.children
    if len(kids) < 2:
        return astral_tree
    sides = [set(c.get_leaf_names()) for c in kids]
    outgroup = min(sides, key=len)
    try:
        og = list(outgroup)
        if len(og) == 1:
            astral_tree.set_outgroup(og[0])
        else:
            astral_tree.set_outgroup(astral_tree.get_common_ancestor(og))
    except Exception:
        pass  # outgroup not monophyletic in ASTRAL -> leave as-is
    return astral_tree


def map_clades_by_jaccard(src_clades, dst_clades):
    """Map each source (ASTRAL) clade to the best-Jaccard destination (true
    backbone) clade. Used to place predictions made on the ASTRAL topology onto
    the true backbone, so we can score 'true backbone + model events' and isolate
    the ASTRAL-backbone cost from the prediction cost."""
    mapped = []
    for c in src_clades:
        best, best_j = c, -1.0
        for d in dst_clades:
            inter = len(c & d)
            if inter == 0:
                continue
            j = inter / len(c | d)
            if j > best_j:
                best, best_j = d, j
        mapped.append(best)
    return mapped


def model_inputs_for(sample, device):
    ei = reorder_edge_index_preorder(sample.species_tree_edge_index)
    sample._edge_index_pre = ei  # used by build_pairwise_feat
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
        partner_pair_feat_dim=int(model_config.get("partner_pair_feat_dim", 2)),
    )
    loaded = False
    for name in ["best_model.pt", "best_partner_model.pt"]:
        p = os.path.join(model_dir, name)
        if os.path.exists(p):
            model.load_state_dict(torch.load(p, map_location=device, weights_only=True))
            loaded = True
            break
    if not loaded:
        raise FileNotFoundError(
            f"No checkpoint (best_model.pt or best_partner_model.pt) found in {model_dir!r}. "
            "Refusing to evaluate a randomly initialized model."
        )
    return model.to(device).eval()


def decode(model, build_tree, clades, wgd_prob, edge_emb, pairwise_feat, threshold,
           event_clades=None):
    """Decode predictions into a MUL-tree at one threshold. Cheap (rows-only).

    `clades` aligns with the model edges (ASTRAL topology) and is used for
    indexing/partner argmax. `event_clades` is what actually goes into the
    WGDEvent and must exist in `build_tree`; it defaults to `clades` (build on
    ASTRAL) but is the Jaccard-mapped true-backbone clades when building on the
    true backbone.
    """
    if event_clades is None:
        event_clades = clades
    n_edges = min(len(clades), wgd_prob.shape[0])
    wgd_edges = [i for i in range(n_edges) if wgd_prob[i].item() >= threshold]
    events = []
    n_auto = n_allo = 0
    if wgd_edges:
        query = torch.tensor(wgd_edges, dtype=torch.long)
        rows = model.compute_partner_scores_rows(edge_emb, query, pairwise_feat)  # [Q, E]
        for q, i in enumerate(wgd_edges):
            j = int(rows[q, :n_edges].argmax())
            if j == i:
                n_auto += 1
            else:
                n_allo += 1
            events.append(WGDEvent(
                wgd_edge_clade=event_clades[i],
                partner_edge_clade=event_clades[j],
                confidence=wgd_prob[i].item(),
            ))
    return build_mul_tree(build_tree, events), n_auto, n_allo


def init_worker(model_dir, model_config, mul_trees_dir, config, replicate,
                thresholds, out_base, edge_feat_dim, backbone):
    _M["model"] = load_model(model_dir, model_config, "cpu")
    _M["mul_trees_dir"] = mul_trees_dir
    _M["config"] = config
    _M["replicate"] = replicate
    _M["thresholds"] = thresholds
    _M["out_base"] = out_base
    _M["edge_feat_dim"] = edge_feat_dim
    _M["backbone"] = backbone


def process_index(idx):
    """Run the model once for sample idx and write a MUL-tree per threshold."""
    idx_str = f"{idx:04d}"
    mul_trees_dir = _M["mul_trees_dir"]; config = _M["config"]
    sample_dir = os.path.join(mul_trees_dir, "training", config, f"sample_{idx_str}")
    gt_path = os.path.join(mul_trees_dir, f"mul_tree_{idx_str}.nex")
    astral_path = os.path.join(
        mul_trees_dir, "simphy", config, idx_str,
        f"replicate_{_M['replicate']}", "astral_species.tre",
    )
    if not (os.path.isdir(sample_dir) and os.path.exists(astral_path)):
        return None

    sample = Gene2NetSample.load(sample_dir)
    ef = sample.species_tree_edge_features
    if ef is None or ef.shape[1] != _M["edge_feat_dim"]:
        return None

    astral_tree = load_astral_tree(astral_path)
    clades = preorder_edge_clades(astral_tree)
    model = _M["model"]

    # Choose the backbone to build on. 'astral' (default) is the real pipeline.
    # 'true' builds the model's predicted events on the true diploid backbone
    # (species_tree_NNNN.nex), mapping ASTRAL clades to true-backbone clades by
    # Jaccard, to isolate the ASTRAL-backbone cost from the prediction cost.
    if _M["backbone"] == "true":
        bb_path = os.path.join(mul_trees_dir, f"species_tree_{idx_str}.nex")
        if not os.path.exists(bb_path):
            return None
        build_tree = load_nexus_tree(bb_path)
        event_clades = map_clades_by_jaccard(clades, preorder_edge_clades(build_tree))
    elif _M["backbone"] == "astral_rerooted":
        # ASTRAL topology, but re-rooted to the true root. Isolates the rooting
        # effect: if this lands near the 'true' score, the gap is just rooting.
        bb_path = os.path.join(mul_trees_dir, f"species_tree_{idx_str}.nex")
        if not os.path.exists(bb_path):
            return None
        true_tree = load_nexus_tree(bb_path)
        build_tree = reroot_to_match(astral_tree.copy(), true_tree)
        event_clades = map_clades_by_jaccard(clades, preorder_edge_clades(build_tree))
    else:
        build_tree = astral_tree
        event_clades = clades

    with torch.no_grad():
        inputs = model_inputs_for(sample, "cpu")
        wgd_logits, edge_emb = model(**inputs)
        wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
        pairwise_feat = build_pairwise_feat(sample)

    gt = open(gt_path).read() if os.path.exists(gt_path) else None
    summary = {}
    for T in _M["thresholds"]:
        mul_tree, n_auto, n_allo = decode(model, build_tree, clades, wgd_prob,
                                          edge_emb, pairwise_feat, T,
                                          event_clades=event_clades)
        case_dir = os.path.join(_M["out_base"], f"t{T}", f"sample_{idx_str}")
        os.makedirs(case_dir, exist_ok=True)
        mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
        if gt is not None:
            with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                f.write(gt)
        summary[T] = (n_auto, n_allo)
    return idx_str, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--mul-trees-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--model-config", default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--n", type=int, default=50)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--thresholds", default="0.5",
                        help="comma-separated, e.g. 0.5,0.7,0.9,0.95")
    parser.add_argument("--out-base", default=None,
                        help="output base dir (default: <model-dir>/mul_trees/<config>)")
    parser.add_argument("--backbone", choices=["astral", "true", "astral_rerooted"],
                        default="astral",
                        help="astral: real pipeline. true: build predicted events on "
                             "the true diploid backbone (species_tree_NNNN.nex) to "
                             "isolate the ASTRAL-backbone cost. astral_rerooted: ASTRAL "
                             "topology re-rooted to the true root, to isolate rooting.")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    with open(cfg_path) as f:
        model_config = yaml.safe_load(f).get("model", {})
    edge_feat_dim = int(model_config.get("edge_feat_dim", 9))

    thresholds = [round(float(t), 4) for t in args.thresholds.split(",")]
    out_base = args.out_base or os.path.join(args.model_dir, "mul_trees", args.config)
    indices = list(range(args.start, args.start + args.n))

    print(f"Reconstructing {len(indices)} samples x {len(thresholds)} thresholds "
          f"{thresholds} on the {args.backbone} backbone with {args.workers} "
          f"workers -> {out_base}")

    init_args = (args.model_dir, model_config, args.mul_trees_dir, args.config,
                 args.replicate, thresholds, out_base, edge_feat_dim, args.backbone)

    done = 0
    if args.workers <= 1:
        init_worker(*init_args)
        for idx in indices:
            r = process_index(idx)
            if r:
                done += 1
    else:
        with ProcessPoolExecutor(max_workers=args.workers,
                                 initializer=init_worker, initargs=init_args) as ex:
            futures = [ex.submit(process_index, idx) for idx in indices]
            for i, fut in enumerate(as_completed(futures), 1):
                if fut.result():
                    done += 1
                if i % 25 == 0:
                    print(f"  {i}/{len(indices)} done")

    print(f"\nDone: {done} reconstructed. Output under {out_base}/t<threshold>/")


if __name__ == "__main__":
    main()
