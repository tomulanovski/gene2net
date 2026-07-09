"""Run the reconstruction model on the 21 benchmark networks (the GRAMPA/
Polyphest test set) and write inferred MUL-trees for scoring.

Uses the SAME inputs as GRAMPA-iter: the clean gene trees + the ASTRAL diploid
species tree. Builds a Gene2NetSample on the fly (from_trees computes the 9 edge
features), runs the model, reconstructs at the given threshold, and writes
output.tre next to the ground-truth network — laid out for score_reconstructions.

Usage (final_project env, needs torch):
    python scripts/benchmark_networks.py \
        --model-dir output/reconstruct_allo --config conf_ils_low_10M \
        --replicate 1 --threshold 0.9 \
        --out-dir output/reconstruct_allo/benchmark/conf_ils_low_10M
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree, WGDEvent
from gene2net_gnn.inference.build_strategies import (
    select_event_edges, build_parent_edge_map, infer_copy_bound,
)
from scripts.reconstruct_mul_tree import (
    load_model, model_inputs_for, preorder_edge_clades, build_pairwise_feat,
)

NETWORKS = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019", "Popp_2005", "Wu_2015",
    "Liu_2023", "Ren_2024", "Marcussen_2011", "Marcussen_2012", "Sessa_2012b", "Zhao_2021",
    "Hori_2014", "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014",
]


def load_inverse_taxa_map(taxa_map_path):
    """Map REPLACEMENT -> ORIGINAL from taxa_map.txt (substring fixes).

    The prep renames a taxon to a substring-safe version (e.g. T -> TX) and
    records ORIGINAL<TAB>REPLACEMENT. We invert it to rename predicted leaves
    back to the original names the ground-truth network uses.
    """
    inv = {}
    if not os.path.exists(taxa_map_path):
        return inv
    with open(taxa_map_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                original, replacement = parts[0].strip(), parts[1].strip()
                inv[replacement] = original
    return inv


def rename_leaves(tree, inv_map):
    """Rename leaves replacement->original in place (no-op if map empty)."""
    if not inv_map:
        return
    for leaf in tree.get_leaves():
        if leaf.name in inv_map:
            leaf.name = inv_map[leaf.name]


def load_gene_trees(path, max_trees=500):
    trees = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                trees.append(Tree(line, format=1))
            except Exception:
                continue
            if len(trees) >= max_trees:
                break
    return trees


def build_for_strategy(model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                       strategy, threshold, parent_edge, copy_bound, polyploid_species=None):
    """Select events for a strategy, resolve partners, build the MUL-tree."""
    event_edges = select_event_edges(strategy, wgd_list, threshold, parent_edge, clades, copy_bound)
    if polyploid_species is not None:
        # Known-ploidy prior: keep only WGD edges whose clade is entirely polyploid
        # species (drops false positives fired on diploid clades).
        event_edges = [i for i in event_edges if set(clades[i]) <= polyploid_species]
    events = []
    n_auto = n_allo = 0
    if event_edges:
        query = torch.tensor(event_edges, dtype=torch.long)
        rows = model.compute_partner_scores_rows(edge_emb, query, pairwise_feat)
        n_edges = len(clades)
        for q, i in enumerate(event_edges):
            j = int(rows[q, :n_edges].argmax())
            # Invalid partner (overlaps the WGD clade) -> treat as autopolyploidy.
            if j != i and (clades[j] & clades[i]):
                j = i
            if j == i:
                n_auto += 1
            else:
                n_allo += 1
            events.append(WGDEvent(
                wgd_edge_clade=clades[i], partner_edge_clade=clades[j],
                confidence=float(wgd_list[i]),
            ))
    mul_tree, n_dropped = build_mul_tree(astral_tree, events, return_dropped=True)
    return mul_tree, n_auto, n_allo, n_dropped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", required=True, help="e.g. conf_ils_low_10M")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--strategies", default="raw,collapse,collapse_cap,bound_driven",
                        help="comma-separated build strategies to generate")
    parser.add_argument("--model-config", default=None)
    parser.add_argument("--sim-base",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
    parser.add_argument("--networks-dir",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks")
    parser.add_argument("--max-gene-trees", type=int, default=500)
    parser.add_argument("--root-species-tree", action="store_true",
                        help="Re-root the ASTRAL species tree (hybrid gene-tree+midpoint) "
                             "before features/build. Use with a model trained on rooted data.")
    parser.add_argument("--ploidy-oracle", action="store_true",
                        help="Known-ploidy prior: mask WGD calls to species that are polyploid "
                             "in the ground truth (measures the ceiling of a known-ploidy input).")
    parser.add_argument("--root-mode", choices=["none", "hybrid", "true"], default="none",
                        help="none: keep ASTRAL's arbitrary root. hybrid: infer the root "
                             "(gene-tree consensus+midpoint, ~80%%). true: root at the KNOWN "
                             "root from the ground-truth network (as real studies do with an "
                             "outgroup); falls back to hybrid if it can't be applied.")
    parser.add_argument("--out-base", default=None,
                        help="output base (default: <model-dir>/benchmark/<config>); "
                             "each strategy writes to <out-base>/<strategy>/<network>/")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    with open(cfg_path) as f:
        model_config = yaml.safe_load(f).get("model", {})
    device = "cpu"
    strategies = [s.strip() for s in args.strategies.split(",")]
    out_base = args.out_base or os.path.join(args.model_dir, "benchmark", args.config)
    # Backward-compat: --root-species-tree is the old name for hybrid rooting.
    root_mode = args.root_mode
    if args.root_species_tree and root_mode == "none":
        root_mode = "hybrid"
    print(f"Root mode: {root_mode}")

    model = load_model(args.model_dir, model_config, device)

    done = skipped = 0
    total_dropped = {}
    for net in NETWORKS:
        rep_dir = os.path.join(args.sim_base, net, "processed", args.config,
                               "grampa_input", f"replicate_{args.replicate}")
        gene_trees_path = os.path.join(rep_dir, "clean_trees.tre")
        astral_path = os.path.join(rep_dir, "species.tre")
        gt_path = os.path.join(args.networks_dir, f"{net}.tre")

        if not (os.path.exists(gene_trees_path) and os.path.exists(astral_path)):
            print(f"  SKIP {net}: missing inputs in {rep_dir}")
            skipped += 1
            continue

        gene_trees = load_gene_trees(gene_trees_path, args.max_gene_trees)
        astral_tree = Tree(open(astral_path).read().strip(), format=1)
        species_list = sorted(set(astral_tree.get_leaf_names()))

        inv_map = load_inverse_taxa_map(os.path.join(rep_dir, "taxa_map.txt"))

        # Root the species tree ONCE here so the sample features, the clades, and the
        # parent-edge map all derive from the same (rooted) tree. from_trees then
        # gets root=False since the tree is already rooted.
        if root_mode == "hybrid":
            from gene2net_gnn.data.rooting import hybrid_root
            astral_tree = hybrid_root(astral_tree, gene_trees)
        elif root_mode == "true":
            # Root at the KNOWN root taken from the ground-truth network (the true
            # root researchers get from an outgroup). inv_map is replacement->original;
            # the ground-truth uses original names, so map original->replacement.
            from gene2net_gnn.data.rooting import root_at_reference, hybrid_root
            forward = {orig: repl for repl, orig in inv_map.items()}
            rooted = False
            if os.path.exists(gt_path):
                try:
                    gt_tree = Tree(open(gt_path).read().strip(), format=1)
                    rooted = root_at_reference(astral_tree, gt_tree, name_map=forward)
                except Exception:
                    rooted = False
            if not rooted:
                print(f"  [{net}] true-root unavailable; falling back to hybrid rooting")
                astral_tree = hybrid_root(astral_tree, gene_trees)

        sample = Gene2NetSample.from_trees(
            species_tree=astral_tree, gene_trees=gene_trees, species_list=species_list,
        )
        clades = preorder_edge_clades(astral_tree)
        parent_edge = build_parent_edge_map(astral_tree)
        copy_bound = infer_copy_bound(gene_trees)

        # One inference pass for the network.
        with torch.no_grad():
            inputs = model_inputs_for(sample, device)
            wgd_logits, edge_emb = model(**inputs)
            wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
            pairwise_feat = build_pairwise_feat(sample)
        wgd_list = wgd_prob.tolist()

        # Known-ploidy prior (optional): species that are polyploid in the ground truth,
        # mapped into the ASTRAL/replacement namespace, to mask WGD calls on diploids.
        polyploid_species = None
        if args.ploidy_oracle and os.path.exists(gt_path):
            try:
                gt_leaves = Tree(open(gt_path).read().strip(), format=1).get_leaf_names()
                poly_orig = {n for n, c in Counter(gt_leaves).items() if c > 1}
                forward = {orig: repl for repl, orig in inv_map.items()}
                polyploid_species = {forward.get(o, o) for o in poly_orig}
            except Exception:
                polyploid_species = None

        # Build + write one MUL-tree per strategy.
        counts = []
        for strat in strategies:
            mul_tree, n_auto, n_allo, n_dropped = build_for_strategy(
                model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                strat, args.threshold, parent_edge, copy_bound, polyploid_species,
            )
            rename_leaves(mul_tree, inv_map)
            case_dir = os.path.join(out_base, strat, net)
            os.makedirs(case_dir, exist_ok=True)
            mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
            if os.path.exists(gt_path):
                with open(gt_path) as f:
                    gt = f.read()
                with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                    f.write(gt)
            total_dropped[strat] = total_dropped.get(strat, 0) + n_dropped
            counts.append(f"{strat}={n_auto + n_allo}" + (f" (dropped {n_dropped})" if n_dropped else ""))

        print(f"[{net}] events per strategy: {', '.join(counts)}")
        done += 1

    print(f"\nDone: {done} networks x {len(strategies)} strategies, {skipped} skipped.")
    if any(total_dropped.values()):
        print("Silently-dropped events (clade unfindable after earlier grafts): "
              + ", ".join(f"{s}={n}" for s, n in total_dropped.items()))
    print(f"Output under {out_base}/<strategy>/")
    print("Score each strategy (gene2net env):")
    for strat in strategies:
        print(f"  python scripts/score_reconstructions.py --recon-dir {out_base}/{strat} --workers 8")


if __name__ == "__main__":
    main()
