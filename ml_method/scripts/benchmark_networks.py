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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.mul_tree_builder import build_mul_tree
from scripts.reconstruct_mul_tree import (
    load_model, model_inputs_for, preorder_edge_clades, decode, build_pairwise_feat,
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", required=True, help="e.g. conf_ils_low_10M")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--model-config", default=None)
    parser.add_argument("--sim-base",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations")
    parser.add_argument("--networks-dir",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/networks")
    parser.add_argument("--max-gene-trees", type=int, default=500)
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base_dir, "configs", "reconstruct.yaml")
    with open(cfg_path) as f:
        model_config = yaml.safe_load(f).get("model", {})
    device = "cpu"
    out_dir = args.out_dir or os.path.join(args.model_dir, "benchmark", args.config)
    os.makedirs(out_dir, exist_ok=True)

    model = load_model(args.model_dir, model_config, device)

    done = skipped = 0
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

        sample = Gene2NetSample.from_trees(
            species_tree=astral_tree, gene_trees=gene_trees, species_list=species_list,
        )
        clades = preorder_edge_clades(astral_tree)
        with torch.no_grad():
            inputs = model_inputs_for(sample, device)
            wgd_logits, edge_emb = model(**inputs)
            wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
            pairwise_feat = build_pairwise_feat(sample)
        mul_tree, n_auto, n_allo = decode(
            model, astral_tree, clades, wgd_prob, edge_emb, pairwise_feat, args.threshold
        )

        # Reverse substring fixes so leaf names match the ground-truth network.
        inv_map = load_inverse_taxa_map(os.path.join(rep_dir, "taxa_map.txt"))
        rename_leaves(mul_tree, inv_map)

        case_dir = os.path.join(out_dir, net)
        os.makedirs(case_dir, exist_ok=True)
        mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
        if os.path.exists(gt_path):
            with open(gt_path) as f:
                gt = f.read()
            with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                f.write(gt)

        print(f"[{net}] events={n_auto + n_allo} (auto={n_auto}, allo={n_allo}) -> {case_dir}/output.tre")
        done += 1

    print(f"\nDone: {done} networks reconstructed, {skipped} skipped. Output under {out_dir}")
    print("Score with: python scripts/score_reconstructions.py "
          f"--recon-dir {out_dir} --workers 8   (in the gene2net env)")


if __name__ == "__main__":
    main()
