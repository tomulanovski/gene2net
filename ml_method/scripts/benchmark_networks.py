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
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetSample
from gene2net_gnn.inference.mul_tree_builder import (
    build_mul_tree, WGDEvent, build_mul_tree_two_parent, TwoParentEvent,
)
from gene2net_gnn.inference.build_strategies import (
    select_event_edges, build_parent_edge_map, infer_copy_bound, infer_copy_bound_kernel,
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


def copy_bound_from_multiset(multiset_path, inv_map):
    """Polyphest's EXACT ploidy input: read its multi_set.txt so our decode uses the
    identical per-species copy numbers Polyphest was fed, instead of our own
    infer_copy_bound. This removes the ploidy-estimator confound from the comparison.

    multi_set.txt lists one species per line, duplicate lines meaning extra copies,
    so Counter(lines) is species -> copies. The multiset uses ORIGINAL names, while
    this benchmark works in the substring-fixed (replacement) names; inv_map is
    replacement->original, so we invert it to rename original->replacement.
    """
    from collections import Counter
    forward = {orig: repl for repl, orig in inv_map.items()}   # original -> replacement
    with open(multiset_path) as f:
        names = [line.strip() for line in f if line.strip()]
    return dict(Counter(forward.get(n, n) for n in names))


def build_for_strategy(model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                       strategy, threshold, parent_edge, copy_bound, polyploid_species=None,
                       build_order="size"):
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
    mul_tree, n_dropped = build_mul_tree(astral_tree, events, return_dropped=True, order=build_order)
    return mul_tree, n_auto, n_allo, n_dropped


def build_for_strategy_two_parent(model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                                  strategy, threshold, parent_edge, copy_bound, polyploid_species=None,
                                  build_order="size"):
    """Two-parent variant: predict BOTH parent edges per WGD edge and build with the
    nested-safe graft build. Auto = both slots pick the WGD edge itself."""
    event_edges = select_event_edges(strategy, wgd_list, threshold, parent_edge, clades, copy_bound)
    if polyploid_species is not None:
        event_edges = [i for i in event_edges if set(clades[i]) <= polyploid_species]
    events = []
    n_auto = n_allo = 0
    if event_edges:
        query = torch.tensor(event_edges, dtype=torch.long)
        rows = model.compute_partner_scores_rows(edge_emb, query, pairwise_feat)  # [Q, E, 2]
        n_edges = len(clades)
        for q, i in enumerate(event_edges):
            pa = int(rows[q, :n_edges, 0].argmax())
            pb = int(rows[q, :n_edges, 1].argmax())
            # A predicted parent that overlaps the target clade is invalid -> self on that slot.
            if pa != i and (clades[pa] & clades[i]):
                pa = i
            if pb != i and (clades[pb] & clades[i]):
                pb = i
            if pa == i and pb == i:
                n_auto += 1
            else:
                n_allo += 1
            events.append(TwoParentEvent(
                target_clade=clades[i], parent_a_clade=clades[pa], parent_b_clade=clades[pb],
                confidence=float(wgd_list[i]),
            ))
    mul_tree, n_dropped = build_mul_tree_two_parent(astral_tree, events, mode="graft",
                                                    return_dropped=True, order=build_order)
    return mul_tree, n_auto, n_allo, n_dropped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--config", required=True, help="e.g. conf_ils_low_10M")
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--build-order", choices=["size", "confidence"], default="size",
                        help="order events are stamped onto the backbone. size: smallest "
                             "clade first (nested-safe, default). confidence: strongest "
                             "prediction first, so conflicts resolve by evidence (cf. Polyphest "
                             "inserting clades by support).")
    parser.add_argument("--strategies", default="bound_driven,cap",
                        help="comma-separated build strategies to generate (bound_driven, cap)")
    parser.add_argument("--copy-bound", choices=["kernel", "infer", "multiset"], default="kernel",
                        help="kernel: the method's own kernel-smoothed ploidy from the gene trees "
                             "(default). infer: majority-consensus infer_copy_bound (ablation). "
                             "multiset: read Polyphest's exact multi_set.txt (bit-exact parity).")
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
    parser.add_argument("--root-mode", choices=["none", "hybrid", "true"], default="hybrid",
                        help="none: keep ASTRAL's arbitrary root. hybrid: infer the root "
                             "(gene-tree consensus+midpoint, ~80%%). true: root at the KNOWN "
                             "root from the ground-truth network (as real studies do with an "
                             "outgroup); falls back to hybrid if it can't be applied.")
    parser.add_argument("--out-base", default=None,
                        help="output base (default: <model-dir>/benchmark/<config>); "
                             "each strategy writes to <out-base>/<strategy>/<network>/")
    parser.add_argument("--profile", action="store_true",
                        help="Time the per-network compute stages (features, forward, build) "
                             "and print a summary. Excludes ASTRAL, which is loaded precomputed.")
    parser.add_argument("--ploidy-check", action="store_true",
                        help="Compare the built MUL-tree's per-species copy counts against "
                             "infer_copy_bound and report how often we fall SHORT of the inferred "
                             "ploidy (headroom for a target-ploidy build). Read-only diagnostic.")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config or os.path.join(base_dir, "configs", "reconstruct_final.yaml")
    with open(cfg_path) as f:
        model_config = yaml.safe_load(f).get("model", {})
    # On-the-fly feature ablations must match how the model was trained.
    from gene2net_gnn.training.trainer_reconstruct import set_feature_opts
    set_feature_opts(
        coclust_condition_on_dup=model_config.get("coclust_condition_on_dup", False),
        use_n_eff=model_config.get("use_n_eff", False),
    )
    device = "cpu"
    strategies = [s.strip() for s in args.strategies.split(",")]
    out_base = args.out_base or os.path.join(args.model_dir, "benchmark", args.config)
    # Backward-compat: --root-species-tree is the old name for hybrid rooting.
    root_mode = args.root_mode
    if args.root_species_tree and root_mode == "none":
        root_mode = "hybrid"
    print(f"Root mode: {root_mode}")

    model = load_model(args.model_dir, model_config, device)
    two_parent = getattr(model, "n_parents", 1) >= 2
    print(f"Decode: {'two-parent graft' if two_parent else 'one-partner'} "
          f"(model n_parents={getattr(model, 'n_parents', 1)})")

    done = skipped = errored = 0
    total_dropped = {}
    prof = {"feat": 0.0, "fwd": 0.0, "build": 0.0, "n": 0, "n_gt": 0}
    # under-shoot accounting over polyploid species (inferred bound >= 2)
    pchk = {"n_poly": 0, "n_short": 0, "n_match": 0, "n_missing": 0, "copies_short": 0}
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

        # One malformed input tree must not kill the whole run. Report it loudly and
        # continue so the other networks still produce output. This is not a silent
        # skip: the network and error are printed and counted in the final summary.
        try:
            gene_trees = load_gene_trees(gene_trees_path, args.max_gene_trees)
            astral_tree = Tree(open(astral_path).read().strip(), format=1)
        except Exception as e:
            print(f"  ERROR {net}: could not parse inputs "
                  f"({type(e).__name__}: {e}) -- skipping this network, continuing.")
            errored += 1
            continue
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
                    gt_tree = Tree(open(gt_path, encoding='utf-8').read().strip(), format=1)
                    rooted = root_at_reference(astral_tree, gt_tree, name_map=forward)
                except Exception:
                    rooted = False
            if not rooted:
                print(f"  [{net}] true-root unavailable; falling back to hybrid rooting")
                astral_tree = hybrid_root(astral_tree, gene_trees)

        _t = time.perf_counter()
        sample = Gene2NetSample.from_trees(
            species_tree=astral_tree, gene_trees=gene_trees, species_list=species_list,
        )
        clades = preorder_edge_clades(astral_tree)
        parent_edge = build_parent_edge_map(astral_tree)
        if args.copy_bound == "multiset":
            ms_path = os.path.join(rep_dir.replace("grampa_input", "polyphest_input"),
                                   "multi_set.txt")
            if not os.path.exists(ms_path):
                print(f"  ERROR {net}: --copy-bound multiset but {ms_path} missing "
                      "-- skipping this network, continuing.")
                errored += 1
                continue
            copy_bound = copy_bound_from_multiset(ms_path, inv_map)
        elif args.copy_bound == "kernel":
            copy_bound = infer_copy_bound_kernel(gene_trees)
        else:
            copy_bound = infer_copy_bound(gene_trees)
        prof["feat"] += time.perf_counter() - _t

        # One inference pass for the network.
        _t = time.perf_counter()
        with torch.no_grad():
            inputs = model_inputs_for(sample, device)
            wgd_logits, edge_emb = model(**inputs)
            wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
            pairwise_feat = build_pairwise_feat(sample)
        wgd_list = wgd_prob.tolist()
        prof["fwd"] += time.perf_counter() - _t
        prof["n"] += 1
        prof["n_gt"] += len(gene_trees)

        # Known-ploidy prior (optional): replace the copy-count-INFERRED bound (which
        # dup/loss inflates -> too many events -> over-prediction, Polyphest's failure
        # mode) with the TRUE ploidy from the ground truth. This is the per-species
        # bound the 'cap' strategy uses to limit events.
        if args.ploidy_oracle and os.path.exists(gt_path):
            try:
                gt_counts = Counter(Tree(open(gt_path, encoding='utf-8').read().strip(), format=1).get_leaf_names())
                true_bound = {s: gt_counts.get(inv_map.get(s, s), 1) for s in copy_bound}
                n_diff = sum(1 for s in copy_bound if true_bound[s] != copy_bound[s])
                n_inflated = sum(1 for s in copy_bound if copy_bound[s] > true_bound[s])
                print(f"  [{net}] ploidy prior: {n_diff} species corrected "
                      f"({n_inflated} were over-estimated)")
                copy_bound = true_bound
            except Exception as e:
                print(f"  [{net}] ploidy prior failed: {e!r}")

        # Build + write one MUL-tree per strategy.
        counts = []
        _t = time.perf_counter()
        for strat in strategies:
            build_fn = build_for_strategy_two_parent if two_parent else build_for_strategy
            mul_tree, n_auto, n_allo, n_dropped = build_fn(
                model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                strat, args.threshold, parent_edge, copy_bound,
                build_order=args.build_order,
            )
            if args.ploidy_check:
                # mul_tree is still in internal names here, matching copy_bound.
                actual = Counter(mul_tree.get_leaf_names())
                for sp, bound in copy_bound.items():
                    if bound < 2:
                        continue  # only polyploids-per-inferred-ploidy
                    a = actual.get(sp, 0)
                    pchk["n_poly"] += 1
                    if a == 0:
                        pchk["n_missing"] += 1
                    if a < bound:
                        pchk["n_short"] += 1
                        pchk["copies_short"] += bound - a
                    elif a == bound:
                        pchk["n_match"] += 1
            rename_leaves(mul_tree, inv_map)
            case_dir = os.path.join(out_base, strat, net)
            os.makedirs(case_dir, exist_ok=True)
            mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
            if os.path.exists(gt_path):
                with open(gt_path, encoding='utf-8') as f:
                    gt = f.read()
                with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                    f.write(gt)
            total_dropped[strat] = total_dropped.get(strat, 0) + n_dropped
            counts.append(f"{strat}={n_auto + n_allo}" + (f" (dropped {n_dropped})" if n_dropped else ""))

        prof["build"] += time.perf_counter() - _t
        print(f"[{net}] events per strategy: {', '.join(counts)}")
        done += 1

    if args.profile and prof["n"]:
        n = prof["n"]
        per = lambda k: prof[k] / n
        print("\n" + "=" * 60)
        print(f"RUNTIME PROFILE  ({n} networks, mean {prof['n_gt']/n:.0f} gene trees each)")
        print(f"  features (from_trees + clades + copy bound): {per('feat'):.3f} s/net")
        print(f"  forward  (model inputs + GNN + pairwise)   : {per('fwd'):.3f} s/net")
        print(f"  build    (decode + MUL-tree per strategy)  : {per('build'):.3f} s/net")
        print(f"  our compute total (excludes ASTRAL)        : "
              f"{per('feat')+per('fwd')+per('build'):.3f} s/net")
        print("  NOTE: ASTRAL is loaded precomputed here. For an end-to-end comparison to")
        print("  Polyphest/GRAMPA, add ASTRAL inference time and use their wall times from logs.")
        print("=" * 60)

    if args.ploidy_check and pchk["n_poly"]:
        p = pchk["n_poly"]
        print("\n" + "=" * 60)
        print(f"PLOIDY UNDER-SHOOT CHECK  ({p} polyploid species, inferred bound >= 2)")
        print(f"  reach the inferred ploidy (actual == bound): {pchk['n_match']}/{p} = {pchk['n_match']/p:.2f}")
        print(f"  fall SHORT (actual < bound)                : {pchk['n_short']}/{p} = {pchk['n_short']/p:.2f}")
        print(f"  missing entirely (0 copies)                : {pchk['n_missing']}/{p} = {pchk['n_missing']/p:.2f}")
        print(f"  total copies short of inferred bound       : {pchk['copies_short']} "
              f"({pchk['copies_short']/p:.2f} per polyploid species)")
        print("  If we fall short a lot, using the inferred ploidy as a TARGET (not just a cap)")
        print("  is real headroom. If most reach the bound, the edit gap is placement, not count.")
        print("=" * 60)

    print(f"\nDone: {done} networks x {len(strategies)} strategies, "
          f"{skipped} skipped (missing inputs), {errored} errored (unparseable inputs).")
    if any(total_dropped.values()):
        print("Silently-dropped events (clade unfindable after earlier grafts): "
              + ", ".join(f"{s}={n}" for s, n in total_dropped.items()))
    print(f"Output under {out_base}/<strategy>/")
    print("Score each strategy (gene2net env):")
    for strat in strategies:
        print(f"  python scripts/score_reconstructions.py --recon-dir {out_base}/{strat} --workers 8")


if __name__ == "__main__":
    main()
