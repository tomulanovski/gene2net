"""Score a trained reconstruction model's WHOLE-NETWORK reconstruction quality
(edit distance / reticulation metrics) on the VALIDATION split of the training
data — the same networks the model was SELECTED on, and never trained on.

Why this script exists
----------------------
Hyperparameters were selected on per-event partner accuracy (a proxy). The
benchmark (benchmark_networks.py) scores whole-network reconstruction with
score_reconstructions.py. This script closes the gap: it runs the model on the
held-out VALIDATION samples and writes inferred MUL-trees laid out exactly like
benchmark_networks.py / oracle_test.py, so score_reconstructions.py can score
them on the SAME metric the benchmark uses. That lets us pick hyperparameters on
the target metric instead of the proxy, WITHOUT leaking: we score only the
validation networks the model did not train on.

=============================================================================
CORRECTNESS-CRITICAL PART 1: reproducing the validation split IDENTICALLY
=============================================================================
train_reconstruct.run_training builds the split like this (verbatim):

    all_samples = []
    for data_dir in data_dirs:                       # dirs in the given order
        dataset = Gene2NetDataset(data_dir, clade_labels=..., away_labels=...)
        for i in range(len(dataset)):
            sample = dataset[i]
            if sample.labels is None:                 # drop no-label
                continue
            ef = sample.species_tree_edge_features
            if ef is None or ef.shape[1] != expected_edge_dim:   # drop wrong dim
                continue
            all_samples.append(sample)
    ...
    val_split = float(train_config.get("val_split", 0.2))
    random.seed(42)
    random.shuffle(all_samples)
    n_val = int(len(all_samples) * val_split)
    val_samples = all_samples[:n_val]

We reproduce this EXACTLY:
  * same `--data-dir` order (must be the six training_rooted/<config> dirs used
    to train the model — see submit_finalize.sh / submit_away_experiment.sh);
  * same label flags (`--away-labels` / `--clade-labels`) — the models here were
    trained with AWAY_LABELS=1, so `--away-labels` is required to match which
    samples pass the `labels is None` filter and therefore which end up in val;
  * same `expected_edge_dim` (model config `edge_feat_dim`, default 9);
  * same `val_split` read from the model config's `training.val_split` (0.2);
  * `random.seed(42); random.shuffle(...)` then the 20% prefix.

The ONLY addition is that we shuffle a parallel list of (example_dir, sample)
tuples instead of bare samples. `random.shuffle` is Fisher-Yates keyed on list
LENGTH only (it swaps by position, independent of element content), so shuffling
the paired list with the same seed produces the identical permutation of
positions. The val prefix therefore selects the identical underlying samples,
now each still carrying its source directory. This is why the split matches
what the model was selected on.

=============================================================================
CORRECTNESS-CRITICAL PART 2: mapping each val sample to mul_tree_{idx}.nex
=============================================================================
Gene2NetSample carries NO index attribute. The index lives in the on-disk
directory name. Gene2NetDataset exposes `self.example_dirs` (sorted), and each
example directory is named `sample_{idx:04d}` (see check_rooted_packaging.py:
`idx = name.replace("sample_", "")`, and reconstruct_mul_tree.py which pairs
`training/<config>/sample_{idx_str}` with `mul_tree_{idx_str}.nex`). So we track
`dataset.example_dirs[i]` alongside each kept sample and parse the index from its
basename:  idx = int(basename.split("_")[-1]);  idx_str = f"{idx:04d}".
The ground truth is then `<mul-trees-dir>/mul_tree_{idx_str}.nex`.

NOTE on the index being GLOBAL: the same `mul_tree_{idx}.nex` is the ground-truth
network for `sample_{idx}` in EVERY config directory (a config only changes the
simulated gene trees under that fixed network). So the same idx can appear in
several of the six pooled config dirs, and both copies can land in val (they are
distinct samples with different gene trees, both legitimately scored against the
same network). To keep them from overwriting each other on disk we name the
output case directory `<config>__sample_{idx_str}` (score_reconstructions.py
simply globs any subdir containing an output.tre and uses the dir name as the
"sample" label, so this layout scores correctly). This is a deliberate deviation
from a bare `sample_{idx}` name, made for correctness.

=============================================================================
CORRECTNESS-CRITICAL PART 3: the build backbone comes FROM the packaged sample
=============================================================================
The training_rooted/ samples are ROOTED (at the true root, via hybrid_root at
packaging time) and their edge ordering/features reflect that rooted tree — which
DIFFERS from the arbitrary-rooted on-disk astral_species.tre (see
check_rooted_packaging.py). To keep the model's per-edge predictions aligned with
the clades we build on, we reconstruct the ete3 species tree DIRECTLY from the
packaged sample's `species_tree_edge_index` + `species_tree_node_names` +
`species_tree_is_leaf`. That guarantees clades[k] (preorder, non-root) aligns
with model edge k (model_inputs_for reorders edges to preorder), exactly as in
benchmark_networks.py where the sample was built from the same tree. We do NOT
load the on-disk astral tree here.

The per-species copy bound needed by `bound_driven` / `infer_copy_bound` is
computed from the packaged gene trees: each packaged gene tree stores its leaf
species (gene_tree_species_ids + gene_tree_leaf_masks + species_list), so we
rebuild a trivial star tree carrying the exact leaf-name multiset per gene tree
and feed those to the unchanged `infer_copy_bound`. infer_copy_bound only counts
per-tree leaf-name multiplicities, so a star tree yields identical bounds to the
original topology — no logic is duplicated.

=============================================================================
Inference + build path
=============================================================================
Identical to benchmark_networks.py's one-partner path (`--parents head`):
  select_event_edges(strategy=bound_driven, ...) -> partner argmax from the
  model's partner head -> build_mul_tree. We import and call
  benchmark_networks.build_for_strategy (one-partner) /
  build_for_strategy_two_parent (>=2-parent heads) unchanged. Both target models
  are one-parter (n_parents: 1), so the one-partner path is used.

Runs in the `final_project` conda env (needs torch). Scoring is a SEPARATE step
in the `gene2net` env; the exact score_reconstructions.py command is printed at
the end.

Usage (final_project env):
    python scripts/score_val_reconstruction.py \
        --data-dir data/mul_trees_2k/training_rooted/ils_low \
                   data/mul_trees_2k/training_rooted/ils_medium \
                   data/mul_trees_2k/training_rooted/ils_high \
                   data/mul_trees_2k/training_rooted/dup_loss_low_ne1M \
                   data/mul_trees_2k/training_rooted/dup_loss_medium_ne1M \
                   data/mul_trees_2k/training_rooted/dup_loss_high_ne1M \
        --model-dir output/reconstruct_final \
        --model-config configs/reconstruct_final.yaml \
        --away-labels \
        --out-base output/reconstruct_final/val_recon
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import torch
import yaml
from ete3 import Tree

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.inference.build_strategies import (
    infer_copy_bound, build_parent_edge_map,
)
from scripts.reconstruct_mul_tree import (
    load_model, model_inputs_for, preorder_edge_clades, build_pairwise_feat,
)
from scripts.benchmark_networks import (
    build_for_strategy, build_for_strategy_two_parent,
)

# Seed for the (documented, deterministic) subsample of the val split taken for
# speed. Separate from the seed-42 split shuffle so it never disturbs the split.
SUBSAMPLE_SEED = 12345


def ete3_from_sample(sample):
    """Rebuild the ROOTED ete3 species tree from the packaged sample so the build
    backbone is exactly the tree the model's edges/features describe.

    Uses species_tree_edge_index (even entries are parent->child, preorder-grouped,
    child order preserved), node_names, and is_leaf. Leaf nodes get their species
    name; internal nodes stay unnamed. Child insertion order is preserved so the
    resulting preorder matches reorder_edge_index_preorder (model edge order).
    """
    ei = sample.species_tree_edge_index
    names = sample.species_tree_node_names
    is_leaf = sample.species_tree_is_leaf
    children = {}
    for k in range(0, ei.shape[1], 2):
        p = int(ei[0, k])
        c = int(ei[1, k])
        children.setdefault(p, []).append(c)
    all_children = {c for cs in children.values() for c in cs}
    root = next(n for n in range(len(names)) if n not in all_children)

    def build(n):
        node = Tree()
        if bool(is_leaf[n]):
            node.name = names[n]
        for c in children.get(n, []):
            node.add_child(build(c))
        return node

    return build(root)


def star_gene_trees_from_sample(sample):
    """Rebuild one star tree per packaged gene tree, carrying the exact leaf-name
    multiset (species copies). infer_copy_bound only counts per-tree leaf-name
    multiplicities, so star topology gives identical bounds to the real gene trees.
    """
    stars = []
    species_list = sample.species_list
    for ids, mask in zip(sample.gene_tree_species_ids, sample.gene_tree_leaf_masks):
        leaf_ids = ids[mask.bool()].tolist()
        t = Tree()
        for i in leaf_ids:
            t.add_child(name=species_list[i])
        stars.append(t)
    return stars


def collect_val_handles(data_dirs, expected_edge_dim, val_split,
                        clade_labels, away_labels):
    """Reproduce run_training's pool -> filter -> seed-42 shuffle -> 20% prefix,
    but store only lightweight (data_dir, index, example_dir) HANDLES, never the
    loaded samples. Holding all ~12k pooled samples (each carrying gene trees) in
    RAM at once can exhaust an interactive session and swap (which looks like a
    hang); we keep only handles here and reload the ~max-samples we actually score.

    Membership is IDENTICAL to run_training: same pool order, same filters (labels
    first, then edge-dim), same seed-42 shuffle of positions, same 20% prefix.
    random.shuffle permutes by position, so shuffling handles yields the same
    permutation as shuffling the samples themselves.
    """
    all_handles = []
    no_label = wrong_dim = 0
    for data_dir in data_dirs:
        print(f"Loading dataset from {data_dir}...", flush=True)
        dataset = Gene2NetDataset(data_dir, clade_labels=clade_labels,
                                  away_labels=away_labels)
        n = len(dataset)
        dir_loaded = 0
        for i in range(n):
            sample = dataset[i]                      # loaded, filtered, then freed
            if sample.labels is None:
                no_label += 1
            elif (sample.species_tree_edge_features is None or
                  sample.species_tree_edge_features.shape[1] != expected_edge_dim):
                wrong_dim += 1
            else:
                all_handles.append((data_dir, i, dataset.example_dirs[i]))
                dir_loaded += 1
            if (i + 1) % 500 == 0:
                print(f"    scanned {i + 1}/{n} ({dir_loaded} kept)", flush=True)
        print(f"  {dir_loaded} samples from {os.path.basename(data_dir)}", flush=True)

    print(f"Total: {len(all_handles)} samples "
          f"({no_label} no-label, {wrong_dim} wrong-edge-dim filtered)")

    # Identical to run_training: seed 42, shuffle the pool, take the 20% prefix.
    random.seed(42)
    random.shuffle(all_handles)
    n_val = int(len(all_handles) * val_split)
    val_handles = all_handles[:n_val]
    print(f"Validation split: {len(val_handles)} samples "
          f"(val_split={val_split}, seed=42 shuffle, 20% prefix)")
    return val_handles


def idx_str_from_dir(example_dir):
    """Parse the network index from a `sample_{idx:04d}` directory basename."""
    base = os.path.basename(os.path.normpath(example_dir))
    idx = int(base.split("_")[-1])
    return f"{idx:04d}"


def config_from_dir(example_dir):
    """The config name is the parent dir of the sample dir (training_rooted/<cfg>)."""
    return os.path.basename(os.path.normpath(os.path.dirname(example_dir)))


def main():
    parser = argparse.ArgumentParser(
        description="Score a model's whole-network reconstruction on the "
                    "held-out validation split (leak-free hyperparameter selection).")
    parser.add_argument("--data-dir", required=True, nargs="+",
                        help="The SAME six training_rooted/<config> dirs, in the "
                             "SAME order, used to train the model. Order matters: "
                             "the val split is reproduced from this exact pooling.")
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--model-config", required=True,
                        help="Model config yaml (e.g. configs/reconstruct_final.yaml). "
                             "Its model.edge_feat_dim and training.val_split are used "
                             "to reproduce the split exactly.")
    parser.add_argument("--mul-trees-dir",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/ml_method/data/mul_trees_2k",
                        help="Directory holding mul_tree_{idx:04d}.nex ground-truth networks.")
    parser.add_argument("--out-base", required=True,
                        help="Output base; each strategy writes to "
                             "<out-base>/<strategy>/<config>__sample_{idx}/.")
    parser.add_argument("--max-samples", type=int, default=300,
                        help="Cap the number of val samples scored (seeded subsample "
                             "for speed). Default 300.")
    parser.add_argument("--strategy", default="bound_driven",
                        help="Event-selection strategy (default bound_driven).")
    parser.add_argument("--parents", choices=["head"], default="head",
                        help="Partner resolution. 'head' = the model's partner head "
                             "(the one-partner path). Only 'head' is supported here: "
                             "coclust would need real gene-tree topology, which is not "
                             "reconstructed from packaged samples.")
    parser.add_argument("--build-order", choices=["size", "confidence"], default="size")
    parser.add_argument("--threshold", type=float, default=0.9,
                        help="Passed through to the build; unused by bound_driven "
                             "(which is threshold-free).")
    parser.add_argument("--clade-labels", action="store_true",
                        help="Reproduce the split with labels_clade.pkl (match a model "
                             "trained with --clade-labels).")
    parser.add_argument("--away-labels", action="store_true",
                        help="Reproduce the split with labels_away.pkl (match a model "
                             "trained with --away-labels). REQUIRED for reconstruct_final "
                             "and the op_away one-parter (both trained AWAY_LABELS=1).")
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()

    base_dir = os.path.join(os.path.dirname(__file__), "..")
    cfg_path = args.model_config
    if not os.path.isabs(cfg_path) and not os.path.exists(cfg_path):
        cfg_path = os.path.join(base_dir, args.model_config)
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    model_config = cfg.get("model", {})
    train_config = cfg.get("training", {})
    expected_edge_dim = int(model_config.get("edge_feat_dim", 9))
    val_split = float(train_config.get("val_split", 0.2))

    # On-the-fly feature ablations must match how the model was trained (same as
    # benchmark_networks.py). Must run before any model forward.
    from gene2net_gnn.training.trainer_reconstruct import set_feature_opts
    set_feature_opts(
        coclust_condition_on_dup=model_config.get("coclust_condition_on_dup", False),
        use_n_eff=model_config.get("use_n_eff", False),
    )

    print("=" * 70)
    print("Validation-split reconstruction scoring")
    print("=" * 70)
    print(f"Model dir:   {args.model_dir}")
    print(f"Model cfg:   {cfg_path}")
    print(f"Strategy:    {args.strategy} | parents: {args.parents} | "
          f"build-order: {args.build_order}")
    print(f"Labels:      {'away' if args.away_labels else ('clade' if args.clade_labels else 'default')}")
    print(f"Out base:    {args.out_base}")
    print("=" * 70)

    # ---- Reproduce the validation split exactly (lightweight handles only). ----
    val_handles = collect_val_handles(
        args.data_dir, expected_edge_dim, val_split,
        clade_labels=args.clade_labels, away_labels=args.away_labels,
    )
    if not val_handles:
        raise RuntimeError("Validation split is empty — check --data-dir and label flags.")

    # ---- Seeded subsample of up to --max-samples for speed. ----
    rng = random.Random(SUBSAMPLE_SEED)
    if len(val_handles) > args.max_samples:
        selected = rng.sample(val_handles, args.max_samples)
        print(f"Subsampled {len(selected)}/{len(val_handles)} val samples "
              f"(seed={SUBSAMPLE_SEED}).")
    else:
        selected = list(val_handles)
        print(f"Scoring all {len(selected)} val samples (<= max-samples).")

    # ---- Load the model (one-partner vs two-parent inferred from checkpoint). ----
    model = load_model(args.model_dir, model_config, args.device)
    two_parent = getattr(model, "n_parents", 1) >= 2
    build_fn = build_for_strategy_two_parent if two_parent else build_for_strategy
    print(f"Decode: {'two-parent graft' if two_parent else 'one-partner'} "
          f"(model n_parents={getattr(model, 'n_parents', 1)})")

    # ---- Reload ONLY the selected samples, grouped by data dir (memory-light). ----
    from collections import defaultdict
    by_dir = defaultdict(list)
    for data_dir, i, example_dir in selected:
        by_dir[data_dir].append((i, example_dir))

    out_strategy_dir = os.path.join(args.out_base, args.strategy)
    done = skipped = 0
    for data_dir, items in by_dir.items():
        dataset = Gene2NetDataset(data_dir, clade_labels=args.clade_labels,
                                  away_labels=args.away_labels)
        for i, example_dir in items:
            sample = dataset[i]
            idx_str = idx_str_from_dir(example_dir)
            config_name = config_from_dir(example_dir)
            case_name = f"{config_name}__sample_{idx_str}"
            mul_path = os.path.join(args.mul_trees_dir, f"mul_tree_{idx_str}.nex")

            # ONLY documented skip: ground truth genuinely absent. Else propagates.
            if not os.path.exists(mul_path):
                print(f"  SKIP {case_name}: ground-truth {mul_path} missing")
                skipped += 1
                continue

            # Backbone + clades + parent map + copy bound, all from the packaged sample.
            astral_tree = ete3_from_sample(sample)
            clades = preorder_edge_clades(astral_tree)
            parent_edge = build_parent_edge_map(astral_tree)
            copy_bound = infer_copy_bound(star_gene_trees_from_sample(sample))

            # One model forward for the sample.
            with torch.no_grad():
                inputs = model_inputs_for(sample, args.device)
                wgd_logits, edge_emb = model(**inputs)
                wgd_prob = torch.softmax(wgd_logits, dim=-1)[:, 1]
                pairwise_feat = build_pairwise_feat(sample)
            wgd_list = wgd_prob.tolist()

            # Fail loud if the model edges and the rebuilt-tree clades disagree.
            if len(clades) != wgd_prob.shape[0]:
                raise RuntimeError(
                    f"{case_name}: clade count {len(clades)} != model edge count "
                    f"{wgd_prob.shape[0]} — backbone/edge alignment broken.")

            mul_tree, n_auto, n_allo, n_dropped = build_fn(
                model, astral_tree, clades, wgd_list, edge_emb, pairwise_feat,
                args.strategy, args.threshold, parent_edge, copy_bound,
                build_order=args.build_order,
            )

            case_dir = os.path.join(out_strategy_dir, case_name)
            os.makedirs(case_dir, exist_ok=True)
            mul_tree.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
            with open(mul_path) as f:
                gt = f.read()
            with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
                f.write(gt)

            done += 1
            if done % 25 == 0:
                print(f"  {done}/{len(selected)} reconstructed "
                      f"(last: {case_name}, events={n_auto + n_allo}"
                      + (f", dropped {n_dropped}" if n_dropped else "") + ")", flush=True)

    print(f"\nDone: {done} val samples reconstructed, {skipped} skipped "
          f"(missing ground truth).")
    print(f"Output under {out_strategy_dir}/")
    print("\nScore the validation reconstructions (gene2net env):")
    print(f"  python scripts/score_reconstructions.py --recon-dir {out_strategy_dir} --workers 8")


if __name__ == "__main__":
    main()
