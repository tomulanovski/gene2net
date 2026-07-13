"""Two-parent construction oracle (Phase 0).

Feed each allopolyploid's TWO parents through the two-parent build primitive, on a
chosen backbone, and score. This separates the CONSTRUCTION question (can a
detach-and-graft-both build represent a two-parent allopolyploid well enough to
beat Polyphest) from the PREDICTION question (can we identify the two parents).

Parent sources (--parents):
  true    : A = target's sibling in the TRUE species tree, B = metadata partner.
            The construction CEILING (perfect placement).
  coclust : A, B = the top-2 co-clustering species from the gene trees (no labels,
            no training). What the construction fix buys with today's signal.
            Applies to single-species allo targets; auto and clade-level allo keep
            their true construction.

Backbones (--backbone):
  astral  : the real pipeline's backbone. The decisive cell.
  true    : faithfulness check; should reproduce ground truth (edit ~0).

Output mirrors oracle_test.py: <out>/sample_NNNN/{output.tre, ground_truth.nex}.
Also writes <out>/mapping_scores.csv (per-event Jaccard snap scores + parent
source) so a low ceiling is not silently a mapping artifact.

Run in the gene2net env, then score with score_reconstructions.py.
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.inference.mul_tree_builder import (
    TwoParentEvent, build_mul_tree_two_parent, _find_node_by_leaf_set,
)
from gene2net_gnn.data.label_extractor import _get_edge_bipartitions, _jaccard
from gene2net_gnn.data.features import compute_clustering_profile
from gene2net_gnn.data.rooting import hybrid_root


def load_nexus_tree(path):
    for line in open(path).read().split("\n"):
        line = line.strip()
        if line.lower().startswith("tree") and "=" in line:
            return Tree(line.split("=", 1)[1].strip(), format=1)
    raise ValueError(f"No tree found in {path}")


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
                pass
            if len(trees) >= max_trees:
                break
    return trees


def sibling_of_clade(true_tree, clade):
    """Leaf set of the clade's sibling in the true tree (the target's OTHER parent)."""
    node = _find_node_by_leaf_set(true_tree, clade)
    if node is None or node.up is None:
        return None
    sib = set()
    for ch in node.up.get_children():
        if ch is not node:
            sib |= set(ch.get_leaf_names())
    sib -= set(clade)
    return frozenset(sib) if sib else None


def snap_to_backbone(clade, backbone):
    """Best-Jaccard node leaf set on the backbone. Guarantees the clade is findable."""
    best_set, best_score = frozenset(), -1.0
    for _, below in _get_edge_bipartitions(backbone):
        s = _jaccard(clade, below)
        if s > best_score:
            best_score, best_set = s, below
    return best_set, best_score


def predicted_parents_coclust(gene_trees, all_species, X):
    """Top co-clustering species for X, most-to-least (the non-learned parent guess)."""
    profile = compute_clustering_profile(gene_trees, X, all_species)
    return sorted(profile, key=profile.get, reverse=True)


def two_parent_events(metadata_events, true_tree, backbone, mode="true",
                      gene_trees=None, all_species=None):
    """Build TwoParentEvents + per-event snap-score records.

    mode="true": parents from ground truth (A = true sibling, B = metadata partner).
    mode="coclust": single-species allo parents from co-cluster top-2; auto and
        clade-level allo fall back to the true construction.
    """
    events, scores = [], []
    for ev in metadata_events:
        etype = ev.get("event_type")
        target = frozenset(ev.get("target_clade") or [])
        if not target:
            continue
        t_snap, t_score = snap_to_backbone(target, backbone)

        if etype == "auto":
            events.append(TwoParentEvent(t_snap, t_snap, t_snap, 1.0))
            scores.append({"type": "auto", "source": "true",
                           "target": t_score, "A": 1.0, "B": 1.0})
            continue
        if etype != "allo":
            continue

        B_true = frozenset(ev.get("partner_clade") or [])
        A_true = sibling_of_clade(true_tree, target)
        if not B_true or A_true is None:
            continue

        use_coclust = (mode == "coclust" and len(target) == 1
                       and gene_trees is not None and all_species is not None)
        if use_coclust:
            X = next(iter(target))
            ranked = predicted_parents_coclust(gene_trees, all_species, X)
            if len(ranked) < 2:
                continue
            a_snap, a_score = snap_to_backbone(frozenset({ranked[0]}), backbone)
            b_snap, b_score = snap_to_backbone(frozenset({ranked[1]}), backbone)
            source = "coclust"
        else:
            a_snap, a_score = snap_to_backbone(A_true, backbone)
            b_snap, b_score = snap_to_backbone(B_true, backbone)
            source = "true" if mode == "true" else "clade_true"

        events.append(TwoParentEvent(t_snap, a_snap, b_snap, 1.0))
        scores.append({"type": "allo", "source": source,
                       "target": t_score, "A": a_score, "B": b_score})
    return events, scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--replicate", type=int, default=1)
    ap.add_argument("--backbone", choices=["astral", "true"], default="astral")
    ap.add_argument("--parents", choices=["true", "coclust"], default="true")
    ap.add_argument("--root-backbone", choices=["none", "hybrid"], default="none",
                    help="hybrid: root the ASTRAL backbone with hybrid_root (matches the real "
                         "pipeline). Ignored for --backbone true (already correctly rooted).")
    ap.add_argument("--max-gene-trees", type=int, default=500)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if args.parents == "coclust" and args.backbone != "astral":
        print("NOTE: coclust parents model the real pipeline; forcing --backbone astral.")
        args.backbone = "astral"

    # Gene trees are needed for co-cluster parents and/or hybrid rooting.
    need_gene_trees = (args.parents == "coclust"
                       or (args.root_backbone == "hybrid" and args.backbone == "astral"))

    os.makedirs(args.out_dir, exist_ok=True)
    all_scores = []
    done = skipped = 0
    for idx in range(args.start, args.start + args.n):
        s = f"{idx:04d}"
        mul_path = os.path.join(args.mul_trees_dir, f"mul_tree_{s}.nex")
        md_path = os.path.join(args.mul_trees_dir, f"metadata_{s}.json")
        true_path = os.path.join(args.mul_trees_dir, f"species_tree_{s}.nex")
        astral_path = os.path.join(args.mul_trees_dir, "simphy", args.config, s,
                                   f"replicate_{args.replicate}", "astral_species.tre")
        gt_path = os.path.join(args.mul_trees_dir, "simphy", args.config, s,
                               f"replicate_{args.replicate}", "gene_trees.tre")
        if not (os.path.exists(mul_path) and os.path.exists(md_path) and os.path.exists(true_path)):
            skipped += 1
            continue
        true_tree = load_nexus_tree(true_path)
        if args.backbone == "true":
            backbone = true_tree.copy()
        else:
            if not os.path.exists(astral_path):
                skipped += 1
                continue
            backbone = Tree(open(astral_path).read().strip(), format=1)

        gene_trees = all_species = None
        if need_gene_trees:
            if not os.path.exists(gt_path):
                skipped += 1
                continue
            gene_trees = load_gene_trees(gt_path, args.max_gene_trees)
            all_species = set()
            for t in gene_trees:
                all_species.update(t.get_leaf_names())

        # Root the ASTRAL backbone to match the real pipeline (rooting-sensitive metric).
        if args.root_backbone == "hybrid" and args.backbone == "astral":
            backbone = hybrid_root(backbone, gene_trees, args.max_gene_trees)

        with open(md_path) as f:
            md_events = json.load(f).get("events", [])
        try:
            events, scores = two_parent_events(md_events, true_tree, backbone,
                                               mode=args.parents,
                                               gene_trees=gene_trees, all_species=all_species)
            oracle = build_mul_tree_two_parent(backbone, events)
        except Exception as e:
            print(f"  SKIP [{s}]: {type(e).__name__}: {e}")
            skipped += 1
            continue

        for sc in scores:
            sc["sample"] = s
            all_scores.append(sc)
        case_dir = os.path.join(args.out_dir, f"sample_{s}")
        os.makedirs(case_dir, exist_ok=True)
        oracle.write(outfile=os.path.join(case_dir, "output.tre"), format=9)
        with open(mul_path) as f:
            gt = f.read()
        with open(os.path.join(case_dir, "ground_truth.nex"), "w") as f:
            f.write(gt)
        done += 1

    if all_scores:
        with open(os.path.join(args.out_dir, "mapping_scores.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["sample", "type", "source", "target", "A", "B"])
            w.writeheader()
            w.writerows(all_scores)

    n_coclust = sum(1 for sc in all_scores if sc["source"] == "coclust")
    print(f"\nDone: {done} oracle trees ({args.backbone} backbone, {args.parents} parents), "
          f"{skipped} skipped. Allo events via coclust: {n_coclust}.")
    print(f"Score: python scripts/score_reconstructions.py --recon-dir {args.out_dir} --workers 8")


if __name__ == "__main__":
    main()
