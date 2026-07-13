"""Rewrite training labels from ground-truth metadata events (clade-level).

Per sample: load metadata events, map them onto the sample's OWN edge
bipartitions, and write labels_clade.pkl (non-destructive). Hard-fails on any
index/dimension mismatch so a misaligned sample can never be silently mislabeled.

Usage:
    python scripts/relabel_from_metadata.py --training-subdir training/ils_low
    python scripts/relabel_from_metadata.py --training-subdir training/ils_low --dry-run
"""
import argparse
import json
import os
import pickle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ete3 import Tree

from gene2net_gnn.data.metadata_labels import (
    labels_from_metadata_for_sample, two_parent_labels_from_metadata,
)


_TREE_CACHE = {}


def load_nexus_tree(path):
    # Cache by file content: the ~21 networks repeat across the 2k samples, so this
    # turns thousands of (slow) ete3 parses into ~21. The tree is only read, never
    # mutated, so sharing the object across samples is safe.
    content = open(path).read()
    cached = _TREE_CACHE.get(content)
    if cached is not None:
        return cached
    tree = None
    for line in content.splitlines():
        s = line.strip()
        if s.lower().startswith("tree") and "=" in s:
            tree = Tree(s.split("=", 1)[1].strip(), format=1)
            break
    if tree is None:
        tree = Tree(content.strip(), format=1)
    _TREE_CACHE[content] = tree
    return tree


def _sample_polyploids_from_copies(sample_dict, min_mode=2):
    """Species whose gene-tree MODE copy count >= min_mode.

    Derived from the data (node copy-count features), not the ASTRAL-mapped labels,
    so it is not corrupted by label mis-mapping. node feature col 2 = mode_copies.
    """
    names = sample_dict["species_tree_node_names"]
    is_leaf = sample_dict["species_tree_is_leaf"]
    feats = sample_dict["species_tree_node_features"]
    species = set(sample_dict["species_list"])
    out = set()
    for i, nm in enumerate(names):
        if bool(is_leaf[i]) and nm in species and float(feats[i, 2]) >= min_mode:
            out.add(nm)
    return out


def relabel_one_sample(sample_dir, metadata, *, dry_run=False, true_tree=None, two_parent=False):
    with open(os.path.join(sample_dir, "sample.pkl"), "rb") as f:
        sample_dict = pickle.load(f)

    # Index-alignment guard 1: species count must match the metadata.
    n_species_md = metadata.get("n_species")
    if n_species_md is not None and len(sample_dict["species_list"]) != n_species_md:
        raise ValueError(
            f"species-count mismatch in {sample_dir}: sample "
            f"{len(sample_dict['species_list'])} != metadata n_species {n_species_md} "
            "(wrong metadata index?)"
        )

    # Index-alignment guard 2: the sample's OBSERVED polyploids (gene-tree copy
    # counts) must be a subset of the metadata's true polyploids. Subset (not
    # equality) is expected because dup/loss can make a true polyploid look diploid
    # in most gene trees; a copy-count polyploid NOT in metadata means the indices
    # don't correspond.
    md_poly = set((metadata.get("polyploid_species") or {}).keys())
    sample_poly = _sample_polyploids_from_copies(sample_dict)
    if sample_poly and not sample_poly.issubset(md_poly):
        raise ValueError(
            f"polyploid mismatch in {sample_dir}: sample copy-count polyploids "
            f"{sorted(sample_poly - md_poly)} not in metadata {sorted(md_poly)} "
            "(wrong metadata index?)"
        )

    if two_parent:
        if true_tree is None:
            raise ValueError("--two-parent needs the true tree (home parent = true-tree sibling)")
        labels = two_parent_labels_from_metadata(metadata["events"], sample_dict, true_tree)
    else:
        labels = labels_from_metadata_for_sample(metadata["events"], sample_dict, true_tree=true_tree)

    # Dimension guard: label edge-count must equal the feature rows.
    n_feat_edges = sample_dict["species_tree_edge_features"].shape[0]
    if labels.n_edges != n_feat_edges:
        raise ValueError(
            f"edge-count mismatch in {sample_dir}: labels {labels.n_edges} != features {n_feat_edges}"
        )

    if not dry_run:
        with open(os.path.join(sample_dir, "labels_clade.pkl"), "wb") as f:
            pickle.dump(labels, f)
    return labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--training-subdir", required=True, help="e.g. training/ils_low")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--away-parent", action="store_true",
                    help="Retarget each allo partner to the true parent NOT at X's ASTRAL "
                         "home (fixes the ~55%% partner==home labelling bug). Needs "
                         "species_tree_<idx>.nex.")
    ap.add_argument("--two-parent", action="store_true",
                    help="Write two-parent labels: partner_edges (B=metadata partner) + "
                         "home_edges (A=target's true-tree sibling). Needs species_tree_<idx>.nex.")
    args = ap.parse_args()

    train_dir = os.path.join(args.data_root, args.training_subdir)
    sample_dirs = sorted(
        d for d in os.listdir(train_dir) if d.startswith("sample_")
        and os.path.isdir(os.path.join(train_dir, d))
    )

    n_done = n_unmappable = total_events = 0
    errors = []
    for name in sample_dirs:
        idx = name.replace("sample_", "")
        md_path = os.path.join(args.data_root, f"metadata_{idx}.json")
        if not os.path.exists(md_path):
            raise FileNotFoundError(f"missing metadata for {name}: {md_path}")
        with open(md_path) as f:
            metadata = json.load(f)
        true_tree = None
        if args.away_parent or args.two_parent:
            tt_path = os.path.join(args.data_root, f"species_tree_{idx}.nex")
            if not os.path.exists(tt_path):
                raise FileNotFoundError(f"--away-parent/--two-parent needs the true tree: {tt_path}")
            true_tree = load_nexus_tree(tt_path)
        try:
            labels = relabel_one_sample(os.path.join(train_dir, name), metadata,
                                        dry_run=args.dry_run, true_tree=true_tree,
                                        two_parent=args.two_parent)
        except ValueError as e:
            errors.append(str(e))
            continue
        n_done += 1
        total_events += len(labels.wgd_edges)
        n_unmappable += labels.n_unmappable

    print(f"Relabeled {n_done} samples in {train_dir}"
          + (" (dry-run, nothing written)" if args.dry_run else ""))
    print(f"  total events: {total_events}, unmappable (clade not in ASTRAL): {n_unmappable} "
          f"({100 * n_unmappable / max(total_events, 1):.1f}%)")

    if errors:
        print(f"\n{len(errors)} samples FAILED consistency checks (NOT relabeled):")
        for e in errors[:10]:
            print("  " + e)
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
