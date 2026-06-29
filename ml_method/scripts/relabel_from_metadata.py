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

from gene2net_gnn.data.metadata_labels import (
    labels_from_metadata_for_sample,
    sample_edge_bipartitions,
)


def _sample_duplicated_species(sample_dict, bip):
    """Species under the sample's existing WGD-positive edges (mappable events)."""
    labels = sample_dict.get("labels")
    if labels is None or not labels.wgd_edges:
        return set()
    bip_map = dict(bip)
    mask = labels.mask or []
    out = set()
    for k, e in enumerate(labels.wgd_edges):
        if k < len(mask) and not mask[k]:
            continue
        out |= set(bip_map.get(e, frozenset()))
    return out


def relabel_one_sample(sample_dir, metadata, *, dry_run=False):
    with open(os.path.join(sample_dir, "sample.pkl"), "rb") as f:
        sample_dict = pickle.load(f)

    bip = sample_edge_bipartitions(sample_dict)

    # Index-alignment guard 1: species count must match the metadata.
    n_species_md = metadata.get("n_species")
    if n_species_md is not None and len(sample_dict["species_list"]) != n_species_md:
        raise ValueError(
            f"species-count mismatch in {sample_dir}: sample "
            f"{len(sample_dict['species_list'])} != metadata n_species {n_species_md} "
            "(wrong metadata index?)"
        )

    # Index-alignment guard 2: the sample's labeled polyploids must be a SUBSET of
    # the metadata's true polyploids. The old pipeline may have dropped some
    # (unmappable events), so subset (not equality) is expected; a species in the
    # sample but NOT in metadata means the indices don't correspond.
    md_poly = set((metadata.get("polyploid_species") or {}).keys())
    sample_poly = _sample_duplicated_species(sample_dict, bip)
    if sample_poly and not sample_poly.issubset(md_poly):
        raise ValueError(
            f"polyploid mismatch in {sample_dir}: sample has "
            f"{sorted(sample_poly - md_poly)} not in metadata {sorted(md_poly)} "
            "(wrong metadata index?)"
        )

    labels = labels_from_metadata_for_sample(metadata["events"], sample_dict)

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
    args = ap.parse_args()

    train_dir = os.path.join(args.data_root, args.training_subdir)
    sample_dirs = sorted(
        d for d in os.listdir(train_dir) if d.startswith("sample_")
        and os.path.isdir(os.path.join(train_dir, d))
    )

    n_done = n_unmappable = total_events = 0
    for name in sample_dirs:
        idx = name.replace("sample_", "")
        md_path = os.path.join(args.data_root, f"metadata_{idx}.json")
        if not os.path.exists(md_path):
            raise FileNotFoundError(f"missing metadata for {name}: {md_path}")
        with open(md_path) as f:
            metadata = json.load(f)
        labels = relabel_one_sample(os.path.join(train_dir, name), metadata, dry_run=args.dry_run)
        n_done += 1
        total_events += len(labels.wgd_edges)
        n_unmappable += labels.n_unmappable

    print(f"Relabeled {n_done} samples in {train_dir}")
    print(f"  total events: {total_events}, unmappable (clade not in ASTRAL): {n_unmappable} "
          f"({100 * n_unmappable / max(total_events, 1):.1f}%)")


if __name__ == "__main__":
    main()
