"""Diagnose whether clade-level WGD events are mislabeled onto the tips.

Hypothesis (from the 100%-single-species finding): a single WGD on an ancestral
clade edge {A,B,C} is labeled in the training data as separate TIP events
({A},{B},{C}), leaving the ancestral edge WGD-negative — even though that edge
carries the clade-duplication features. We never put the event where it happened.

This measures, on the packaged samples (features + labels, no retraining), how
often an internal edge is:
  - WGD-negative in the labels,
  - has a high clade-duplication feature (frac_clade_duplicated),
  - yet >=2 of its descendant single-species tips ARE WGD-positive.
That pattern is the signature of a clade-level event fragmented onto the tips.

Run in the final_project env. Usage:
  python scripts/diagnose_clade_label_mismatch.py \
      --data-dir data/mul_trees_2k/training/ils_low --max-samples 500
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset
from gene2net_gnn.data.tree_io import reorder_edge_index_preorder
from gene2net_gnn.data.features import edge_clades_species

FRAC_DUP_IDX = 8  # frac_clade_duplicated, per dataset.py edge-feature order


def wgd_positive_edges(labels):
    """Set of preorder edge indices labeled WGD (mappable events only)."""
    pos = set()
    if labels is None or not labels.wgd_edges:
        return pos
    mask = labels.mask or []
    for k, e in enumerate(labels.wgd_edges):
        if k < len(mask) and not mask[k]:
            continue
        pos.add(e)
    return pos


def edge_to_partner_map(labels):
    """wgd_edge -> partner_edge for mappable events."""
    m = {}
    if labels is None or not labels.wgd_edges:
        return m
    mask = labels.mask or []
    for k, e in enumerate(labels.wgd_edges):
        if k < len(mask) and not mask[k]:
            continue
        if k < len(labels.partner_edges):
            m[e] = labels.partner_edges[k]
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--max-samples", type=int, default=500)
    ap.add_argument("--frac-dup-threshold", type=float, default=0.5)
    args = ap.parse_args()

    ds = Gene2NetDataset(args.data_dir)
    n = min(args.max_samples, len(ds))

    size_hist = Counter()        # clade size of every WGD-positive edge
    n_tip_pos = n_internal_pos = 0
    n_sig_edges = 0
    n_sig_edges_inward = 0       # signature edges whose tips point at each other (true fragmentation)
    n_samples_sig = 0
    examined = 0

    for i in range(n):
        s = ds[i]
        if s.labels is None or not s.labels.wgd_edges:
            del s
            continue
        examined += 1

        ei = reorder_edge_index_preorder(s.species_tree_edge_index)
        sp_to_idx = {sp: j for j, sp in enumerate(s.species_list)}
        node_species = [sp_to_idx.get(nm, -1) for nm in s.species_tree_node_names]
        clades = edge_clades_species(ei, s.species_tree_is_leaf, node_species)
        E = len(clades)
        ef = s.species_tree_edge_features

        pos = {e for e in wgd_positive_edges(s.labels) if 0 <= e < E}
        partner_of = edge_to_partner_map(s.labels)

        for e in pos:
            sz = len(clades[e])
            size_hist[sz] += 1
            if sz == 1:
                n_tip_pos += 1
            else:
                n_internal_pos += 1

        # species index -> its tip edge (clade size 1)
        sp_to_tip = {}
        for e in range(E):
            if len(clades[e]) == 1:
                sp_to_tip[next(iter(clades[e]))] = e

        sample_sig = False
        for e in range(E):
            if len(clades[e]) < 2 or e in pos:
                continue  # want internal, WGD-negative edges
            fdup = float(ef[e, FRAC_DUP_IDX]) if (ef is not None and e < ef.shape[0]) else 0.0
            if fdup < args.frac_dup_threshold:
                continue
            desc_tips = {sp_to_tip[sp] for sp in clades[e] if sp in sp_to_tip}
            desc_pos_tips = [t for t in desc_tips if t in pos]
            if len(desc_pos_tips) >= 2:
                n_sig_edges += 1
                sample_sig = True
                # Do the positive tips point at each other (partner inside the clade)?
                inward = sum(1 for t in desc_pos_tips if partner_of.get(t) in desc_tips)
                if inward >= 2 and inward > len(desc_pos_tips) / 2:
                    n_sig_edges_inward += 1
        if sample_sig:
            n_samples_sig += 1
        del s

    print(f"\nExamined {examined} samples with WGD events from {args.data_dir}\n")
    print("WGD-positive edges by clade size:")
    for sz in sorted(size_hist):
        print(f"  size {sz}: {size_hist[sz]}")
    print(f"\nWGD-positive tips (size 1):  {n_tip_pos}")
    print(f"WGD-positive internal (>=2): {n_internal_pos}")
    print("\nClade-fragmentation signature")
    print("  (internal edge: WGD-negative, frac_clade_duplicated >= "
          f"{args.frac_dup_threshold}, >=2 descendant tips WGD-positive)")
    print(f"  signature edges:            {n_sig_edges}")
    print(f"  ...of which 'inward' (tips point at each other = true fragmentation): "
          f"{n_sig_edges_inward} ({100 * n_sig_edges_inward / max(n_sig_edges, 1):.1f}%)")
    print(f"  samples with >=1 signature: {n_samples_sig}/{examined} "
          f"({100 * n_samples_sig / max(examined, 1):.1f}%)")
    print("\nReading:")
    print("  high 'inward' share -> genuine clade-level events fragmented onto reciprocal tips -> label fix is high-leverage")
    print("  low 'inward' share  -> mostly independent co-located tip events -> fragmentation is rarer than the raw signature suggests")


if __name__ == "__main__":
    main()
