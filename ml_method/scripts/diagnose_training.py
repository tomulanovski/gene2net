"""Diagnose training data quality: label mapping, class balance, feature-label correlation.

Run on cluster:
    python scripts/diagnose_training.py --dir /path/to/mul_trees_2k/training/ils_low

Saves report to the training directory as diagnostics.txt
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetSample, Gene2NetDataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", required=True, help="Training data directory (with sample_NNNN/ folders)")
    parser.add_argument("--max-samples", type=int, default=500, help="Max samples to analyze")
    args = parser.parse_args()

    dataset = Gene2NetDataset(args.dir)
    n_total = min(len(dataset), args.max_samples)
    print(f"Analyzing {n_total} / {len(dataset)} samples...")

    # Collectors
    all_n_edges = []
    all_n_wgd_edges = []       # edges with wgd_count > 0
    all_n_events = []          # total events per sample
    all_n_unmappable = []
    all_n_mappable = []
    all_wgd_fractions = []     # fraction of edges that are WGD per sample

    # Feature-label correlation
    # For each sample, collect (feature_value, is_wgd) pairs for edges
    concordance_wgd = []       # (concordance_factor, is_wgd)
    concordance_no_wgd = []

    # Node features for polyploid vs non-polyploid species
    copy_feat_wgd_nodes = []     # node features adjacent to WGD edges
    copy_feat_normal_nodes = []  # node features adjacent to non-WGD edges

    # WGD count distribution
    wgd_count_dist = {}  # count -> number of edges

    skipped = 0

    for i in range(n_total):
        try:
            sample = dataset[i]
        except Exception as e:
            skipped += 1
            continue

        labels = sample.labels
        if labels is None:
            skipped += 1
            continue

        n_edges = labels.n_edges
        wgd_counts = labels.wgd_counts
        n_unmappable = labels.n_unmappable
        n_events_total = len(labels.mask)
        n_mappable = sum(1 for m in labels.mask if m)

        if isinstance(wgd_counts, list):
            wgd_counts_arr = np.array(wgd_counts)
        else:
            wgd_counts_arr = np.array(wgd_counts)

        n_wgd_edges = int((wgd_counts_arr > 0).sum())
        wgd_frac = n_wgd_edges / max(n_edges, 1)

        all_n_edges.append(n_edges)
        all_n_wgd_edges.append(n_wgd_edges)
        all_n_events.append(n_events_total)
        all_n_unmappable.append(n_unmappable)
        all_n_mappable.append(n_mappable)
        all_wgd_fractions.append(wgd_frac)

        # WGD count distribution
        for c in wgd_counts_arr:
            c_int = int(c)
            wgd_count_dist[c_int] = wgd_count_dist.get(c_int, 0) + 1

        # Edge features correlation
        edge_feats = sample.species_tree_edge_features  # [n_edges, 4]
        if edge_feats is not None and edge_feats.shape[0] == n_edges:
            for e_idx in range(n_edges):
                cf = float(edge_feats[e_idx, 0])  # concordance factor
                is_wgd = int(wgd_counts_arr[e_idx]) > 0
                if is_wgd:
                    concordance_wgd.append(cf)
                else:
                    concordance_no_wgd.append(cf)

        # Node features: check if nodes adjacent to WGD edges have different copy patterns
        node_feats = sample.species_tree_node_features  # [n_nodes, 13]
        edge_index = sample.species_tree_edge_index     # [2, 2*n_edges]
        if node_feats is not None and edge_index is not None:
            wgd_edge_set = set()
            for e_idx in range(n_edges):
                if wgd_counts_arr[e_idx] > 0:
                    wgd_edge_set.add(e_idx)

            # Each undirected edge i corresponds to directed edges 2*i and 2*i+1
            for e_idx in range(n_edges):
                child_node = int(edge_index[1, 2 * e_idx])  # child in parent->child direction
                feat = node_feats[child_node].numpy()
                if e_idx in wgd_edge_set:
                    copy_feat_wgd_nodes.append(feat)
                else:
                    copy_feat_normal_nodes.append(feat)

        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{n_total}")

    # =========================================================================
    # BUILD REPORT
    # =========================================================================
    lines = []

    def p(msg=""):
        lines.append(msg)

    analyzed = n_total - skipped
    p(f"{'='*70}")
    p(f"TRAINING DATA DIAGNOSTICS ({analyzed} samples analyzed)")
    p(f"{'='*70}")

    # --- Label mapping quality ---
    p(f"\n{'='*70}")
    p(f"1. LABEL MAPPING QUALITY")
    p(f"{'='*70}")
    total_events = sum(all_n_events)
    total_unmappable = sum(all_n_unmappable)
    total_mappable = sum(all_n_mappable)
    p(f"  Total WGD events across all samples: {total_events}")
    p(f"  Mappable (Jaccard >= 0.5):   {total_mappable} ({100*total_mappable/max(total_events,1):.1f}%)")
    p(f"  Unmappable (Jaccard < 0.5):  {total_unmappable} ({100*total_unmappable/max(total_events,1):.1f}%)")

    samples_with_unmappable = sum(1 for u in all_n_unmappable if u > 0)
    p(f"  Samples with >=1 unmappable event: {samples_with_unmappable}/{analyzed} ({100*samples_with_unmappable/max(analyzed,1):.1f}%)")

    zero_event_samples = sum(1 for e in all_n_events if e == 0)
    p(f"  Samples with 0 events (negatives): {zero_event_samples}/{analyzed} ({100*zero_event_samples/max(analyzed,1):.1f}%)")

    all_wgd_zero = sum(1 for w in all_n_wgd_edges if w == 0)
    p(f"  Samples with 0 labeled WGD edges (after mapping): {all_wgd_zero}/{analyzed} ({100*all_wgd_zero/max(analyzed,1):.1f}%)")

    # --- Class balance ---
    p(f"\n{'='*70}")
    p(f"2. EDGE-LEVEL CLASS DISTRIBUTION")
    p(f"{'='*70}")
    total_edge_count = sum(wgd_count_dist.values())
    for c in sorted(wgd_count_dist.keys()):
        cnt = wgd_count_dist[c]
        pct = 100 * cnt / max(total_edge_count, 1)
        bar = "#" * int(pct / 2)
        p(f"  Class {c} (={c} WGD events): {cnt:>8} edges ({pct:5.1f}%) {bar}")

    total_wgd_edges = sum(cnt for c, cnt in wgd_count_dist.items() if c > 0)
    total_no_wgd = wgd_count_dist.get(0, 0)
    ratio = total_no_wgd / max(total_wgd_edges, 1)
    p(f"\n  Non-WGD : WGD ratio = {ratio:.0f} : 1")

    p(f"\n  WGD edges per sample: min={min(all_n_wgd_edges)}, max={max(all_n_wgd_edges)}, "
      f"mean={np.mean(all_n_wgd_edges):.1f}")
    p(f"  WGD fraction per sample: min={min(all_wgd_fractions):.3f}, max={max(all_wgd_fractions):.3f}, "
      f"mean={np.mean(all_wgd_fractions):.3f}")

    # --- Feature correlation ---
    p(f"\n{'='*70}")
    p(f"3. CONCORDANCE FACTOR vs WGD")
    p(f"{'='*70}")
    if concordance_wgd and concordance_no_wgd:
        p(f"  WGD edges     — concordance: mean={np.mean(concordance_wgd):.4f}, "
          f"std={np.std(concordance_wgd):.4f}, median={np.median(concordance_wgd):.4f}")
        p(f"  Non-WGD edges — concordance: mean={np.mean(concordance_no_wgd):.4f}, "
          f"std={np.std(concordance_no_wgd):.4f}, median={np.median(concordance_no_wgd):.4f}")
        diff = abs(np.mean(concordance_wgd) - np.mean(concordance_no_wgd))
        p(f"  Difference in means: {diff:.4f}")
        if diff < 0.02:
            p(f"  >>> WARNING: Concordance shows NO meaningful difference between WGD and non-WGD edges")
        elif diff < 0.05:
            p(f"  >>> WEAK signal: small difference in concordance")
        else:
            p(f"  >>> GOOD signal: meaningful difference in concordance")
    else:
        p(f"  Insufficient data for concordance analysis")

    # --- Node features correlation ---
    p(f"\n{'='*70}")
    p(f"4. NODE FEATURES (child of edge) vs WGD")
    p(f"{'='*70}")
    feat_names = ["mean_copies", "var_copies", "mode_copies", "p_absent",
                  "p_1_copy", "p_2_copies", "p_3plus_copies", "max_copies",
                  "clust_mean", "clust_std", "clust_max", "clust_min", "clust_median"]
    if copy_feat_wgd_nodes and copy_feat_normal_nodes:
        wgd_arr = np.array(copy_feat_wgd_nodes)
        normal_arr = np.array(copy_feat_normal_nodes)
        p(f"  {'Feature':<16} | {'WGD mean':>10} | {'Non-WGD mean':>12} | {'Diff':>8} | Signal")
        p(f"  {'-'*16}-+-{'-'*10}-+-{'-'*12}-+-{'-'*8}-+-------")
        for j, fname in enumerate(feat_names):
            w_mean = np.mean(wgd_arr[:, j])
            n_mean = np.mean(normal_arr[:, j])
            diff = abs(w_mean - n_mean)
            # Relative signal strength
            denom = max(abs(n_mean), 0.001)
            rel_diff = diff / denom
            signal = "STRONG" if rel_diff > 0.3 else ("weak" if rel_diff > 0.1 else "none")
            p(f"  {fname:<16} | {w_mean:>10.4f} | {n_mean:>12.4f} | {diff:>8.4f} | {signal}")
    else:
        p(f"  Insufficient data for node feature analysis")

    # --- Summary ---
    p(f"\n{'='*70}")
    p(f"5. SUMMARY & RECOMMENDATIONS")
    p(f"{'='*70}")
    if total_events > 0:
        unmappable_pct = 100 * total_unmappable / total_events
        if unmappable_pct > 30:
            p(f"  [CRITICAL] {unmappable_pct:.0f}% of events are unmappable — ASTRAL topology mismatch is severe")
            p(f"             Consider: node-level prediction (polyploid species) instead of edge-level")
        elif unmappable_pct > 15:
            p(f"  [WARNING] {unmappable_pct:.0f}% of events are unmappable — significant topology mismatch")
        else:
            p(f"  [OK] Only {unmappable_pct:.0f}% unmappable — topology mapping is reasonable")

    if ratio > 20:
        p(f"  [CRITICAL] Class ratio {ratio:.0f}:1 — extreme imbalance, current weights may be insufficient")
    elif ratio > 10:
        p(f"  [WARNING] Class ratio {ratio:.0f}:1 — consider stronger class weights or oversampling")

    # Print and save
    output = "\n".join(lines)
    print(output)

    out_path = os.path.join(args.dir, "diagnostics.txt")
    with open(out_path, "w") as f:
        f.write(output + "\n")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
