"""Common-subset comparison for the method chapter.

The per-method means elsewhere are each over the networks that method completed, so the
GNN is averaged over 21 while Polyphest is averaged over the 17 to 20 it finishes. Those
subsets differ, and the difference flatters whichever method completes less, because the
networks it skips are the harder ones. This script reports the unbiased view instead: for
each configuration it restricts every method to the networks ALL of them completed, and
reports the means there, on all three measures.

It also reports what is lost by that restriction, namely the GNN's score on the networks
outside the common subset, so the chapter can say how the method does on the networks the
competitors never returned.

Reuses the loaders in head_to_head.py so the competitor side is parsed and the Polyphest
percentiles merged exactly as they are everywhere else. Unlike head_to_head's CLI it does
NOT drop Bendiksby_2011 and Marcussen_2011: that default predates the re-simulation that
fixed their duplication and loss rates, so excluding them now discards good data.

Run from ml_method/:
    python scripts/common_subset_tables.py
    python scripts/common_subset_tables.py --strategy detect_only --metric ret_leaf_jaccard
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from head_to_head import load_gnn_scores, load_competitor_scores  # noqa: E402

CONFIGS = [
    "conf_ils_low_10M", "conf_ils_medium_10M", "conf_ils_high_10M",
    "conf_dup_loss_low_10M", "conf_dup_loss_medium_10M", "conf_dup_loss_high_10M",
    "conf_dup_loss_low_10M_ne1M", "conf_dup_loss_medium_10M_ne1M", "conf_dup_loss_high_10M_ne1M",
    "conf_dup_loss_low_10M_ne2M", "conf_dup_loss_medium_10M_ne2M", "conf_dup_loss_high_10M_ne2M",
    "conf_dup_loss_medium_10M_ne1M_fix025", "conf_dup_loss_medium_10M_ne1M_fix050",
    "conf_dup_loss_medium_10M_ne1M_fix075",
]
METRICS = ["mu_distance", "ret_leaf_jaccard", "ret_sisters_jaccard"]
# Methods that must all have a value for a network to count as common.
COMMON_OVER = ["gnn", "polyphest", "grampa_iter", "grampa_iter_prior"]
SHOW = ["gnn", "polyphest", "grampa_iter", "grampa_iter_prior"]


def per_network(config, strategy, base):
    """Wide frame: index network, columns method, one column block per metric."""
    gnn_dir = os.path.join(base, "final", config, "agg", strategy)
    if not os.path.isdir(gnn_dir):
        raise SystemExit(f"missing GNN scores directory: {gnn_dir}")
    gnn = load_gnn_scores(gnn_dir, config, METRICS)
    comp = load_competitor_scores(None, config, METRICS)
    return pd.concat([gnn, comp], ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strategy", default="bound_driven", choices=["bound_driven", "detect_only"])
    ap.add_argument("--base", default="output/reconstruct_final")
    ap.add_argument("--out", default="docs/common_subset_tables.md")
    args = ap.parse_args()

    lines = [f"# Common-subset comparison ({args.strategy})", "",
             "Every method restricted to the networks all of them completed. `n` is the size of",
             "that common subset. `gnn_extra` is the GNN mean on the networks outside it, with",
             "the count in `n_extra`, which are the networks at least one competitor failed on.",
             ""]

    for metric in METRICS:
        lines += [f"## {metric}", "",
                  "| Configuration | " + " | ".join(SHOW) + " | n | gnn_extra | n_extra |",
                  "| --- " * (len(SHOW) + 4) + "|"]
        for config in CONFIGS:
            tidy = per_network(config, args.strategy, args.base)
            pm = tidy[tidy["metric"] == metric]
            if pm.empty:
                raise SystemExit(f"no rows for metric {metric} in {config}")
            wide = pm.pivot_table(index="network", columns="method", values="value", aggfunc="mean")
            for m in COMMON_OVER:
                if m not in wide.columns:
                    raise SystemExit(f"{config}: method {m!r} absent; cannot form a common subset")
            common = wide[COMMON_OVER].notna().all(axis=1)
            cells = []
            for m in SHOW:
                vals = wide.loc[common, m].dropna() if m in wide.columns else pd.Series(dtype=float)
                cells.append(f"{vals.mean():.3f}" if len(vals) else "-")
            extra = wide.loc[~common, "gnn"].dropna()
            extra_s = f"{extra.mean():.3f}" if len(extra) else "-"
            lines.append(f"| {config} | " + " | ".join(cells)
                         + f" | {int(common.sum())} | {extra_s} | {len(extra)} |")
        lines.append("")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    text = "\n".join(lines)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
