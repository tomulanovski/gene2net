"""Emit the benchmark (4a) and fractionation (4e) head-to-head tables as markdown,
per config, on the mu-distance headline metric.

It recomputes the GNN-vs-competitor join using head_to_head.py's tested loaders,
rather than reading the h2h_*.csv the jobs wrote, so the tables are authoritative
and independent of any network exclusion baked into those files.

NO networks are excluded. Bendiksby_2011 / Marcussen_2011 were only ever wrong-rate
in conf_dup_loss_high_10M (Ne 200k), which is NOT in the method benchmark, and that
config is fixed now anyway (check_dup_loss_rates.sh all OK, 2026-08-19). The old
DEFAULT_EXCLUDE in head_to_head therefore dropped two good networks from every
config; here we pass all 21.

For each config it shows both GNN decodes as separate rows -- ploidy-informed
(bound_driven) and ploidy-free (detect_only) -- against the competitors, which are
identical across the two GNN decodes. Values are the mean over the benchmark
networks each method scored; `n` is how many networks that method has a finite
mu-distance on (mu-unscorable networks, mu_scored=0, do not count).

Run from ml_method/ in the gene2net env (needs the competitor comparisons_raw.csv):
    python scripts/make_results_tables.py
    python scripts/make_results_tables.py --base output/reconstruct_final \
        --out docs/generated_tables_mu.md
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from head_to_head import load_gnn_scores, load_competitor_scores  # noqa: E402

METRICS = ["mu_distance", "ret_leaf_jaccard", "ret_sisters_jaccard"]

# Config groups for the two results sections. ne2M is included in the benchmark
# group; it is dropped automatically (with a printed note) if not yet scored.
BENCHMARK_CONFIGS = [
    "conf_ils_low_10M", "conf_ils_medium_10M", "conf_ils_high_10M",
    "conf_dup_loss_low_10M_ne1M", "conf_dup_loss_medium_10M_ne1M", "conf_dup_loss_high_10M_ne1M",
    "conf_dup_loss_low_10M_ne2M", "conf_dup_loss_medium_10M_ne2M", "conf_dup_loss_high_10M_ne2M",
]
FRACTIONATION_CONFIGS = [
    "conf_dup_loss_medium_10M_ne1M_fix025",
    "conf_dup_loss_medium_10M_ne1M_fix050",
    "conf_dup_loss_medium_10M_ne1M_fix075",
]

STRAT_LABEL = {
    "bound_driven": "GNN (ploidy-informed)",
    "detect_only": "GNN (ploidy-free)",
}
# head_to_head canonical competitor names -> thesis display names.
COMPETITOR_LABEL = {
    "polyphest": "Polyphest",
    "grampa_iter": "GRAMPA-iter",
    "grampa_iter_prior": "GRAMPA-iter + prior",
    "grampa": "GRAMPA",
}
ROW_ORDER = [
    "GNN (ploidy-informed)", "GNN (ploidy-free)",
    "Polyphest", "GRAMPA-iter", "GRAMPA-iter + prior", "GRAMPA",
]

CONFIG_TITLE = {
    "conf_ils_low_10M": "ILS low",
    "conf_ils_medium_10M": "ILS medium",
    "conf_ils_high_10M": "ILS high",
    "conf_dup_loss_low_10M_ne1M": "dup/loss low (Ne 1M)",
    "conf_dup_loss_medium_10M_ne1M": "dup/loss medium (Ne 1M)",
    "conf_dup_loss_high_10M_ne1M": "dup/loss high (Ne 1M)",
    "conf_dup_loss_low_10M_ne2M": "dup/loss low (Ne 2M)",
    "conf_dup_loss_medium_10M_ne2M": "dup/loss medium (Ne 2M)",
    "conf_dup_loss_high_10M_ne2M": "dup/loss high (Ne 2M)",
    "conf_dup_loss_medium_10M_ne1M_fix025": "retention 0.25 (severe)",
    "conf_dup_loss_medium_10M_ne1M_fix050": "retention 0.50 (medium)",
    "conf_dup_loss_medium_10M_ne1M_fix075": "retention 0.75 (mild)",
}


def agg_csv(base, config, strat):
    return os.path.join(base, "final", config, "agg", strat, "reconstruction_scores.csv")


def per_method_means(combined):
    """combined: tidy [network, config, method, metric, value] -> per-method table
    with one column per metric plus n (networks with a finite value)."""
    means = (combined.groupby(["method", "metric"])["value"].mean()
                     .unstack("metric"))
    n = (combined.dropna(subset=["value"]).groupby("method")["network"].nunique())
    means["n"] = n
    return means


def build_config_table(base, config, competitor_path):
    """Load competitors once and both GNN decodes, return the per-method mean table,
    or None if the competitor side is missing (printed loudly)."""
    try:
        comp = load_competitor_scores(competitor_path, config, METRICS)
    except SystemExit as exc:
        print(f"  SKIP {config}: competitor scores unavailable ({exc})")
        return None
    comp = comp.copy()
    comp["method"] = comp["method"].map(lambda m: COMPETITOR_LABEL.get(m, m))

    frames = [comp]
    for strat, label in STRAT_LABEL.items():
        csv = agg_csv(base, config, strat)
        if not os.path.exists(csv):
            print(f"  NOTE {config}: GNN {strat} not scored yet ({csv})")
            continue
        g = load_gnn_scores(csv, config, METRICS)
        g = g.copy()
        g["method"] = label
        frames.append(g)

    combined = pd.concat(frames, ignore_index=True)
    return per_method_means(combined)


def table_md(means, title):
    def f(x):
        return f"{x:.3f}" if pd.notna(x) else "-"
    lines = [
        f"**{title}**", "",
        "| method | mu | ret_leaf | ret_sisters | n |",
        "| --- | --- | --- | --- | --- |",
    ]
    for m in ROW_ORDER:
        if m not in means.index:
            continue
        r = means.loc[m]
        n = int(r["n"]) if pd.notna(r.get("n")) else 0
        lines.append(
            f"| {m} | {f(r.get('mu_distance'))} | {f(r.get('ret_leaf_jaccard'))} "
            f"| {f(r.get('ret_sisters_jaccard'))} | {n} |"
        )
    lines.append("")
    return "\n".join(lines)


def section(base, configs, heading, competitor_path):
    out = [f"## {heading}", ""]
    for config in configs:
        means = build_config_table(base, config, competitor_path)
        if means is None:
            continue
        out.append(table_md(means, CONFIG_TITLE.get(config, config)))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default="output/reconstruct_final",
                    help="GNN output base holding final/<config>/agg/<strat>/reconstruction_scores.csv")
    ap.add_argument("--competitor-scores", default=None,
                    help="Competitor comparisons_raw.csv or a dir under which it lives. "
                         "Default: simulations/analysis/summary/<config>/comparisons_raw.csv")
    ap.add_argument("--out", default="docs/generated_tables_mu.md",
                    help="Markdown output file.")
    args = ap.parse_args()

    parts = [
        "# Results tables (mu-distance)", "",
        "Lower is better for every column. `mu` is the normalized mu-distance "
        "(headline overall-quality metric). ret_leaf and ret_sisters are the "
        "reticulation Jaccard distances. `n` is the number of benchmark networks "
        "the method has a finite mu-distance on (out of 21). All 21 networks are "
        "used; no exclusions.", "",
        section(args.base, BENCHMARK_CONFIGS, "Benchmark (4a)", args.competitor_scores),
        section(args.base, FRACTIONATION_CONFIGS, "Fractionation (4e)", args.competitor_scores),
    ]
    doc = "\n".join(parts)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(doc)
    print(doc)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
