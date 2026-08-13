"""Join the GNN benchmark scores with the competitor (Polyphest / GRAMPA /
GRAMPA-iter) scores into ONE head-to-head table, per config, on the same metric.

The GNN benchmark (scripts/benchmark_networks.py -> scripts/score_reconstructions.py)
emits ONLY the GNN's own scores, one CSV per build strategy:
    <out-base>/<strategy>/reconstruction_scores.csv
with columns:
    sample, edit_distance_multree, rf_distance, num_rets_diff,
    ret_leaf_jaccard, ret_sisters_jaccard, ploidy_diff
Here `sample` IS the network name (benchmark_networks.py writes one dir per
network), and there is no `config` column (the config is fixed by the run).

The competitor numbers are produced on the `simulations/` side, NOT by anything
in ml_method. The real competitor score file is (CONFIRMED to exist locally):
    simulations/analysis/summary/<config>/comparisons_raw.csv
a LONG-format table with columns:
    network, config, method, replicate, metric, value, jaccard_mode, status
and methods:
    grampa, grandma_split, mpsugar, padre, alloppnet,
    polyphest_p50, polyphest_p70, polyphest_p90
Confirmed method meanings (from simulations/scripts/create_polyphest_vs_grampaiter.py,
METHOD_DISPLAY):
    grampa         -> GRAMPA (single pass)
    grandma_split  -> GRAMPA-iter  (displayed GRAMPA^Iter)
    polyphest_p50/p70/p90 -> Polyphest at the 50/70/90 percentile cutoff; these
                             are merged into one `polyphest` the same way the
                             simulations figures do: per (network, replicate)
                             take the LOWEST percentile that produced output
                             (see simulations/scripts/create_analysis_figures.py
                             merge_polyphest_comparisons).

Metric-name reconciliation (CONFIRMED): the GNN CSV stores the Jaccard/ploidy
*distance* under the bare column name (score_reconstructions.py extracts
comp[metric]['dist']), whereas the competitor long table stores it as
'<metric>.dist' (e.g. ret_leaf_jaccard.dist) alongside .TP/.FP/.FN rows. Scalar
metrics (edit_distance_multree, rf_distance, num_rets_diff) share the same name
on both sides. The `jaccard_mode` column (partial vs standard) tags which
scoring convention was used; in the observed data each (network, method,
replicate, metric) appears exactly once, so no de-duplication across modes is
needed, but we keep only status == SUCCESS rows.

Join key: NETWORK (config is fixed by --config; the config column is used only to
filter, so a merged multi-config competitor file also works). Each output row is
one network with the GNN and each competitor's value side by side.

Read-only w.r.t. inputs; writes only the combined CSV given by --out.

Usage:
    python scripts/head_to_head.py \
        --gnn-scores output/reconstruct_allo/benchmark/conf_ils_low_10M/collapse \
        --config conf_ils_low_10M \
        --out output/reconstruct_allo/benchmark/conf_ils_low_10M/head_to_head.csv
    # --competitor-scores defaults to
    #   ../simulations/analysis/summary/<config>/comparisons_raw.csv
"""
import argparse
import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Constants describing the competitor schema (see module docstring).
# ---------------------------------------------------------------------------
POLYPHEST_THRESHOLDS = ["polyphest_p50", "polyphest_p70", "polyphest_p90"]

# competitor method label -> canonical name used in the output columns.
# grandma_split is the iterative GRAMPA variant (GRAMPA^Iter in the figures).
COMPETITOR_METHOD_MAP = {
    "grampa": "grampa",
    "grandma_split": "grampa_iter",
    "grandma_split_prior": "grampa_iter_prior",  # GRAMPA-iter given the inferred ploidy prior
    "polyphest": "polyphest",  # after the p50/p70/p90 merge below
}

# Order of methods in the output (gnn first).
OUTPUT_METHODS = ["gnn", "polyphest", "grampa", "grampa_iter", "grampa_iter_prior"]

# The two wrong-rate networks the clean-config restriction drops by default.
DEFAULT_EXCLUDE = ["Bendiksby_2011", "Marcussen_2011"]

# Metrics stored as '<metric>.dist' in the competitor long table but under the
# bare name in the GNN CSV. Used only to translate the competitor side.
JACCARD_METRICS = {"ret_leaf_jaccard", "ret_sisters_jaccard", "ploidy_diff"}

# Extra metrics to carry alongside the primary --metric when both sides have
# them (kept short and thesis-relevant).
EXTRA_METRICS = ["ret_leaf_jaccard", "ret_sisters_jaccard"]

# Friendly short suffix for the primary edit-distance metric in column names.
METRIC_SHORT = {"edit_distance_multree": "edit"}


def short_name(metric):
    return METRIC_SHORT.get(metric, metric)


def competitor_metric_name(metric):
    """GNN column name -> competitor long-table metric name."""
    return f"{metric}.dist" if metric in JACCARD_METRICS else metric


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
def _repo_root():
    # this file lives at <repo>/ml_method/scripts/head_to_head.py
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def default_competitor_path(config):
    return os.path.join(_repo_root(), "simulations", "analysis", "summary",
                        config, "comparisons_raw.csv")


def resolve_gnn_path(path):
    """Accept the reconstruction_scores.csv file OR a dir containing it."""
    if os.path.isdir(path):
        cand = os.path.join(path, "reconstruction_scores.csv")
        if os.path.exists(cand):
            return cand
        raise SystemExit(
            f"--gnn-scores is a directory but has no reconstruction_scores.csv: {path}\n"
            f"  Point it at the strategy dir (e.g. .../benchmark/<config>/collapse) "
            f"or directly at the CSV."
        )
    return path


def resolve_competitor_path(path, config):
    """Accept the comparisons_raw.csv file, or a dir under which it lives."""
    if path and os.path.isfile(path):
        return path
    if path and os.path.isdir(path):
        for cand in (
            os.path.join(path, "comparisons_raw.csv"),
            os.path.join(path, config, "comparisons_raw.csv"),
            os.path.join(path, "summary", config, "comparisons_raw.csv"),
        ):
            if os.path.exists(cand):
                return cand
    # fall through to the default location
    return path or default_competitor_path(config)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_gnn_scores(path, config, metrics):
    """Load the GNN reconstruction_scores.csv -> tidy [network, method, metric, value].

    method is always 'gnn'. `sample` is renamed to `network`. Only the requested
    metric columns that exist are kept.
    """
    csv_path = resolve_gnn_path(path)
    if not os.path.exists(csv_path):
        raise SystemExit(f"GNN scores not found: {csv_path}")
    df = pd.read_csv(csv_path)
    if "sample" not in df.columns:
        raise SystemExit(
            f"GNN scores {csv_path} has no 'sample' column (found: {list(df.columns)}). "
            f"Expected score_reconstructions.py output."
        )
    df = df.rename(columns={"sample": "network"})

    present = [m for m in metrics if m in df.columns]
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        print(f"  NOTE: GNN scores missing metric column(s): {missing}")
    if not present:
        raise SystemExit(f"None of the requested metrics {metrics} are in {csv_path}")

    tidy = df.melt(id_vars=["network"], value_vars=present,
                   var_name="metric", value_name="value")
    tidy["value"] = pd.to_numeric(tidy["value"], errors="coerce")
    tidy["method"] = "gnn"
    tidy["config"] = config
    return tidy[["network", "config", "method", "metric", "value"]]


def merge_polyphest(long_df):
    """Collapse polyphest_p50/p70/p90 into a single 'polyphest'.

    Mirrors simulations/scripts/create_analysis_figures.py::merge_polyphest_comparisons:
    per (config, network, replicate) keep the LOWEST percentile that has any
    SUCCESS row, and relabel it 'polyphest'. Non-polyphest rows pass through.
    """
    non_poly = long_df[~long_df["method"].isin(POLYPHEST_THRESHOLDS)].copy()
    poly = long_df[long_df["method"].isin(POLYPHEST_THRESHOLDS)].copy()
    if poly.empty:
        return non_poly

    group_cols = [c for c in ("config", "network", "replicate") if c in poly.columns]
    kept = []
    for _, group in poly.groupby(group_cols):
        chosen = None
        for method in POLYPHEST_THRESHOLDS:
            rows = group[group["method"] == method]
            if len(rows) > 0 and (rows["status"] == "SUCCESS").any():
                chosen = method
                break
        if chosen is None:
            continue
        sel = group[group["method"] == chosen].copy()
        sel["method"] = "polyphest"
        kept.append(sel)

    if kept:
        return pd.concat([non_poly] + kept, ignore_index=True)
    return non_poly


def load_competitor_scores(path, config, metrics):
    """Load competitor comparisons_raw.csv -> tidy [network, method, metric, value].

    Filters to the config and status == SUCCESS, merges polyphest percentiles,
    renames methods to canonical names, and translates each requested GNN metric
    to its competitor '<metric>.dist' form where applicable. `value` is averaged
    over replicates per (network, method, metric).
    """
    csv_path = resolve_competitor_path(path, config)
    if not os.path.exists(csv_path):
        # This is the ONE place the join can silently degrade. Fail loudly with
        # the exact expected schema so the user can point us at the real file.
        raise SystemExit(
            "\n*** COMPETITOR SCORES NOT FOUND ***\n"
            f"  Looked for: {csv_path}\n"
            "  TODO: pass --competitor-scores pointing at the competitor score file.\n"
            "  Expected: simulations/analysis/summary/<config>/comparisons_raw.csv\n"
            "  Expected columns (long format):\n"
            "    network, config, method, replicate, metric, value, jaccard_mode, status\n"
            "  Expected method labels: grampa, grandma_split (=GRAMPA-iter),\n"
            "    polyphest_p50/p70/p90 (merged to 'polyphest').\n"
        )

    df = pd.read_csv(csv_path)
    required = {"network", "method", "metric", "value", "status"}
    if not required.issubset(df.columns):
        raise SystemExit(
            f"Competitor file {csv_path} missing columns "
            f"{required - set(df.columns)} (found: {list(df.columns)})."
        )
    if "config" in df.columns:
        df = df[df["config"] == config].copy()
    else:
        df["config"] = config
    if "replicate" not in df.columns:
        df["replicate"] = 1

    df = df[df["status"] == "SUCCESS"].copy()
    df = merge_polyphest(df)

    # Keep only the competitor metric names we care about, then translate back to
    # the GNN's bare metric names for a clean join.
    comp_to_gnn = {competitor_metric_name(m): m for m in metrics}
    df = df[df["metric"].isin(comp_to_gnn.keys())].copy()
    df["metric"] = df["metric"].map(comp_to_gnn)

    # Canonical method names; drop methods we are not comparing (mpsugar, padre,
    # alloppnet) so they do not clutter the head-to-head.
    df = df[df["method"].isin(COMPETITOR_METHOD_MAP.keys())].copy()
    df["method"] = df["method"].map(COMPETITOR_METHOD_MAP)

    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    # Average over replicates -> one value per (network, method, metric).
    tidy = (df.groupby(["network", "config", "method", "metric"], as_index=False)["value"]
              .mean())
    return tidy[["network", "config", "method", "metric", "value"]]


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def build_wide(tidy, primary_metric):
    """Pivot the combined tidy frame to one row per network, columns {method}_{metric}."""
    tidy = tidy.copy()
    tidy["col"] = tidy.apply(
        lambda r: f"{r['method']}_{short_name(r['metric'])}", axis=1
    )
    wide = tidy.pivot_table(index="network", columns="col", values="value",
                            aggfunc="mean")

    # Column order: for each metric (primary first), each method in OUTPUT_METHODS.
    metrics_in_order = [primary_metric] + [m for m in EXTRA_METRICS if m != primary_metric]
    ordered = []
    for metric in metrics_in_order:
        for method in OUTPUT_METHODS:
            col = f"{method}_{short_name(metric)}"
            if col in wide.columns:
                ordered.append(col)
    # Append any stragglers not covered by the ordering, deterministically.
    ordered += [c for c in wide.columns if c not in ordered]
    wide = wide.reindex(columns=ordered)
    wide = wide.reset_index()
    return wide


def print_summary(tidy, config, primary_metric):
    """Per-config mean per method for each metric, plus network coverage counts."""
    print(f"\n{'='*72}")
    print(f"Head-to-head summary  |  config = {config}")
    print(f"{'='*72}")

    metrics_in_order = [primary_metric] + [m for m in EXTRA_METRICS if m != primary_metric]
    metrics_in_order = [m for m in metrics_in_order if m in tidy["metric"].unique()]

    mean_tbl = (tidy.groupby(["method", "metric"])["value"].mean()
                    .unstack("metric").reindex(index=OUTPUT_METHODS, columns=metrics_in_order))
    n_tbl = (tidy.dropna(subset=["value"]).groupby("method")["network"].nunique()
                 .reindex(OUTPUT_METHODS))

    out = mean_tbl.copy()
    out.insert(0, "n_networks", n_tbl)
    out = out.dropna(how="all", subset=metrics_in_order)
    print("\nMean per method (lower is better for distances):")
    print(out.to_string(float_format=lambda v: f"{v:.4f}"))
    print()


def print_common_subset(combined, config, primary_metric):
    """FAIR comparison: means on the networks BOTH the GNN and Polyphest completed,
    plus completion counts, plus the GNN's mean on Polyphest-completed vs -failed
    networks (to check whether Polyphest skipped the harder ones). The per-method
    means over *different* subsets (GNN on 21, Polyphest on the 15-18 it finishes)
    are biased against the method that completes more; this is the unbiased view."""
    pm = combined[combined["metric"] == primary_metric]
    wide = pm.pivot_table(index="network", columns="method", values="value", aggfunc="mean")
    if "gnn" not in wide.columns or "polyphest" not in wide.columns:
        print("  (common-subset comparison needs both gnn and polyphest columns)")
        return

    print(f"{'-'*72}")
    print(f"FAIR comparison on {primary_metric}  |  config = {config}")
    print(f"{'-'*72}")
    gnn_ok = wide["gnn"].notna()
    poly_ok = wide["polyphest"].notna()
    common = gnn_ok & poly_ok
    print(f"Completion: total={wide.shape[0]} | gnn={int(gnn_ok.sum())} | "
          f"polyphest={int(poly_ok.sum())} | both={int(common.sum())}")

    print("\nMean on the COMMON subset (networks BOTH completed):")
    for m in OUTPUT_METHODS:
        if m in wide.columns:
            vals = wide.loc[common, m].dropna()
            if len(vals):
                print(f"  {m:<12} {vals.mean():.4f}  (n={len(vals)})")

    poly_failed = gnn_ok & ~poly_ok
    if poly_failed.any():
        print(f"\nGNN {primary_metric}: where Polyphest FAILED vs COMPLETED "
              "(higher 'failed' => Polyphest skipped the harder networks):")
        print(f"  polyphest-completed : {wide.loc[common, 'gnn'].mean():.4f} "
              f"(n={int(common.sum())})")
        print(f"  polyphest-failed    : {wide.loc[poly_failed, 'gnn'].mean():.4f} "
              f"(n={int(poly_failed.sum())})")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Join GNN benchmark scores with Polyphest/GRAMPA/GRAMPA-iter "
                    "into one per-network head-to-head table.")
    parser.add_argument("--gnn-scores", required=True,
                        help="GNN reconstruction_scores.csv or the strategy dir "
                             "containing it (e.g. .../benchmark/<config>/collapse).")
    parser.add_argument("--competitor-scores", default=None,
                        help="Competitor comparisons_raw.csv (or a dir under which "
                             "it lives). Default: "
                             "simulations/analysis/summary/<config>/comparisons_raw.csv")
    parser.add_argument("--config", required=True, help="e.g. conf_ils_low_10M")
    parser.add_argument("--out", required=True, help="output combined per-network CSV")
    parser.add_argument("--metric", default="edit_distance_multree",
                        help="primary comparison metric (default: edit_distance_multree)")
    parser.add_argument("--exclude-networks", default=",".join(DEFAULT_EXCLUDE),
                        help="comma-separated networks to drop (clean-config "
                             f"restriction). Default: {','.join(DEFAULT_EXCLUDE)}")
    args = parser.parse_args()

    exclude = {n.strip() for n in args.exclude_networks.split(",") if n.strip()}

    # Metrics to pull: the primary metric plus the extra jaccard metrics.
    metrics = [args.metric] + [m for m in EXTRA_METRICS if m != args.metric]

    print(f"Config: {args.config}")
    print(f"Primary metric: {args.metric}")
    print(f"Excluding networks: {sorted(exclude) if exclude else '(none)'}")

    gnn = load_gnn_scores(args.gnn_scores, args.config, metrics)
    comp = load_competitor_scores(args.competitor_scores, args.config, metrics)

    combined = pd.concat([gnn, comp], ignore_index=True)

    # Clean-config restriction: drop the wrong-rate networks from BOTH sides.
    if exclude:
        combined = combined[~combined["network"].isin(exclude)].copy()

    if combined.empty:
        raise SystemExit("No rows after loading/filtering — check paths, config, "
                         "and metric name.")

    # Report which methods actually contributed (helps spot a missing side).
    have = sorted(combined["method"].unique())
    print(f"Methods present: {have}")
    for m in OUTPUT_METHODS:
        if m not in have:
            print(f"  WARNING: no rows for '{m}' — column(s) will be empty.")

    wide = build_wide(combined, args.metric)
    wide.insert(1, "config", args.config)
    wide = wide.sort_values("network").reset_index(drop=True)

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    wide.to_csv(args.out, index=False)
    print(f"\nWrote combined per-network table ({len(wide)} networks) to: {args.out}")

    print_summary(combined, args.config, args.metric)
    print_common_subset(combined, args.config, args.metric)


if __name__ == "__main__":
    main()
