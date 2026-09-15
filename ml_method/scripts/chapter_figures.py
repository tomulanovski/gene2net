#!/usr/bin/env python3
"""Chapter figures: degradation lines for the benchmark and for the fractionation series.

Two figures, in the style of simulations/scripts/create_polyphest_vs_grampaiter.py:

  discordance_degradation     3 x 4 grid. Rows are the three measures, columns the four
                              configuration families, x the Low, Medium, High level.
  fractionation_degradation   1 x 3 panels over retention 1.00, 0.75, 0.50, 0.25. The 1.00
                              point is the unfractionated medium dup/loss, medium ILS config.

Every point is restricted to the networks that every method completed in that configuration,
as in common_subset_tables.py, so the figures show the same means as the appendix tables.

Data sources:
  default         per-network scores on the cluster, loaded through
                  common_subset_tables.per_network. Error bars are the standard error across
                  the common-subset networks. The means are then checked against the appendix
                  tables, and any disagreement stops the run.
  --from-tables   the labeled tables in docs/chapter_appendix_tables_draft.md. Means only, no
                  error bars. A local layout preview for when the per-network scores are absent.

Run from ml_method/:
    python scripts/chapter_figures.py                # cluster, final figures
    python scripts/chapter_figures.py --from-tables  # local preview
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
APPENDIX = "chapter_appendix_tables_draft.md"

MEASURES = [  # (metric key in the score files, y-axis label); row order of the figure
    ("ret_leaf_jaccard", "Ret. descendants distance"),
    ("ret_sisters_jaccard", "Ret. sister distance"),
    ("mu_distance", r"$\mu$-distance"),
]
LEVELS = ["Low", "Medium", "High"]
FAMILIES = [  # (column title, x-axis label, level -> config)
    ("ILS only", "ILS level",
     {"Low": "conf_ils_low_10M", "Medium": "conf_ils_medium_10M", "High": "conf_ils_high_10M"}),
    ("Dup/loss, low ILS", "Dup/loss rate",
     {"Low": "conf_dup_loss_low_10M", "Medium": "conf_dup_loss_medium_10M",
      "High": "conf_dup_loss_high_10M"}),
    ("Dup/loss, medium ILS", "Dup/loss rate",
     {"Low": "conf_dup_loss_low_10M_ne1M", "Medium": "conf_dup_loss_medium_10M_ne1M",
      "High": "conf_dup_loss_high_10M_ne1M"}),
    ("Dup/loss, high ILS", "Dup/loss rate",
     {"Low": "conf_dup_loss_low_10M_ne2M", "Medium": "conf_dup_loss_medium_10M_ne2M",
      "High": "conf_dup_loss_high_10M_ne2M"}),
]
FRACTIONATION = [  # (retention label, config); lower retention means more fractionation
    ("1.00", "conf_dup_loss_medium_10M_ne1M"),
    ("0.75", "conf_dup_loss_medium_10M_ne1M_fix075"),
    ("0.50", "conf_dup_loss_medium_10M_ne1M_fix050"),
    ("0.25", "conf_dup_loss_medium_10M_ne1M_fix025"),
]

# Hue is the method. Within a method, the dashed line with an open marker is the variant
# that takes ploidy information. Polyphest orange and GRAMPA-iter sky blue are the thesis
# colors. PlaceNet purple was chosen against them with the dataviz palette validator, all
# pairs, light mode: worst colorblind separation dE 22.2, normal vision dE 26.6.
SERIES = [  # (key, legend label, color, linestyle)
    ("free", "PlaceNet ploidy-free", "#7A5195", "-"),
    ("informed", "PlaceNet ploidy-informed", "#7A5195", "--"),
    ("polyphest", "Polyphest", "#DE8F05", "-"),
    ("grampa_iter", r"GRAMPA$^{Iter}$", "#56B4E9", "-"),
    ("grampa_iter_prior", r"GRAMPA$^{Iter}$ with prior", "#56B4E9", "--"),
]
KEYS = [s[0] for s in SERIES]

# Sized for a Word page, so the text stays legible when the figure is placed at 6.3 inches.
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.linewidth": 0.8,
    "savefig.dpi": 300, "lines.linewidth": 1.5, "lines.markersize": 4.5,
})

# ---------------------------------------------------------------- appendix tables
TABLE_METRIC = {"tab:descendants": "ret_leaf_jaccard", "tab:sisters": "ret_sisters_jaccard",
                "tab:mu": "mu_distance"}
TABLE_FRAC = {"tab:frac075": "conf_dup_loss_medium_10M_ne1M_fix075",
              "tab:frac050": "conf_dup_loss_medium_10M_ne1M_fix050",
              "tab:frac025": "conf_dup_loss_medium_10M_ne1M_fix025"}
ROW_CONFIG = {"ILS low": "conf_ils_low_10M", "ILS medium": "conf_ils_medium_10M",
              "ILS high": "conf_ils_high_10M"}
for _rate in ("low", "medium", "high"):
    ROW_CONFIG[f"dup/loss {_rate}, Ne 200k"] = f"conf_dup_loss_{_rate}_10M"
    ROW_CONFIG[f"dup/loss {_rate}, Ne 1M"] = f"conf_dup_loss_{_rate}_10M_ne1M"
    ROW_CONFIG[f"dup/loss {_rate}, Ne 2M"] = f"conf_dup_loss_{_rate}_10M_ne2M"
BENCH_COLUMNS = ["informed", "free", "polyphest", "grampa_iter", "grampa_iter_prior"]
FRAC_ROW = {"PlaceNet ploidy-informed": "informed", "PlaceNet ploidy-free": "free",
            "Polyphest": "polyphest", "GRAMPA-iter": "grampa_iter",
            "GRAMPA-iter with prior": "grampa_iter_prior"}


def table_rows(text, label):
    """Data rows of the pipe table under the caption carrying {#label}, or None if absent."""
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if line.strip().startswith("Table:") and "{#" + label + "}" in line:
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            rows = []
            while j < len(lines) and lines[j].strip().startswith("|"):
                s = lines[j].strip()
                if not set(s) <= set("|-: "):
                    rows.append([c.strip() for c in s.strip("|").split("|")])
                j += 1
            return rows[1:]
    return None


def load_tables(docs):
    """{(config, metric): {series: (mean, None, n)}} from the appendix's labeled tables."""
    path = os.path.join(docs, APPENDIX)
    if not os.path.exists(path):
        raise SystemExit(f"appendix tables not found: {path}")
    text = open(path, encoding="utf-8").read()
    out = {}
    for label, metric in TABLE_METRIC.items():
        rows = table_rows(text, label)
        if not rows:
            raise SystemExit(f"table {{#{label}}} not found in {APPENDIX}")
        for r in rows:
            if r[0] not in ROW_CONFIG:
                raise SystemExit(f"{label}: unrecognized configuration row {r[0]!r}")
            width = len(BENCH_COLUMNS)
            n = int(r[width + 1])
            out[(ROW_CONFIG[r[0]], metric)] = {
                k: (float(v), None, n) for k, v in zip(BENCH_COLUMNS, r[1:width + 1])}
    metrics = [m for m, _ in MEASURES]   # fractionation columns are in measure order
    for label, config in TABLE_FRAC.items():
        rows = table_rows(text, label)
        if not rows:
            raise SystemExit(f"table {{#{label}}} not found in {APPENDIX}")
        for r in rows:
            if r[0] not in FRAC_ROW:
                raise SystemExit(f"{label}: unrecognized method row {r[0]!r}")
            for metric, v in zip(metrics, r[1:4]):
                out.setdefault((config, metric), {})[FRAC_ROW[r[0]]] = (float(v), None, None)
    return out


# ---------------------------------------------------------------- cluster scores
def load_cluster(base, configs):
    """{(config, metric): {series: (mean, sem, n)}} on each configuration's common subset."""
    from common_subset_tables import per_network
    out = {}
    for config in configs:
        informed = per_network(config, "bound_driven", base)
        free = per_network(config, "detect_only", base)
        tidy = pd.concat([
            informed[informed["method"] == "gnn"].assign(method="informed"),
            free[free["method"] == "gnn"].assign(method="free"),
            informed[informed["method"] != "gnn"],
        ], ignore_index=True)
        for metric, _ in MEASURES:
            wide = tidy[tidy["metric"] == metric].pivot_table(
                index="network", columns="method", values="value", aggfunc="mean")
            missing = [k for k in KEYS if k not in wide.columns]
            if missing:
                raise SystemExit(f"{config} {metric}: no scores for {missing}")
            sub = wide.loc[wide[KEYS].notna().all(axis=1), KEYS]
            if sub.empty:
                raise SystemExit(f"{config} {metric}: the common subset is empty")
            n = len(sub)
            out[(config, metric)] = {
                k: (sub[k].mean(), sub[k].std(ddof=1) / np.sqrt(n) if n > 1 else 0.0, n)
                for k in KEYS}
    return out


def cross_check(computed, tables):
    """Stop if any figure mean or n differs from the appendix tables at their printed precision."""
    bad = []
    for key, points in tables.items():
        for k, (mean, _, n) in points.items():
            cm, _, cn = computed[key][k]
            if f"{cm:.3f}" != f"{mean:.3f}" or (n is not None and cn != n):
                bad.append(f"{key[0]} {key[1]} {k}: figure {cm:.3f} n={cn}, table {mean:.3f} n={n}")
    if bad:
        raise SystemExit("Figure means disagree with the appendix tables:\n  " + "\n  ".join(bad))
    print("cross-check: every point matches the appendix tables")


# ---------------------------------------------------------------- drawing
DODGE = np.linspace(-0.12, 0.12, len(SERIES))  # small x offsets so error bars do not overlap


def draw(ax, points):
    xs = np.arange(len(points))
    for (key, label, color, ls), dx in zip(SERIES, DODGE):
        y = [p[key][0] for p in points]
        errs = [p[key][1] for p in points]
        ax.errorbar(xs + dx, y, yerr=None if None in errs else errs, color=color, linestyle=ls,
                    marker="o", markerfacecolor=color if ls == "-" else "white",
                    markeredgewidth=1.2, capsize=2, label=label)
    ax.set_xticks(xs)
    ax.grid(True, alpha=0.25, linestyle="--")


# Legend columns fill top to bottom, so this order puts one method in each column:
# the two PlaceNet modes, the two GRAMPA-iter variants, then Polyphest.
LEGEND_ORDER = ["free", "informed", "grampa_iter", "grampa_iter_prior", "polyphest"]


def legend(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    order = [KEYS.index(k) for k in LEGEND_ORDER]
    fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper center",
               ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.join(out_dir, name)}.png and .pdf")


def discordance_figure(data, out_dir):
    fig, axes = plt.subplots(len(MEASURES), len(FAMILIES), figsize=(7.2, 6.2),
                             sharey="row", squeeze=False)
    for r, (metric, ylabel) in enumerate(MEASURES):
        for c, (title, xlabel, cfgs) in enumerate(FAMILIES):
            ax = axes[r, c]
            draw(ax, [data[(cfgs[level], metric)] for level in LEVELS])
            ax.set_xticklabels(LEVELS)
            if r == 0:
                ax.set_title(title, fontsize=8.5, fontweight="bold")
            if r == len(MEASURES) - 1:
                ax.set_xlabel(xlabel)
            if c == 0:
                ax.set_ylabel(ylabel)
    legend(fig, axes[0, 0])
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, out_dir, "discordance_degradation")


def fractionation_figure(data, out_dir):
    fig, axes = plt.subplots(1, len(MEASURES), figsize=(7.2, 2.7), squeeze=False)
    for ax, (metric, ylabel) in zip(axes[0], MEASURES):
        draw(ax, [data[(config, metric)] for _, config in FRACTIONATION])
        ax.set_xticklabels([label for label, _ in FRACTIONATION])
        ax.set_xlabel("Retention rate")
        ax.set_ylabel(ylabel)
    legend(fig, axes[0, 0])
    fig.tight_layout(rect=[0, 0, 1, 0.84])
    save(fig, out_dir, "fractionation_degradation")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-tables", action="store_true",
                    help="preview from the appendix tables, means only, no error bars")
    ap.add_argument("--no-check", action="store_true",
                    help="cluster mode: skip the cross-check against the appendix tables")
    ap.add_argument("--base", default="output/reconstruct_final")
    ap.add_argument("--docs", default=os.path.join(HERE, "..", "docs"))
    ap.add_argument("--out-dir", default=os.path.join(HERE, "..", "docs", "figures"))
    args = ap.parse_args()

    configs = []
    for _, _, cfgs in FAMILIES:
        configs += [cfgs[level] for level in LEVELS]
    configs += [c for _, c in FRACTIONATION if c not in configs]

    if args.from_tables:
        data = load_tables(args.docs)
        print("source: appendix tables, means only, no error bars")
    else:
        data = load_cluster(args.base, configs)
        print("source: per-network scores, error bars are standard errors on the common subset")
        if not args.no_check:
            cross_check(data, load_tables(args.docs))

    gaps = [f"{c} {m}" for c in configs for m, _ in MEASURES
            if (c, m) not in data or any(k not in data[(c, m)] for k in KEYS)]
    if gaps:
        raise SystemExit("no value for:\n  " + "\n  ".join(gaps))

    discordance_figure(data, args.out_dir)
    fractionation_figure(data, args.out_dir)


if __name__ == "__main__":
    main()
