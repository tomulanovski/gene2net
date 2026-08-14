"""Decompose the ASTRAL-vs-true backbone edit-distance gap by event type.

backbone_error_analysis.py shows the aggregate gap (astral edit - true edit) and
whether it tracks ASTRAL's topology error; localize_faithfulness.py groups the
true-backbone residual by event composition but only PRINTS. This script joins
the same two per-sample score files (astral vs true), computes the per-sample
gap the same way, then splits it by the sample's event composition (auto-only /
single-species-allo / clade-allo, exactly as localize_faithfulness classifies)
and WRITES the numbers to a CSV so the thesis diagnostic section can cite them.

The two inputs are the per-sample score files from the backbone 2x2 (same
predicted events, different backbone):
    <astral-dir>/reconstruction_scores.csv   (built on ASTRAL)
    <true-dir>/reconstruction_scores.csv     (built on the true backbone)

Event composition is read from the mul-tree metadata JSON, keyed by sample index
(sample_0007 -> metadata_0007.json), identical to localize_faithfulness.py.

Buckets follow localize_faithfulness.classify_sample, so a sample contributes to
EVERY tag it matches (a clade-allo sample that also has an auto event lands in
both has_clade_allo and has_auto). The "ALL" row is one row per sample.

Run in the final_project env (needs pandas). Read-only on all inputs.

Usage:
    python scripts/decompose_edit_gap.py \
        --astral-dir output/reconstruct_allo/backbone_exp/ils_low_astral/t0.9 \
        --true-dir   output/reconstruct_allo/backbone_exp/ils_low_true/t0.9 \
        --mul-trees-dir data/mul_trees_2k --config ils_low \
        --out output/backbone_exp/ils_low_edit_gap_by_event.csv
"""
import argparse
import json
import os

import pandas as pd


def classify_sample(events):
    """Tags describing a sample's event composition.

    Verbatim from localize_faithfulness.py so the buckets match exactly:
      - auto event                 -> has_auto
      - allo event, |target| == 1  -> has_single_allo
      - allo event, |target| >  1  -> has_clade_allo
      - no allo events at all       -> auto_only
    A sample can carry several tags at once.
    """
    if not events:
        return {"no_events"}
    tags = set()
    n_allo = 0
    for ev in events:
        etype = ev.get("event_type")
        tgt = ev.get("target_clade") or []
        if etype == "auto":
            tags.add("has_auto")
        elif etype == "allo":
            n_allo += 1
            if len(tgt) == 1:
                tags.add("has_single_allo")
            else:
                tags.add("has_clade_allo")
    if n_allo == 0:
        tags.add("auto_only")
    return tags


def _sample_idx(name):
    # "sample_0007" -> "0007"  (matches localize_faithfulness._sample_idx)
    return name.split("_")[-1]


def load_events(mul_trees_dir, sample_name):
    """Event list from the sample's metadata JSON, or None if it is missing."""
    idx = _sample_idx(sample_name)
    md_path = os.path.join(mul_trees_dir, f"metadata_{idx}.json")
    if not os.path.exists(md_path):
        return None
    with open(md_path) as f:
        return json.load(f).get("events", [])


def summarize(bucket, records):
    """One tidy row for a bucket. records = list of (astral, true, gap)."""
    astral = pd.Series([r[0] for r in records], dtype="float64")
    true = pd.Series([r[1] for r in records], dtype="float64")
    gap = pd.Series([r[2] for r in records], dtype="float64")
    return {
        "bucket": bucket,
        "n": len(records),
        "mean_astral": astral.mean(),
        "median_astral": astral.median(),
        "mean_true": true.mean(),
        "median_true": true.median(),
        "mean_gap": gap.mean(),
        "median_gap": gap.median(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--astral-dir", required=True,
                        help="dir with reconstruction_scores.csv built on ASTRAL")
    parser.add_argument("--true-dir", required=True,
                        help="dir with reconstruction_scores.csv built on the true backbone")
    parser.add_argument("--mul-trees-dir", required=True,
                        help="dir with metadata_NNNN.json per sample")
    parser.add_argument("--config", required=True,
                        help="config label, recorded in the output for provenance")
    parser.add_argument("--metric", default="edit_distance_multree")
    parser.add_argument("--out", required=True, help="output CSV path")
    args = parser.parse_args()

    a = pd.read_csv(os.path.join(args.astral_dir, "reconstruction_scores.csv"))
    t = pd.read_csv(os.path.join(args.true_dir, "reconstruction_scores.csv"))
    m = args.metric

    # Same join + gap as backbone_error_analysis.py.
    df = a[["sample", m]].merge(t[["sample", m]], on="sample", suffixes=("_astral", "_true"))
    ca, ct = f"{m}_astral", f"{m}_true"
    df = df.dropna(subset=[ca, ct])
    df["gap"] = df[ca] - df[ct]

    if df.empty:
        print("No samples with both astral and true scores — check the score CSVs.")
        return

    buckets = {}          # tag -> list of (astral, true, gap)
    all_records = []      # every joined sample (one per sample)
    missing_md = 0
    for _, r in df.iterrows():
        rec = (float(r[ca]), float(r[ct]), float(r["gap"]))
        all_records.append(rec)
        events = load_events(args.mul_trees_dir, r["sample"])
        if events is None:
            missing_md += 1
            continue
        for tag in classify_sample(events):
            buckets.setdefault(tag, []).append(rec)

    rows = [summarize("ALL", all_records)]
    for tag in sorted(buckets):
        rows.append(summarize(tag, buckets[tag]))

    out_df = pd.DataFrame(rows)
    out_df.insert(0, "config", args.config)
    out_df.insert(1, "metric", m)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out_df.to_csv(args.out, index=False)

    print(f"\nEdit-gap decomposition ({args.config}, metric={m})")
    print(f"joined samples: {len(all_records)}   "
          f"(metadata missing for {missing_md}, excluded from event buckets)\n")
    hdr = (f"{'bucket':<16} {'n':>5} {'mean_astral':>12} {'mean_true':>10} "
           f"{'mean_gap':>10} {'median_gap':>11}")
    print(hdr)
    print("-" * len(hdr))
    for row in rows:
        print(f"{row['bucket']:<16} {row['n']:>5} {row['mean_astral']:>12.4f} "
              f"{row['mean_true']:>10.4f} {row['mean_gap']:>10.4f} "
              f"{row['median_gap']:>11.4f}")
    print(f"\nSaved to {args.out}")


if __name__ == "__main__":
    main()
