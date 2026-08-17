"""Average per-network reconstruction scores across replicates.

benchmark_networks writes one reconstruction_scores.csv per replicate (scored under
<base>/rep<r>/<strategy>/). This averages the per-network metrics across those
replicate CSVs into a single reconstruction_scores.csv that head_to_head consumes,
matching how the competitor side is aggregated over replicates by run_full_summary.

The 'sample' column is the network name (score_reconstructions uses the output dir
name, which benchmark_networks names by network). We group on it and take the mean
of every numeric metric over the replicates.

Usage:
    python scripts/aggregate_replicates.py \
        --rep-dirs output/.../rep1/detect_only output/.../rep2/detect_only [...] \
        --out output/.../agg/detect_only/reconstruction_scores.csv
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rep-dirs", nargs="+", required=True,
                    help="One dir per replicate, each holding a reconstruction_scores.csv "
                         "(or the csv path directly).")
    ap.add_argument("--out", required=True, help="Output aggregated reconstruction_scores.csv")
    ap.add_argument("--key", default="sample", help="Per-network key column (default: sample).")
    args = ap.parse_args()

    frames = []
    for d in args.rep_dirs:
        path = os.path.join(d, "reconstruction_scores.csv") if os.path.isdir(d) else d
        if not os.path.exists(path):
            raise FileNotFoundError(f"missing replicate scores: {path}")
        frames.append(pd.read_csv(path))

    df = pd.concat(frames, ignore_index=True)
    if args.key not in df.columns:
        raise ValueError(f"key column {args.key!r} not in the scores "
                         f"(columns: {list(df.columns)})")
    num_cols = [c for c in df.select_dtypes("number").columns if c != args.key]
    agg = df.groupby(args.key, as_index=False)[num_cols].mean()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    agg.to_csv(args.out, index=False)
    print(f"Aggregated {len(frames)} replicate files over {len(agg)} networks -> {args.out}")


if __name__ == "__main__":
    main()
