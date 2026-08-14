"""Summarize an oracle run: mean scores + allo snap-score means (the artifact guard)."""
import argparse
import csv


def _mean(vals):
    vals = [float(v) for v in vals if v not in ("", None)]
    return sum(vals) / len(vals) if vals else float("nan")


def summarize(scores_csv, mapping_csv=None):
    out = {}
    with open(scores_csv) as f:
        rows = list(csv.DictReader(f))
    for m in ("edit_distance_multree", "ret_leaf_jaccard", "ret_sisters_jaccard"):
        out[m] = _mean([r.get(m) for r in rows])
    if mapping_csv:
        with open(mapping_csv) as f:
            mrows = [r for r in csv.DictReader(f) if r.get("type") == "allo"]
        out["allo_A_snap"] = _mean([r.get("A") for r in mrows])
        out["allo_B_snap"] = _mean([r.get("B") for r in mrows])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="reconstruction_scores.csv")
    ap.add_argument("--mapping", default=None, help="mapping_scores.csv")
    args = ap.parse_args()
    for k, v in summarize(args.scores, args.mapping).items():
        print(f"{k}: {v:.4f}")


if __name__ == "__main__":
    main()
