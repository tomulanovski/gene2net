"""Localize the true-backbone faithfulness residual (edit ~0.159, not ~0).

Groups the per-sample edit distances of an oracle run by the sample's event
composition, so we can tell WHERE the build fails to round-trip the ground truth:
auto-only samples isolate the auto duplication; clade-level-allo samples isolate
the 32% clade path; single-allo samples isolate two-parent tip placement.

Run in the gene2net env, pointed at a --backbone true --parents true run:
  python scripts/localize_faithfulness.py \
      --scores output/two_parent_oracle/ils_low_true_bb/reconstruction_scores.csv \
      --mul-trees-dir data/mul_trees_2k
"""
import argparse
import csv
import json
import os


def classify_sample(events):
    """Tags describing a sample's event composition."""
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
    # "sample_0007" -> "0007"
    return name.split("_")[-1]


def localize(scores_csv, mul_trees_dir):
    with open(scores_csv) as f:
        rows = list(csv.DictReader(f))
    groups = {}  # tag -> list of edit values
    overall = []
    for r in rows:
        edit = r.get("edit_distance_multree")
        if edit in ("", None):
            continue
        edit = float(edit)
        overall.append(edit)
        idx = _sample_idx(r.get("sample", ""))
        md_path = os.path.join(mul_trees_dir, f"metadata_{idx}.json")
        if not os.path.exists(md_path):
            continue
        with open(md_path) as f:
            events = json.load(f).get("events", [])
        for tag in classify_sample(events):
            groups.setdefault(tag, []).append(edit)

    out = {"overall": (len(overall), sum(overall) / len(overall) if overall else float("nan"))}
    for tag, vals in sorted(groups.items()):
        out[tag] = (len(vals), sum(vals) / len(vals))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="reconstruction_scores.csv of a true_bb run")
    ap.add_argument("--mul-trees-dir", required=True)
    args = ap.parse_args()
    res = localize(args.scores, args.mul_trees_dir)
    print(f"{'group':<18} {'n':>5} {'mean_edit':>10}")
    for tag, (n, mean) in res.items():
        print(f"{tag:<18} {n:>5} {mean:>10.4f}")


if __name__ == "__main__":
    main()
