"""Count auto/allo and single-species/clade-level events across metadata.

Resolves the memory contradiction (100% single-species tips vs 131/196) and
reports the max target multiplicity, which decides whether the two-parent
primitive's clade path and >2-parent handling matter in practice.

Run in the gene2net env on the cluster:
  python scripts/metadata_event_stats.py --mul-trees-dir data/mul_trees_2k \
      --config ils_low --n 2000
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def event_stats(metadata_events):
    n_auto = n_allo = n_allo_single = n_allo_clade = 0
    max_target = 0
    for ev in metadata_events:
        tgt = ev.get("target_clade") or []
        max_target = max(max_target, len(tgt))
        if ev.get("event_type") == "auto":
            n_auto += 1
        elif ev.get("event_type") == "allo":
            n_allo += 1
            if len(tgt) == 1:
                n_allo_single += 1
            else:
                n_allo_clade += 1
    return {
        "n_auto": n_auto, "n_allo": n_allo,
        "n_allo_single": n_allo_single, "n_allo_clade": n_allo_clade,
        "max_target_copies": max_target,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul-trees-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()

    agg = {"n_auto": 0, "n_allo": 0, "n_allo_single": 0, "n_allo_clade": 0, "max_target_copies": 0}
    n_samples = 0
    for idx in range(args.start, args.start + args.n):
        md_path = os.path.join(args.mul_trees_dir, f"metadata_{idx:04d}.json")
        if not os.path.exists(md_path):
            continue
        with open(md_path) as f:
            events = json.load(f).get("events", [])
        s = event_stats(events)
        for k in ("n_auto", "n_allo", "n_allo_single", "n_allo_clade"):
            agg[k] += s[k]
        agg["max_target_copies"] = max(agg["max_target_copies"], s["max_target_copies"])
        n_samples += 1

    print(f"config={args.config} samples={n_samples}")
    for k, v in agg.items():
        print(f"  {k}: {v}")
    if agg["n_allo"]:
        pct = 100 * agg["n_allo_clade"] / agg["n_allo"]
        print(f"  clade-level allo fraction: {pct:.1f}%  "
              f"(single={agg['n_allo_single']}, clade={agg['n_allo_clade']})")


if __name__ == "__main__":
    main()
