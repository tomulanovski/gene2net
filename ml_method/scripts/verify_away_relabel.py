"""Pre-flight check for the away-relabel experiment, on real data.

Assumes both sidecars already exist for the config (run relabel_from_metadata.py
once without --away-parent and once with it first). Verifies, per sample:
  1. detection is UNCHANGED  -- wgd_edges and wgd_counts identical between
     labels_clade.pkl and labels_away.pkl (the relabel must touch only the partner);
  2. the partner DID change for a substantial share of allo events (~the audit's
     ~56% partner==home rate) -- otherwise the relabel is a no-op;
  3. the away labels stay mappable -- the retargeted partner clade still matches an
     ASTRAL edge, so we are not silently dropping allo training examples;
  4. the new dataset path loads labels_away.pkl (away_labels=True).
If all four pass, the experiment is sound to run.

Read-only except it loads the two sidecars. gene2net or final_project env.
Usage:
  python scripts/verify_away_relabel.py --config training_rooted/ils_low --n 300
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from gene2net_gnn.data.dataset import Gene2NetSample


def load_sidecar(sample_dir, name):
    p = os.path.join(sample_dir, name)
    if not os.path.exists(p):
        return None
    with open(p, "rb") as f:
        return pickle.load(f)


def is_allo(lab, i):
    # allo event: partner edge differs from the wgd edge (auto is self-partner)
    return i < len(lab.partner_edges) and lab.partner_edges[i] != lab.wgd_edges[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data/mul_trees_2k")
    ap.add_argument("--config", required=True, help="e.g. training_rooted/ils_low")
    ap.add_argument("--n", type=int, default=300)
    args = ap.parse_args()

    root = os.path.join(args.data_root, args.config)
    dirs = sorted(d for d in os.listdir(root)
                  if d.startswith("sample_") and os.path.isdir(os.path.join(root, d)))[:args.n]

    n_samples = 0
    det_mismatch = 0
    allo_clade = 0
    partner_changed = 0
    away_allo = 0
    away_allo_mappable = 0
    load_ok = 0
    load_fail = 0

    for d in dirs:
        sd = os.path.join(root, d)
        clade = load_sidecar(sd, "labels_clade.pkl")
        away = load_sidecar(sd, "labels_away.pkl")
        if clade is None or away is None:
            continue
        n_samples += 1

        # 1. detection unchanged
        if list(clade.wgd_edges) != list(away.wgd_edges) or list(clade.wgd_counts) != list(away.wgd_counts):
            det_mismatch += 1

        # 2. partner changed on allo events
        for i in range(min(len(clade.partner_edges), len(away.partner_edges))):
            if is_allo(clade, i):
                allo_clade += 1
                if clade.partner_edges[i] != away.partner_edges[i]:
                    partner_changed += 1

        # 3. away allo mappability
        amask = away.mask or [True] * len(away.wgd_edges)
        for i in range(len(away.wgd_edges)):
            if is_allo(away, i):
                away_allo += 1
                if i < len(amask) and amask[i]:
                    away_allo_mappable += 1

        # 4. dataset load path
        try:
            s = Gene2NetSample.load(sd, away_labels=True)
            if s.labels is not None and list(s.labels.partner_edges) == list(away.partner_edges):
                load_ok += 1
            else:
                load_fail += 1
        except Exception:
            load_fail += 1

    print(f"\nAway-relabel verification: {args.config}  ({n_samples} samples with both sidecars)\n")
    if n_samples == 0:
        print("No samples had BOTH sidecars. Relabel first:")
        print(f"  python scripts/relabel_from_metadata.py --training-subdir {args.config}")
        print(f"  python scripts/relabel_from_metadata.py --training-subdir {args.config} --away-parent")
        return
    print(f"1. detection unchanged        : {n_samples - det_mismatch}/{n_samples} samples "
          f"({'PASS' if det_mismatch == 0 else 'FAIL — relabel altered detection!'})")
    frac = partner_changed / allo_clade if allo_clade else 0.0
    print(f"2. allo partner changed       : {partner_changed}/{allo_clade} = {frac:.2f} "
          f"({'PASS — matches ~0.56 audit' if 0.3 <= frac <= 0.8 else 'CHECK — off from audit'})")
    mfrac = away_allo_mappable / away_allo if away_allo else 0.0
    print(f"3. away allo still mappable    : {away_allo_mappable}/{away_allo} = {mfrac:.2f} "
          f"({'PASS' if mfrac >= 0.9 else 'CHECK — losing allo examples'})")
    print(f"4. dataset loads labels_away   : {load_ok}/{load_ok + load_fail} "
          f"({'PASS' if load_fail == 0 and load_ok > 0 else 'FAIL'})")

    ok = det_mismatch == 0 and 0.3 <= frac <= 0.8 and mfrac >= 0.9 and load_fail == 0 and load_ok > 0
    print(f"\nOVERALL: {'PASS — experiment is sound to run' if ok else 'REVIEW the flags above before running'}")


if __name__ == "__main__":
    main()
