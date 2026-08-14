"""How many true WGD events get DROPPED by the Jaccard label mapping?

Each true event is matched to ASTRAL edges (wgd clade + partner clade) by Jaccard.
If either side scores below 0.5 the event is unmappable and excluded from training
and from our in-distribution metrics (labels.mask[i] = False). This sizes the
'events outside our per-edge framing' limitation: small fraction = footnote, large
fraction = a real limitation to state to the PI.

Reads packaged sample labels only (no model). Reports, per config:
  - total true events, dropped events, drop rate
  - fraction of samples that lose at least one event

Usage:
    python scripts/event_drop_rate.py --data-dir data/mul_trees_2k/training/ils_low \
        data/mul_trees_2k/training/dup_loss_high data/mul_trees_2k/training/ils_high
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gene2net_gnn.data.dataset import Gene2NetDataset


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, nargs="+")
    parser.add_argument("--max-samples", type=int, default=2000)
    args = parser.parse_args()

    for d in args.data_dir:
        config = os.path.basename(os.path.normpath(d))
        try:
            ds = Gene2NetDataset(d)
        except Exception as e:
            print(f"{config}: could not load ({e})")
            continue

        n_samples = 0
        total_events = 0
        dropped = 0
        samples_with_drop = 0

        for i in range(min(len(ds), args.max_samples)):
            try:
                s = ds[i]
            except Exception:
                continue
            if s.labels is None or s.labels.mask is None:
                continue
            mask = s.labels.mask
            if len(mask) == 0:
                continue
            n_samples += 1
            n_ev = len(mask)
            n_drop = sum(1 for m in mask if not m)
            total_events += n_ev
            dropped += n_drop
            if n_drop > 0:
                samples_with_drop += 1

        rate = 100 * dropped / max(total_events, 1)
        samp_rate = 100 * samples_with_drop / max(n_samples, 1)
        print(f"\n{config} ({n_samples} samples with events):")
        print(f"  true events: {total_events}, dropped: {dropped} ({rate:.1f}%)")
        print(f"  samples losing >=1 event: {samples_with_drop}/{n_samples} ({samp_rate:.1f}%)")

    print("\n<5% dropped = footnote. >20% = a real limitation to state to the PI.")


if __name__ == "__main__":
    main()
