"""Aggregate statistics from training data metadata files."""
import argparse
import json
import os
import glob
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description="Training data statistics")
    parser.add_argument("--dir", required=True, help="Directory with metadata_NNNN.json files")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, "metadata_*.json")))
    if not files:
        print(f"No metadata files found in {args.dir}")
        return

    n_species_list = []
    n_events_requested = []
    n_events_applied = []
    n_leaves_mul = []
    n_polyploid_species = []
    event_types = Counter()  # auto vs allo
    events_per_tree = Counter()  # distribution of n_events_applied
    species_bins = Counter()  # binned species counts
    zero_event_count = 0

    for f in files:
        with open(f) as fh:
            meta = json.load(fh)

        n_sp = meta["n_species"]
        n_req = meta["n_events_requested"]
        n_app = meta["n_events_applied"]

        n_species_list.append(n_sp)
        n_events_requested.append(n_req)
        n_events_applied.append(n_app)
        n_leaves_mul.append(meta.get("n_leaves_mul_tree", n_sp))
        n_polyploid_species.append(meta.get("n_polyploid_species", 0))

        events_per_tree[n_app] += 1
        if n_app == 0:
            zero_event_count += 1

        for evt in meta.get("events", []):
            event_types[evt["event_type"]] += 1

        # Bin species counts
        if n_sp <= 10:
            species_bins["5-10"] += 1
        elif n_sp <= 20:
            species_bins["11-20"] += 1
        elif n_sp <= 40:
            species_bins["21-40"] += 1
        elif n_sp <= 60:
            species_bins["41-60"] += 1
        else:
            species_bins["61-80"] += 1

    total = len(files)
    print(f"{'='*60}")
    print(f"Training Data Statistics ({total} trees)")
    print(f"{'='*60}")

    print(f"\n--- Species count ---")
    print(f"  Min: {min(n_species_list)}, Max: {max(n_species_list)}, "
          f"Mean: {sum(n_species_list)/total:.1f}, Median: {sorted(n_species_list)[total//2]}")
    print(f"  Distribution:")
    for bin_name in ["5-10", "11-20", "21-40", "41-60", "61-80"]:
        count = species_bins.get(bin_name, 0)
        pct = 100 * count / total
        bar = "#" * int(pct / 2)
        print(f"    {bin_name:>5}: {count:>5} ({pct:5.1f}%) {bar}")

    print(f"\n--- WGD events applied ---")
    print(f"  Min: {min(n_events_applied)}, Max: {max(n_events_applied)}, "
          f"Mean: {sum(n_events_applied)/total:.1f}")
    print(f"  Zero-event trees: {zero_event_count} ({100*zero_event_count/total:.1f}%)")
    print(f"  Distribution:")
    for n_evt in sorted(events_per_tree.keys()):
        count = events_per_tree[n_evt]
        pct = 100 * count / total
        bar = "#" * int(pct / 2)
        print(f"    {n_evt:>2} events: {count:>5} ({pct:5.1f}%) {bar}")

    print(f"\n--- Events requested vs applied ---")
    total_req = sum(n_events_requested)
    total_app = sum(n_events_applied)
    print(f"  Total requested: {total_req}, Total applied: {total_app} "
          f"({100*total_app/max(total_req,1):.1f}% success rate)")

    print(f"\n--- Event types ---")
    for etype, count in event_types.most_common():
        pct = 100 * count / max(sum(event_types.values()), 1)
        print(f"  {etype}: {count} ({pct:.1f}%)")

    print(f"\n--- MUL-tree leaves ---")
    print(f"  Min: {min(n_leaves_mul)}, Max: {max(n_leaves_mul)}, "
          f"Mean: {sum(n_leaves_mul)/total:.1f}")

    print(f"\n--- Polyploid species per tree ---")
    print(f"  Min: {min(n_polyploid_species)}, Max: {max(n_polyploid_species)}, "
          f"Mean: {sum(n_polyploid_species)/total:.1f}")


if __name__ == "__main__":
    main()
