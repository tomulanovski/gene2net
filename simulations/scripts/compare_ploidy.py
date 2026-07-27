#!/usr/bin/env python3
"""Compare Polyphest ploidy inference between two configs (e.g. original vs diploidized).

Polyphest writes a per-replicate ``distribution.tsv`` with a
``RepresentativeCopyNumber`` (its inferred ploidy) per species. This script pairs
that file between two configs and reports how the inferred ploidy changed, with
the headline metric being how often a species Polyphest called polyploid
(copy number >= 2) in config A is still called polyploid in config B.

Layout:
  {base}/{network}/processed/{config}/polyphest_input/replicate_*/distribution.tsv

Usage:
  python compare_ploidy.py conf_dup_loss_medium_10M_ne1M conf_dup_loss_medium_10M_ne1M_fix060
  python compare_ploidy.py ORIG DIP --csv ploidy_changes.csv
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

REL = "*/processed/{config}/polyphest_input/replicate_*/distribution.tsv"


def parse_distribution(text):
    """Parse a distribution.tsv into {species: representative_copy_number}."""
    result = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("Species"):
            continue
        parts = line.split()
        result[parts[0]] = int(parts[1])
    return result


def _replicate_of(path):
    for part in path.parts:
        if part.startswith("replicate_"):
            return part.split("_", 1)[1]
    return None


def collect_pairs(base_dir, config_a, config_b):
    """Pair distribution.tsv files between the two configs by network/replicate/species."""
    base = Path(base_dir)
    records = []
    for path_a in sorted(base.glob(REL.format(config=config_a))):
        network = path_a.relative_to(base).parts[0]
        replicate = _replicate_of(path_a)
        path_b = (base / network / "processed" / config_b / "polyphest_input"
                  / f"replicate_{replicate}" / "distribution.tsv")
        if not path_b.exists():
            print(f"WARNING: no matching file for {network} replicate {replicate} in {config_b}")
            continue
        da = parse_distribution(path_a.read_text())
        db = parse_distribution(path_b.read_text())
        for sp, a in da.items():
            if sp in db:
                records.append({"network": network, "replicate": replicate,
                                "species": sp, "a": a, "b": db[sp]})
    return records


def aggregate(records):
    """Per-network and overall polyploid-retention metrics.

    A species is 'polyploid' in a config when its representative copy number >= 2.
    """
    per_net = defaultdict(lambda: Counter())
    for r in records:
        net = r["network"]
        if r["a"] >= 2:
            per_net[net]["n_polyploid_orig"] += 1
            if r["b"] >= 2:
                per_net[net]["n_still_polyploid"] += 1
            else:
                per_net[net]["n_lost"] += 1
        if r["b"] != r["a"]:
            per_net[net]["n_changed"] += 1

    def finalize(c):
        orig = c["n_polyploid_orig"]
        return {
            "n_polyploid_orig": orig,
            "n_still_polyploid": c["n_still_polyploid"],
            "n_lost": c["n_lost"],
            "n_changed": c["n_changed"],
            "retained_frac": round(c["n_still_polyploid"] / orig, 4) if orig else None,
        }

    per_network = {net: finalize(c) for net, c in sorted(per_net.items())}
    overall = Counter()
    for c in per_net.values():
        overall += c
    return {"per_network": per_network, "overall": finalize(overall)}


def _print_table(agg):
    header = f"{'network':30s} {'poly_orig':>9s} {'still_poly':>10s} " \
             f"{'lost':>5s} {'retained':>8s}"
    print(header)
    print("-" * len(header))
    for net, m in agg["per_network"].items():
        ret = "NA" if m["retained_frac"] is None else f"{m['retained_frac']:.3f}"
        print(f"{net:30s} {m['n_polyploid_orig']:>9d} {m['n_still_polyploid']:>10d} "
              f"{m['n_lost']:>5d} {ret:>8s}")
    o = agg["overall"]
    ret = "NA" if o["retained_frac"] is None else f"{o['retained_frac']:.3f}"
    print("-" * len(header))
    print(f"{'OVERALL':30s} {o['n_polyploid_orig']:>9d} {o['n_still_polyploid']:>10d} "
          f"{o['n_lost']:>5d} {ret:>8s}")


def _write_csv(records, path):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "replicate", "species", "ploidy_a", "ploidy_b", "changed"])
        for r in records:
            w.writerow([r["network"], r["replicate"], r["species"],
                        r["a"], r["b"], int(r["a"] != r["b"])])
    print(f"\nWrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config_a", help="baseline config (e.g. original)")
    parser.add_argument("config_b", help="comparison config (e.g. diploidized)")
    parser.add_argument("--base-dir", default=None,
                        help="dir with network subdirs (default: simulations/simulations)")
    parser.add_argument("--csv", default=None, help="write per-species changes to this CSV")
    args = parser.parse_args(argv)

    base_dir = (args.base_dir if args.base_dir
                else str(Path(__file__).resolve().parent.parent / "simulations"))
    records = collect_pairs(base_dir, args.config_a, args.config_b)
    if not records:
        raise SystemExit("no paired distribution.tsv files found")

    print(f"Comparing ploidy: {args.config_a}  ->  {args.config_b}")
    print(f"(polyploid = RepresentativeCopyNumber >= 2; 'retained' = still polyploid in B)\n")
    _print_table(aggregate(records))
    if args.csv:
        _write_csv(records, args.csv)


if __name__ == "__main__":
    main()
