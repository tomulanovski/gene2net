#!/usr/bin/env python3
"""Quantify how well Algorithm 2 (kernel smoothing) classifies polyploids, per config.

Polyphest is given a ploidy list inferred from the gene trees by kernel smoothing
(``copies_smoothing_with_multiset.py``), written per replicate as
``distribution.tsv`` (a ``RepresentativeCopyNumber`` per species). Fractionation
returns duplicated genes to single copy, so the inferred ploidy of true polyploids
degrades as the retention rate q falls. This script measures that degradation as
**recall on true polyploids**: of the species that are genuinely polyploid (from
the ground-truth ``mul_tree_final_stats.csv`` ``Polyploid_Names`` column), the
fraction whose inferred ``RepresentativeCopyNumber`` is >= 2.

The ground-truth names and the ``distribution.tsv`` species share one namespace:
prep_polyphest collapses each gene tip to the part before its first underscore, so
both polyploid copies ``X_1``/``X_2`` become the base accession ``X`` (which is the
name listed in ``Polyploid_Names``).

Layout (same as compare_ploidy.py):
  {base}/{network}/processed/{config}/polyphest_input/replicate_*/distribution.tsv

Usage:
  # one row per config -> the recall-vs-q curve
  python ploidy_accuracy.py conf_dup_loss_medium_10M_ne1M \
      conf_dup_loss_medium_10M_ne1M_fix075 \
      conf_dup_loss_medium_10M_ne1M_fix050 \
      conf_dup_loss_medium_10M_ne1M_fix025 \
      conf_dup_loss_medium_10M_ne1M_fix010 \
      conf_dup_loss_medium_10M_ne1M_fix000
"""
import argparse
import csv
import io
import re
from collections import Counter, defaultdict
from pathlib import Path

REL = "*/processed/{config}/polyphest_input/replicate_*/distribution.tsv"


def parse_true_polyploids(text):
    """Parse mul_tree_final_stats.csv text into {network: set(polyploid names)}.

    ``network`` is the ``Filename`` with the ``.tre`` suffix stripped. A blank
    ``Polyploid_Names`` field (zero polyploids) yields an empty set.
    """
    result = {}
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        network = row["Filename"]
        if network.endswith(".tre"):
            network = network[: -len(".tre")]
        names = {n.strip() for n in (row["Polyploid_Names"] or "").split(",") if n.strip()}
        result[network] = names
    return result


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


def records_for_distribution(network, replicate, true_polyploids, dist):
    """One record per true polyploid of ``network``, scored against ``dist``.

    A record is ``present`` when the species has a row in the distribution (it
    always should, since fractionation keeps >= 1 copy) and ``correct`` when its
    inferred copy number is >= 2 (still called polyploid).
    """
    records = []
    for species in sorted(true_polyploids):
        present = species in dist
        correct = present and dist[species] >= 2
        records.append({
            "network": network, "replicate": replicate, "species": species,
            "inferred": dist.get(species), "present": present, "correct": correct,
        })
    return records


def _replicate_of(path):
    for part in path.parts:
        if part.startswith("replicate_"):
            return part.split("_", 1)[1]
    return None


def collect(base_dir, config, true_polyploids, exclude):
    """Score every replicate distribution.tsv of ``config`` against ground truth."""
    base = Path(base_dir)
    records = []
    for path in sorted(base.glob(REL.format(config=config))):
        network = path.relative_to(base).parts[0]
        if network in exclude:
            continue
        if network not in true_polyploids:
            raise ValueError(f"network {network!r} not in stats CSV")
        replicate = _replicate_of(path)
        dist = parse_distribution(path.read_text())
        records.extend(records_for_distribution(
            network, replicate, true_polyploids[network], dist))
    return records


def aggregate(records):
    """Per-network and overall recall on true polyploids."""
    per_net = defaultdict(Counter)
    for r in records:
        c = per_net[r["network"]]
        if not r["present"]:
            c["n_absent"] += 1
            continue
        c["n_true"] += 1
        if r["correct"]:
            c["n_correct"] += 1

    def finalize(c):
        n_true = c["n_true"]
        return {
            "n_true": n_true,
            "n_correct": c["n_correct"],
            "n_absent": c["n_absent"],
            "recall": round(c["n_correct"] / n_true, 4) if n_true else None,
        }

    per_network = {net: finalize(c) for net, c in sorted(per_net.items())}
    overall = Counter()
    for c in per_net.values():
        overall += c
    return {"per_network": per_network, "overall": finalize(overall)}


def q_from_config(config):
    """Recover the retention rate q from a config name (``_fixNNN`` -> NNN/100).

    A config without a ``_fixNNN`` suffix is the un-fractionated baseline, q = 1.
    """
    m = re.search(r"_fix(\d{3})$", config)
    return int(m.group(1)) / 100 if m else 1.0


def _print_curve(rows):
    header = f"{'q':>5s} {'config':45s} {'n_true':>7s} {'correct':>8s} " \
             f"{'recall':>7s} {'absent':>7s}"
    print(header)
    print("-" * len(header))
    for q, config, agg in rows:
        o = agg["overall"]
        rec = "NA" if o["recall"] is None else f"{o['recall'] * 100:.1f}%"
        print(f"{q:>5.2f} {config:45s} {o['n_true']:>7d} {o['n_correct']:>8d} "
              f"{rec:>7s} {o['n_absent']:>7d}")


def _print_per_network(config, agg):
    print(f"\nPer-network recall for {config}:")
    hdr = f"  {'network':30s} {'n_true':>7s} {'correct':>8s} {'recall':>7s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for net, m in agg["per_network"].items():
        rec = "NA" if m["recall"] is None else f"{m['recall'] * 100:.1f}%"
        print(f"  {net:30s} {m['n_true']:>7d} {m['n_correct']:>8d} {rec:>7s}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("configs", nargs="+",
                        help="config names; a bare name (no _fixNNN) is q=1 baseline")
    parser.add_argument("--stats", default=None,
                        help="mul_tree_final_stats.csv (default: ../networks/...)")
    parser.add_argument("--base-dir", default=None,
                        help="dir with network subdirs (default: ../simulations)")
    parser.add_argument("--exclude", nargs="*", default=["Morales-Briones_2021"],
                        help="networks to skip (default: Morales-Briones_2021)")
    parser.add_argument("--per-network", action="store_true",
                        help="also print per-network recall for each config")
    parser.add_argument("--csv", default=None, help="write per-species records here")
    args = parser.parse_args(argv)

    here = Path(__file__).resolve().parent
    stats_path = Path(args.stats) if args.stats else here.parent / "networks" / "mul_tree_final_stats.csv"
    base_dir = Path(args.base_dir) if args.base_dir else here.parent / "simulations"
    true_polyploids = parse_true_polyploids(Path(stats_path).read_text())

    print(f"Ground truth: {stats_path}")
    print(f"Base dir:     {base_dir}")
    print(f"Excluding:    {', '.join(args.exclude) or '(none)'}")
    print(f"Metric: recall on true polyploids "
          f"(inferred RepresentativeCopyNumber >= 2)\n")

    rows = []
    all_records = []
    for config in args.configs:
        records = collect(base_dir, config, true_polyploids, set(args.exclude))
        if not records:
            raise SystemExit(f"no distribution.tsv found for config {config!r}")
        agg = aggregate(records)
        rows.append((q_from_config(config), config, agg))
        all_records.extend((config, r) for r in records)

    rows.sort(key=lambda x: x[0], reverse=True)  # q=1 first, down to q=0
    _print_curve(rows)
    if args.per_network:
        for _, config, agg in rows:
            _print_per_network(config, agg)

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["config", "q", "network", "replicate", "species",
                        "inferred_copy", "present", "correct"])
            for config, r in all_records:
                w.writerow([config, q_from_config(config), r["network"],
                            r["replicate"], r["species"], r["inferred"],
                            int(r["present"]), int(r["correct"])])
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
