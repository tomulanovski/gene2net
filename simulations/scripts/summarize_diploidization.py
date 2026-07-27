#!/usr/bin/env python3
"""Summarize diploidization across all networks/replicates of a config.

Reads every ``diploidization_summary.json`` written by ``apply_diploidization.py``
for a config and reports, per network and overall: the retention rate used, how
many replicates and polyploid observations, the mean copy number of polyploids
after diploidization, and the fraction collapsed to a single copy.

Note: on dup_loss configs the copy numbers reflect SimPhy loss AND diploidization
combined; on ILS configs they reflect diploidization alone.

Usage:
  python summarize_diploidization.py conf_ils_low_10M_dip050
  python summarize_diploidization.py conf_dup_loss_medium_10M_ne1M_fix060 --csv stats.csv
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def _polyploid_species(summary):
    """Species that belong to a WGD event (the polyploids)."""
    species = set()
    for event in summary.get("events", []):
        species.update(event.get("species", []))
    return species


def aggregate(records):
    """Aggregate a list of ``{"network", "summary"}`` records.

    Returns ``{"per_network": {net: metrics}, "overall": metrics}`` where metrics
    are computed over *polyploid* (species, gene) observations only.
    """
    net_counts = defaultdict(Counter)   # network -> Counter(copy_count -> n_obs)
    net_replicates = Counter()
    net_polyploids = defaultdict(set)
    net_qs = defaultdict(list)

    for rec in records:
        net, summary = rec["network"], rec["summary"]
        net_replicates[net] += 1
        polys = _polyploid_species(summary)
        net_polyploids[net] |= polys
        net_qs[net].extend(e["q"] for e in summary.get("events", []))
        for sp in polys:
            for copies, n in summary.get("copy_numbers", {}).get(sp, {}).items():
                net_counts[net][int(copies)] += n

    def metrics(counter):
        n_obs = sum(counter.values())
        mean = sum(c * n for c, n in counter.items()) / n_obs if n_obs else 0.0
        single = counter.get(1, 0) / n_obs if n_obs else 0.0
        return n_obs, mean, single

    per_network = {}
    overall_counter = Counter()
    overall_qs = []
    for net in sorted(net_counts):
        n_obs, mean, single = metrics(net_counts[net])
        qs = net_qs[net]
        per_network[net] = {
            "q": round(sum(qs) / len(qs), 4) if qs else None,
            "n_replicates": net_replicates[net],
            "n_polyploids": len(net_polyploids[net]),
            "n_obs": n_obs,
            "mean_copies": mean,
            "single_copy_frac": single,
        }
        overall_counter += net_counts[net]
        overall_qs.extend(qs)

    n_obs, mean, single = metrics(overall_counter)
    overall = {
        "q": round(sum(overall_qs) / len(overall_qs), 4) if overall_qs else None,
        "n_networks": len(per_network),
        "n_obs": n_obs,
        "mean_copies": mean,
        "single_copy_frac": single,
    }
    return {"per_network": per_network, "overall": overall}


def collect_summaries(base_dir, config):
    """Load every diploidization_summary.json for ``config`` under ``base_dir``."""
    base = Path(base_dir)
    pattern = f"*/data/{config}/replicate_*/diploidization_summary.json"
    records = []
    for path in sorted(base.glob(pattern)):
        network = path.relative_to(base).parts[0]
        records.append({"network": network, "summary": json.loads(path.read_text())})
    return records


def _print_table(agg):
    header = f"{'network':30s} {'q':>5s} {'reps':>4s} {'polyp':>5s} " \
             f"{'mean_copies':>11s} {'single_frac':>11s}"
    print(header)
    print("-" * len(header))
    for net, m in agg["per_network"].items():
        print(f"{net:30s} {m['q']:>5.2f} {m['n_replicates']:>4d} "
              f"{m['n_polyploids']:>5d} {m['mean_copies']:>11.3f} "
              f"{m['single_copy_frac']:>11.3f}")
    o = agg["overall"]
    print("-" * len(header))
    print(f"{'OVERALL ('+str(o['n_networks'])+' networks)':30s} {o['q']:>5.2f} "
          f"{'':>4s} {'':>5s} {o['mean_copies']:>11.3f} {o['single_copy_frac']:>11.3f}")


def _write_csv(agg, path):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["network", "q", "n_replicates", "n_polyploids",
                    "n_obs", "mean_copies", "single_copy_frac"])
        for net, m in agg["per_network"].items():
            w.writerow([net, m["q"], m["n_replicates"], m["n_polyploids"],
                        m["n_obs"], f"{m['mean_copies']:.4f}",
                        f"{m['single_copy_frac']:.4f}"])
        o = agg["overall"]
        w.writerow(["OVERALL", o["q"], "", "", o["n_obs"],
                    f"{o['mean_copies']:.4f}", f"{o['single_copy_frac']:.4f}"])
    print(f"\nWrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="diploidized config name (e.g. ..._fix060)")
    parser.add_argument("--base-dir", default=None,
                        help="dir with network subdirs (default: simulations/simulations)")
    parser.add_argument("--csv", default=None, help="write per-network stats to this CSV")
    args = parser.parse_args(argv)

    base_dir = (args.base_dir if args.base_dir
                else str(Path(__file__).resolve().parent.parent / "simulations"))
    records = collect_summaries(base_dir, args.config)
    if not records:
        raise SystemExit(f"no diploidization summaries found for {args.config} under {base_dir}")

    agg = aggregate(records)
    _print_table(agg)
    if args.csv:
        _write_csv(agg, args.csv)


if __name__ == "__main__":
    main()
