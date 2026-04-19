#!/usr/bin/env python3
"""
Analyze method runtimes across all configurations from SLURM log files.

Parses log files for each method, extracts Duration (seconds) from completed jobs,
and produces a summary CSV + console report.

Usage (on cluster):
    python analyze_runtime.py

Output:
    simulations/analysis/summary/cross_config/runtime/
        runtime_summary.csv          - per (method, config): mean, median, std, min, max, n
        runtime_per_job.csv          - every successful job with its duration
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict
import csv
import statistics

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
SIM_ROOT = SCRIPT_DIR.parent  # simulations/

# Detect cluster vs local
CLUSTER_BASE = Path("/groups/itay_mayrose/tomulanovski/gene2net/simulations")
if CLUSTER_BASE.exists():
    LOG_DIR = CLUSTER_BASE / "logs"
    BASE_DIR = CLUSTER_BASE / "simulations"
    print(f"Running on cluster. Log dir: {LOG_DIR}")
else:
    # Local fallback — check if logs were synced
    LOG_DIR = SIM_ROOT / "logs"
    BASE_DIR = SIM_ROOT / "simulations"
    print(f"Running locally. Log dir: {LOG_DIR}")

OUT_DIR = SIM_ROOT / "analysis" / "summary" / "cross_config" / "runtime"

# ── Configurations ─────────────────────────────────────────────────────────────
CONFIGS = [
    "conf_ils_low_10M",
    "conf_ils_medium_10M",
    "conf_ils_high_10M",
    "conf_dup_loss_low_10M",
    "conf_dup_loss_medium_10M",
    "conf_dup_loss_high_10M",
    "conf_dup_loss_low_10M_ne1M",
    "conf_dup_loss_medium_10M_ne1M",
    "conf_dup_loss_high_10M_ne1M",
]

# ── Methods and their log patterns ────────────────────────────────────────────
# Log filename patterns:
#   grampa:         run_grampa_<config>_<jobid>_<taskid>.out
#   polyphest:      run_polyphest_<config>_<jobid>_<taskid>.out
#   padre:          run_padre_<config>_<jobid>_<taskid>.out
#   mpsugar:        run_mpsugar_<config>_i<iter>_c<chains>_<jobid>_<taskid>.out
#   grandma_split:  run_grandma_split_<config>_<taskid>.out  (uses %a only, no %A)
#   alloppnet:      alloppnet_run_<config>_rep<rep>_<jobid>_<taskid>.out
#                   OR alloppnet_<config>_rep<rep>_<jobid>_<taskid>.out

NETWORKS_21 = [
    "Bendiksby_2011", "Koenen_2020", "Brysting_2007", "Lawrence_2016",
    "Diaz-Perez_2018", "Wisecaver_2023", "Ding_2023", "Liang_2019",
    "Popp_2005", "Wu_2015", "Liu_2023", "Ren_2024", "Marcussen_2011",
    "Marcussen_2012", "Sessa_2012b", "Zhao_2021", "Hori_2014",
    "Marcussen_2015", "Shahrestani_2015", "Morales-Briones_2021", "Soza_2014",
]

TETRA_NETWORKS = [
    "Bendiksby_2011", "Ding_2023", "Koenen_2020", "Liu_2023",
    "Shahrestani_2015", "Wisecaver_2023", "Wu_2015", "Zhao_2021",
]

# Display names for thesis
METHOD_DISPLAY = {
    "grampa": "GRAMPA",
    "grandma_split": "GRAMPAIter",
    "polyphest": "Polyphest",
    "padre": "PADRE",
    "mpsugar": "MPAllopp",
    "alloppnet": "AlloppNET",
}


def task_to_network_replicate(task_id, networks=NETWORKS_21):
    """Convert 1-based task_id to (network, replicate)."""
    idx = (task_id - 1) // 5
    rep = ((task_id - 1) % 5) + 1
    if idx < len(networks):
        return networks[idx], rep
    return None, None


def parse_duration_from_log(filepath):
    """
    Extract duration in seconds from a log file.
    Looks for 'Duration: X seconds' and 'COMPLETED SUCCESSFULLY'.
    Returns duration (int) or None if not found / not completed.
    """
    try:
        with open(filepath, 'r', errors='replace') as f:
            content = f.read()
    except (IOError, OSError):
        return None

    # Check completion
    if not re.search(r'COMPLETED SUCCESSFULLY', content, re.IGNORECASE):
        return None

    # Extract duration
    m = re.search(r'Duration:\s*(\d+)\s*seconds', content)
    if m:
        dur = int(m.group(1))
        return dur if dur > 0 else None

    # Fallback: look for elapsed time patterns like "Elapsed: 01:23:45"
    m = re.search(r'Elapsed.*?(\d+):(\d+):(\d+)', content)
    if m:
        dur = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
        return dur if dur > 0 else None

    return None


def find_method_logs(method, config):
    """
    Find all log files for a given method and config.
    Returns list of (filepath, job_id, task_id) tuples.
    """
    if not LOG_DIR.exists():
        return []

    results = []

    if method == "alloppnet":
        # Pattern: alloppnet_run_<config>_rep<rep>_<jobid>_<taskid>.out
        #      or: alloppnet_<config>_rep<rep>_<jobid>_<taskid>.out
        for f in LOG_DIR.glob(f"alloppnet*_{config}_rep*_*.out"):
            name = f.stem  # remove .out
            # Extract rep, jobid, taskid from the end
            # alloppnet_run_conf_xxx_rep1_12345_3.out
            # or alloppnet_conf_xxx_rep1_12345_3.out
            parts = name.split('_')
            if len(parts) < 4:
                continue
            try:
                task_id = int(parts[-1])
                job_id = int(parts[-2])
                # Find rep number
                rep = None
                for p in parts:
                    if p.startswith('rep'):
                        rep = int(p[3:])
                        break
                if rep is not None and 1 <= task_id <= 8:
                    # For alloppnet, task_id is 1-8 (8 tetra networks)
                    # Convert to a unique key: (rep-1)*8 + task_id
                    composite_id = (rep - 1) * 8 + task_id
                    results.append((f, job_id, composite_id, rep, task_id))
            except (ValueError, IndexError):
                continue

    elif method == "grandma_split":
        # Pattern: run_grandma_split_<config>_<taskid>.out
        #      or: run_grandma_split_<config>_<jobid>_<taskid>.out
        for f in LOG_DIR.glob(f"run_grandma_split_{config}_*.out"):
            name = f.stem
            suffix = name.replace(f"run_grandma_split_{config}_", "")
            parts = suffix.split('_')
            try:
                if len(parts) == 1:
                    # Just task_id
                    task_id = int(parts[0])
                    job_id = 0
                elif len(parts) == 2:
                    job_id = int(parts[0])
                    task_id = int(parts[1])
                else:
                    continue
                results.append((f, job_id, task_id, None, None))
            except ValueError:
                continue

    elif method == "mpsugar":
        # Pattern: run_mpsugar_<config>_i<iter>_c<chains>_<jobid>_<taskid>.out
        for f in LOG_DIR.glob(f"run_mpsugar_{config}_i*_c*_*.out"):
            name = f.stem
            parts = name.split('_')
            try:
                task_id = int(parts[-1])
                job_id = int(parts[-2])
                results.append((f, job_id, task_id, None, None))
            except (ValueError, IndexError):
                continue

    else:
        # grampa, polyphest, padre: run_<method>_<config>_<jobid>_<taskid>.out
        for f in LOG_DIR.glob(f"run_{method}_{config}_*.out"):
            name = f.stem
            suffix = name.replace(f"run_{method}_{config}_", "")
            parts = suffix.split('_')
            try:
                if len(parts) == 2:
                    job_id = int(parts[0])
                    task_id = int(parts[1])
                elif len(parts) == 1:
                    task_id = int(parts[0])
                    job_id = 0
                else:
                    continue
                results.append((f, job_id, task_id, None, None))
            except ValueError:
                continue

    return results


def collect_durations(method, config):
    """
    Collect durations for a method+config, keeping only the latest job_id per task.
    Returns list of dicts with task info and duration.
    """
    logs = find_method_logs(method, config)
    if not logs:
        return []

    # Group by composite task key, keep highest job_id
    best = {}  # task_key -> (filepath, job_id, ...)
    for entry in logs:
        filepath, job_id, composite_id = entry[0], entry[1], entry[2]
        if composite_id not in best or job_id > best[composite_id][1]:
            best[composite_id] = entry

    results = []
    for composite_id, entry in best.items():
        filepath = entry[0]
        duration = parse_duration_from_log(filepath)
        if duration is None:
            continue

        if method == "alloppnet":
            rep, task_id = entry[3], entry[4]
            if 1 <= task_id <= len(TETRA_NETWORKS):
                network = TETRA_NETWORKS[task_id - 1]
            else:
                network = f"unknown_task{task_id}"
        else:
            network, rep = task_to_network_replicate(composite_id)
            if network is None:
                continue

        results.append({
            "method": METHOD_DISPLAY.get(method, method),
            "config": config,
            "network": network,
            "replicate": rep,
            "duration_sec": duration,
            "duration_hours": round(duration / 3600, 3),
        })

    return results


def format_time(seconds):
    """Format seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}min"
    else:
        return f"{seconds / 3600:.1f}h"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_methods = ["grampa", "grandma_split", "polyphest", "padre", "mpsugar", "alloppnet"]
    all_jobs = []
    summary_rows = []

    print(f"\n{'='*80}")
    print("METHOD RUNTIME ANALYSIS — ALL CONFIGS")
    print(f"{'='*80}\n")

    for method in all_methods:
        display = METHOD_DISPLAY.get(method, method)
        print(f"\n--- {display} ---")

        for config in CONFIGS:
            # AlloppNET runs on all configs (but may have low completion on dup_loss)

            durations = collect_durations(method, config)
            all_jobs.extend(durations)

            if not durations:
                print(f"  {config}: no completed jobs found")
                continue

            times = [d["duration_sec"] for d in durations]
            n = len(times)
            mean_t = statistics.mean(times)
            median_t = statistics.median(times)
            std_t = statistics.stdev(times) if n > 1 else 0
            min_t = min(times)
            max_t = max(times)

            # Expected tasks
            if method == "alloppnet":
                expected = 8 * 5  # 8 networks × 5 replicates
            else:
                expected = 21 * 5  # 21 networks × 5 replicates

            summary_rows.append({
                "method": display,
                "config": config,
                "n_completed": n,
                "n_expected": expected,
                "completion_pct": round(100 * n / expected, 1),
                "mean_sec": round(mean_t, 1),
                "median_sec": round(median_t, 1),
                "std_sec": round(std_t, 1),
                "min_sec": min_t,
                "max_sec": max_t,
                "mean_hours": round(mean_t / 3600, 3),
                "median_hours": round(median_t / 3600, 3),
            })

            print(f"  {config}: n={n}/{expected}, "
                  f"mean={format_time(mean_t)}, median={format_time(median_t)}, "
                  f"range=[{format_time(min_t)}, {format_time(max_t)}]")

    # ── Write per-job CSV ──────────────────────────────────────────────────────
    per_job_path = OUT_DIR / "runtime_per_job.csv"
    if all_jobs:
        with open(per_job_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_jobs[0].keys())
            writer.writeheader()
            writer.writerows(all_jobs)
        print(f"\nPer-job data: {per_job_path} ({len(all_jobs)} rows)")

    # ── Write summary CSV ──────────────────────────────────────────────────────
    summary_path = OUT_DIR / "runtime_summary.csv"
    if summary_rows:
        with open(summary_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"Summary data:  {summary_path} ({len(summary_rows)} rows)")

    # ── Print grand summary table ──────────────────────────────────────────────
    if summary_rows:
        print(f"\n{'='*80}")
        print("GRAND SUMMARY — Mean Runtime per Method (across all configs)")
        print(f"{'='*80}")
        print(f"{'Method':<15} {'Configs':>8} {'Total Jobs':>11} "
              f"{'Mean':>10} {'Median':>10} {'Min':>10} {'Max':>10}")
        print("-" * 80)

        method_groups = defaultdict(list)
        for row in summary_rows:
            method_groups[row["method"]].append(row)

        for method in METHOD_DISPLAY.values():
            if method not in method_groups:
                continue
            rows = method_groups[method]
            total_jobs = sum(r["n_completed"] for r in rows)
            # Weighted mean across configs
            all_times = []
            for r in rows:
                # Approximate: use mean * n to get sum
                all_times.extend([r["mean_sec"]] * r["n_completed"])
            if all_times:
                grand_mean = statistics.mean(all_times)
                grand_median = statistics.median(all_times)
                grand_min = min(r["min_sec"] for r in rows)
                grand_max = max(r["max_sec"] for r in rows)
                print(f"{method:<15} {len(rows):>8} {total_jobs:>11} "
                      f"{format_time(grand_mean):>10} {format_time(grand_median):>10} "
                      f"{format_time(grand_min):>10} {format_time(grand_max):>10}")
            else:
                print(f"{method:<15} {len(rows):>8} {'0':>11} {'—':>10} {'—':>10} {'—':>10} {'—':>10}")

    print(f"\n{'='*80}")
    print("Done.")
    if not all_jobs:
        print("\nNo jobs found. Make sure you're running this on the cluster")
        print(f"where log files exist at: {LOG_DIR}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
