#!/usr/bin/env python3
"""
diagnose_failed_runs.py - Why did GRAMPA^Iter / Polyphest runs not finish?

For every failed run in the summary inventories (inferred_exists == False), finds the
matching SLURM log by parsing the log header (Network / Replicate / Configuration /
Percentile) and classifies the failure from the log tails:

    TIMEOUT        slurm "CANCELLED ... DUE TO TIME LIMIT"
    OOM            out-of-memory kill
    CANCELLED      cancelled for another reason (e.g. scancel, node failure)
    PYTHON_ERROR   Python traceback (last line shown)
    NO_OUTPUT      job reported success but the output file is missing
    NO_LOG         no log file found for this run (e.g. overwritten)
    UNKNOWN        log found but none of the above matched (last stderr line shown)

Note: run_grandma_split.sh logs are named run_grandma_split_<task>.out with no job id,
so each config's submission overwrites the previous one. Failed runs from configs whose
logs were overwritten show up as NO_LOG (a log only matches if its header config matches).

Usage:
    python diagnose_failed_runs.py
    python diagnose_failed_runs.py --configs conf_ils_low_10M conf_ils_high_10M
    python diagnose_failed_runs.py --csv failed_runs.csv
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_BASE = SCRIPT_DIR.parent / "analysis" / "summary"
LOG_DIR = SCRIPT_DIR.parent / "logs"

POLYPHEST_THRESHOLDS = ['polyphest_p50', 'polyphest_p70', 'polyphest_p90']

# method -> glob for its log files (the digit keeps prior / real variants out)
LOG_GLOBS = {
    'grandma_split': 'run_grandma_split_[0-9]*.out',
    'polyphest': 'run_polyphest_[0-9]*.out',
}

HEADER_RE = {
    'network': re.compile(r'^Network:\s*(\S+)'),
    'replicate': re.compile(r'^Replicate:\s*(\d+)'),
    'config': re.compile(r'^Configuration:\s*(\S+)'),
    'percentile': re.compile(r'^\s*Percentile:\s*(\d+)'),
}

TAIL_BYTES = 64 * 1024


def read_head(path, n_lines=40):
    lines = []
    with open(path, errors='replace') as fh:
        for _ in range(n_lines):
            line = fh.readline()
            if not line:
                break
            lines.append(line.rstrip('\n'))
    return lines


def read_tail(path):
    if not path.exists():
        return ''
    with open(path, 'rb') as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - TAIL_BYTES))
        return fh.read().decode(errors='replace')


def index_logs(method):
    """Map (config, network, replicate, percentile) -> list of .out log paths."""
    index = defaultdict(list)
    for out in LOG_DIR.glob(LOG_GLOBS[method]):
        header = {}
        for line in read_head(out):
            for key, rx in HEADER_RE.items():
                if key not in header:
                    m = rx.match(line)
                    if m:
                        header[key] = m.group(1)
        if not {'network', 'replicate', 'config'} <= header.keys():
            continue
        pct = header.get('percentile') if method == 'polyphest' else None
        key = (header['config'], header['network'], int(header['replicate']), pct)
        index[key].append(out)
    return index


def classify(out_path):
    err_path = out_path.with_suffix('.err')
    text = read_tail(out_path) + '\n' + read_tail(err_path)
    if re.search(r'DUE TO TIME LIMIT', text):
        return 'TIMEOUT', ''
    if re.search(r'oom[-_ ]kill|Out Of Memory|OUT_OF_MEMORY', text, re.IGNORECASE):
        return 'OOM', ''
    if re.search(r'CANCELLED', text):
        m = re.search(r'.*CANCELLED.*', text)
        return 'CANCELLED', m.group(0).strip()
    if 'Traceback (most recent call last)' in text:
        tb = text[text.rfind('Traceback (most recent call last)'):]
        last = [l for l in tb.splitlines() if l.strip()][-1]
        return 'PYTHON_ERROR', last.strip()
    if 'COMPLETED SUCCESSFULLY' in text:
        return 'NO_OUTPUT', ''
    err_lines = [l for l in read_tail(err_path).splitlines() if l.strip()]
    return 'UNKNOWN', err_lines[-1].strip() if err_lines else '(empty stderr)'


def failed_runs(config):
    inv_file = SUMMARY_BASE / config / "inventory.csv"
    inv = pd.read_csv(inv_file)
    inv['inferred_exists'] = inv['inferred_exists'].astype(bool)
    rows = []

    g = inv[(inv['method'] == 'grandma_split') & ~inv['inferred_exists']]
    for _, r in g.iterrows():
        rows.append(('grandma_split', r['network'], int(r['replicate']), None))

    # Polyphest counts as failed only if no threshold finished; report every threshold.
    poly = inv[inv['method'].isin(POLYPHEST_THRESHOLDS)]
    for (net, rep), grp in poly.groupby(['network', 'replicate']):
        if grp['inferred_exists'].any():
            continue
        for method in grp['method']:
            rows.append(('polyphest', net, int(rep), method.split('_p')[-1]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--configs', nargs='+', default=None,
                        help='Configs to check (default: every summary dir with inventory.csv)')
    parser.add_argument('--csv', default=None, help='Also write the per-run table to this CSV')
    args = parser.parse_args()

    if not LOG_DIR.is_dir():
        sys.exit(f"ERROR: log directory not found: {LOG_DIR}")

    if args.configs:
        configs = args.configs
        for c in configs:
            if not (SUMMARY_BASE / c / "inventory.csv").exists():
                sys.exit(f"ERROR: no inventory.csv for {c} under {SUMMARY_BASE}")
    else:
        configs = sorted(p.parent.name for p in SUMMARY_BASE.glob("conf_*/inventory.csv"))

    print(f"Indexing logs in {LOG_DIR} ...")
    log_index = {m: index_logs(m) for m in LOG_GLOBS}
    for m, idx in log_index.items():
        print(f"  {m}: {sum(len(v) for v in idx.values())} logs with a parseable header")

    records = []
    for config in configs:
        for method, net, rep, pct in failed_runs(config):
            logs = log_index[method].get((config, net, rep, pct), [])
            if not logs:
                reason, detail, log_name = 'NO_LOG', '', ''
            else:
                newest = max(logs, key=lambda p: p.stat().st_mtime)
                reason, detail = classify(newest)
                log_name = newest.name + (f' (+{len(logs) - 1} older)' if len(logs) > 1 else '')
            records.append({
                'config': config, 'method': method,
                'percentile': pct or '', 'network': net, 'replicate': rep,
                'reason': reason, 'detail': detail, 'log': log_name,
            })

    if not records:
        print("\nNo failed GRAMPA^Iter or Polyphest runs found.")
        return

    df = pd.DataFrame(records)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 90)
    pd.set_option('display.max_rows', None)

    print("\n=== Reason counts ===")
    print(df.groupby(['method', 'reason']).size().to_string())

    print("\n=== Reason by network ===")
    print(df.groupby(['method', 'network', 'reason']).size().to_string())

    print("\n=== All failed runs ===")
    print(df.to_string(index=False))

    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nWrote {args.csv}")


if __name__ == '__main__':
    main()
