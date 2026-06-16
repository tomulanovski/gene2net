"""Score inferred MUL-trees against ground truth, reusing the existing
comparison core (reticulate_tree + compare_reticulations) so the numbers are
directly comparable to how GRAMPA/Polyphest were scored.

Expects the layout produced by reconstruct_mul_tree.py:
    <recon-dir>/sample_NNNN/output.tre        (inferred MUL-tree)
    <recon-dir>/sample_NNNN/ground_truth.nex  (true MUL-tree)

Each sample is scored independently, so this parallelises across CPU cores
(--workers). The per-sample cost is the graph edit distance / network
conversion, which dominates; running many samples at once is the speedup.

Run in the env where the comparison code imports work (the `gene2net` env),
and on a compute node (not the login node) when using many workers.

Usage:
    python scripts/score_reconstructions.py \
        --recon-dir output/reconstruct_allo/mul_trees/ils_low_t09 --workers 8
"""
import argparse
import multiprocessing as mp
import os
import sys

import pandas as pd

METRICS = ["edit_distance", "ret_leaf_jaccard", "ret_sisters_jaccard",
           "num_rets_diff", "ploidy_diff"]

# Set by each worker (and the parent) so the comparison imports resolve.
_SIM_SCRIPTS = None


def _ensure_imports(sim_scripts):
    """Import the comparison core, adding sim_scripts to sys.path once."""
    global _SIM_SCRIPTS
    if sim_scripts and sim_scripts not in sys.path:
        sys.path.insert(0, sim_scripts)
    _SIM_SCRIPTS = sim_scripts
    from reticulate_tree import ReticulateTree
    from compare_reticulations import pairwise_comparison
    return ReticulateTree, pairwise_comparison


def newick_from_file(path):
    """Read a Newick string, extracting it from NEXUS if needed."""
    txt = open(path).read().strip()
    low = txt.lower()
    if low.startswith("#nexus") or "begin trees" in low:
        for line in txt.splitlines():
            s = line.strip()
            if s.lower().startswith("tree") and "=" in s:
                return s.split("=", 1)[1].strip()
    return txt


def build_rt(tree_str, ReticulateTree):
    """Mirror run_comparison_analysis: apply MUL-tree conversion when duplicated."""
    rt_temp = ReticulateTree(tree_str)
    if rt_temp.check_duplicated():
        return ReticulateTree(tree_str, is_multree=True)
    return rt_temp


def score_one(task):
    """Score a single sample dir. task = (name, sample_dir, sim_scripts)."""
    name, sdir, sim_scripts = task
    inf_path = os.path.join(sdir, "output.tre")
    gt_path = os.path.join(sdir, "ground_truth.nex")
    if not (os.path.exists(inf_path) and os.path.exists(gt_path)):
        return None
    try:
        ReticulateTree, pairwise_comparison = _ensure_imports(sim_scripts)
        inf = newick_from_file(inf_path)
        gt = newick_from_file(gt_path)
        rt_gt = build_rt(gt, ReticulateTree)
        rt_inf = build_rt(inf, ReticulateTree)
        df = pd.DataFrame([
            {"name": "ground_truth", "object": rt_gt, **rt_gt.measure(printout=False)},
            {"name": "gene2net", "object": rt_inf, **rt_inf.measure(printout=False)},
        ]).set_index("name")
        comp = pairwise_comparison(df, debug=False)
        row = {"sample": name}
        for m in METRICS:
            try:
                row[m] = float(comp[m].loc["ground_truth", "gene2net"])
            except Exception:
                row[m] = float("nan")
        return row
    except Exception as e:
        return {"sample": name, "error": f"{type(e).__name__}: {e}"}


def _score_worker(task, q):
    q.put(score_one(task))


def score_isolated(task, timeout):
    """Score one sample in a disposable process with a timeout.

    graph_edit_distance can hang or OOM on malformed networks (>2-parent
    reticulations from over-ploidy); isolating each sample means one blow-up is
    skipped instead of killing the whole run.
    """
    q = mp.Queue()
    p = mp.Process(target=_score_worker, args=(task, q))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        return {"sample": task[0], "error": "timeout"}
    try:
        return q.get_nowait()
    except Exception:
        return {"sample": task[0], "error": "crashed (worker died)"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recon-dir", required=True)
    parser.add_argument("--sim-scripts",
                        default="/groups/itay_mayrose/tomulanovski/gene2net/simulations/scripts")
    parser.add_argument("--workers", type=int, default=4)  # kept for CLI compat; scoring is isolated/sequential
    parser.add_argument("--timeout", type=int, default=120,
                        help="per-sample seconds before skipping a hung/blown-up comparison")
    parser.add_argument("--out-csv", default=None)
    args = parser.parse_args()

    # Any subdir that contains an output.tre (works for sample_NNNN and the
    # network-named dirs from benchmark_networks.py).
    sample_dirs = sorted(
        d for d in os.listdir(args.recon_dir)
        if os.path.isdir(os.path.join(args.recon_dir, d))
        and os.path.exists(os.path.join(args.recon_dir, d, "output.tre"))
    )
    tasks = [(d, os.path.join(args.recon_dir, d), args.sim_scripts) for d in sample_dirs]
    print(f"Scoring {len(tasks)} samples (isolated, timeout {args.timeout}s)...")

    rows = []
    for i, t in enumerate(tasks, 1):
        r = score_isolated(t, args.timeout)
        if r:
            rows.append(r)
        if i % 10 == 0:
            print(f"  {i}/{len(tasks)} done")

    ok = [r for r in rows if "error" not in r]
    errs = [r for r in rows if "error" in r]
    if errs:
        print(f"\n{len(errs)} samples errored (showing up to 5):")
        for r in errs[:5]:
            print(f"  {r['sample']}: {r['error']}")

    if not ok:
        print("No samples scored — check the recon-dir layout.")
        return

    df = pd.DataFrame(ok)
    print(f"\nScored {len(df)} samples from {args.recon_dir}\n")
    print("Per-metric summary (distances — lower is better):")
    print(df[METRICS].describe().loc[["mean", "50%", "min", "max"]].to_string())

    out_csv = args.out_csv or os.path.join(args.recon_dir, "reconstruction_scores.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nSaved per-sample scores to {out_csv}")


if __name__ == "__main__":
    main()
