# Publishability checklist + the backbone go/no-go

Two tracks. Part A is the one shot at beating Polyphest (the backbone). Part B is the
consolidation/publishability checklist that is worth doing regardless of A's outcome —
it fills both the thesis results chapter and the paper.

Prose in the thesis follows the style rule (no semicolons, non-math parens, em-dashes);
this doc is working notes.

---

## Part A — the backbone go/no-go (spend the optimism here first)

The true-backbone floor (~0.11 sim) is well below Polyphest (~0.42), so a correct backbone
WOULD beat it. The only lever is placing polyploids right on the (good) diploid scaffold.
The earlier diploid oracle (0.642) was a BROKEN assembly (per-species grafting flattened
clades), now fixed to graft polyploid clades as whole subtrees (`build_diploid_skeleton_backbone`
v2, clade-preserving).

**Run the clean oracle** (true parents = the ceiling of a diploid-scaffold backbone):
```
python scripts/two_parent_oracle.py --mul-trees-dir data/mul_trees_2k --config ils_low --n 200 \
    --backbone diploid --parents true --root-backbone hybrid \
    --out-dir output/two_parent_oracle/ils_low_diploid_v2
python scripts/score_reconstructions.py --recon-dir output/two_parent_oracle/ils_low_diploid_v2
python scripts/summarize_oracle.py --scores output/two_parent_oracle/ils_low_diploid_v2/reconstruction_scores.csv
```
Smoke-test with `--n 20` first (tree surgery, may need iteration).

**Read against:** true backbone floor 0.114, ASTRAL ceiling 0.543 (both ils_low, same metric).
- Near the floor (~0.11-0.25) → a realistic diploid scaffold + correct placement recovers the
  edit. GO: build the realistic version (diploid-only ASTRAL + co-cluster parents, `--parents
  coclust`) and benchmark it vs Polyphest across configs. This is the score lever.
- Stuck near 0.54 → the diploid scaffold alone is not enough; NO-GO for a cheap fix, phasing
  is the (RA) project.
- In between → partial; likely narrows the gap, may win at hard configs where Polyphest's
  ploidy heuristic fails, unlikely to beat it at low ILS.

If GO: the realistic variant is the "poor man's phasing" pilot and a strong thesis result
either way.

---

## Part B — publishability checklist (do regardless of A)

Ordered by value-per-cost. Items also populate the thesis results/ablation chapter.

### B0 — free (reads existing logs)
- **Overfitting check:** parse the training slurm logs, plot train vs val loss/F1 over epochs.
  Confirms the dropout 0.2 + weight-decay 1e-4 regularization is adequate (small model, ~82k
  params, ~14k samples -> expected fine, but must show it).

### B1 — eval-only (fast, headline results)
- **Head-to-head vs Polyphest + GRAMPA-iter, all 6 clean configs.** The benchmark table.
  ```
  for c in conf_ils_low_10M conf_ils_medium_10M conf_ils_high_10M \
           conf_dup_loss_low_10M_ne1M conf_dup_loss_medium_10M_ne1M conf_dup_loss_high_10M_ne1M; do
    python scripts/benchmark_networks.py --model-dir output/reconstruct_cladelabels_rooted \
        --config $c --threshold 0.9 --strategies cap --root-mode true --build-order confidence \
        --out-base output/bench_final/$c    # (final env for benchmark; score in gene2net)
    python scripts/score_reconstructions.py --recon-dir output/bench_final/$c/cap
    python scripts/head_to_head.py --gnn-scores output/bench_final/$c/cap --config $c \
        --out output/bench_final/$c/head_to_head.csv
  done
  ```
- **Runtime comparison:** time the end-to-end pipeline (ASTRAL + GNN + build) vs Polyphest /
  GRAMPA-iter on a few networks. The likely unclaimed win; motivates the fast-prior-free framing.

### B2 — retrain batch (the ablations; run together)
Needs two small code additions first (below):
- **Learning curve (data sufficiency):** retrain on 250 / 500 / 1000 / 2000 networks, plot val
  F1 / partner-acc vs size. Answers "did 2000 converge or need more/less." The #1 reviewer ask.
- **GAT-depth ablation:** train n_gat_layers 1 / 2 / 3 / 4, compare F1/partner. Closes the
  "why 3 layers" hole (currently never ablated).
- **Feature ablation:** drop the detection-dead features (`branch_length`, `depth`,
  `dup_synchrony`); pairwise-off (`configs/reconstruct_nopair.yaml`, exists); add `n_eff`
  (effective-parent-count node feature) and test vs the 5 co-clustering moments.
- **Threshold on a held-out split** (not the benchmark) for an honest operating point.

### Code additions needed for B2 (small, write before the retrain batch)
1. `--max-train-samples N` flag on `train_reconstruct.py` (subsample training set for the learning curve).
2. `n_eff` node feature in `features.py` = effective number of parents from the co-clustering row
   = (sum p)^2 / sum(p^2); add to the node feature vector (and a config to toggle it for the ablation).

---

## Publishability verdict

Not "new SOTA" (does not beat Polyphest broadly). IS a solid methods+analysis paper framed as:
a fast, prior-free GNN for polyploid networks + a rigorous decomposition proving the single-copy
backbone (not prediction) is the fundamental limit of detect-then-place, motivating phasing.
Beats the fair prior-free peer (GRAMPA-iter); honest negatives (two-parent, build fix, metric
audit) are a strength. Target: Bioinformatics / BMC / Syst. Biol. methods. Part A, if GO, raises
the ceiling toward a stronger venue.
