# Thesis experiment plan — GNN polyploid method chapter

Execution roadmap and provenance record for the "detect-then-place" GNN method chapter.
Every number in the chapter must trace to one run listed here, against one checkpoint,
on the clean config set, with one metric definition.

This chapter is a **new method** added to an existing benchmarking thesis. No intro,
no related work, no polyploidy background. Sections: method, diagnostic, ablations,
benchmark, limitations/future work.

---

## Canonical setup (confirm before running anything)

- **Model:** `SpeciesTreeGNNv2`, trained via `scripts/train_reconstruct.py`.
- **Checkpoint (current, rooted):** `output/reconstruct_cladelabels_rooted/best_model.pt`.
- **Config — MUST CONFIRM:** diagnostics load it against `configs/reconstruct_base.yaml`
  (`partner_pair_feat_dim: 2`, co-clustering only, cluster-support OFF). The train
  submit-script default disagrees (`configs/reconstruct.yaml`, dim 4). Read the training
  log for `reconstruct_cladelabels_rooted` and record which config produced it. This
  decides whether "cluster-support off" is a property of the shipped model or an ablation.
- **Shared hyperparameters:** node 13, edge 9, hidden 64, GAT layers 3, heads 4,
  dropout 0.2, lr 1e-3, wd 1e-4, batch 8 (grad-accum), max epochs 200, patience 30,
  pos_class_weight 12, focal alpha 0.25 gamma 2.0, val split 0.2 (seed 42).
- **Metric of record:** `edit_distance_multree` (plus Jaccard ret_leaf / ret_sisters),
  scored by `score_reconstructions.py` exactly as the other benchmark methods were.

### Clean config set (use these ONLY for reported results)
- ILS sweep: `ils_low`, `ils_medium`, `ils_high` — all clean.
- Dup/loss sweep: `dup_loss_low_ne1M`, `dup_loss_medium_ne1M`, `dup_loss_high_ne1M`.
- **Exclude from reported results:** `conf_dup_loss_medium_10M_ne2M` and
  `conf_dup_loss_high_10M_ne2M` (fully wrong-rate), and 2 networks
  (Bendiksby_2011, Marcussen_2011) in `conf_dup_loss_high_10M`. See `project_dup_loss_rate_bug`.

---

## Step 0 — Provenance checks (do first, quick)

1. **Confirm the shipped config** (base dim-2 vs reconstruct dim-4) from the training log.
2. **Verify the TRAINING data rates** are correct — this was assumed, never checked.
   The documented rate bug is about the *benchmark* configs, not necessarily
   `data/mul_trees_2k/training/dup_loss_*`. Measure leaves-per-tree spread or read the
   training data's sim config. If clean, reuse the existing checkpoint for all eval-only work.
3. **Decision:** if training data is clean → reuse existing base checkpoint for Steps 1-2.
   Retraining is only needed for Step 3 regardless.

---

## Step 1 — Eval suite on the existing checkpoint (no training)

All eval-only, point `--model-dir output/reconstruct_cladelabels_rooted`, pass the rooted
clean-config data dirs. Feeds §method / §diagnostic / §ablations.

| Experiment | Command | Output | Chapter |
|---|---|---|---|
| WGD F1 / threshold sweep | `python scripts/tune_threshold.py --data-dir <clean rooted dirs> --model-dir <ckpt> --config configs/reconstruct_base.yaml` | P/R/F1 + `threshold_tuning.txt` | detection |
| Single-feature baseline | `python scripts/feature_baseline.py --data-dir <dirs> --feature frac_clade_dup` | best single-feature F1 | detection |
| Permutation importance | `python scripts/permutation_importance.py --data-dir <dirs> --model-dir <ckpt> --features both --threshold 0.9` | per-feature F1 drop | detection/partner |
| False-positive analysis | `python scripts/fp_analysis.py --data-dir <dirs> --model-dir <ckpt> --threshold 0.9` | FP clade-size histogram | detection |
| Auto/allo partner accuracy | `python scripts/partner_diagnostic.py --data-dir <dirs> --model-dir <ckpt>` | auto/allo acc, near/far split | partner |
| Node-feature permutation (partner) | add `--permute-node-features` to above | allo acc drop | partner |
| Per-config error budget | `python scripts/error_budget.py --model-dir <ckpt> --model-config configs/reconstruct_base.yaml --configs ils_low ils_medium ils_high dup_loss_low_ne1M dup_loss_medium_ne1M dup_loss_high_ne1M` | per-config prec/rec/auto/allo | results |
| Event-selection strategies | `python scripts/benchmark_networks.py --model-dir <ckpt> --config conf_ils_low_10M --strategies raw,collapse,cap,collapse_cap,bound_driven --threshold 0.9` then `score_reconstructions.py` | per-strategy edit | build |
| Rooting accuracy vs known root | `python scripts/root_from_gene_trees.py --mul-trees-dir data/mul_trees_2k --config ils_low --n 200` | recovery rates | build |
| Event-drop rate | `python scripts/event_drop_rate.py --data-dir <dirs>` | per-config drop rate | build |
| Backbone RF localization (full vs pruned) | `python scripts/backbone_polyploid_localization.py --config ils_low --max-samples 300` | **record the actual full/pruned RF means, not just the 86% headline** | diagnostic |
| Polyploid placement accuracy | `python scripts/backbone_polyploid_placement_accuracy.py --config ils_low --threshold 0.5` | one-parent-hit rate | diagnostic |
| Gene-trees-per-sample saturation | subsample gene trees at inference (100/250/500/1000), eval — cheap, eval-only | saturation curve | diagnostic |

---

## Step 2 — The decisive oracle (single run, gates phasing)

Does perfect placement on the ASTRAL backbone reach the true-backbone floor?

```
sbatch jobs/two_parent_oracle_job.sh          # runs the array incl. astral_ceiling_rooted
# or directly:
python scripts/two_parent_oracle.py --backbone astral --parents true --root-backbone hybrid --build graft ...
python scripts/score_reconstructions.py --recon-dir <oracle out> ; python scripts/summarize_oracle.py
```

Also run the true-backbone cell (`--backbone true --parents true`) and the edit-gap
decomposition (new script `decompose_edit_gap.py`, below) to get the numeric split.

- If ASTRAL+true-parents approaches the true-backbone floor (~0.30) → bounded phasing viable → Step 5.
- If it stays ~0.68 → diploid skeleton is the wall → phasing is future work, documented with proof.

---

## Step 3 — Retrain batch (only compute-heavy part)

Batch these trainings. All use `jobs/submit_train_reconstruct.sh ... --rooted --clade-labels`
restricted to the clean config set.

- **Cluster-support ablation:** `reconstruct.yaml` (dim 4, on) vs `reconstruct_base.yaml` (dim 2, off).
- **Pairwise-off ablation:** new `configs/reconstruct_nopair.yaml` (`partner_pair_feat_dim: 0`, below).
- **Hyperparameter tuning:** sweep lr / hidden / GAT layers / dropout / class weight; pick by val.
  Lock final numbers to the tuned config.
- **Training-size saturation:** retrain on 250 / 500 / 1000 / 2000 samples; plot val vs size.
- **Leave-one-config-out CV:** train on all clean configs but one, eval on the held-out one.
  Light version: leave out only high-ILS.
- **(Optional) two-parent phaser re-bench:** `reconstruct_two_parent.yaml`, warm-start from the
  one-parter; prereq `relabel_from_metadata.py --two-parent`.

The three "rule out the boring explanation" checks: ablations (not underfit),
saturation (not data-starved), LOCO-CV (not regime-overfit). Together they isolate the backbone.

---

## Step 4 — Writing (parallel, starts now)

Draft order: method (stable now) → diagnostic §2 (as Step-1 numbers land) → ablations §3 →
benchmark §4 → limitations/future-work §5. Prose style: no semicolons, no non-math
parentheses, no em-dashes.

---

## Step 5 — Bounded phasing (conditional on Step 2)

Build only if Step 2 says a phased/diploid backbone can beat the wall. It slots in at the
backbone stage: build the diploid-only backbone, phase each polyploid into subgenomes,
attach each to its two parents. Detection stays with the GNN. Otherwise: future-work section
backed by the Step-2 diagnostic.

---

## Missing scripts (being written)

1. `scripts/decompose_edit_gap.py` — writes the numeric ASTRAL-vs-true edit-gap decomposition
   by event type (auto / single-allo / clade-allo). Adapts `localize_faithfulness.py`.
2. `scripts/backbone_displacement.py` — per-polyploid topological displacement (ASTRAL vs true
   position) classified local (NNI) vs long-range. Adapts `backbone_polyploid_placement_accuracy.py`.
3. `scripts/head_to_head.py` — joins the GNN's `reconstruction_scores.csv` with the
   Polyphest / GRAMPA scored outputs into one comparison table per config.
4. `configs/reconstruct_nopair.yaml` — `partner_pair_feat_dim: 0` for the pairwise-off ablation.

---

## What feeds which chapter section

- **§method:** architecture, features, training (stable, write now).
- **§diagnostic:** Step-1 backbone localization + Step-2 oracle + `decompose_edit_gap.py`
  + `backbone_displacement.py`. The three-walls story with corrected numbers.
- **§ablations:** Step-1 permutation + Step-3 feature ablations + saturation + LOCO-CV.
- **§benchmark:** Step-1 benchmark + `head_to_head.py` vs Polyphest / GRAMPA-iter, clean configs only.
- **§limitations/future-work:** one-parent limitation, backbone ceiling, bounded phasing.
