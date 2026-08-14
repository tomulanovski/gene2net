# ml_method codebase audit & cleanup plan

Purpose: state the ONE canonical pipeline, fix two latent bugs, and propose what to remove so
the repo is thesis-grade. Nothing is deleted until each section below is approved. The method is:
**one-parter GNN, away labels, always hybrid-rooted ASTRAL, bound_driven decode (cap kept as a
secondary knob), final config `reconstruct_final.yaml`.**

---

## 1. The canonical pipeline — the files that ARE the method

Data prep → train → select → decode → score → compare. Keep all of these.

**Library (`gene2net_gnn/`)**
- `data/tree_io.py` — Newick/ETE3 ↔ edge_index.
- `data/dataset.py` — `Gene2NetSample` / `Gene2NetDataset`; loads away/clade label sidecars.
- `data/features.py` — the 13 node + 9 edge + 4 pairwise features.
- `data/label_extractor.py` — `decompose_mul_tree` → events on ASTRAL edges. **(carries the sp08 faithfulness bug — see §4)**
- `data/metadata_labels.py` — clade labels + away-relabel.
- `data/rooting.py` — `hybrid_root` (the rooting we always use).
- `data/mul_tree_generator.py` — synthetic network generation (data gen only).
- `model/species_gnn_v2.py` — `SpeciesTreeGNNv2` (the model).
- `training/trainer_reconstruct.py` — the trainer.
- `training/trainer_phase1.py` — dependency (provides `prepare_sample`, `focal_loss`).
- `inference/mul_tree_builder.py` — `build_mul_tree` (grafts events → MUL-tree). **(sp08 bug lives here too)**
- `inference/build_strategies.py` — `bound_driven` + `infer_copy_bound`.

**Scripts**
- Data gen: `generate_mul_trees.py`, `generate_one_example.py`, `package_training_data.py`, `augment_edge_features.py`, `relabel_from_metadata.py`.
- Train/select: `train_reconstruct.py`, `optuna_reconstruct.py`, `score_val_reconstruction.py`.
- Decode/eval: `benchmark_networks.py`, `reconstruct_mul_tree.py`.
- Score/compare: `score_reconstructions.py`, `head_to_head.py`.
- Active analysis (keep): `oracle_test.py`, `compare_oracle_faithfulness.py`, `partner_diagnostic.py`, `permutation_importance.py`, `partner_permutation_importance.py`, `summarize_oracle.py`, `localize_faithfulness.py`.

**Configs**: `reconstruct_final.yaml` (final), `reconstruct_oneparter.yaml` (HPO base). 
**Jobs**: the data-gen chain, `train_reconstruct_job.sh` / `submit_finalize.sh`, `relabel_away_job.sh`, `optuna_array.sh`, `run_rebenchmark_final.sh`, `score_final.sh`.

---

## 2. Two latent bugs to fix first (not cleanup — real footguns)

1. **Default config is the wrong model.** `train_reconstruct.py`, `benchmark_networks.py`, `reconstruct_mul_tree.py` all default `--config`/`--model-config` to `reconstruct.yaml` = hidden-64 / depth-3 / **implicit `n_parents=2` (two-parent)**. Any run without an explicit config silently uses the wrong architecture.
   - Fix: default to `reconstruct_final.yaml` (and delete `reconstruct.yaml`), or make the config required.
2. **Default rooting is `none`.** Both decode scripts default `--root-mode none` (arbitrary ASTRAL root) despite "always hybrid." This is the exact bug that cost a day on the benchmark.
   - Fix: default `--root-mode hybrid`; drop the redundant `--root-species-tree` alias.

---

## 3. Cleanup plan (proposal — approve per group before removal)

### 3a. Delete — dead / deprecated (safe, high confidence)
- `scripts/rescore_canonical_edit.py` — self-labeled DEPRECATED, superseded by `score_reconstructions.py`.
- `gene2net_gnn/data/collate.py` — only importer is its own test.

### 3b. Delete — superseded whole branches (confirm each is abandoned)
Each is a self-contained cluster reachable only from its own entry point; removing it also removes its tests/config/job.
- **Legacy `Gene2NetModel` lineage:** `model/gene2net_model.py`, `model/species_gnn.py`, `model/heads.py`, `model/aggregator.py`, `model/tree_encoder.py`, `inference/predict.py`, `training/trainer.py`, `training/loss.py`, `training/metrics.py`; scripts `train.py`, `evaluate.py`, `infer.py`; job `train_job.sh`, `submit_train.sh`; config `default.yaml`; tests `test_model_forward/test_aggregator/test_tree_encoder/test_loss`.
- **Phase-1 detection branch:** `scripts/train_phase1.py`, `compute_phase2_features.py`, `compute_phase3_features.py`; `data/gene_tree_features.py`, `data/gene_tree_topology.py`; configs `phase1*.yaml`, `phase2.yaml`, `phase3.yaml`; jobs `train_phase1_job.sh`, `submit_train_phase1/2/3.sh`, `compute_phase2_job.sh`, `submit_compute_phase2.sh`.
- **WGD-detector branch:** `model/wgd_detector.py`, `model/gene_tree_gin.py`, `training/trainer_wgd.py`; `scripts/train_wgd.py`, `tune_threshold.py`; config `wgd_detector.yaml`; jobs `train_wgd_job.sh`, `submit_train_wgd.sh`.

### 3c. Refactor core files — strip the dead two-parent subsystem (careful; do with tests)
The one-parter is final; the two-parent path is dead but is *interwoven* in core files, so this is edits, not deletes:
- `species_gnn_v2.py`: drop `n_parents≥2` head slots and the unused full-matrix `compute_partner_scores` (only `compute_partner_scores_rows` is used).
- `trainer_reconstruct.py`: drop the `two_parent_loss` branch and its eval.
- `metadata_labels.py`: drop `home_edges` / `two_parent_labels_from_metadata` / `home_edge_for_event`.
- `mul_tree_builder.py`: drop `build_mul_tree_two_parent`, `TwoParentEvent`, both graft modes.
- `benchmark_networks.py`: drop `--parents coclust` and the two-parent build path.
- Remove tests `test_two_parent_*`, `test_load_model_infers_head`; keep `conv_types`, `pick_objective`, `train_subset`, `optuna_config`, `away_parent`, `mul_tree_builder`, `rooting`, etc.
- Also delete `conv_type` gin/gcn branches only if you don't want them reported; otherwise keep (they were part of the honest HPO). `reconstruct_two_parent.yaml`, `reconstruct_neff/dupcond/dupcond_neff/nopair/base.yaml` and the `abl_*` configs go.

### 3d. Trim decode strategies
Keep `bound_driven` + `cap`. Remove `raw`, `collapse`, `collapse_cap`, `dedup`, `dedup_cap`, `clade_collapse`, `clade_collapse_cap` from `build_strategies.py`, and set `benchmark_networks.py` default `--strategies` to `bound_driven`.

### 3e. Archive (don't delete) the diagnostic one-offs
~30 scripts (`backbone_*`, `partner_confuser/contrast/ks_signal/top2`, `phaser_baseline`, `debug_*`, `check_*`, `diagnose_*`, `error_*`, `*_audit`, `event_drop_rate`, `decompose_edit_gap`, `feature_baseline`, `baseline_test`, `smoke_new_features`, `profile_features`, `root_from_gene_trees`, `verify_away_relabel`, `training_data_stats`, `metadata_event_stats`, `two_parent_oracle`). Each documents a concluded experiment behind a MEMORY.md finding. Move to `ml_method/scripts/archive/` so the repo is clean but the record survives.

### 3f. Trim configs
From 20 → 2 (`reconstruct_final.yaml`, `reconstruct_oneparter.yaml`). The rest are ablation/branch leftovers removed by 3b/3c.

---

## 4. Error map — where the error is and what's improvable

| source | evidence | improvable? |
|---|---|---|
| **Backbone (ASTRAL)** | oracle-on-ASTRAL ceiling; RF error at polyploids | Only via **phasing** — a different method. Future work. |
| **Build/label faithfulness** | oracle with TRUE events on TRUE backbone still scores edit 0.33 / ret_sister 0.34; sp08 reticulation gets right descendants, wrong sisters | **Yes — the one real bug-fix lever.** `decompose_mul_tree` + `build_mul_tree` grafting convention. Fixing it cleans ret_sister and un-corrupts allo labels. This is the focused look. |
| **Copy number (ploidy)** | ploidy_diff ~0.08–0.11; edit is copy-number-heavy | Inferred-bound heuristic; force-fill rejected. Mostly inherent to the decode. |
| **Detection** | F1 ~0.87; over-produces events (num_rets_diff) | Marginal; decode-side. |
| **Placement (partner)** | proven near-exhausted: 2× val partner → ~2% reconstruction | Not a lever. |

**Conclusion:** the single genuine improvement lever is **build/label faithfulness** (§4 row 2), which is also what corrupts the ret_sister metric. Everything else is inherent (copy number), structural (backbone → phasing), or exhausted (placement).

---

## 5. Suggested execution order

1. Fix the two footguns (§2) — they're bugs.
2. Delete 3a (dead/deprecated).
3. Delete 3b branches you confirm abandoned.
4. Archive 3e diagnostics.
5. Trim 3d strategies + 3f configs.
6. Refactor 3c (two-parent strip) — carefully, run the kept tests after.
7. Then the `build_mul_tree` / `decompose_mul_tree` faithfulness fix (§4) — the real improvement.

Each step is its own reviewable commit; nothing is removed until you approve that step.
