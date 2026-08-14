# Archived branches

These are earlier approaches that were superseded by the current method and moved
out of the working tree on 2026-08-14. Nothing in the live codebase imports any of
this (verified before archiving). Files keep their original relative paths, so a
branch can be revived by moving its subtree back under `gene2net_gnn/` and `scripts/`.

The live method is the reconstruction pipeline: `model/species_gnn_v2.py` +
`training/trainer_reconstruct.py`, config `configs/reconstruct_final.yaml`, decoded
through `scripts/benchmark_networks.py`. None of the branches below are part of it.

---

## A. Legacy Gene2NetModel (gene-tree encoding)

The first design. Instead of hand-crafting features on the ASTRAL species tree, it
**encoded the gene trees directly**: a per-tree `GeneTreeEncoder` fed a cross-tree
attention `CrossTreeAggregator`, combined with a species-tree GNN, into WGD/partner
heads. Superseded by the features-on-species-tree approach (`species_gnn_v2`), which
was simpler and scored better.

Kept because we may revisit learned gene-tree encoding instead of hand-crafted
features. This subtree is the record of what that first attempt looked like.

Files:
- `gene2net_gnn/model/`: `gene2net_model.py` (top-level model), `tree_encoder.py`
  (GeneTreeEncoder), `aggregator.py` (CrossTreeAggregator), `species_gnn.py` (v1
  species GNN), `heads.py` (WGDHead, PartnerHead)
- `gene2net_gnn/training/`: `trainer.py` (Trainer), `loss.py` (Gene2NetLoss,
  focal_loss), `metrics.py` (wgd_f1, partner_accuracy, event_count_error)
- `gene2net_gnn/inference/predict.py` (predict_mul_tree; also held the two-parent
  decode helper `events_from_two_parent_scores`)
- `scripts/`: `train.py`, `evaluate.py`, `infer.py`
- `configs/default.yaml`
- `jobs/train_job.sh`, `jobs/submit_train.sh`
- `tests/`: `test_model_forward.py`, `test_tree_encoder.py`, `test_aggregator.py`,
  `test_species_gnn.py`, `test_loss.py`, `test_two_parent_decode.py`

To revive: move `gene2net_gnn/*` and `scripts/*` back to their original locations.
The internal imports (`from gene2net_gnn.model.tree_encoder import ...`) will then
resolve again.

## B. Phase-1/2/3 detection

A stepping stone before the joint reconstruction model: first solve WGD detection
alone (Phase 1), then experiments piling on more gene-tree-derived features
(Phase 2/3). Superseded by `trainer_reconstruct`, which does detection and partner
prediction jointly.

Note: `gene2net_gnn/training/trainer_phase1.py` did NOT move — the live
reconstruction trainer still borrows `prepare_sample` and `focal_loss` from it.

Files:
- `gene2net_gnn/data/`: `gene_tree_features.py`, `gene_tree_topology.py`
- `scripts/`: `train_phase1.py`, `compute_phase2_features.py`,
  `compute_phase3_features.py`
- `configs/`: `phase1.yaml`, `phase1_cw6.yaml`, `phase2.yaml`, `phase3.yaml`
- `jobs/`: `train_phase1_job.sh`, `submit_train_phase1.sh`, `compute_phase2_job.sh`,
  `submit_compute_phase2.sh`, `submit_train_phase2.sh`, `submit_train_phase3.sh`

## C. WGD-detector

An alternative detector architecture: a separate GIN (`gene_tree_gin.py`) + GAT
model (`wgd_detector.py`) for the detection task. Replaced by the detection head
inside `species_gnn_v2`.

Files:
- `gene2net_gnn/model/`: `wgd_detector.py`, `gene_tree_gin.py`
- `gene2net_gnn/training/trainer_wgd.py`
- `scripts/`: `train_wgd.py`, `tune_threshold.py`
- `configs/wgd_detector.yaml`
- `jobs/train_wgd_job.sh`, `jobs/submit_train_wgd.sh`
