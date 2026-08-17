# Phaser prototype — spec

**Question this answers:** can a polyploid's second parent be recovered from gene-tree
*coherence* better than the current partner head's **0.44** (allo accuracy)? That single
number decides whether the subgenome-phasing direction is worth building.

## Background (what "second parent" means)

Metadata stores each allopolyploid event as:
- `target_clade` — the polyploid itself (usually a single species X),
- `partner_clade` — its **second parent** lineage.

The current model predicts `partner_clade` and gets **0.44** on allo events. The home
(first) parent is ASTRAL's placement of X (~75% correct, separate).

The coherence signal: X's two homeolog copies cluster with its two parents across gene
trees. `compute_clustering_profile(gene_trees, X, all_species)` gives, for X, how often a
copy of X sits sister to each other species — so the true parents should be near the top.

## Stage 0 — non-learned baseline (fast, ~1–2 days) — `scripts/phaser_baseline.py`

No training. For each allo event (single-species target X):
1. Compute X's co-clustering profile.
2. Rank other species by co-clustering frequency.
3. Check whether `partner_clade` appears in the top-1 / top-2 / top-3.

**Report:** partner-recovery rate at top-1/2/3, per config, vs the 0.44 baseline.
- One of X's two peaks is the home parent, so the partner often sits at rank 2 — hence
  top-2 is the fair comparison.

**Interpretation:**
- top-2 recovery **≫ 0.44** → the signal is rich and the current head *under-uses* it →
  a learned phaser will help → **build Stage 1**.
- top-2 recovery **≈ 0.44** → the signal is genuinely limited → phasing won't beat it →
  reconsider (consolidate / different angle).

## Stage 1 — learned phaser (only if Stage 0 is promising, ~1–2 weeks)

Extend the GNN: for each polyploid X, score every backbone edge as a candidate parent
(edge embedding + X's co-clustering with that edge's clade), predict the **top-2** parents
with a permutation-invariant loss against the true two parents (home + partner).
Train on the sim data (true parents from metadata). Same metric; compare to Stage 0 and 0.44.

**Success criterion for the whole direction:** Stage 0 or Stage 1 clearly beats 0.44 →
proceed to the full phasing-based method. Otherwise → the signal is the wall, not the method.

## Metrics to log
- partner-in-top-k (k=1,2,3), overall and per config.
- fraction of allo events that are single-species (sanity on the "mostly tips" assumption).
- (Stage 1 also) wall-clock vs the current pipeline, for the speed story.
