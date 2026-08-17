# Polyploid Network Reconstruction — Design Rationale & Discussion

Reference document for thesis defense / advisor discussion. For every stage of the
pipeline: **what** we do, **why**, the **alternatives**, the **tradeoffs**, and the
**open debates**. Numbers are from the rooted model (`reconstruct_cladelabels_rooted`)
and the diagnostics of July 2026.

---

## 0. Problem statement

Given a set of **gene trees** from species that include **polyploids**, reconstruct the
**evolutionary network** — a phylogeny with **reticulation nodes** where two lineages
merged (allopolyploidy) or a lineage duplicated itself (autopolyploidy).

- **Allopolyploid** X = hybrid of two *different* diploid lineages P1, P2. In gene trees X
  appears as two **homeolog** copies, one clustering with P1, one with P2. In the network
  X is a reticulation node with **two parents**.
- **Autopolyploid** X = whole-genome duplication *within* one lineage. Its two copies are
  paralogs that diverged at the duplication; in folding these simplify away (H_strict = 0).
- Homeologs are therefore **ortholog-like** to their progenitors (they diverged at the
  *parental speciation*, not at a duplication) — a fact that matters for backbone choice (§2.1).

**Output object:** a MUL-tree (multi-labelled tree; a polyploid appears as multiple leaves)
that folds into the reticulate network. Scored against the true network.

---

## 1. Pipeline overview

```
gene trees → ASTRAL backbone → hybrid rooting → GNN (detection + partner) → build MUL-tree → fold → score
```

The framing is **detect-then-place**: find where WGD happened, then place the reticulation.
Its virtue is interpretability and speed; its ceiling (§5) is that it builds on a *single-copy*
backbone that collapses each allopolyploid into one tip.

---

## 2. Stage-by-stage design decisions

### 2.1 Backbone: species tree from gene trees

**What:** ASTRAL-IV ("ASTRAL 4") summarises the gene trees into a species tree; that tree is
the scaffold the GNN reasons over and the build stamps onto.

**Why ASTRAL:** fast, statistically consistent under ILS, standard in the field.

**Alternatives & tradeoffs:**
- **ASTRAL-IV (current):** treats multiple copies of a species as *multiple individuals*
  (population sampling). Cheap and it runs, but for homeologs this is the wrong model — it
  forces the polyploid to *one* location and averages over its two parents.
- **ASTRAL-Pro:** models copies as *paralogs* (tags gene-tree nodes as duplication vs
  speciation, discounts duplication-driven quartets). **But** for allopolyploids the homeolog
  split *is* a speciation (the parental split); ASTRAL-Pro would tag it a duplication and
  **discount the very signal we need.** So it is not obviously right either — neither tool is
  designed for polyploidy. *(Open debate: worth an empirical A/B, but low priority — see §5.)*
- **True backbone:** only available in simulation; the oracle ceiling (§5).

**Rooting:** ASTRAL outputs an *unrooted* tree; MUL-tree edit distance is rooting-sensitive.
We root via `hybrid_root`: use the **gene-tree consensus outgroup** if it forms a clean clade,
else **midpoint**. Validated accuracy ~82/82/64% (low/med/high ILS) vs ~4% for the arbitrary
root. A **bug** was fixed here (§4). Rooting at the *known* root changed nothing (0.680→0.682),
so rooting is **exhausted** as a lever.

**Debate for PI:** homeologs are ortholog-like, so the "right" backbone tool for this problem
is arguably neither ASTRAL-IV nor ASTRAL-Pro. This motivates *not* relying on a single-copy
backbone at all (§7).

### 2.2 Features

**What:** per-edge features derived from the gene trees, fed to the GNN.
- **Bipartition concordance** — fraction of gene trees supporting each backbone edge
  (grouping all copies of a species to the species' side). *The biggest detection win:
  lifted WGD-detection F1 from 0.51 to 0.79.*
- **Copy-count statistics** — how many copies of each species (the ploidy signal).
- **Co-clustering** — for each species pair, how often one species' copy sits sister to the
  other across gene trees (the partner signal).
- **Cluster-support** (copy-aware) — attempted partner feature. **Came out flat** — no
  improvement in or out of distribution. Documented as a negative result.

**Tradeoff:** features are hand-designed and species-level; they capture concordance well but
the *partner* signal they expose is weak (hence allo partner 0.44, §4).

### 2.3 Model

**What:** a Graph Attention Network (GAT) over the rooted backbone produces edge embeddings,
with two heads:
1. **Detection head** — per-edge binary: is there a WGD at this edge?
2. **Partner head** — for a WGD edge, a pairwise score over other edges → the partner.

~82k parameters. **Why a GNN:** the input is a tree (graph); message passing lets each edge's
prediction use its topological neighbourhood. **Alternatives:** set/transformer over gene trees
(Phyloformer-style) — more flexible, less inductive bias, harder to train (§7 direction C).

### 2.4 Labels

**What:** training labels = the true WGD events mapped onto each sample's *own* backbone edges,
built from ground-truth **metadata** (not by decomposing the MUL-tree).

**Why:** the naive route, `decompose_mul_tree`, **fragments clade-level allopolyploidy into
reciprocal tip labels** (a confirmed bug, §4). Metadata labels avoid this and guarantee
alignment between features, labels, and model edges.

**Tradeoff:** metadata is only available in simulation; on real data we rely on the model.

### 2.5 Detection & over-prediction

**What:** per-edge sigmoid, thresholded (0.5 by default; swept on the benchmark).

**Result:** F1 ~0.80, **precision ~0.68, recall ~0.985** → the model over-predicts,
**num_rets_diff mean ≈ 6–7** (median only 1–3; a subset of high-reticulation networks blows up).

**Why it over-predicts — three compounding causes:**
1. **No count control.** Each edge is scored independently; nothing constrains the *number* of
   events or says "these adjacent edges are one event."
2. **Recall-favoring operating point.** 0.985 recall necessarily admits ~1/3 false positives.
3. **Signal smearing.** A real WGD elevates probability on *adjacent* edges too, so one event
   fires several edges.

**Tradeoffs / alternatives:**
- Raise the threshold → fewer false positives, but hurts recall, and **does not fix edit
  distance** (backbone-bound, §5).
- **Structured / set decoding** — predict the *set* of events jointly with count control. Or
  **per-species ploidy** prediction (predict how many copies each species has) → count falls
  out naturally. Tradeoff: per-species handling of *clade-level* events (a whole clade
  duplicated) is less direct than per-edge.

**Debate for PI:** is over-prediction worth fixing? It helps ret metrics but not edit; the real
fix is architectural (§7).

### 2.6 Partner prediction & the one-parent limitation

**What:** for a WGD edge, the partner head picks one other edge = the second parent.

**Result:** **auto 0.97, allo 0.44.** Autopolyploidy (partner = self) is nearly solved;
**allopolyploid partner (which two lineages hybridised) is the weak axis.**

**The structural limitation (central to the thesis).** The WGD edge is the polyploid's *home*
branch in the backbone, which doubles as **parent 1**; the partner head predicts **parent 2**.
But parent 1 is **not chosen** — it is fixed to the polyploid's existing ASTRAL position. The
binary classifier can only *flag* an existing edge; it cannot **relocate** the polyploid to a
better parent 1. So when ASTRAL mis-places a hybrid, parent 1 is wrong and unfixable.

*Concrete:* allopolyploid **sp39**, true parents sp22 & sp42. ASTRAL buries sp39 in a big clade
(parent 1 = wrong); partner head predicts sp42 (parent 2 = right). Result: one copy in the wrong
clade, one at sp42 — **sp22 never recovered.** Even the *oracle* (true events) reproduces this,
so it is a construction limitation, not a prediction error.

**Why allo is hard:** choosing the correct second parent among ~2N edges from a weak co-cluster
signal, on top of a backbone that may already have parent 1 wrong. It is hard for everyone
(Polyphest ret_sisters rises with dup; GRAMPA-Iter 0.61–0.76).

**The fix (§7): predict both parents.** Free the model to choose parent 1 *and* parent 2 over
all edges, and have the build detach the polyploid and attach a copy at each.

### 2.7 Build / decode (MUL-tree construction)

Not straightforward; this is where a key limitation lives. Four sub-decisions:

**(a) Event selection** from per-edge scores. Strategies:
- `raw` — every edge above threshold.
- `collapse` — merge adjacent firing edges into one event (fights signal smearing).
- `cap` (used) — select events bounded by `copy_bound` (the max copies observed), preventing
  runaway over-prediction.
- `bound_driven` — driven primarily by the copy-count bound.
Tradeoff: aggressive collapsing/capping lowers false positives but can drop true events.

**(b) Stamping** (`_apply_wgd_event`): find the node with the WGD clade, **duplicate its
subtree**, graft the copy next to the partner. Auto (partner = self) creates identical sibling
subtrees; allo grafts at a distinct partner.

**(c) Ordering:** events applied bottom-up (smallest clades first) so nested events compose.

**(d) Folding:** MUL-tree → network via `ReticulateTree` (folds duplicated leaves into
reticulation nodes), then scored.

**Known failure modes:**
- **Silent event-drops** — an event whose clade was consumed by an earlier graft can't be found
  and is dropped. Now *counted*; measured tiny (6/137 on the hardest config).
- **Two-parent misplacement** — the one-partner encoding + backbone-home loses parent 1 (§2.6).
- **Whole-tree guard** — a whole-tree clade would alias the root and hang; now skipped.

### 2.8 Metrics

Scored with the same core as the baselines (`reticulate_tree` + `compare_reticulations`), so
numbers are comparable to Polyphest/GRAMPA-Iter.
- **edit_distance_multree** — edit distance on the MUL-trees. Rooting-sensitive. *(The folded-
  network graph edit distance is NP-hard and hangs on highly-reticulate nets, so it is dropped.)*
- **rf_distance** — Robinson-Foulds on the MUL-trees.
- **num_rets_diff** — |inferred − true| reticulation count.
- **ret_leaf_jaccard** — reticulation *descendants*: do we identify which species participate.
- **ret_sisters_jaccard** — reticulation *sisters*: do we identify the parental lineages
  (placement). This is the hard one.
- **ploidy_diff** — per-species copy-count error.
Lower is better throughout. **Strict vs partial match:** strict penalises unmatched
reticulations; partial normalises over matched pairs only (lenient, for count-limited methods).

---

## 3. Evaluation design

- **12 configurations:** 3 ILS levels (Ne = 200k / 1M / 2M) × duplication/loss rate
  (none / low / medium / high).
- **Data:** SimPhy gene trees down 21 published allopolyploid networks (the HomeoSorter
  dataset). Both simulation sweeps and the real network structures.
- **Splits:** train 14.3k / val 3.6k samples (rooted); the 21 networks are the held-out benchmark.
- **DATA CAVEAT:** a duplication-rate simulation bug means **only the Ne = 1M family is a valid
  low→med→high dup sweep**; the `_10M` and `ne2M` dup configs were simulated at the wrong rates.
  So the reportable grid is 6 configs: ils {low, med, high} + dup ne1M {low, med, high}.
- **Threshold note:** the operating threshold was selected on the benchmark (optimistic); the
  key finding is threshold-independent (§5), so it does not affect conclusions — but state it plainly.

---

## 4. Bugs found & fixed (methodological rigor)

1. **Label fragmentation** — `decompose_mul_tree` split clade-level allopolyploidy into
   reciprocal tip labels. Fixed with metadata labels. With the rooting fix, moved ret_leaf
   0.38→0.24 and ret_sisters 0.62→0.46 at ils_low.
2. **Silent rooting fallback** — `set_outgroup` failed ~42% of the time when the outgroup
   spanned the arbitrary root and quietly returned it; `robust_set_outgroup` (two-step reroot)
   fixed it (RF=0 reroot success 58%→100%).
3. **Detection mask clobbering positives** — an unmappable event could switch off a genuine
   positive edge from the loss; now guarded.
4. **Silent event-drops in build** — now counted and reported.
5. **Whole-tree build hang** — guarded.

All covered by a 90-test regression suite.

---

## 5. Where the error comes from (the core analysis)

Oracle decomposition (ils_low, n=200): feed the *true* events through the build on different
backbones.

| condition | edit | ret_leaf | ret_sisters |
|---|---|---|---|
| A · true backbone + true events (floor) | 0.302 | 0.178 | 0.351 |
| B · ASTRAL backbone + true events | 0.740 | 0.178 | 0.371 |
| actual · ASTRAL + predicted events | 0.680 | 0.245 | 0.459 |

**Reading:**
- **Edit distance is backbone-bound.** ASTRAL + *perfect* events (0.74, ≈0.68 rooted) equals the
  current predicted pipeline (0.68). Prediction is not the bottleneck; the backbone is. And the
  true-backbone floor (0.30) **beats Polyphest (0.42)**. → No threshold/tuning crosses this.
- **Ret metrics are prediction-bound.** A ≈ B on ret_leaf/ret_sisters — the ASTRAL backbone
  barely hurts reticulation ID — so the gap to `actual` (0.245, 0.459) is *prediction*, improvable.

**Three convergent confirmations that the single-copy backbone is the ceiling:**
1. Rooting exhausted (fixed; true ≈ hybrid).
2. Edit backbone-bound (above).
3. Placement construction-bound: the one-partner event can't hold two parents (sp39).

All three are the same object failing — a single-copy backbone + one-partner event cannot
represent a two-parent allopolyploid.

**Caveat / open test:** condition B used the *one-partner* encoding, so it conflates backbone
topology error with the one-partner limitation. The **decisive untested experiment** is a
**two-parent oracle** (place both true parents on ASTRAL, correct build). If edit → ~0.35, the
ASTRAL diploid skeleton is fine and the redesign is *bounded* (a placement head); if it stays
~0.6, the backbone itself is the wall. The diploid skeleton is ~86% topologically correct, which
makes the bounded outcome plausible. (§7)

---

## 6. Related methods

| method | approach | strengths | weaknesses |
|---|---|---|---|
| **Polyphest** | cluster copies → star-tree assembly; ploidy from copy-count distribution | near-perfect at low noise; accurate count | ploidy heuristic breaks at high dup → over-predicts |
| **GRAMPA-Iter** | MUL-tree reconciliation search | principled | worst throughout; systematically under-predicts |
| **HomeoSorter** (2024) | permutation phasing of homeologs → MUL-tree | learned-free phasing; same data | permutation cost; not yet compared against (must add) |
| **homologizer** | Bayesian phasing into subgenomes | principled | slow / doesn't scale |
| **ASTRAL-Pro** | paralogy-aware species tree | handles multi-copy | models homeologs as paralogs (wrong for allo) |
| **Ours** | GNN detect-then-place | fast, interpretable, noise-robust; beats GRAMPA-Iter throughout | single-copy backbone ceiling; allo placement |

**Positioning:** we beat GRAMPA-Iter on every metric/config; reach/pass Polyphest only at high
dup/loss; Polyphest dominates easy settings and placement. **HomeoSorter must be added to the
comparison** (it shares our data and is the closest method).

---

## 7. The redesign: predict both parents

**Current vs proposed:**

| stage | current | two-parent |
|---|---|---|
| detection | per-edge: WGD on this edge? | per-species: is X polyploid, how many copies? (count control) |
| placement | 1 partner; parent 1 = ASTRAL position (fixed) | **2 parent edges**, both predicted over all edges (permutation-invariant) |
| build | keep in place + graft 1 copy | **detach** polyploid, attach a copy at *each* predicted parent |

**Why it should help:** the diploid skeleton is ~86% right; the *only* broken thing is placing
allopolyploids at their two parents. Predicting both parents (and not inheriting parent 1 from
ASTRAL) targets exactly that — and could reach the true-backbone floor (0.30 < Polyphest 0.42)
*on the ASTRAL backbone*, without a full rebuild.

**Design choices / debates:**
- **Emit two edges:** two independent softmaxes (simple) vs joint top-2 vs pairwise scorer.
- **Symmetry:** parents are an *unordered* pair → permutation-invariant loss.
- **Auto vs allo:** auto = both parents identical; allo = distinct. The head must represent both.
- **Higher ploidy:** >2 copies → >2 parents, or cap.
- **Detaching risk:** placing both copies purely from prediction removes the ASTRAL-collapse
  ceiling but also removes the ASTRAL fallback — full exposure to model error.
- **Per-species vs per-edge detection:** per-species gives count control (fixes over-prediction)
  but must handle clade-level events.

**Decisive next experiment:** implement the two-parent build and run the two-parent oracle
(above). It costs one build change and tells you whether the redesign is bounded (placement head,
keep ASTRAL) or big (better backbone / phasing). It also gives the thesis its build/decode section.

**Bigger alternative (if the oracle says the backbone is the wall):** a phasing-based
reconstruction — cluster copies into subgenomes, place each subgenome by its parent — which is
HomeoSorter's paradigm done with a learned, robust, fast model. Higher risk, higher novelty.

---

## 8. Open questions for the advisor

1. **Bar:** is "beats GRAMPA-Iter, matches Polyphest at high dup, plus a principled ceiling
   analysis" a sufficient thesis, or is beating Polyphest broadly required (→ redesign)?
2. **Scope:** bounded two-parent placement head vs full phasing reconstruction — decided by the
   two-parent oracle.
3. **Publication:** conditional on the redesign clearing SOTA + adding HomeoSorter + runtime.
   Which venue (bioinformatics methods vs ML-for-biology workshop)?
4. **Framing of "predict both parents":** unordered-pair prediction, permutation-invariant — is
   there a cleaner formulation (e.g. edge-pair scoring, or a subgenome-assignment head)?
5. **Backbone:** accept ASTRAL-IV's limitation, or test ASTRAL-Pro / a phasing-derived backbone?

---

*Companion artifacts: the advisor summary report (published), and the reproducible diagnostic
scripts `oracle_test.py`, `compare_oracle_faithfulness.py`, `debug_allo_build.py`.*
