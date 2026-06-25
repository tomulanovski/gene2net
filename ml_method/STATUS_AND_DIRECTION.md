# Gene2Net-GNN: Status and Direction

A consolidated summary of the method, what we measured, the rooting investigation and
its outcome, where we stand against Polyphest and GRAMPA-iter, and the proposed next direction.

---

## 1. The method, in one breath

A two-stage learned approach that is deliberately intuitive:

1. Infer a diploid species tree from the gene trees with ASTRAL.
2. For each edge of that species tree, a features-plus-GNN model predicts whether a whole
   genome duplication happened on it (detection).
3. For each predicted WGD edge, a partner head predicts which edge the duplicated copy
   attaches to. Partner equal to the edge itself means autopolyploidy, a different edge means
   allopolyploidy.
4. Build a multi-labeled tree by stamping those events onto the species tree, then fold it to
   a network with the Huber algorithm (the same folding is applied to every method, so it is
   not a source of difference).

The whole thing is a forward pass plus deterministic feature extraction, no combinatorial search.

---

## 2. Data

- Training and validation: fully simulated. SimPhy birth-death species trees with 0 to 15
  random WGD events, gene trees simulated down them. About 17,900 packaged samples across 9
  configurations spanning low/medium/high ILS and low/medium/high duplication-and-loss. Split
  80/20 train/val by whole sample.
- Test: the 21 empirical polyploid networks, evaluated across 12 evolutionary configurations.
  The model never trains on these. They get the same inputs GRAMPA-iter receives, so the
  comparison is apples to apples.

---

## 3. Detection: strong, and well characterized

- Tuned binary WGD-edge F1 is about 0.80 in distribution. Adding the 5 detection features took
  it from about 0.51 to about 0.80.
- The errors are false positives, not false negatives. Recall is saturated near 1.0 (the model
  almost never misses a real duplication), precision is about 0.77 at the operating threshold of
  0.90. We over-flag rather than miss.
- The false positives are characterized: 99.7 percent are single-species tips, and 87 percent
  sit inside a truly-duplicated clade (redundant with the real ancestral event). On simulated
  data a containment dedup removes them, but on the real networks the flagged edges are separate
  tips, not nested, so dedup is a no-op there.

### Feature importance (permutation, all 9 configs)

Detection rides on two edge features:

| feature | F1 drop when permuted |
| --- | --- |
| frac_clade_dup | 0.625 |
| clade_size | 0.239 |
| copy_pair_div_mean | 0.067 |
| mirrored_sister_frac | 0.043 |
| concordance | 0.042 |
| everything else (incl. all 13 node features) | about 0 or negative |

Implication: for detection the node features (copy-count stats and co-clustering summaries)
contribute nothing. They could be pruned for speed, pending a check that the partner head does
not need them. The GNN earns its place through localization (clade_size plus the attention
picking the right edge in the high-frac ancestor chain), lifting precision over a single-feature
threshold baseline of about 0.62.

---

## 4. Partner prediction: good, and it already proves the key idea

- Auto partner accuracy about 0.99 (solved by the mirrored-sister feature).
- Allo partner accuracy about 0.80 to 0.86.
- The lever that got allo from 0.52 to 0.80 was a pairwise co-clustering feature, namely how
  often clade i and clade j species appear as sisters across the gene trees. This is a thin
  slice of gene-tree cluster support, and it is the single most important precedent for the
  direction in section 9.
- Conditioned on the true WGD edges, the remaining allo errors are about 70 percent far misses
  (the predicted partner is a disjoint clade), which means the residual is a signal limitation,
  not a tie-breaking one.

---

## 5. Head-to-head comparison

All numbers are MUL-tree edit distance, lower is better, scored in the same strict mode as the
thesis (this is the metric and the scoring function the thesis uses for every method). All 21
networks score in every configuration.

Edit distance:

| condition | Polyphest | GRAMPA-iter (no prior) | GRAMPA-iter (with prior) | Ours |
| --- | --- | --- | --- | --- |
| low ILS, no dup loss | 0.42 | 0.80 | 0.73 | 0.79 |
| medium ILS, no dup loss | 0.44 | 0.83 | - | 0.76 |
| high ILS, no dup loss | 0.56 | 0.91 | - | 0.79 |
| low ILS, high dup loss | 0.75 | 0.86 | - | 0.82 |
| medium ILS, high dup loss | 0.77 | 0.90 | - | 0.81 |
| high ILS, high dup loss | 0.53 | 0.88 | - | 0.78 |

Reticulation descendants (lower is better): Polyphest 0.02 to 0.34, GRAMPA-iter no prior 0.39
to 0.61, ours 0.37 to 0.56. Reticulation sister: GRAMPA-iter and ours both about 0.62 to 0.76.

What it says:
- We are flat across all 12 configurations (edit 0.76 to 0.82). Both competitors degrade with
  ILS and dup/loss. That noise robustness is a genuine strength.
- Versus prior-free GRAMPA-iter we win on edit distance and reticulation descendants, and the
  margin grows under noise.
- The with-prior GRAMPA-iter is given the true polyploid set, which we are not. It is stronger
  at the easy end.
- Polyphest is the most accurate on edit distance and reticulation metrics everywhere.

---

## 6. Build strategies we tried

Several MUL-tree build strategies from the same predictions: raw, collapse, cap, collapse_cap,
bound_driven, dedup. Cap (limit duplications per species to the inferred consensus copy number)
is the well-formed default we benchmark with. None of them moved edit distance, which the oracle
in section 7 explains.

---

## 7. The backbone investigation (the main finding)

We chased why we are behind Polyphest on edit distance.

### Oracle: edit distance is backbone-driven, reticulation metrics are not

Feeding the model's predicted events onto two different backbones (same events, simulated
ils_low):

| backbone | edit (MUL) | ret_leaf | ret_sisters |
| --- | --- | --- | --- |
| ASTRAL (our pipeline) | 0.762 | 0.255 | 0.431 |
| true diploid backbone | 0.348 | 0.250 | 0.409 |

Two conclusions that shaped everything after:
- Edit distance is set by the backbone. Predicted events on the true backbone reach 0.348,
  which beats Polyphest at 0.42. So our predictions are good enough, the backbone is the cost.
- Reticulation metrics are backbone-independent. ret_leaf and ret_sisters are essentially
  identical on both backbones, so they are driven by the predictions and where reticulations
  attach, not by the topology.

### Why the backbone hurts: ASTRAL is unrooted

The gap was not ASTRAL's topology. Even on samples where ASTRAL recovers the exact unrooted
topology, the build still scored badly. The cause is rooting. ASTRAL estimates an unrooted
species tree (confirmed in its own documentation), and our pipeline runs it with no outgroup
(run_astral.sh, no --root), so the species tree gets an arbitrary root. The MUL-tree edit
distance is computed on a rooted (directed) graph, so a correct topology with the wrong root
still scores far from the truth. Direct inspection confirmed it: on exact-topology samples the
true-backbone build matched the ground truth exactly while the ASTRAL build was 0.7 to 0.9 off,
differing only in the root.

This also explains the pattern across methods. We and GRAMPA-iter both build on the ASTRAL tree
and both inherit the bad rooting, so both have inflated edit distance. Polyphest never uses the
ASTRAL tree (it builds structure from clusters), so it is unaffected.

### A rooting method that works

We validated several ways to root the ASTRAL tree against the known true root:

| method | ils_low | dup_loss_high | ils_high |
| --- | --- | --- | --- |
| arbitrary ASTRAL (current) | 4% | 6% | 2% |
| midpoint | 72% | 64% | 34% |
| gene-tree consensus (mode split) | 48% | 62% | 44% |
| hybrid (consensus if clean, else midpoint) | 82% | 82% | 64% |

The hybrid is best and robust. Notably, under high ILS midpoint collapses (the clock breaks) but
the gene-tree consensus holds, so the hybrid keeps most of the signal where each individual method
fails. We wired this hybrid rooting into the feature pipeline.

### What we tested re-rooting on, and whether it worked

We tested on conf_ils_low_10M (the 21 real networks at low ILS), three ways:

| setup | edit distance |
| --- | --- |
| unrooted (original model) | 0.790 |
| root + infer (original model on rooted trees) | 0.741 |
| root + retrain (model retrained on rooted ils_low) | 0.725 |

It worked, but only a little. Retraining restored the reticulation metrics that root-plus-infer
had nudged the wrong way, but the total edit-distance gain is about 0.05 to 0.065, nowhere near
the 0.42 needed to match Polyphest, and part of even that is confounded (the rooted model was
trained on a single config).

### Why it only helped a little

Rooting fixes only one of the two errors we inherit from ASTRAL. The root is corrected, but the
topology errors remain, and on a 91 percent-correct backbone those errors still amplify once
clades are duplicated. The true-backbone 0.348 came from correct topology AND root, which is not
reachable by fixing only the root. On the real networks, which are larger and more reticulate,
the gain is smaller still.

Decision: the rooting code is wired and is a legitimate small improvement and a clean thesis
subsection (we identified and fixed the unrooted-ASTRAL problem), but the full 9-config
re-package and retrain is not worth the marginal gain. It is shelved as a polish, not a headline.

---

## 8. Where we stand, honestly

- Detection of WGD edges and polyploid lineages: strong, noise robust.
- Reticulation identification: we beat prior-free GRAMPA-iter, behind Polyphest.
- Edit distance: behind Polyphest, ahead of prior-free GRAMPA-iter, and we understand exactly why
  (Polyphest builds structure jointly from clusters, we build on a fixed backbone with rooting
  and topology errors).
- Speed: very likely fastest, since ours is deterministic feature extraction plus a forward pass
  while Polyphest solves an ILP and GRAMPA does reconciliation search, but this is not yet
  measured end to end.
- Intuitiveness: a real asset. Our method explains as find the duplication, find its partner,
  build the network. Polyphest explains as enumerate clusters, weight them, build a compatibility
  graph, solve a max-weight integer program, greedily assemble. For adoption, explainability and
  trust matter.

---

## 9. What to do now

Goal restated: beat Polyphest on edit distance, reticulation descendants, and reticulation
sister, and become the method researchers actually reach for. Keep the intuitive framing, which
is a strategic asset, and take inspiration from Polyphest where useful.

Key analytical constraint from section 7: edit distance is a backbone/structure problem, while
the two reticulation metrics are a prediction/placement problem. They need different levers.

The aligned direction: keep our intuitive detect-then-place method, and borrow the one ingredient
that makes Polyphest accurate, namely gene-tree cluster support, to drive reticulation placement.
We already proved this works in miniature: the pairwise co-clustering feature (a thin slice of
cluster support) is what lifted allo placement from 0.52 to 0.80. The plan is to feed the partner
and placement step the full weighted cluster support that Polyphest exploits, instead of only
pairwise sisters. This targets reticulation descendants and reticulation sister directly, the
metrics the backbone cannot help, and it stays a learned forward pass.

Proposed sequencing:
1. First win: enrich the cluster-support signal feeding reticulation placement, measured against
   Polyphest's reticulation metrics. Contained, testable, builds on a proven result.
2. Decide later whether to extend the cluster signal to detect reticulations the species tree
   misses entirely, which is more powerful but edges toward Polyphest's full cluster approach.
3. Separately and cheaply: measure end-to-end speed to lock in that likely win.
4. Edit distance via the backbone stays a later, harder lever (cluster-informed or co-inferred
   topology), not the immediate focus.

Open scouting tasks that make the above concrete and de-risk it:
- Read exactly how Polyphest turns clusters into the MUL-tree (we have its code).
- Per-network comparison of where and why Polyphest beats us, to confirm the lever.
