# Diagnostic: where the edit distance comes from

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. All edit distances here use the canonical
multi-labeled-tree edit distance, which orders children by a topology-determined canonical form
before scoring so the result does not depend on the arbitrary Newick child order. Numbers are
the oracle on the validation split, n = 200 per configuration.

## The question

A reconstruction can lose edit distance for three reasons. The build and scoring themselves
leave a residual even on a perfect input. The backbone can be wrong, because the events are
placed on the ASTRAL species tree and ASTRAL misplaces lineages. And the events can be wrong,
namely the wrong branches are flagged, the wrong copy number is produced, or the wrong parent is
chosen. This section separates the three so we know which one to work on.

## The oracle experiment

We remove event error by feeding the TRUE events through the build, and we vary only the
backbone. In the first setting the true events are stamped onto the true species tree, so
everything is correct and only the build residual remains. We call this the floor. In the second
setting the same true events are stamped onto the ASTRAL backbone, so the only added error is the
backbone. We call this the ceiling, because it is the best any predictor could reach on the
ASTRAL backbone. The difference between the two is the backbone's contribution, and whatever a
real method loses beyond the ceiling is its event-prediction error.

## Results

The table reports mean edit distance across the validation networks. Lower is better.

| Configuration | Floor (true backbone) | Ceiling (ASTRAL backbone) | Backbone gap |
| --- | --- | --- | --- |
| ils_low | 0.109 | 0.251 | 0.142 |
| ils_medium | 0.109 | 0.310 | 0.201 |
| ils_high | 0.109 | 0.391 | 0.282 |
| dup_loss_low | 0.109 | 0.311 | 0.202 |
| dup_loss_medium | 0.109 | 0.341 | 0.232 |

Two structural facts stand out. The floor is exactly 0.109 in every configuration. This is
expected, because the floor uses the true species tree and the true events, which are the same
underlying networks across configurations, and only the gene trees change. So 0.109 is the
irreducible residual of the build and the scoring convention, not a modelling error. The ceiling
rises from 0.251 at low incomplete lineage sorting to 0.391 at high, and sits near 0.31 to 0.34
under gene duplication and loss. The ceiling depends on the gene trees only through ASTRAL, so
its rise is the degradation of ASTRAL as discordance grows.

## Interpretation

The backbone is a moderate contributor, not a wall. Placing true events on the ASTRAL backbone
still reconstructs well, with a ceiling between 0.25 and 0.39 depending on conditions, against a
floor of 0.109. The backbone's contribution is the gap between them, from 0.14 at low incomplete
lineage sorting to 0.28 at high. It grows with discordance, as expected, but it does not dominate
the edit distance, and on the easier conditions half of the networks reach the floor exactly even
on the ASTRAL backbone.

The larger source of error is event prediction. A real method must also detect the right
branches, produce the right copy number, and choose the right parents, and whatever it loses
beyond the ceiling is that event error. Our method sits well above the ceiling, so event
prediction, not the backbone, is where most of its edit distance is lost. This is encouraging,
because event prediction is improvable within the present design. Two concrete gains already
follow from it. Filling each species to its inferred copy number rather than capping at a
detection threshold recovers events the threshold suppressed and improves copy number and
reticulation together. Better placement of the second parent is the remaining lever, and it is
the partner problem analysed separately.

## Consequence for the method

The diagnostic redirects the improvement effort. The backbone is worth improving mainly at high
incomplete lineage sorting, where its contribution grows, and a rebuilt backbone through phasing
is a genuine but conditional gain rather than the only path. The dominant and more tractable
lever is event prediction, which the reconstruction controls directly. This reverses an earlier
reading of this diagnostic, discussed next, that had attributed most of the edit distance to the
backbone.

## A note on the metric

The numbers above required correcting the edit-distance measure. The multi-labeled-tree edit
distance is computed by a greedy graph-edit-distance search, and that search visits nodes in the
order the tree was written, so its result depended on the phylogenetically meaningless Newick
child order. Two identical trees could score far from zero when their child orders happened to
misalign. Ordering both trees canonically before the search removes this dependence. Under the
uncorrected measure the ASTRAL ceiling appeared to be about 0.54 and the backbone looked like the
dominant wall. Under the corrected canonical measure the ceiling is 0.25 to 0.39, which is why
the conclusion changed. All edit distances in this work use the canonical measure.
