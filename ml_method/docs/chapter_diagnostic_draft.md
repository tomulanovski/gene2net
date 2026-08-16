# Diagnostic: where the edit distance comes from

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. All edit distances here use the canonical
multi-labelled-tree edit distance, which orders children by a topology-determined canonical form
before scoring so the result does not depend on the arbitrary Newick child order. Numbers are
computed on the held-out validation split, n = 300 networks pooled across the six training
configurations, and the oracle events are read from the ground-truth simulation metadata in
clade-level form rather than decomposed from the MUL-tree.

## The question

A reconstruction can lose edit distance for three reasons. The build and the scoring themselves
can leave a residual even on a perfect input. The backbone can be wrong, because the events are
placed on the ASTRAL species tree and ASTRAL can misplace lineages. And the events can be wrong,
namely the wrong branches are flagged, the wrong copy number is produced, or the wrong parent is
chosen. This section separates the three so we know which one to work on.

## The oracle experiment

We remove event error by feeding the true events through the build, and we vary only the
backbone. The true events are the ground-truth whole genome duplication events recorded in the
simulation metadata, in clade-level form. In the first setting these true events are stamped onto
the true species tree, so everything is correct and only the build residual remains. We call this
the floor. In the second setting the same true events are stamped onto the ASTRAL backbone, so
the only added error is the backbone. We call this the ceiling, because it is the best any
predictor could reach on the ASTRAL backbone. The difference between the two is the backbone
contribution, and whatever a real method loses beyond the ceiling is its event-prediction error.

## Results

Feeding the true events through the build reconstructs the ground-truth network exactly for most
networks. On the true species tree, 79 percent of the validation networks reconstruct perfectly,
at edit distance zero and reticulation sister distance zero, and the mean floor is 0.091 with a
median of zero. The residual is not spread across all networks. It comes from a minority of about
21 percent that contain nested or overlapping reticulations, where one event is placed inside a
region another event has already reshaped and the build and fold convention cannot always recover
the exact sister structure of the second event. The reticulation leaf sets still match on these
networks, which is why the leaf distance stays near zero while the sister distance carries the
residual.

The table reports the mean edit distance per configuration on the same validation networks. The
floor uses the true species tree, the ceiling uses the ASTRAL backbone, and the model column is
our trained method with its predicted events on the ASTRAL backbone. Lower is better.

| Configuration | Floor, true backbone | Ceiling, ASTRAL backbone | Model | Backbone gap |
| --- | --- | --- | --- | --- |
| ils_low | 0.128 | 0.326 | 0.335 | 0.198 |
| ils_medium | 0.108 | 0.322 | 0.314 | 0.214 |
| ils_high | 0.071 | 0.421 | 0.413 | 0.350 |
| dup_loss_low | 0.054 | 0.352 | 0.374 | 0.298 |
| dup_loss_medium | 0.101 | 0.362 | 0.334 | 0.261 |
| dup_loss_high | 0.088 | 0.416 | 0.375 | 0.328 |
| overall | 0.091 | 0.366 | 0.357 | 0.275 |

Two facts stand out. The backbone gap, from floor to ceiling, is large, from about 0.20 to 0.35
depending on the configuration, and it grows with discordance, reaching its largest value at high
incomplete lineage sorting. The gap is the cost of placing correct events on an imperfect species
tree, and it is the dominant term. The Robinson-Foulds distance of the reconstruction confirms
the source, rising from 0.007 on the true tree to 0.132 on the ASTRAL tree, which is the ASTRAL
backbone error appearing directly.

The second fact is the more important one. The model column sits at the ceiling in every
configuration, sometimes a little below it and sometimes a little above, with an overall mean of
0.357 against a ceiling of 0.366. Feeding the true events onto the ASTRAL backbone does no better
than the model's predicted events. The model can even edge below the ceiling, because it was
trained on the ASTRAL backbone and learned to place events in a way that partly compensates for
its errors, whereas the true events are stamped on without that adjustment.

## Interpretation

The backbone is the wall and event prediction is saturated. On the validation networks the model
already extracts essentially everything the ASTRAL backbone allows, so more work on detection or
partner prediction cannot lower the edit distance further. This reverses an earlier reading of
this diagnostic, which had compared the benchmark model against the validation ceiling and
concluded that event prediction was the dominant source of error. That comparison mixed two
different test sets. Measured correctly, with the model and the ceiling on the same validation
networks, the model is at the ceiling and the backbone carries the error.

The consequence for the method is that the remaining lever for in-distribution quality is the
species tree, not the events. A rebuilt backbone, for example through phasing before tree
inference, would move the ceiling down and is the only path that can lower the edit distance on
these networks. Improving the event head cannot, because it is already at the ceiling.

## The benchmark is a separate question

The numbers above are on the validation split, which is drawn from the same simulator as the
training data. On the empirical benchmark of 21 published networks the model reaches a higher edit
distance, near 0.6, and that gap is not explained by the decomposition above. It is distribution
shift. The benchmark topologies differ from the simulated training topologies, and the model
generalises to them imperfectly. This is a distinct source of error from the backbone and the
events, and it is analysed in the results chapter. The point here is only that the
in-distribution decomposition and the benchmark gap must be kept apart, because conflating them
is what produced the earlier incorrect conclusion.

## A note on the floor

The floor is a property of the networks and the build, not of the discordance. The oracle on the
true species tree does not use the gene trees at all, so its residual does not depend on the
incomplete lineage sorting or the duplication and loss rate. The per-configuration variation in
the floor column reflects only the different networks that fall into each configuration's share
of the validation split, not a real effect of the condition. The single meaningful floor figure
is the overall 0.091, and the honest reading of it is that the build is exact for the typical
network and leaves a residual only under overlapping reticulations.

## A note on the metric

The numbers above required correcting the edit-distance measure. The multi-labelled-tree edit
distance is computed by a greedy graph-edit-distance search, and that search visits nodes in the
order the tree was written, so its result depended on the phylogenetically meaningless Newick
child order. Two identical trees could score far from zero when their child orders happened to
misalign. Ordering both trees canonically before the search removes this dependence. All edit
distances in this work use the canonical measure.
