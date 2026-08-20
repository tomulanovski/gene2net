# Diagnostic: where the reconstruction error comes from

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. All distances are the normalized mu-distance.
Numbers are computed on the held-out validation split, n = 300 networks pooled across the six
training configurations, and the oracle events are read from the ground-truth simulation metadata
in clade-level form rather than decomposed from the MUL-tree.

## The question

A reconstruction can lose mu-distance for three reasons. The build and the scoring themselves can
leave a residual even on a perfect input. The backbone can be wrong, because the events are placed
on the ASTRAL species tree and ASTRAL can misplace lineages. And the events can be wrong, namely
the wrong branches are flagged, the wrong copy number is produced, or the wrong parent is chosen.
This section separates the three so we know which one to work on.

## The oracle experiment

We remove event error by feeding the true events through the build, and we vary only the backbone.
The true events are the ground-truth whole genome duplication events recorded in the simulation
metadata, in clade-level form. In the first setting these true events are stamped onto the true
species tree, so everything is correct and only the build residual remains. We call this the floor.
In the second setting the same true events are stamped onto the ASTRAL backbone, so the only added
error is the backbone. We call this the ceiling, because it is the best any predictor could reach
on the ASTRAL backbone. The difference between the two is the backbone contribution, and whatever a
real method loses beyond the ceiling is its event-prediction error.

## Results

The three points on the pooled validation split are as follows.

| Setting | mu-distance |
| --- | --- |
| Floor, true events on the true backbone | 0.028 |
| Ceiling, true events on the ASTRAL backbone | 0.099 |
| Model, predicted events on the ASTRAL backbone | 0.113 |

The model's error of 0.113 decomposes into three parts. The build residual is 0.028, the ASTRAL
backbone adds 0.071 on top of that to reach the ceiling of 0.099, and the model's own event
prediction adds a further 0.014 to reach 0.113. So the backbone is 0.071 of the total, the build
residual is 0.028, and event prediction is only 0.014.

The build residual is small and concentrated. Feeding the true events onto the true backbone
reconstructs the ground-truth network exactly for most networks, with a median mu-distance of
zero, so more than half of the validation networks reconstruct perfectly. The mean floor of 0.028
comes from a minority that contain nested or overlapping reticulations, where one event is placed
inside a region another event has already reshaped and the build and fold convention cannot always
recover the exact sister structure of the second event. On these networks the reticulation leaf
sets still match, with a mean leaf distance of 0.076 against a sister distance of 0.116, which is
why the residual sits in the sister structure and not in which lineages are reticulate.

## Interpretation

The backbone is the wall and event prediction is nearly saturated. The gap from the ceiling to the
model is 0.014, while the gap from the floor to the ceiling is 0.071, so the backbone contributes
about five times as much error as the event head. On the validation networks the model already
extracts almost everything the ASTRAL backbone allows, so more work on detection or partner
prediction can recover at most 0.014 of mu-distance, whereas the backbone accounts for 0.071. This
holds under the mu-distance exactly as it did under the earlier edit-distance reading. It also
corrects a still earlier reading of this diagnostic, which had compared the benchmark model against
the validation ceiling and concluded that event prediction was the dominant source of error. That
comparison mixed two different test sets. Measured correctly, with the model and the ceiling on the
same validation networks, the model is essentially at the ceiling and the backbone carries the
error.

The consequence for the method is that the remaining lever for in-distribution quality is the
species tree, not the events. A rebuilt backbone, for example through phasing before tree inference,
would move the ceiling down and is the only path that can lower the mu-distance materially on these
networks. Improving the event head cannot, because it is already within 0.014 of the ceiling.

## The benchmark is a separate question

The numbers above are on the validation split, which is drawn from the same simulator as the
training data. On the empirical benchmark of 21 published networks the model reaches a higher
mu-distance, around 0.27 to 0.30, and that gap from the validation 0.113 is not explained by the
decomposition above. It is distribution shift. The benchmark topologies differ from the simulated
training topologies, and the model generalises to them imperfectly. This is a distinct source of
error from the backbone and the events, and it is analysed in the results chapter. The point here
is only that the in-distribution decomposition and the benchmark gap must be kept apart, because
conflating them is what produced the earlier incorrect conclusion.

## A note on the floor

The floor is a property of the networks and the build, not of the discordance. The oracle on the
true species tree does not use the gene trees at all, so its residual does not depend on the
incomplete lineage sorting or the duplication and loss rate. The single meaningful floor figure is
the overall 0.028, and the honest reading of it is that the build is exact for the typical network
and leaves a residual only under overlapping reticulations, concentrated in the sister structure.

## A note on the metric

The mu-distance is order-invariant by construction, because it compares two networks by the
multiset of per-node path-count vectors rather than by a search over the tree structure. It
therefore does not depend on the order in which children are written, which is a phylogenetically
meaningless choice, and it needs no canonical-ordering correction. This is one reason it is
preferred here over a graph-edit-distance measure, whose greedy search visited nodes in the written
order and could score two identical networks far from zero when their child orders misaligned.
