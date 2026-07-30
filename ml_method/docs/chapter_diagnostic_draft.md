# Diagnostic: locating the error

DRAFT for the thesis diagnostic section. Prose style follows the thesis convention of no
semicolons, no non-mathematical parentheses, and no em-dashes. Numbers are the measured rooted
two-parent oracle on the validation split, unless stated otherwise.

## The question

The reconstructed network can be wrong for two reasons. Either the predicted events are wrong,
namely the method flags the wrong branches or assigns the wrong parents, or the backbone the
events are placed on is wrong, namely ASTRAL misplaces lineages before any event is stamped. The
two failure modes call for very different fixes. If the events are the problem, a better detector
or a better placement head would help. If the backbone is the problem, no amount of better
prediction can help, because the method never rearranges the backbone. This section isolates
which one dominates.

## The oracle experiment

We measure edit distance under an oracle that removes prediction error, so that only the backbone
and the build remain. We reconstruct each validation network twice. In the first setting the
events are placed on the ASTRAL backbone. In the second they are placed on the true species tree.
In both settings the events are the true events with their true parents, so detection and
placement are perfect. The difference between the two settings is therefore attributable to the
backbone alone. The first setting is the best that any predictor could achieve on the ASTRAL
backbone, so we call it the ceiling. The second is the residual left by the build itself on a
correct backbone, so we call it the floor.

## Results

The table reports mean edit distance across the validation networks of each configuration. Lower
is better.

| Configuration | ASTRAL backbone, perfect placement | True backbone, perfect placement | Backbone gap |
| --- | --- | --- | --- |
| ils_low | 0.543 | 0.114 | 0.429 |
| ils_medium | 0.582 | 0.114 | 0.468 |
| ils_high | 0.599 | 0.114 | 0.485 |
| dup_loss_low | 0.566 | 0.114 | 0.452 |
| dup_loss_medium | 0.566 | 0.114 | 0.452 |
| dup_loss_high | 0.618 | 0.114 | 0.503 |

Two structural facts stand out. The floor is constant at 0.114 across every configuration. This
is expected, because the floor uses the true species tree and the true events, and these are the
same underlying networks in every configuration. Only the gene trees change between
configurations, and the floor does not depend on the gene trees. The ceiling, in contrast, rises
as the conditions harden, from 0.543 at low incomplete lineage sorting to 0.599 at high, and from
0.452 to 0.503 as gene duplication and loss increase. The ceiling depends on the gene trees only
through ASTRAL, so its degradation is the degradation of ASTRAL.

## Interpretation

The backbone is the dominant source of error. On the ASTRAL backbone, perfect placement still
leaves an edit distance of 0.543 or more, whereas the true backbone reaches 0.114. The gap
between the two, between 0.43 and 0.50 depending on the configuration, is entirely the backbone,
because everything else in the two settings is identical. No detector and no placement head can
close this gap, since the method places events on the backbone it is given and never rebuilds it.

The floor of 0.114 is the irreducible residual of the build and the scoring convention on a
correct backbone. It is small and constant, so it is not where the method loses.

Placement is close to solved on the ASTRAL backbone. When the two parents are chosen by a simple
co-clustering heuristic rather than given as an oracle, the edit distance at low incomplete
lineage sorting is 0.558, against 0.543 for the oracle. The heuristic is within 0.015 of perfect,
so the signal needed to place parents on the ASTRAL backbone is already present in the gene trees
and is easy to extract. This further isolates the backbone as the bottleneck, because even the
placement the method does control is near its ceiling.

There is a second and smaller lever. Giving the build two parents rather than one improves edit
distance by roughly 0.19 on both backbones. This is a real gain, but it is far smaller than the
backbone gap, and the learned two-parent head did not realize it in practice. The head reached an
allopolyploid set accuracy of only about 0.29 in training and did not improve the reconstruction
on the benchmark. So the two-parent construction helps in principle, but the prediction of two
parents under gene tree noise remains hard, and the practical method uses a single parent.

## Consequence for the method

The diagnostic reframes the improvement problem. The events are not the lever, because perfect
events on the ASTRAL backbone are still capped near 0.55. The lever is the backbone. Any method
that keeps the ASTRAL backbone fixed and only decorates it with events inherits a floor near 0.55
under perfect placement, which is why better detection, better thresholds, and better event
selection do not move edit distance. The path to the 0.114 floor runs through a better backbone,
which motivates the backbone repair and phasing direction developed in the next section.
