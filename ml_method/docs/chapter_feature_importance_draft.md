## Feature Importance

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes.

We ranked the features by permutation importance on the validation split. Each feature was shuffled
across the validation examples, which breaks its link to the labels while preserving its
distribution, and the resulting drop in a target measure gives the reliance of the trained model on
that feature. The procedure was applied to both heads on the shipped model with no features removed,
with detection F1 as the target measure for one head and allopolyploid partner accuracy for the
other.

Detection rests on one feature. The clade-duplicated fraction gives an F1 drop of 0.284, more than
three times that of the next feature, clade size at 0.079, while the mean copy number, the
mirrored-sister fraction and the copy-pair divergence each lower F1 by less than 0.05. These are the
features designed to detect polyploidization events, which supports that design. The node
co-clustering summaries contribute nothing to detection, at or below 0.01, as expected of a signal
intended for placement.

Placement draws on many features. The pairwise co-clustering mean leads with an accuracy drop of
0.40. Clade size, depth and mean copy number follow at 0.31 to 0.36, the pairwise co-clustering
maximum is at 0.27, and the copy-number distribution, the node co-clustering maximum, the
concordance and the copy-pair divergence each contribute between 0.17 and 0.25. Cluster support adds
0.17, through its total rather than its peak. Depth is notable, since it does nothing for detection
yet is among the strongest inputs to placement, the model apparently using the position of a clade
in the tree to constrain where its partner can lie.

Six of the 26 features contribute nothing. The branch length and the duplication synchrony lower
neither head by more than 0.005. Of the five node co-clustering summaries only the maximum matters,
through a placement drop of 0.23, while the mean, minimum, standard deviation and median are near
zero in both heads and redundant with the pairwise feature computed from the same co-clustering
fractions. These six could be pruned, but doing so changes the input dimensions and would require
retraining and re-evaluating the whole pipeline for a negligible expected gain, so the shipped model
retains the full feature set.
