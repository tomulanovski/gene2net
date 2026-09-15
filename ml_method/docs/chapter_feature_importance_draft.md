# Feature importance

DRAFT for the thesis, to sit inside the ablations section. Prose style follows the thesis
convention of no semicolons, no non-mathematical parentheses, and no em-dashes.

We rank features by permutation importance on the validation split. Each feature is shuffled across
the validation examples, which breaks its link to the labels while keeping its distribution, and the
drop in a target metric measures how much the trained model relies on it. We do this for both heads,
with detection F1 as the target for one and allopolyploid partner accuracy for the other, on the
shipped model with no features removed.

Detection rests on one feature. The fraction of the clade duplicated per gene tree gives an F1 drop
of 0.284, more than three times the next feature, clade size at 0.079, and mean copy number,
mirrored-sister fraction and copy-pair divergence each lower it by less than 0.05. These are the
features designed for whole genome duplication detection, which supports that design. The
co-clustering node summaries contribute nothing to detection, at or below 0.01, as expected for a
signal meant for partner prediction.

Partner prediction draws on many features. The pairwise co-clustering mean leads with an accuracy
drop of 0.40. Clade size, depth and mean copy number follow at 0.31 to 0.36, the pairwise
co-clustering maximum is at 0.27, and the copy-number distribution, the node co-clustering maximum,
concordance and copy-pair divergence each contribute 0.17 to 0.25. Cluster support adds 0.17 through
its summed channel. Depth stands out, since it does nothing for detection but is among the strongest
inputs to partner prediction.

Only two groups of features contribute nothing. Branch length and duplication synchrony lower
neither head by more than 0.005. Among the co-clustering node summaries only the maximum matters,
through its partner drop of 0.23, while the mean, minimum, standard deviation and median are near
zero in both heads and redundant with the pairwise feature. These could be pruned, but that would
change the input dimensions and require retraining and re-evaluating the whole pipeline for a
negligible expected gain, so the shipped model keeps the full feature set.
