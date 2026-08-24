# Feature importance

DRAFT for the thesis, to sit inside the ablations section. Prose style follows the thesis
convention of no semicolons, no non-mathematical parentheses, and no em-dashes.

We rank features by permutation importance. For each feature we shuffle its values across the
validation examples, breaking its association with the labels while keeping its marginal
distribution, rerun the trained model, and measure the drop in a target metric. A large drop
means the feature matters to the trained model. We do this twice with two targets, detection F1
and allopolyploid partner accuracy, because the method has two heads with different needs. The
model is the shipped one, and features are not removed, so this ranks the features of the
deployed model rather than proposing a new one.

## Detection

Detection is carried by a single feature and a short tail. The fraction of the clade duplicated
per gene tree dominates, with an F1 drop of 0.284, more than three times the next feature. The
rest of the signal is clade size and copy number, namely clade size at 0.079 and then mean copy
number, mirrored-sister fraction, and copy-pair divergence, each with a drop below 0.05. Two
things are notable. First, the features designed for whole genome duplication detection are
exactly the important ones, which is direct evidence for that design choice and is consistent with
the F1 gain those features produced when they were added. Second, the co-clustering node
summaries contribute essentially nothing to detection, with drops at or below 0.01. That is
expected, because co-clustering is a partner signal, not a detection signal.

## Partner

Partner prediction draws on a much broader set of features, and the ranking is different. The
pairwise co-clustering feature is the single most important input, with partner-accuracy drops of
0.40 for its mean and 0.27 for its maximum, which validates the pairwise design. But unlike
detection, partner does not rest on one feature. Clade size, depth, and mean copy number follow
immediately, each with a drop between 0.31 and 0.36, and the copy-number distribution, the node
co-clustering maximum, concordance, and copy-pair divergence each contribute between 0.17 and
0.25. The copy-aware cluster-support pairwise feature contributes 0.17 through its intensity
channel, so it does help partner even though it was flat for the earlier two-parent head. Partner
is therefore a genuinely multi-feature decision, reading the pairwise co-clustering together with
where the clade sits in the tree and how many copies its species carry.

## What earns its place, and what is prunable

Reading the two rankings together gives a clean account of every feature.

- The detection features earn their place through detection, led by the clade-duplicated fraction.
- The pairwise co-clustering features earn their place through partner, as the dominant partner
  signal, and the cluster-support pairwise feature earns a smaller place there through its
  intensity channel.
- The copy-number and structural features earn their place mainly through partner, where clade
  size, depth, and mean copy number are near the top of the ranking, and secondarily through
  detection.
- Branch length and duplication synchrony contribute nothing to either head, with drops at or
  below 0.005, and are prunable.
- Among the co-clustering node summaries, only the maximum earns its place, through a partner drop
  of 0.23. The mean, minimum, standard deviation, and median are near-zero in both heads, so they
  are largely redundant with the pairwise co-clustering feature and are candidates for pruning.
- One feature is worth singling out. Depth is inert for detection but contributes a partner drop
  of 0.32, so it earns its place through the partner head rather than through detection, and it
  should be kept.

We report these rankings but retain the full feature set in the deployed model. Pruning the
low-importance features would change the input dimensions and require retraining and re-evaluating
the whole pipeline, and because those features contribute little by construction, the expected
change in accuracy is negligible. So the value here is the analysis, namely knowing which features
drive each head, rather than a leaner model.
