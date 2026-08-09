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
per gene tree dominates, with an F1 drop of 0.38, roughly seven times the next feature. The rest
of the signal is copy number and clade size, namely mean copy number, clade size,
mirrored-sister fraction, and copy-pair divergence, each with a drop below 0.05. Two things are
notable. First, the features designed for whole genome duplication detection are exactly the
important ones, which is direct evidence for that design choice and is consistent with the F1
gain those features produced when they were added. Second, the co-clustering node summaries
contribute essentially nothing to detection, with drops at or below 0.001. That is expected,
because co-clustering is a partner signal, not a detection signal.

## Partner

Partner prediction draws on a broader set of features, and the ranking is different. The pairwise
co-clustering feature dominates, with partner-accuracy drops of 0.33 for its mean channel and
0.23 for its max channel. This is the direct signal the partner head consumes and it is the
single most important input, which validates the pairwise design. Below it sit two groups. The
copy-aware cluster-support pairwise feature contributes a drop of 0.15, so it does help partner
even though it was flat for the two-parent head. And the structural and copy features act through
the shared trunk and the edge embeddings the partner head reads, so clade size, mean copy number,
concordance, depth, and copy-pair divergence each contribute drops between 0.11 and 0.19.

## What earns its place, and what is prunable

Reading the two rankings together gives a clean account of every feature.

- The detection features earn their place through detection, led by the clade-duplicated fraction.
- The pairwise co-clustering features earn their place through partner, as the dominant partner
  signal, and the cluster-support pairwise feature earns a smaller place there as well.
- The copy-number and structural features earn their place through both heads, primarily by
  informing where events sit.
- Branch length and duplication synchrony contribute nothing to either head and are prunable.
- The co-clustering node summaries, in contrast to the pairwise co-clustering feature, contribute
  little to either head, with only the maximum channel showing a modest partner drop of 0.08. So
  the node co-clustering summaries are largely redundant with the pairwise co-clustering feature
  and are candidates for pruning.
- One feature is worth singling out. Depth is inert for detection but contributes a partner drop
  of 0.13, so it earns its place through the partner head rather than through detection, and it
  should be kept.

We report these rankings but retain the full feature set in the deployed model. Pruning the
low-importance features would change the input dimensions and require retraining and re-evaluating
the whole pipeline, and because those features contribute little by construction, the expected
change in accuracy is negligible. So the value here is the analysis, namely knowing which features
drive each head, rather than a leaner model.
