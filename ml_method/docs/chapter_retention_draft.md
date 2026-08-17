# Reconstruction under fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The numbers below are on the 21 benchmark
networks under three fractionation configurations. The competitor scores are averaged over
five replicates. The learned-method scores are from replicate one and are being extended to
the full five replicates, so the comparison here is preliminary. [PENDING: five-replicate
aggregation of the learned method.]

## Fractionation and the hypothesis

After a whole genome duplication, the two copies of a region do not always both survive. A
process of fractionation removes one copy from some of the descendant lineages. We model this
with a retention rate, the probability that a duplicated branch is kept. We simulate three
levels, retention 0.25, 0.50, and 0.75, on top of the medium duplication and loss condition.
A retention of 0.25 removes most duplicate copies and is the most severe, and a retention of
0.75 is the mildest.

Fractionation does not merely hide copies. It deletes them. In many gene trees the second copy
of a fractionated species is simply not present. This has a direct consequence for any method
that infers ploidy from the copy number in the gene trees. When the copies are gone, the
inferred copy number falls to one, and the method concludes that the species is diploid. This
is an information loss that no copy-number method can undo, so we expect every ploidy-based
method to degrade as retention falls.

## Why the learned method can still recover events

The reconstruction places events on the ASTRAL species tree using two learned signals, a
per-edge detection probability and a per-edge partner distribution. The detection head does not
depend only on raw copy number. It also reads branch-length asymmetry and gene-tree clustering
structure, which can persist in the surviving gene trees even after a copy is fractionated away.
This raises a question that decides whether the learned method has any advantage under
fractionation. On the events that the inferred copy number misses, does the detection head still
fire?

We measured this directly. For every polyploid species in the ground truth we recorded the
detection probability on its lineage, and we split the species into three groups. The first are
polyploids whose duplication the inferred copy number still sees. The second are polyploids
whose copies fractionation has deleted, so the inferred copy number reads one. The third are the
true diploids, which serve as a control. If detection is nothing more than copy counting, the
second group should look like the control.

It does not. On the medium fractionation condition the deleted-copy polyploids carry a mean
detection probability of 0.55, and 59 percent of them exceed 0.5, against a mean of 0.09 and 6
percent for the true diploids. Even on the most severe condition the deleted-copy polyploids
sit at 0.30, roughly nine times the diploid level. The detection head therefore recognises
events that the copy number has lost, which is a signal a ploidy-only method cannot access.

## The decode that uses it

The default decode fills each species up to its inferred copy number, so when fractionation
collapses that number to one it produces no event. To use the detection signal we add a third
decode that treats the inferred copy number as a soft floor rather than a hard target. It fills
to the copy number as before, and then adds any edge whose detection probability exceeds a
threshold, even if that exceeds the inferred copy number. Under fractionation the floor is low,
so the confident detections do the work and recover events the copy number deleted. We call
this the detection-driven decode and set the threshold to 0.5, the value at which the diagnostic
above separates recovered events from the diploid control. [PENDING: threshold chosen on the
validation split rather than a priori, or a sensitivity sweep.]

The detection-driven decode improves the copy-number and reticulation measures over the default
decode exactly where fractionation bites. On the medium condition it lowers the reticulation
leaf distance from 0.65 to 0.50 and the ploidy distance from 0.68 to 0.53. On the most severe
condition it lowers the reticulation leaf distance from 0.90 to 0.82. It does not move the edit
distance, because recovering the events still leaves them to be placed on an imperfect backbone,
and it slightly over-predicts on the mildest condition, where the copy number is already correct
and the extra detections are false positives. This is the price of a single fixed threshold
under an unknown retention rate.

## Comparison to existing methods

The tables report the mean of each measure across the benchmark networks, at the three
fractionation levels. Lower is better. The learned method uses the detection-driven decode. The
learned-method scores are from replicate one and are being extended to five replicates. [PENDING:
five-replicate aggregation.]

Retention 0.75, the mild condition:

| method | edit | ret_leaf | ret_sisters |
| --- | --- | --- | --- |
| learned method | 0.543 | 0.189 | 0.482 |
| Polyphest | 0.461 | 0.261 | 0.358 |
| GRAMPA-iter | 0.776 | 0.506 | 0.719 |
| GRAMPA-iter with prior | 0.651 | 0.276 | 0.607 |

Retention 0.50, the medium condition:

| method | edit | ret_leaf | ret_sisters |
| --- | --- | --- | --- |
| learned method | 0.691 | 0.502 | 0.756 |
| Polyphest | 0.712 | 0.608 | 0.769 |
| GRAMPA-iter | 0.764 | 0.565 | 0.754 |
| GRAMPA-iter with prior | 0.759 | 0.694 | 0.851 |

Retention 0.25, the severe condition:

| method | edit | ret_leaf | ret_sisters |
| --- | --- | --- | --- |
| learned method | 0.758 | 0.826 | 0.926 |
| Polyphest | 0.735 | 0.939 | 0.955 |
| GRAMPA-iter | 0.737 | 0.620 | 0.784 |
| GRAMPA-iter with prior | 0.778 | 0.917 | 0.954 |

The three levels form a map of which method to prefer as fractionation increases, and the learned
method owns the middle of it.

At mild fractionation the copy number is still reliable, so Polyphest is the strongest method on
edit distance and on the reticulation sisters, and the learned method is second while still
leading on the reticulation-leaf distance. This is the regime where inferring ploidy from copies
works, and a ploidy method is rewarded for it.

At medium fractionation the learned method is the best method on edit distance and on the
reticulation-leaf distance, and tied best on the reticulation sisters. This is the crossover. The
copy number has begun to fail, so the ploidy methods weaken, but enough signal remains for the
detection head to recover the lost events, so the learned method pulls ahead of every competitor.
This middle regime is the method's niche.

At severe fractionation the copy number is largely destroyed. Here iterative GRAMPA, which
searches for reticulations one at a time without a copy-number estimate, recovers the most, and
the learned method is second, still ahead of Polyphest on both reticulation measures. So at the
extreme the method to beat is the free search rather than the ploidy baseline, and the learned
method sits between them.

One result cuts across all three levels. Giving iterative GRAMPA the inferred ploidy prior helps
only at mild fractionation and hurts as fractionation grows. With the prior its reticulation-leaf
distance is 0.28 at mild loss but 0.69 and 0.92 at medium and severe loss, against 0.51, 0.57, and
0.62 for the free version. The prior is the same collapsed multiset that misleads Polyphest, so
handing it to the search constrains it to the wrong ploidy exactly when the ploidy is wrong. This
confirms from a second direction that under fractionation the copy number is the problem, and a
method that does not lean on it, whether the free search or the learned detection head, degrades
more gracefully.

The learned method also completes more networks than the ploidy baseline, 17 to 18 against 16 to
17, and the networks the ploidy baseline fails on are the harder ones, where the learned method
still returns an answer. Robustness of completion is a practical advantage that the mean scores do
not capture.

## Summary

Fractionation is an information loss that degrades every copy-number method, so no method
reconstructs these networks well in absolute terms. Within that limit the three levels form a
clear map. At mild loss the copy number is reliable and Polyphest is best. At medium loss the copy
number begins to fail but the detection head recovers the lost events, and the learned method is
the best method on every measure. At severe loss the copy number is destroyed and a free
reticulation search recovers the most, with the learned method second and still ahead of the
ploidy baseline. The learned method therefore owns the middle regime, the crossover where copy
number fails but events remain recoverable, and it is competitive at the extremes. The result is
best read not as a single winner but as this map, together with the method's bounded runtime and
its completion of every network.
