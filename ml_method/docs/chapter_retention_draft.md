# Reconstruction under fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The numbers below are on the 21 benchmark
networks under three fractionation configurations, five replicates. The headline metric is the
normalized mu-distance, and the reticulation-leaf and reticulation-sister Jaccard distances
measure whether the reticulations involve the right lineages. Lower is better throughout.

## Fractionation and the hypothesis

After a whole genome duplication the two copies of a region do not always both survive. A process
of fractionation removes one copy from some of the descendant lineages. We model this with a
retention rate, the probability that a duplicated branch is kept. We simulate three levels,
retention 0.25, 0.50, and 0.75, on top of the medium duplication and loss condition. A retention
of 0.25 removes most duplicate copies and is the most severe, and a retention of 0.75 is the
mildest.

Fractionation does not merely hide copies. It deletes them. In many gene trees the second copy of a
fractionated species is simply not present. This has a direct consequence for any method that
infers ploidy from the copy number in the gene trees. When the copies are gone the inferred copy
number falls to one, and the method concludes that the species is diploid. This is an information
loss that no copy-number method can undo, so every ploidy-based method is expected to degrade as
retention falls.

## Why PlaceNet can still recover events

The reconstruction places events on the ASTRAL species tree using two learned signals, a per-edge
detection probability and a per-edge partner distribution. The detection head does not depend only
on raw copy number. It also reads branch-length asymmetry and gene-tree clustering structure, which
can persist in the surviving gene trees even after a copy is fractionated away. This raises the
question that decides whether PlaceNet has any advantage under fractionation. On the
events that the inferred copy number misses, does the detection head still fire?

We measured this directly. For every polyploid species in the ground truth we recorded the
detection probability on its lineage, and we split the species into three groups. The first are
polyploids whose duplication the inferred copy number still sees. The second are polyploids whose
copies fractionation has deleted, so the inferred copy number reads one. The third are the true
diploids, which serve as a control. If detection were nothing more than copy counting, the second
group would look like the control.

It does not. On the medium fractionation condition the deleted-copy polyploids carry a mean
detection probability of 0.55, and 59 percent of them exceed 0.5, against a mean of 0.09 and 6
percent for the true diploids. Even on the most severe condition the deleted-copy polyploids sit at
0.30, roughly nine times the diploid level. The detection head therefore recognises events that the
copy number has lost, which is a signal a ploidy-only method cannot access.

## The decode that uses it

This is where the two decode modes matter. The ploidy-informed mode fills each species up to its
inferred copy number, so when fractionation collapses that number to one it produces no event. The
ploidy-free mode ignores copy number entirely and keeps every edge the detection head is confident
about. Under fractionation the copy number is exactly the signal that has been destroyed, so the
ploidy-free mode is the one that can recover the deleted events, and the results below bear this
out. The threshold of the ploidy-free mode is fixed at the value used throughout the benchmark.

## Comparison to existing methods

The tables report the mean of each measure across the benchmark networks at the three fractionation
levels. Lower is better. Both decode modes of PlaceNet are shown. The reticulation-count
difference is the mean number of reticulations by which a method misses the truth, the
reticulation-leaf Jaccard penalizes a method for reticulations it never finds, and the matched
reticulation-leaf Jaccard scores only the reticulations a method commits to.

Retention 0.75, the mild condition:

| method | mu | num_rets | ret_leaf | ret_leaf matched |
| --- | --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.271 | 3.79 | 0.170 | 0.067 |
| PlaceNet ploidy-free | 0.285 | 4.24 | 0.242 | 0.105 |
| Polyphest | 0.161 | 2.46 | 0.233 | 0.170 |
| GRAMPA-iter | 0.345 | 4.42 | 0.499 | 0.272 |
| GRAMPA-iter with prior | 0.310 | 4.40 | 0.266 | 0.129 |

Retention 0.50, the medium condition:

| method | mu | num_rets | ret_leaf | ret_leaf matched |
| --- | --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.346 | 7.51 | 0.670 | 0.340 |
| PlaceNet ploidy-free | 0.322 | 6.07 | 0.506 | 0.227 |
| Polyphest | 0.254 | 5.94 | 0.594 | 0.299 |
| GRAMPA-iter | 0.343 | 5.11 | 0.556 | 0.350 |
| GRAMPA-iter with prior | 0.341 | 7.33 | 0.677 | 0.326 |

Retention 0.25, the severe condition:

| method | mu | num_rets | ret_leaf | ret_leaf matched |
| --- | --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.329 | 9.58 | 0.903 | 0.617 |
| PlaceNet ploidy-free | 0.326 | 8.51 | 0.832 | 0.561 |
| Polyphest | 0.277 | 8.58 | 0.922 | 0.670 |
| GRAMPA-iter | 0.354 | 5.40 | 0.615 | 0.402 |
| GRAMPA-iter with prior | 0.354 | 8.93 | 0.902 | 0.645 |

Two readings, one for each family of metric.

On the mu-distance Polyphest is the most accurate method at all three levels, 0.161, 0.254, and
0.277 as retention falls. This is the same pattern as the rest of the benchmark, and for the same
reason, namely that the mu-distance rewards the copy-number structure that folding a multiset
recovers. So on the overall metric no method beats Polyphest here.

On reticulation recovery, which is the question fractionation actually poses, the three levels form
a map of which method to prefer, and the decode modes move across it. The reticulation-leaf figures
here are the penalized Jaccard, which charges a method for reticulations it never finds. At mild
retention the copy number is still informative, and the ploidy-informed mode has the lowest
reticulation-leaf distance of any method, 0.170 against Polyphest's 0.233. At medium retention the
copy number has begun to fail, the ploidy-informed mode degrades to 0.670, and the ploidy-free mode
takes over at 0.506, which is the lowest of any method, below Polyphest's 0.594 and GRAMPA-iter's
0.556. The matched Jaccard, which scores only the reticulations each method commits to, tells the
same story at the crossover, with the ploidy-free mode at 0.227 against Polyphest's 0.299, so the
ploidy-free advantage at medium retention is in placement and not only in count. At severe retention
the copy number is largely destroyed, and iterative GRAMPA, which searches for reticulations one at
a time without a copy-number estimate, recovers the most at 0.615, with the ploidy-free mode second
at 0.832 and still ahead of Polyphest's 0.922. So the ploidy-free mode owns the middle of the map,
the crossover where copy number fails but detection can still recover the events, and it is
competitive at the extremes.

One result cuts across all three levels. Giving iterative GRAMPA the inferred ploidy prior helps
only at mild fractionation and hurts as fractionation grows. With the prior its reticulation-leaf
distance is 0.266 at mild loss but 0.677 and 0.902 at medium and severe loss, against 0.499, 0.556,
and 0.615 for the free version. The prior is the same collapsed multiset that limits Polyphest, so
handing it to the search constrains it to the wrong ploidy exactly when the ploidy is wrong. This
confirms from a second direction that under fractionation the copy number is the problem, and a
method that does not lean on it, whether the free search or the learned detection head in its
ploidy-free mode, degrades more gracefully.

The PlaceNet also completes every network at all three levels, 21 of 21, against Polyphest's
18, 19, and 19. Robustness of completion is a practical advantage that the mean scores do not
capture.

## Summary

Fractionation is an information loss that degrades every copy-number method, so no method
reconstructs these networks well in absolute terms, and Polyphest remains the most accurate on the
overall mu-distance because that metric rewards the copy-number structure. The contribution here is
on reticulation recovery and on the decode. The three retention levels form a clear map. At mild
loss the copy number is reliable and the ploidy-informed mode recovers reticulate lineages best. At
medium loss the copy number begins to fail but the detection head recovers the deleted events, and
the ploidy-free mode is the best method on the reticulation-leaf distance. At severe loss the copy
number is destroyed and a free reticulation search recovers the most, with the ploidy-free mode
second and still ahead of the ploidy baseline. The PlaceNet therefore owns the crossover
regime on reticulation recovery, degrades gracefully because its ploidy-free mode does not lean on
the copy number, and completes every network throughout.
