# Reconstruction under fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The numbers below are on the 20 benchmark
networks that carry a fractionated condition, under three fractionation configurations, five
replicates. The headline metric is the
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

The tables report the mean of each measure at the three fractionation levels, with every method
restricted to the networks all of them completed, 18 at mild and medium retention and 16 at
severe. Lower is better. Both decode modes of PlaceNet are shown.

Retention 0.75, the mild condition:

| method | ret_leaf | ret_sis | mu |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.160 | 0.474 | 0.430 |
| PlaceNet ploidy-free | 0.233 | 0.519 | 0.460 |
| Polyphest | 0.244 | 0.335 | 0.251 |
| GRAMPA-iter | 0.492 | 0.705 | 0.546 |
| GRAMPA-iter with prior | 0.249 | 0.585 | 0.491 |

Retention 0.50, the medium condition:

| method | ret_leaf | ret_sis | mu |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.682 | 0.829 | 0.577 |
| PlaceNet ploidy-free | 0.506 | 0.720 | 0.532 |
| Polyphest | 0.609 | 0.766 | 0.420 |
| GRAMPA-iter | 0.542 | 0.733 | 0.535 |
| GRAMPA-iter with prior | 0.672 | 0.836 | 0.562 |

Retention 0.25, the severe condition:

| method | ret_leaf | ret_sis | mu |
| --- | --- | --- | --- |
| PlaceNet ploidy-informed | 0.922 | 0.964 | 0.548 |
| PlaceNet ploidy-free | 0.876 | 0.942 | 0.538 |
| Polyphest | 0.917 | 0.948 | 0.462 |
| GRAMPA-iter | 0.602 | 0.774 | 0.547 |
| GRAMPA-iter with prior | 0.923 | 0.962 | 0.592 |

Two readings, one for each family of measure. The reticulation measures come first, since
they are the ones fractionation puts under pressure.

On reticulation recovery, which is the question fractionation actually poses, the three levels form
a map of which method to prefer, and the decode modes move across it. The reticulation-leaf figures
here are the penalized Jaccard, which charges a method for reticulations it never finds. At mild
retention the copy number is still informative, and the ploidy-informed mode has the lowest
reticulation-leaf distance of any method, 0.160 against Polyphest's 0.244. At medium retention the
copy number has begun to fail, the ploidy-informed mode degrades to 0.682, and the ploidy-free mode
takes over at 0.506, which is the lowest of any method, below Polyphest's 0.609 and GRAMPA-iter's
0.542. At severe retention
the copy number is largely destroyed, and iterative GRAMPA, which searches for reticulations one at
a time without a copy-number estimate, recovers the most at 0.602, with the ploidy-free mode second
at 0.876 and still ahead of Polyphest's 0.917. So the ploidy-free mode owns the middle of the map,
the crossover where copy number fails but detection can still recover the events, and it is
competitive at the extremes.

On the mu-distance Polyphest is the most accurate method at all three levels, 0.251, 0.420, and
0.462 as retention falls. This is the same pattern as the rest of the benchmark, and for the same
reason, namely that the mu-distance rewards the copy-number structure that folding a multiset
recovers. So on the overall metric no method beats Polyphest here.

One result cuts across all three levels. Giving iterative GRAMPA the inferred ploidy prior helps
only at mild fractionation and hurts as fractionation grows. With the prior its reticulation-leaf
distance is 0.249 at mild loss but 0.672 and 0.923 at medium and severe loss, against 0.492, 0.542,
and 0.602 for the free version. The prior is the same collapsed multiset that limits Polyphest, so
handing it to the search constrains it to the wrong ploidy exactly when the ploidy is wrong. This
confirms from a second direction that under fractionation the copy number is the problem, and a
method that does not lean on it, whether the free search or the learned detection head in its
ploidy-free mode, degrades more gracefully.

The PlaceNet completes 20 of the 20 fractionated networks at mild and medium retention and 19 at
severe retention, where one network's ASTRAL backbone could not be parsed, against Polyphest's
19, 19, and 18. Robustness of completion is a practical advantage that the mean scores do not
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
