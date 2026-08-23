# Results: comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The headline metric is the normalized mu-distance.
It is a distance on the phylogenetic network, so for scoring, the multi-labeled-tree output of every
method, including the present one, is folded to a network by the shared Holm-algorithm folding of
the evaluation, and the mu-distance is computed on the result. This folding is part of the scoring,
applied identically to every method, and not part of any method. Lower is better for every distance.
The PlaceNet appears in its two decode modes, ploidy-informed and ploidy-free. Every method
infers its own ploidy from the gene trees, so Polyphest here is the inferred-ploidy variant and the
comparison is prior-free on both sides. Scores are the mean over five replicates. The PlaceNet
is scored on all 21 networks. Each competitor mean is over the networks that competitor completes,
which is fewer, so the two sides are also compared on the common subset.

## The benchmark

We compare against the GRAMPA family and Polyphest on the 21 benchmark networks. The six
discordance configurations are shown first, namely three levels of incomplete lineage sorting and
three rates of duplication and loss at a fixed effective population size of one million. The two
other effective population sizes and the fractionation configurations follow and confirm the same
pattern. The table reports mean mu-distance.

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter |
| --- | --- | --- | --- | --- |
| ILS low | 0.268 | 0.270 | 0.021 | 0.321 |
| ILS medium | 0.290 | 0.300 | 0.036 | 0.335 |
| ILS high | 0.288 | 0.303 | 0.090 | 0.375 |
| dup/loss low | 0.284 | 0.290 | 0.043 | 0.350 |
| dup/loss medium | 0.287 | 0.293 | 0.059 | 0.340 |
| dup/loss high | 0.339 | 0.301 | 0.216 | 0.359 |

Five findings.

First, PlaceNet beats the GRAMPA family on mu in every configuration. Against
GRAMPA-iter, the fair peer that also searches for reticulations without a supplied ploidy, the
ploidy-free mode is lower in all six configurations, by about 0.05 to 0.09. It is also below
single-pass GRAMPA everywhere, and below iterative GRAMPA given the inferred ploidy prior, though
that last margin is narrow, for example 0.301 against 0.319 at high duplication and loss. So among
methods that reconstruct a full network without being told the ploidy, PlaceNet is the
most accurate.

Second, Polyphest is more accurate than PlaceNet on mu in every configuration, and the
gap is large at low difficulty. This traces to the metric and the method together. The mu-distance
rewards recovering each species' copy number and the reticulation structure it implies, and folding
Polyphest's inferred multiset recovers that structure closely on clean data. The advantage narrows as duplication and loss rise, from a
gap of about 0.25 at low sorting to about 0.09 at high duplication and loss, but it does not
reverse. This corrects an earlier reading under a different metric. Under the mu-distance, Polyphest
is the more accurate reconstruction wherever it completes.

Third, the common subset does not change this ordering. Restricting both methods to the networks
each completes leaves Polyphest ahead everywhere, for example 0.257 against 0.021 at low sorting
and 0.288 against 0.211 at high duplication and loss. So the ordering is not an artifact of the two
sides being scored on different subsets of networks.

Fourth, completion separates the methods in the other direction. The PlaceNet reconstructs
all 21 networks in every configuration. Polyphest completes 17 to 20 of the 21 at these six
configurations, and as few as 5 at the hardest configuration in the wider sweep below. On the
networks Polyphest does not complete, PlaceNet still returns an answer, and on the three
such networks at high duplication and loss its mu is 0.245, which is lower than its own mean on the
networks Polyphest does complete. So the networks Polyphest fails on are not the ones PlaceNet
finds hardest.

Fifth, the two decode modes trade off with copy-number reliability, and this is visible in the
table. On the clean and moderate configurations the two modes are close, with ploidy-informed
slightly ahead. At high duplication and loss the copy number is corrupted, and the ploidy-free mode
overtakes the ploidy-informed mode, from 0.339 to 0.301 on mu and far more on the reticulation-leaf
distance, from 0.782 to 0.298. This is the decode principle developed in the decode section, and it
is why the method offers two modes rather than one.

Sixth, decomposing the reticulation error into count and placement makes both the crossover and the
comparison to Polyphest precise. Alongside the reticulation-leaf Jaccard, which penalizes a method
for reticulations it never finds, we report a matched reticulation-leaf Jaccard that scores only the
reticulations a method commits to, and the reticulation-count difference that shows how many it
commits. Two facts emerge. The PlaceNet over-predicts reticulations relative to Polyphest,
with a count difference around 5 against Polyphest's 1 to 2, which is why its penalized
reticulation-leaf distance is higher than its matched one. But conditional on the reticulations it
commits to, its placement is strong. In the corrupted regimes the ploidy-free mode has the best
matched reticulation-leaf distance of any method, below Polyphest at every high duplication and loss
configuration, for example 0.118 against 0.180 at high duplication and loss with an effective
population size of one million, and 0.124 against 0.131 and 0.119 against 0.173 at the other two
population sizes. So the ploidy-free mode both places reticulations better than Polyphest where copy
number is corrupted and, from the count difference, over-predicts their number, and the two effects
are separated cleanly rather than conflated in one distance. The failure of the ploidy-informed mode
at high duplication and loss shows in both terms at once, a count difference near 10 and a matched
distance near 0.45, so the corrupted copy bound breaks it on number and placement together.

## The effective-population-size sweep

The six findings hold unchanged at effective population sizes of two hundred thousand and two
million. Polyphest is more accurate on mu at every configuration and population size. PlaceNet
completes all 21 networks throughout, while Polyphest's completion falls as difficulty
rises, reaching 5 of 21 at high duplication and loss at the largest population size. The ploidy-free
mode overtakes the ploidy-informed mode at every high duplication and loss configuration, with the
reticulation-leaf distance dropping from about 0.77 to about 0.30 in each case. So none of the
conclusions is specific to one effective population size, and the reticulation-recovery advantage of
the ploidy-free mode over Polyphest strengthens as sorting increases.

## Runtime

The PlaceNet reconstructs a network in seconds. Its per-network compute beyond the ASTRAL step
is about eight seconds on the benchmark, of which the graph neural network forward pass is a fraction
of a second and the remainder is feature extraction from the gene trees. ASTRAL adds a few seconds.
Polyphest, by contrast, has a runtime that ranges over four orders of magnitude and often does not
terminate at all. On one configuration's twenty-one networks its wall time ran from about four
minutes on the easiest network to more than three days on the hardest network that completed, and a
majority did not complete within reasonable limits. Five networks reached a five-day time limit
without finishing, nine exhausted memory, and several of those that did complete took between one
and three and a half days. So the method offers a bounded runtime of seconds per network against a
search whose runtime is variable, can reach days, and frequently ends in a timeout or an
out-of-memory failure rather than an answer.

## Summary

Read on a single axis the result is clear and not in the method's favor. On the mu-distance,
Polyphest is the more accurate reconstruction in every configuration, even on the common subset,
because the metric rewards the exact copy-number structure that folding an inferred multiset
recovers on clean data. Read across axes the method has a defensible place. It completes every
network where Polyphest completes 5 to 20, it runs in seconds where Polyphest can take days, it
beats GRAMPA-iter, the peer that also reconstructs without a ploidy input, on mu in every
configuration, and it recovers reticulate lineages better than Polyphest exactly where copy number
is corrupted. The contribution is therefore not accuracy supremacy on the mu-distance but a fast
method that completes every network and, through a decode that does not depend on the copy count,
recovers reticulate lineages robustly where fractionation corrupts it.
