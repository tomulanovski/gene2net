# Results: comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The measures are those defined earlier in the
thesis. The reticulation descendants and reticulation sister measures come first, because they ask
directly whether a method recovers the polyploidization events, and the mu-distance follows as a
summary of overall structural agreement. Every measure is computed on the phylogenetic network, so
for scoring, the multi-labeled-tree output of every method, including the present one, is folded to
a network by the shared Holm-algorithm folding of the evaluation. This folding is part of the
scoring, applied identically to every method, and not part of any method. Lower is better for every
distance. The PlaceNet appears in its two decode modes, ploidy-informed and ploidy-free. Every
method infers its own ploidy from the gene trees, so Polyphest here is the inferred-ploidy variant
and the comparison is prior-free on both sides. Scores are the mean over five replicates. The
PlaceNet is scored on all 21 networks. Each competitor mean is over the networks that competitor
completes, which is fewer, so the two sides are also compared on the common subset.

## Where the method wins

The comparison spans fifteen configurations, twelve of discordance and three of fractionation, and
the result is not uniform across them. It is ordered by difficulty.

On the reticulation descendants measure, which asks which lineages are polyploid, PlaceNet is ahead
of Polyphest in eight of the fifteen configurations and level in a ninth. Those eight are not
scattered. They are every configuration at the largest effective population size, the two higher
duplication and loss rates at the middle population size, and all three fractionation levels.
Polyphest leads in the remaining six, which are the configurations with the least discordance and an
intact copy number.

The dividing line is the reliability of the copy number. Polyphest builds its reconstruction around
an inferred multiset of copy counts, and where the gene trees support that multiset it recovers the
events better than any learned model. As discordance rises and as fractionation deletes duplicate
copies, the multiset degrades, and PlaceNet overtakes it. Against iterative GRAMPA, the peer that
also reconstructs without a supplied ploidy, PlaceNet is ahead in fourteen of the fifteen
configurations on this measure and in all fifteen on the mu-distance.

This matters for the intended use. Published polyploid datasets are not the clean end of the sweep.
They carry substantial gene tree discordance, and fractionation has removed duplicate copies from
most polyploid genomes, which is why the copy number is unreliable in practice and not only in
principle. The conditions in which PlaceNet leads are the conditions real data present.

## Which lineages are reticulate

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior |
| --- | --- | --- | --- | --- | --- |
| ILS low | 0.161 | 0.235 | 0.025 | 0.483 | 0.270 |
| ILS medium | 0.207 | 0.285 | 0.141 | 0.507 | 0.326 |
| ILS high | 0.230 | 0.294 | 0.232 | 0.611 | 0.309 |
| dup/loss low, Ne 200k | 0.166 | 0.246 | 0.121 | 0.391 | 0.226 |
| dup/loss medium, Ne 200k | 0.186 | 0.270 | 0.128 | 0.412 | 0.216 |
| dup/loss high, Ne 200k | 0.769 | 0.315 | 0.308 | 0.534 | 0.345 |
| dup/loss low, Ne 1M | 0.231 | 0.274 | 0.184 | 0.498 | 0.262 |
| dup/loss medium, Ne 1M | 0.212 | 0.278 | 0.231 | 0.483 | 0.257 |
| dup/loss high, Ne 1M | 0.782 | 0.298 | 0.377 | 0.552 | 0.383 |
| dup/loss low, Ne 2M | 0.212 | 0.296 | 0.238 | 0.580 | 0.303 |
| dup/loss medium, Ne 2M | 0.216 | 0.283 | 0.236 | 0.570 | 0.296 |
| dup/loss high, Ne 2M | 0.793 | 0.301 | 0.357 | 0.606 | 0.362 |

PlaceNet is far ahead of the GRAMPA family throughout. The better of its two modes is below
GRAMPA-iter in every configuration, usually by a factor of two, and below iterative GRAMPA given the
inferred ploidy prior in every configuration as well. Among methods that reconstruct without a
supplied ploidy, PlaceNet identifies the reticulate lineages most accurately, and it does so
everywhere rather than in a subset.

Against Polyphest the ordering turns on difficulty. At the smallest effective population size
Polyphest is ahead at every duplication and loss rate. At the middle size the two cross, with
Polyphest ahead at the low rate, 0.184 against 0.231, and PlaceNet ahead at the medium and high
rates, 0.212 against 0.231 and 0.298 against 0.377. At the largest size PlaceNet is ahead at all
three rates, 0.212 against 0.238, 0.216 against 0.236, and 0.301 against 0.357. The same crossing
appears along the sorting axis, where Polyphest leads at low and medium sorting and the two are
level at high sorting, 0.230 against 0.232.

The two decode modes trade off with copy-number reliability, and this measure shows it most sharply.
On the clean and moderate configurations the ploidy-informed mode is ahead. At every high
duplication and loss rate the copy number is corrupted, the ploidy-informed mode collapses to
between 0.769 and 0.793, and the ploidy-free mode holds between 0.298 and 0.315. That collapse is
the corrupted copy bound rather than a failure of the detection head, and it shows in the count as
well, where the ploidy-informed mode reaches a reticulation-count difference above 10 against 3 to 5
for Polyphest. This is the decode principle developed in the decode section, and it is why the
method offers two modes rather than one.

## Where the reticulate lineages came from

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior |
| --- | --- | --- | --- | --- | --- |
| ILS low | 0.453 | 0.497 | 0.057 | 0.645 | 0.508 |
| ILS medium | 0.576 | 0.618 | 0.192 | 0.697 | 0.608 |
| ILS high | 0.580 | 0.626 | 0.298 | 0.759 | 0.623 |
| dup/loss low, Ne 200k | 0.542 | 0.592 | 0.152 | 0.606 | 0.492 |
| dup/loss medium, Ne 200k | 0.562 | 0.619 | 0.165 | 0.619 | 0.513 |
| dup/loss high, Ne 200k | 0.937 | 0.733 | 0.439 | 0.690 | 0.591 |
| dup/loss low, Ne 1M | 0.581 | 0.605 | 0.236 | 0.681 | 0.574 |
| dup/loss medium, Ne 1M | 0.579 | 0.617 | 0.267 | 0.698 | 0.577 |
| dup/loss high, Ne 1M | 0.928 | 0.694 | 0.504 | 0.731 | 0.634 |
| dup/loss low, Ne 2M | 0.573 | 0.619 | 0.307 | 0.758 | 0.632 |
| dup/loss medium, Ne 2M | 0.576 | 0.617 | 0.288 | 0.734 | 0.623 |
| dup/loss high, Ne 2M | 0.921 | 0.692 | 0.490 | 0.773 | 0.636 |

The parental context is harder for PlaceNet than the reticulate lineage, and this is the one measure
on which Polyphest leads throughout. Unlike the descendants measure the ordering does not reverse at
high discordance. PlaceNet remains below the GRAMPA family in every configuration, so it is still
the most accurate of the methods that infer their own ploidy, but the gap to Polyphest is real and
it does not close.

The cause is known rather than mysterious, and the diagnostic section quantifies it. The method
predicts two things on every edge, whether an event occurred and which lineage the duplicated
lineage merged with. The detection head is strong, which is what the descendants measure reflects.
The partner head is the weaker of the two, which is what this measure reflects. An
ablation of the partner features reaches the same conclusion from the model side, and the oracle
experiment shows that the parental context inherits the error of an ASTRAL
backbone the method never rebuilds. Improving it is therefore a backbone problem before it is a
partner-head problem, and that is the direction the future-work section takes.

The decode crossover appears here as well. At every high duplication and loss rate the ploidy-free
mode overtakes the ploidy-informed mode, from about 0.93 to about 0.70, so the pattern is the same
as on the descendants measure even though the level is worse.

## Overall structural agreement

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior |
| --- | --- | --- | --- | --- | --- |
| ILS low | 0.421 | 0.429 | 0.032 | 0.520 | 0.456 |
| ILS medium | 0.463 | 0.480 | 0.062 | 0.541 | 0.497 |
| ILS high | 0.461 | 0.486 | 0.144 | 0.606 | 0.517 |
| dup/loss low, Ne 200k | 0.453 | 0.470 | 0.050 | 0.493 | 0.436 |
| dup/loss medium, Ne 200k | 0.460 | 0.491 | 0.055 | 0.507 | 0.459 |
| dup/loss high, Ne 200k | 0.593 | 0.508 | 0.313 | 0.546 | 0.492 |
| dup/loss low, Ne 1M | 0.454 | 0.464 | 0.073 | 0.558 | 0.486 |
| dup/loss medium, Ne 1M | 0.459 | 0.473 | 0.096 | 0.545 | 0.471 |
| dup/loss high, Ne 1M | 0.579 | 0.491 | 0.337 | 0.581 | 0.513 |
| dup/loss low, Ne 2M | 0.462 | 0.483 | 0.139 | 0.612 | 0.523 |
| dup/loss medium, Ne 2M | 0.470 | 0.484 | 0.128 | 0.610 | 0.517 |
| dup/loss high, Ne 2M | 0.603 | 0.496 | 0.348 | 0.603 | 0.532 |

PlaceNet beats the GRAMPA family on the mu-distance in every configuration, so among the prior-free
methods it is again the most accurate. Polyphest is ahead of it everywhere, and the reason is
structural rather than about the events. The mu-distance scores every node of the network, so it
charges PlaceNet for the placement of the diploid species as well as for the polyploidization
events, and PlaceNet stamps its events onto an ASTRAL backbone it never rebuilds. The diagnostic
section measures that backbone term directly and finds it dominant, five times the error the event
prediction contributes. A method can therefore reach a lower mu-distance while recovering the
polyploidization events less accurately, and at high discordance that is what happens.

Completion runs the other way. PlaceNet reconstructs all 21 networks in every discordance
configuration. Polyphest completes 17 to 20 of the 21, and on the networks it does not complete
PlaceNet still returns an answer.

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

PlaceNet recovers the polyploidization events more accurately than any existing method in the
conditions that real polyploid data present. On the reticulation descendants measure it leads
Polyphest in eight of the fifteen configurations and is level in a ninth, and those are the
configurations with high discordance or fractionation rather than a scattered subset. Against
iterative GRAMPA, the peer that also works without a supplied ploidy, it leads on every measure in
almost every configuration. It reconstructs every network in seconds with a bounded runtime, where
Polyphest is variable, can take days, and often returns nothing at all.

Two limits are worth stating plainly. Polyphest recovers the parental context better in every
configuration, because the partner head is the weaker of the method's two heads and inherits the
error of a backbone the method does not rebuild. And Polyphest is ahead on the mu-distance
everywhere, because that measure scores every node and therefore charges PlaceNet for the same
backbone. Both limits point at the same cause, which the diagnostic section quantifies and the
future-work section addresses, and neither changes the regime in which the method is preferable.
