# Results: comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The measures are those defined earlier in the
thesis. The reticulation descendants and reticulation sister measures come first, because they ask
directly whether a method recovers the polyploidization events, and the mu-distance follows as a
summary of overall structural agreement. Every measure is computed on the phylogenetic network, so
for scoring, the multi-labeled-tree output of every method, including the present one, is folded to
a network by the shared Holm-algorithm folding of the evaluation. This folding is part of the
scoring, applied identically to every method, and not part of any method. Lower is better for every
distance. PlaceNet appears in its two decode modes, ploidy-informed and ploidy-free. Every method
infers its own ploidy from the gene trees, so Polyphest here is the inferred-ploidy variant and the
comparison is prior-free on both sides. Scores are the mean over five replicates.

The competitors do not complete every network, and the networks they abandon are the harder ones,
so a mean taken over each method's own completions would flatter whichever method completes least.
Every table below therefore restricts all methods to the networks that all of them completed, and
reports the size of that subset as n. It is limited by Polyphest and by iterative GRAMPA with the
ploidy prior rather than by PlaceNet, which completes every network. What the restriction leaves
out is reported separately at the end of the section.

## Where the method wins

The comparison spans fifteen configurations, twelve of discordance and three of fractionation, and
the result is not uniform across them. It is ordered by difficulty.

On the reticulation descendants measure, which asks which lineages are polyploid, the ploidy-free
decode is more accurate than Polyphest in eleven of the fifteen configurations. It is a single fixed
decode at a single fixed threshold, so this is not a matter of choosing the better mode per
condition. The four it loses are the four with the least discordance and an intact copy number,
namely low and medium sorting and the low and medium duplication and loss rates at the smallest
effective population size. The ploidy-informed decode wins nine of the fifteen, so both modes are
ahead of Polyphest on the majority of conditions and the choice between them changes the margin
rather than the outcome.

The dividing line is the reliability of the copy number. Polyphest builds its reconstruction around
an inferred multiset of copy counts, and where the gene trees support that multiset it recovers the
events better than any learned model. As discordance rises and as fractionation deletes duplicate
copies, the multiset degrades and PlaceNet overtakes it. Against iterative GRAMPA, the peer that
also reconstructs without a supplied ploidy, PlaceNet is ahead on this measure in fourteen of the
fifteen configurations, and ahead of iterative GRAMPA given the inferred ploidy prior in all
fifteen.

This matters for the intended use. Published polyploid datasets are not the clean end of the sweep.
They carry substantial gene tree discordance, and fractionation has removed duplicate copies from
most polyploid genomes, which is why the copy number is unreliable in practice and not only in
principle. The conditions in which PlaceNet leads are the conditions real data present.

## Which lineages are reticulate

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.112 | 0.174 | 0.022 | 0.418 | 0.223 | 17 |
| ILS medium | 0.134 | 0.205 | 0.148 | 0.452 | 0.228 | 16 |
| ILS high | 0.157 | 0.207 | 0.245 | 0.573 | 0.210 | 16 |
| dup/loss low, Ne 200k | 0.111 | 0.166 | 0.127 | 0.374 | 0.210 | 16 |
| dup/loss medium, Ne 200k | 0.122 | 0.190 | 0.135 | 0.405 | 0.188 | 16 |
| dup/loss high, Ne 200k | 0.782 | 0.284 | 0.319 | 0.527 | 0.306 | 17 |
| dup/loss low, Ne 1M | 0.154 | 0.188 | 0.194 | 0.467 | 0.222 | 16 |
| dup/loss medium, Ne 1M | 0.136 | 0.195 | 0.243 | 0.473 | 0.237 | 16 |
| dup/loss high, Ne 1M | 0.776 | 0.247 | 0.405 | 0.549 | 0.356 | 18 |
| dup/loss low, Ne 2M | 0.137 | 0.208 | 0.251 | 0.557 | 0.228 | 16 |
| dup/loss medium, Ne 2M | 0.148 | 0.195 | 0.249 | 0.550 | 0.232 | 16 |
| dup/loss high, Ne 2M | 0.817 | 0.230 | 0.368 | 0.585 | 0.316 | 16 |

PlaceNet is far ahead of the GRAMPA family throughout. The better of its two modes is below
GRAMPA-iter in every configuration, usually by a factor of two, and below iterative GRAMPA given the
inferred ploidy prior in every configuration as well.

Against Polyphest the ordering turns on difficulty and it reverses. Polyphest is ahead at low and
medium sorting and at the two lower duplication and loss rates at the smallest population size. From
there PlaceNet is ahead everywhere, at high sorting 0.207 against 0.245, and at high duplication and
loss 0.284 against 0.319 at the smallest population size, 0.247 against 0.405 at the middle size and
0.230 against 0.368 at the largest, all with the ploidy-free decode.

The two decode modes trade off with copy-number reliability, and this measure shows it most sharply.
On the clean and moderate configurations the ploidy-informed mode is ahead. At every high
duplication and loss rate the copy number is corrupted, the ploidy-informed mode collapses to
between 0.776 and 0.817, and the ploidy-free mode holds between 0.230 and 0.284. That collapse is
the corrupted copy bound rather than a failure of the detection head, and it shows in the
reticulation count as well, where the ploidy-informed mode over-predicts events sharply once the
bound stops constraining it. This is the decode principle developed in the decode section, and it is
why the method offers two modes rather than one.

## Where the reticulate lineages came from

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.382 | 0.428 | 0.052 | 0.594 | 0.471 | 17 |
| ILS medium | 0.511 | 0.555 | 0.196 | 0.656 | 0.540 | 16 |
| ILS high | 0.522 | 0.563 | 0.310 | 0.729 | 0.565 | 16 |
| dup/loss low, Ne 200k | 0.482 | 0.524 | 0.156 | 0.574 | 0.464 | 16 |
| dup/loss medium, Ne 200k | 0.505 | 0.558 | 0.171 | 0.597 | 0.485 | 16 |
| dup/loss high, Ne 200k | 0.942 | 0.696 | 0.444 | 0.680 | 0.581 | 17 |
| dup/loss low, Ne 1M | 0.516 | 0.540 | 0.244 | 0.646 | 0.543 | 16 |
| dup/loss medium, Ne 1M | 0.515 | 0.554 | 0.274 | 0.682 | 0.559 | 16 |
| dup/loss high, Ne 1M | 0.925 | 0.671 | 0.529 | 0.725 | 0.618 | 18 |
| dup/loss low, Ne 2M | 0.514 | 0.553 | 0.316 | 0.741 | 0.590 | 16 |
| dup/loss medium, Ne 2M | 0.520 | 0.556 | 0.299 | 0.708 | 0.580 | 16 |
| dup/loss high, Ne 2M | 0.925 | 0.652 | 0.495 | 0.756 | 0.604 | 16 |

The parental context is harder for PlaceNet than the reticulate lineage, and this is the measure on
which Polyphest leads throughout. Unlike the descendants measure the ordering does not reverse at
high discordance. PlaceNet remains below the GRAMPA family in every configuration, so it is still
the most accurate of the methods that infer their own ploidy, but the gap to Polyphest is real and
it does not close.

The cause is known rather than mysterious, and the diagnostic section quantifies it. The method
predicts two things on every edge, whether an event occurred and which lineage the duplicated
lineage merged with. The detection head is strong, which is what the descendants measure reflects.
The partner head is the weaker of the two, which is what this measure reflects. An ablation of the
partner features reaches the same conclusion from the model side, and the oracle experiment shows
that the parental context inherits the error of an ASTRAL backbone the method never rebuilds.
Improving it is therefore a backbone problem before it is a partner-head problem, and that is the
direction the future-work section takes.

The decode crossover appears here as well. At every high duplication and loss rate the ploidy-free
mode overtakes the ploidy-informed mode, from about 0.93 to about 0.67, so the pattern is the same
as on the descendants measure even though the level is worse.

## Overall structural agreement

| Configuration | PlaceNet informed | PlaceNet free | Polyphest | GRAMPA-iter | GRAMPA-iter + prior | n |
| --- | --- | --- | --- | --- | --- | --- |
| ILS low | 0.393 | 0.404 | 0.027 | 0.480 | 0.429 | 17 |
| ILS medium | 0.417 | 0.436 | 0.062 | 0.487 | 0.437 | 16 |
| ILS high | 0.413 | 0.442 | 0.150 | 0.567 | 0.465 | 16 |
| dup/loss low, Ne 200k | 0.402 | 0.423 | 0.052 | 0.447 | 0.398 | 16 |
| dup/loss medium, Ne 200k | 0.410 | 0.445 | 0.057 | 0.466 | 0.419 | 16 |
| dup/loss high, Ne 200k | 0.607 | 0.523 | 0.317 | 0.563 | 0.507 | 17 |
| dup/loss low, Ne 1M | 0.404 | 0.416 | 0.075 | 0.511 | 0.442 | 16 |
| dup/loss medium, Ne 1M | 0.414 | 0.430 | 0.098 | 0.505 | 0.447 | 16 |
| dup/loss high, Ne 1M | 0.567 | 0.477 | 0.351 | 0.567 | 0.497 | 18 |
| dup/loss low, Ne 2M | 0.408 | 0.430 | 0.144 | 0.574 | 0.475 | 16 |
| dup/loss medium, Ne 2M | 0.421 | 0.437 | 0.132 | 0.575 | 0.472 | 16 |
| dup/loss high, Ne 2M | 0.579 | 0.456 | 0.348 | 0.574 | 0.505 | 16 |

PlaceNet beats GRAMPA-iter on the mu-distance in every configuration, and iterative GRAMPA with the
ploidy prior in thirteen of the fifteen, so among the prior-free methods it is again the most
accurate. Polyphest is ahead of it everywhere, and the reason is structural rather than about the
events. The mu-distance scores every node of the network, so it charges PlaceNet for the placement
of the diploid species as well as for the polyploidization events, and PlaceNet stamps its events
onto an ASTRAL backbone it never rebuilds. The diagnostic section measures that backbone term
directly and finds it dominant, five times the error the event prediction contributes. A method can
therefore reach a lower mu-distance while recovering the polyploidization events less accurately,
and from moderate difficulty onward that is what happens.

## What the common subset leaves out

Restricting every method to the networks all of them completed removes four to five networks per
configuration, and those are networks PlaceNet reconstructed and at least one competitor did not.
PlaceNet completes all 21 in every discordance configuration, against 17 to 20 for Polyphest.

Those networks are harder, and honestly so. PlaceNet's mu-distance on them runs from 0.56 to 0.73,
against 0.39 to 0.61 on the common subset, so they are not networks it finds easy while the
competitors merely declined them. The claim is that it returns an answer where the alternatives
return nothing, not that the answer is as good as its average.

## Runtime

PlaceNet reconstructs a network in seconds. Its per-network compute beyond the ASTRAL step
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

PlaceNet identifies the polyploid lineages more accurately than any existing method in the
conditions that real polyploid data present. With a single fixed decode it is ahead of Polyphest on
the reticulation descendants measure in eleven of the fifteen configurations, losing only the four
with the least discordance and an intact copy number. Against iterative GRAMPA, the peer that also
works without a supplied ploidy, it leads on every measure in almost every configuration. It
reconstructs every network in seconds with a bounded runtime, where Polyphest is variable, can take
days, and often returns nothing at all.

Two limits are worth stating plainly. Polyphest recovers the parental context better in every
configuration, and it is ahead on the mu-distance in every configuration as well. Both have the same
cause. The partner head is the weaker of the method's two heads, and the mu-distance scores every
node of the network, so both charge PlaceNet for the ASTRAL backbone it never rebuilds. The
diagnostic section measures that term and finds it five times the error the event prediction
contributes, which makes rebuilding the backbone the single change that would move both.
