# Results: comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The headline metric is the normalized
mu-distance, a distance on the folded network that does not depend on Newick child order and
equals zero exactly when two networks are isomorphic. Lower is better for every distance. The
learned method appears in its two decode modes, ploidy-informed and ploidy-free. Every method
infers its own ploidy from the gene trees, so Polyphest here is the inferred-ploidy variant and
the comparison is prior-free on both sides. Scores are the mean over five replicates. The learned
method is scored on all 21 networks. Each competitor mean is over the networks that competitor
completes, which is fewer, so the two sides are also compared on the common subset.

## The benchmark

We compare against the GRAMPA family and Polyphest on the 21 benchmark networks. The six
discordance configurations are shown first, namely three levels of incomplete lineage sorting and
three rates of duplication and loss at a fixed effective population size of one million. The two
other effective population sizes and the fractionation configurations follow and confirm the same
pattern. The table reports mean mu-distance.

| Configuration | GNN informed | GNN free | Polyphest | GRAMPA-iter |
| --- | --- | --- | --- | --- |
| ILS low | 0.268 | 0.270 | 0.021 | 0.321 |
| ILS medium | 0.290 | 0.300 | 0.036 | 0.335 |
| ILS high | 0.288 | 0.303 | 0.090 | 0.375 |
| dup/loss low | 0.284 | 0.290 | 0.043 | 0.350 |
| dup/loss medium | 0.287 | 0.293 | 0.059 | 0.340 |
| dup/loss high | 0.339 | 0.301 | 0.216 | 0.359 |

Five findings.

First, the learned method beats the GRAMPA family on mu in every configuration. Against
GRAMPA-iter, the fair prior-free peer that also searches for reticulations without a ploidy input,
the ploidy-free mode is lower in all six configurations, by about 0.05 to 0.09. It is also below
single-pass GRAMPA everywhere, and below iterative GRAMPA given the inferred ploidy prior, though
that last margin is narrow, for example 0.301 against 0.319 at high duplication and loss. So among
methods that reconstruct a full network without being told the ploidy, the learned method is the
most accurate.

Second, Polyphest is more accurate than the learned method on mu in every configuration, and the
gap is large at low difficulty. This traces to the metric and the method together. The mu-distance
is a distance on the folded network, where recovering each species' copy number and the reticulation
structure it implies dominates the count, and Polyphest's folding of its inferred multiset recovers
that structure closely on clean data. The advantage narrows as duplication and loss rise, from a
gap of about 0.25 at low sorting to about 0.09 at high duplication and loss, but it does not
reverse. This corrects an earlier reading under a different metric. Under the mu-distance, Polyphest
is the more accurate reconstruction wherever it completes.

Third, the common subset does not change this ordering. Restricting both methods to the networks
each completes leaves Polyphest ahead everywhere, for example 0.257 against 0.021 at low sorting
and 0.288 against 0.211 at high duplication and loss. So the ordering is not an artifact of the two
sides being scored on different subsets of networks.

Fourth, completion separates the methods in the other direction. The learned method reconstructs
all 21 networks in every configuration. Polyphest completes 17 to 20 of the 21 at these six
configurations, and as few as 5 at the hardest configuration in the wider sweep below. On the
networks Polyphest does not complete, the learned method still returns an answer, and on the three
such networks at high duplication and loss its mu is 0.245, which is lower than its own mean on the
networks Polyphest does complete. So the networks Polyphest fails on are not the ones the learned
method finds hardest.

Fifth, the two decode modes trade off with copy-number reliability, and this is visible in the
table. On the clean and moderate configurations the two modes are close, with ploidy-informed
slightly ahead. At high duplication and loss the copy number is corrupted, and the ploidy-free mode
overtakes the ploidy-informed mode, from 0.339 to 0.301 on mu and far more on the reticulation-leaf
distance, from 0.782 to 0.298. This is the decode principle developed in the decode section, and it
is why the method offers two modes rather than one. The same reticulation-leaf result places the
ploidy-free mode ahead of Polyphest on reticulation recovery at high duplication and loss, 0.298
against 0.377, even though Polyphest's overall mu is lower there.

## The effective-population-size sweep

The six findings hold unchanged at effective population sizes of two hundred thousand and two
million. Polyphest is more accurate on mu at every configuration and population size. The learned
method completes all 21 networks throughout, while Polyphest's completion falls as difficulty
rises, reaching 5 of 21 at high duplication and loss at the largest population size. The ploidy-free
mode overtakes the ploidy-informed mode at every high duplication and loss configuration, with the
reticulation-leaf distance dropping from about 0.77 to about 0.30 in each case. So none of the
conclusions is specific to one effective population size, and the reticulation-recovery advantage of
the ploidy-free mode over Polyphest strengthens as sorting increases.

## Runtime

The learned method reconstructs a network in seconds. Its per-network compute beyond the ASTRAL step
is about eight seconds on the benchmark, of which the graph neural network forward pass is a fraction
of a second and the remainder is feature extraction from the gene trees. ASTRAL adds a few seconds.
Polyphest ranges from seconds on easy networks to days on hard ones, and does not complete the
networks it is missing from the tables, with individual runs still unfinished after several days on
the hardest configurations. [TODO insert the concrete Polyphest wall-time range from the run logs.]
So the method offers a bounded runtime of seconds per network against a search whose runtime is
variable and can be prohibitive.

## Summary

Read on a single axis the result is clear and not in the method's favor. On the mu-distance,
Polyphest is the more accurate reconstruction in every configuration, even on the common subset,
because the metric rewards the exact copy-number structure that folding an inferred multiset
recovers on clean data. Read across axes the method has a defensible place. It completes every
network where Polyphest completes 5 to 20, it runs in seconds where Polyphest can take days, it
beats the fair prior-free baseline GRAMPA-iter on mu in every configuration, and it recovers
reticulate lineages better than Polyphest exactly where copy number is corrupted. The contribution
is therefore not accuracy supremacy on the mu-distance but a fast, complete, prior-free method with
a clear niche, together with a decode that adapts to whether copy number can be trusted.
