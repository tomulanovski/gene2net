# Conclusion

This chapter developed a learned detect-then-place method for polyploid network reconstruction and
evaluated it against the existing methods on fifteen simulated configurations under the mu-distance.
The honest summary is a map of where each method is preferable rather than a single winner.

On the overall mu-distance Polyphest is the most accurate method wherever it completes, including on
the networks both methods finish, because the metric rewards the copy-number structure that folding
an inferred multiset recovers. The PlaceNet does not overturn that. What it offers instead is
a distinct and practical profile. It reconstructs every network in seconds, with a bounded and
predictable runtime, where Polyphest searches for days and often does not finish. It completes all
twenty-one benchmark networks in every configuration, including the harder ones that the
copy-number methods leave unfinished. It is the most accurate of the methods that infer their own
ploidy, beating iterative GRAMPA on the mu-distance in every configuration. And through its two
decode modes it adapts to whether the copy number can be trusted, so that in the corrupted regimes,
high duplication and loss and moderate fractionation, its ploidy-free mode recovers the reticulate
lineages more faithfully than the ploidy baseline, which is the regime that fractionation creates in
real polyploid data.

The diagnostic locates the method's remaining error precisely. On in-distribution networks the model
already predicts events almost as well as an oracle allows, and the dominant source of error is the
ASTRAL backbone, which the method never rebuilds. This turns the main limitation into a concrete
next step, namely rebuilding the backbone, for example through phasing, which the diagnostic
identifies as the only lever that can lower the in-distribution distance materially. Training on the
target distributions, in particular on fractionated data, is the second direction the results point
to, because the gap to the ploidy baseline is smallest exactly where the method already leads on
reticulation recovery and where it currently pays a distribution-shift cost.

The method is therefore best read not as a replacement for the accurate but slow search, but as a
fast, self-contained, and robust complement to it, together with a clear account of which method to
prefer under which conditions.
