# Introduction

The preceding chapters established the problem of reconstructing polyploid phylogenetic networks
from gene trees and benchmarked the existing methods. Two limitations of those methods motivate the
approach developed here. The first is runtime. The most accurate method, Polyphest, searches over
network space and can take days on a single network, and on the hardest networks it does not finish
at all. The second is a dependence on copy number. The leading methods infer or are given a per
species ploidy and build their reconstruction around it, so when the copy number in the gene trees
is corrupted, by high duplication and loss or by fractionation after a whole genome duplication,
their reconstruction degrades with it.

This chapter develops a learned alternative. The method is a graph neural network that follows a
detect-then-place strategy. It builds a single-copy species tree backbone with ASTRAL, reads that
backbone together with features summarising the gene trees, and predicts on every backbone edge
whether a whole genome duplication occurred there and, if so, which lineage is the second parent.
The predicted events are stamped onto the backbone and folded into a network. The method infers its
own per species copy number from the gene trees rather than being supplied a ploidy, so it is
directly comparable to the prior-free iterative GRAMPA baseline.

The contribution is not accuracy supremacy on the overall distance measure. Polyphest, given its
inferred ploidy, reconstructs more accurately on the mu-distance across the benchmark. The
contribution is a method with a different and practical profile. It runs in seconds where Polyphest
runs for days. It completes every benchmark network where Polyphest completes between a quarter and
all of them. It beats the fair prior-free baseline, iterative GRAMPA, on the mu-distance in every
configuration. And it recovers the reticulate lineages better than Polyphest exactly in the regimes
where copy number is corrupted, through a decode that can be run in a ploidy-free mode that does not
depend on the copy count at all. The chapter presents the method, the benchmark comparison across
fifteen simulated configurations, a diagnostic that locates the method's error in the species tree
backbone rather than in its event prediction, and an analysis of its behaviour under fractionation.
