# Reconstruction under fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. The numbers below are on the 20 benchmark
networks that carry a fractionated condition, under three fractionation configurations, five
replicates. The headline metric is the
normalized mu-distance, and the reticulation-leaf and reticulation-sister Jaccard distances
measure whether the reticulations involve the right lineages. Lower is better throughout.

After a whole genome duplication, one copy of a region is often lost from some descendant lineages,
a process called fractionation. We simulate it with a retention rate, the probability that a
duplicated branch is kept, at 0.75, 0.50 and 0.25 on top of the medium duplication and loss, medium
ILS condition. Fractionation deletes copies, so a method that infers ploidy from the copy counts in
the gene trees sees fewer copies, concludes that polyploids are diploid, and is expected to degrade
as retention falls.

PlaceNet's detection head reads more than copy counts, including branch-length asymmetry and gene
tree clustering structure, which can survive after a copy is lost. We tested whether it still fires
on the events that copy counting misses. On the medium fractionation condition, polyploids whose
copies were deleted carry a mean detection probability of 0.55, and 59 percent of them exceed 0.5,
against 0.09 and 6 percent for true diploids. Even at severe fractionation they sit at 0.30, roughly
nine times the diploid level. The ploidy-informed decode cannot use this signal, because it fills
each species only up to an inferred copy number that fractionation has reduced to one. The
ploidy-free decode keeps every confident edge and can recover these events.

@fig:fractionation shows the three measures as retention falls from 1.00, the unfractionated
condition, to 0.25, and @tab:frac075, @tab:frac050 and @tab:frac025 in the appendix give the exact
values. Every method is restricted to the networks all of them completed, 18 at mild and medium
retention and 16 at severe.

Figure: {#fig:fractionation} figures/fractionation_degradation.png | Reconstruction accuracy as retention falls, on the medium duplication and loss, medium ILS condition, where retention 1.00 is the unfractionated condition. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants that use ploidy information. Lower is better.

On the reticulation descendants measure the decode modes trade places as retention falls. At mild
retention the copy number is still informative and the ploidy-informed mode is the most accurate
method, 0.160 against Polyphest's 0.244. At medium retention the ploidy-free mode takes over at
0.506, the lowest of any method. At severe retention iterative GRAMPA recovers the most at 0.602,
and the ploidy-free mode is second at 0.876, still ahead of Polyphest's 0.917. On the reticulation
sister measure the ploidy-free mode is ahead of Polyphest at medium and severe retention.

Giving iterative GRAMPA the ploidy prior helps only at mild fractionation. With the prior its
reticulation descendants distance is 0.249, 0.672 and 0.923 as retention falls, against 0.492, 0.542
and 0.602 without it, because the prior is the same collapsed copy number that fractionation
corrupts. On the mu-distance Polyphest is the most accurate method at every level. PlaceNet
completes 20 of the 20 fractionated networks at mild and medium retention and 19 at severe
retention, where one ASTRAL backbone could not be parsed, against 19, 19 and 18 for Polyphest.
