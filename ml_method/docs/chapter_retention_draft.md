## Reconstruction under fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Numbers are on the 20 benchmark networks that carry
a fractionated condition, at three retention levels, five replicates. Lower is better throughout.

We simulate fractionation at retention rates of 0.75, 0.50 and 0.25 on top of the medium duplication
and loss, medium ILS condition. Since it deletes copies, any method that reads ploidy from the copy
counts in the gene trees comes to see a fractionated polyploid as a diploid.

@fig:fractionation shows the three measures as retention falls from 1.00, the unfractionated
condition, to 0.25, and @tab:frac075, @tab:frac050 and @tab:frac025 in the appendix give the exact
values. Every method is restricted to the networks all of them completed, 18 at 0.75 and 0.50
retention and 16 at 0.25.

Figure: {#fig:fractionation} figures/fractionation_degradation.png | Reconstruction accuracy as retention falls, on the medium duplication and loss, medium ILS condition, where retention 1.00 is the unfractionated condition. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants that use ploidy information. Lower is better.

On the reticulation descendants measure the decode modes trade places as retention falls. At 0.75
the copy number is still informative and the ploidy-informed mode is the most accurate method, 0.160
against Polyphest's 0.244. At 0.50 the ploidy-free mode takes over at 0.506, the lowest of any
method. At 0.25 iterative GRAMPA recovers the most at 0.602, and the ploidy-free mode is second at
0.876, still ahead of Polyphest's 0.917. On the reticulation sister measure the ploidy-free mode is
ahead of Polyphest at 0.50 and 0.25.

Giving iterative GRAMPA the ploidy prior helps only at 0.75. With the prior its reticulation
descendants distance is 0.249, 0.672 and 0.923 as retention falls, against 0.492, 0.542 and 0.602
without it, because the prior is the same collapsed copy number that fractionation corrupts. On the
mu-distance Polyphest is the most accurate method at every level. PlaceNet completes 20 of the 20
fractionated networks at 0.75 and 0.50 retention and 19 at 0.25, where one ASTRAL backbone could not
be parsed, against 19, 19 and 18 for Polyphest.
