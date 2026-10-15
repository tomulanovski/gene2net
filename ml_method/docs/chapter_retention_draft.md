## Reconstruction Under Fractionation

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes.

Fractionation was simulated at retention rates of 0.75, 0.50 and 0.25 on top of the medium
duplication and loss, medium ILS configuration, on the 20 benchmark networks carrying a fractionated
condition. Because fractionation deletes gene copies, any method reading ploidy from the copy
numbers in the gene trees comes to see a fractionated polyploid as a diploid.

@fig:fractionation shows the three measures as retention falls from 1.00, the unfractionated
condition, to 0.25, and @tab:frac075, @tab:frac050 and @tab:frac025 in the appendix give the exact
values. Every method is again restricted to the networks all of them completed, 18 at retention 0.75
and 0.50 and 16 at 0.25.

Figure: {#fig:fractionation} figures/fractionation_degradation.png | Reconstruction accuracy as retention falls, on the medium duplication and loss, medium ILS configuration, where retention 1.00 is the unfractionated condition. Each point is the mean on the networks that every method completed, with error bars showing the standard error across those networks. Dashed lines with open markers are the variants using ploidy information. Lower is better.

Alt text: Three line plots side by side, for the reticulation descendants distance, the reticulation sister distance and the μ-distance, each plotting the five methods against four retention levels from 1.00 down to 0.25. All lines rise as retention falls. In the left panel the two PlaceNet lines cross the Polyphest line between retention 0.75 and 0.50.

The two decode modes trade places on the reticulation descendants measure as retention falls. At
0.75 the copy number remains informative and the ploidy-informed mode is the most accurate of the
five methods, at 0.160 against 0.244 for Polyphest. At 0.50 the ploidy-free mode takes over at 0.506,
the lowest of any method. At 0.25 GRAMPAIter recovers the most at 0.602 and the ploidy-free mode is
second at 0.876, still ahead of Polyphest at 0.917. On the reticulation sister measure the
ploidy-free mode is ahead of Polyphest at retention 0.50 and 0.25, the only two configurations of
the fifteen where it is.

Supplying GRAMPAIter with the inferred ploidy list helps only at retention 0.75. With the prior its
reticulation descendants distance is 0.249, 0.672 and 0.923 as retention falls, against 0.492, 0.542
and 0.602 without it, because the prior is the same collapsed copy number that fractionation
corrupts. On the μ-distance Polyphest is the most accurate method at every retention level. PlaceNet
reconstructed 20 of the 20 fractionated networks at retention 0.75 and 0.50 and 19 at 0.25, the one
exception being a network whose ASTRAL backbone could not be parsed, against 19, 19 and 18 for
Polyphest.
