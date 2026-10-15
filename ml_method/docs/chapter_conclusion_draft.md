NOTE: superseded and no longer built. Its summary now opens chapter_limitations_draft.md as the
Discussion, and its closing paragraph ends chapter_future_work_draft.md. Edits made here will not
appear in the document.

# Conclusion

This chapter introduced PlaceNet, a graph neural network reconstructing polyploid networks by
detecting polyploidization events on a single-copy species tree and then placing them, and compared
it against Polyphest and GRAMPAIter on fifteen simulated configurations.

PlaceNet identifies the polyploid lineages more accurately than Polyphest in all but the mildest
conditions. Its ploidy-free decode leads on the reticulation descendants measure in eleven of the
fifteen configurations at a single fixed threshold, and against GRAMPAIter, the other method
requiring no supplied ploidy, it leads in at least thirteen of the fifteen on each of the three
measures. It reconstructed every network of the discordance benchmark in approximately eight seconds
beyond the ASTRAL step, where Polyphest can take days or fail to return at all. Its two decode modes
allow it to use a trusted ploidy list or to work without one, and at retention 0.50, where copy
counting has begun to fail, the ploidy-free mode is the most accurate of the five methods on both
event measures.

Polyphest remains more accurate on the reticulation sister measure outside the two more severe
fractionation levels, and on the μ-distance throughout. The diagnostic analysis traces both to the
ASTRAL backbone that PlaceNet never rebuilds, which contributes approximately five times the error
of its event prediction. Loosening that dependence is therefore the clearest next step, whether by
phasing the subgenomes before the backbone is built or by reconstructing on several candidate
backbones and selecting among the results.

PlaceNet is therefore best read as a fast and robust complement to the slower search methods rather
than as a replacement for them.
