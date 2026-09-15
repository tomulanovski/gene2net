# Conclusion

This chapter introduced PlaceNet, a learned detect-then-place method for polyploid network
reconstruction, and compared it with Polyphest and iterative GRAMPA on fifteen simulated
configurations.

PlaceNet identifies the polyploid lineages more accurately than Polyphest in all but the mildest
conditions. Its ploidy-free decode leads on the reticulation descendants measure in eleven of the
fifteen configurations at a single fixed threshold, and against iterative GRAMPA, the other method
that needs no supplied ploidy, it leads on every measure in at least thirteen of the fifteen. It
reconstructs every discordance-benchmark network in seconds, where Polyphest can take days or not
finish. Its two decode modes let it use a trusted ploidy list or work without one, and the
ploidy-free mode is the most accurate method at medium fractionation, where copy counting fails.

Polyphest remains ahead on the reticulation sister measure outside the two more severe fractionation
levels, and on the mu-distance everywhere. The diagnostic traces both to the ASTRAL backbone that
PlaceNet never rebuilds, which contributes five times the error of its event prediction, so
rebuilding the backbone, for example by phasing, is the clearest next step.

PlaceNet is best read as a fast and robust complement to the slower search methods rather than a
replacement for them.
