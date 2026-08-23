# Introduction

Up to this point we benchmarked the existing methods for polyploid phylogenetic network inference.
Building on that comparison, this chapter develops and evaluates a new method of our own, which we
call PlaceNet.

The benchmark exposed two practical limitations that a new method could try to address. The first is
runtime. Polyphest for example solves an expensive integer linear program and can run for days on a
single network, and on the hardest networks it does not finish at all, reaching a time limit or
exhausting memory. The second concerns the copy number in the gene trees. Some methods lean on a
supplied per species copy number and build their reconstruction around it. Copy number is a reliable
signal on clean data, but it is degraded by fractionation, the gradual loss of one duplicate copy
after a WGD, and real polyploids fractionate extensively. A method whose reconstruction does not
depend on the copy number could therefore be more robust in the regime that real data occupies.

PlaceNet is a graph neural network (GNN) that follows a detect-then-place strategy. It builds a
single-copy species tree backbone with ASTRAL-4, reads that backbone together with features
summarising the gene trees, and predicts on every backbone edge whether a WGD occurred there and, if
so, which lineage is the second parent. The predicted events are grafted onto the backbone to
produce a MUL tree, which is the output of the method.

On the overall distance measure, the mu-distance, Polyphest reconstructs more accurately across the
benchmark. What PlaceNet offers is a different and practical profile. It runs in seconds for every
input where Polyphest can run for days. It completes every benchmark network, whereas the completion
rate of the other methods varies, and in our implementation only GRAMPA-Iter also reached full
completion. It recovers the reticulate lineages more faithfully than any other method where the copy
number is corrupted, through a decode that can be run in a mode that does not use the copy count at
all, which is the regime that fractionation creates in real data. The chapter presents the method,
its comparison on the same benchmark used in the previous chapters, a diagnostic that locates the
method's error in the species tree backbone rather than in its event prediction, and its behaviour
under fractionation.
