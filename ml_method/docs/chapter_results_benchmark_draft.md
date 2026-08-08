# Results: detection and comparison to baselines

DRAFT for the thesis. Prose style follows the thesis convention of no semicolons, no
non-mathematical parentheses, and no em-dashes. Edit distance is the canonical multi-labeled-tree
edit distance, so it does not depend on Newick child order. Lower is better for every distance.
The method here is the final one, namely the one-partner model with the away-parent labels and
the fill-to-inferred-ploidy build. Items marked TODO are to confirm before submission. GNN scores
are over all 21 networks. The competitor columns are over the 15 to 18 networks each competitor
completes, so the competitor means are taken on the subset each can finish.

## Detection

The whole genome duplication detector reaches an F1 of about 0.80 on the held-out validation
split. [TODO confirm precision and recall from the training log.] Detection is deliberately
recall-oriented, because the fill-to-inferred-ploidy build reaches each species' inferred copy
number from the ranked candidate events rather than from a hard probability threshold, so a
moderately confident true event is still used.

## Reconstruction benchmark

We compare against three baselines on the 21 benchmark networks across all six configurations.
All four methods infer ploidy from the gene trees, so the comparison is between prior-free
methods on equal inputs. The table reports mean edit distance.

| Configuration | GNN | Polyphest | GRAMPA | GRAMPA-iter |
| --- | --- | --- | --- | --- |
| ils_low | 0.573 | 0.065 | 0.682 | 0.681 |
| ils_medium | 0.571 | 0.112 | 0.695 | 0.720 |
| ils_high | 0.554 | 0.209 | 0.709 | 0.825 |
| dup_loss_low | 0.570 | 0.196 | 0.700 | 0.729 |
| dup_loss_medium | 0.564 | 0.177 | 0.707 | 0.756 |
| dup_loss_high | 0.603 | 0.636 | 0.716 | 0.812 |

Four findings.

First, the method beats both GRAMPA and iterative GRAMPA in every configuration, by roughly 0.1
to 0.25 in edit distance. Among the prior-free methods that reconstruct a full network from gene
trees, the method is clearly better than the GRAMPA family across all conditions.

Second, Polyphest is the strongest method at low and moderate difficulty, and by a large margin.
This traces to copy number. The edit distance on multi-labeled trees is dominated by getting each
species' copy number right, and Polyphest builds its start tree from the inferred consensus
multiset, which on clean simulations is essentially the true ploidy. So on the easy conditions
Polyphest is effectively handed the copy-number structure and only has to place events, which is
why it reaches an edit distance as low as 0.065.

Third, that advantage erodes as duplication and loss rise and reverses at the hardest condition.
Under high duplication and loss the copy-number signal is corrupted, the inferred multiset stops
matching the truth, and Polyphest's start tree is wrong. There the method reaches 0.603 against
Polyphest's 0.636, so it is the better method, and the reticulation-leaf distance agrees at that
condition, with the method at 0.264 against Polyphest's 0.412. The method reaches each species'
inferred ploidy from its ranked events rather than from a fixed start tree, so it degrades more
gracefully when ploidy inference is hard.

Fourth, the comparison understates the method's robustness. Polyphest is scored on 15 to 18 of
the networks per configuration because it does not complete the rest, while the method completes
all 21. So Polyphest's means are taken over the subset it can finish, and the method always
produces an answer.

On the reticulation distances the ordering matches the edit distance. Polyphest leads on the
networks it completes at low and moderate difficulty, the method beats iterative GRAMPA
everywhere and GRAMPA on the reticulation-leaf distance in most configurations, and the method
overtakes Polyphest at high duplication and loss.

## Runtime

The method reconstructs a network in a few seconds. Its per-network compute, excluding the ASTRAL
step, is about 8 seconds on the benchmark, of which the graph neural network forward pass is a
fraction of a second and the rest is feature extraction from the gene trees. ASTRAL adds a few
seconds. Polyphest, by contrast, ranges from seconds on easy networks to days on hard ones, and
does not complete the networks it is missing from the table above. [TODO insert a concrete
Polyphest wall time or range from the run logs.] So the method offers a bounded, predictable
runtime of seconds per network, against a search-based competitor whose runtime is variable and
can be prohibitive.

## Summary

The method is a fast, prior-free reconstruction that beats the GRAMPA family across all
conditions and runs in seconds where Polyphest can take days and sometimes does not finish.
Polyphest is more accurate at low and moderate difficulty, where inferring ploidy from clean
simulations is easy and the copy-number-dominated edit distance rewards it, but the method
overtakes Polyphest at high duplication and loss, where ploidy inference breaks. The remaining
gap at the easier conditions is copy number and placement, which the diagnostic traces to event
prediction rather than to the backbone.
