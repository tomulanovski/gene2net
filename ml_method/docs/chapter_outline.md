# Method chapter — outline and writing checklist

This maps the current, corrected section drafts into one coherent chapter, in order,
and flags what each section still needs. The chapter is added to the existing
benchmarking thesis, so it does not need its own field-wide intro or related work.
The spine is: a fast learned method, and where it is best.

Framing to hold throughout: the contribution is a fast, self-contained, learned
detect-then-place method that is competitive on accuracy, orders of magnitude faster,
always completes, and robust where copy-number methods break. Not accuracy supremacy.

## Section order

1. **Introduction (short, ~1 page). [WRITE NEW]**
   - The problem is already set up by the thesis. Here: existing methods are slow
     (Polyphest runs for days and fails on hard networks) and copy-number-dependent
     (they break under fractionation). We introduce a fast learned alternative.
   - State the four contributions: learned method, runtime, completion robustness,
     fractionation robustness with a ploidy-independent decode.

2. **Method.** Source: `chapter_method_draft` (archived — salvage the architecture prose)
   + `chapter_setup_draft` (data/metrics/baselines).
   - Architecture: ASTRAL backbone, GAT trunk, detection head + one-partner partner head.
   - Ploidy: the method's own kernel-smoothed copy-number inference. [UPDATE: describe
     kernel smoothing, not the old majority consensus.]
   - Decode: two modes on a ploidy-dependence axis — ploidy-informed (`bound_driven`,
     copy number as a ceiling) and ploidy-free (`detect_only`, no ploidy at all), selected
     by whether a trusted copy-number list exists. [DRAFTED: `chapter_decode_draft.md`.
     Still needs the calibrated ploidy-free threshold from the val sweep.]
   - Training: away-parent labels, the one-parter head. [confirm HPO-final config.]

3. **Experimental setup.** Source: `chapter_setup_draft` (CURRENT).
   - Six discordance/dup-loss configs plus three fractionation configs, 21 benchmark
     networks, five replicates. Metrics: mu-distance, ret_leaf, ret_sisters, plus
     completion and runtime. Baselines: Polyphest, GRAMPA-iter, GRAMPA-iter+prior.

4. **Results.**
   - 4a **Benchmark comparison.** Source: `chapter_results_benchmark_draft` (fixed).
     Beats the GRAMPA family everywhere; overtakes Polyphest at high dup/loss; loses to
     Polyphest on easy conditions (copy-number-dominated). [FILL: finalized 5-replicate,
     kernel-ploidy numbers from the sweep, with error bars.]
   - 4b **Runtime.** Source: `chapter_results_benchmark_draft` runtime section.
     Bounded seconds vs Polyphest days-to-timeout-to-OOM; GRAMPA-iter hours. [FILL: the
     sacct wall-time table you pulled + `--profile` seconds/network.] This is the headline.
   - 4c **Completion robustness.** All 21 vs competitors' 15-18; competitors fail on the
     harder networks. Woven into 4a or its own short subsection.
   - 4d **Diagnostic: where the error comes from.** Source: `chapter_diagnostic_draft`
     (CURRENT). Backbone is the wall, event prediction saturated. Separate the benchmark's
     distribution-shift component.
   - 4e **Fractionation.** Source: `chapter_retention_draft` (CURRENT). The three-level
     regime map, the ploidy-independent `detect_only` win at the crossover, the
     prior-backfires finding. [FILL: finalized 5-replicate numbers.]
   - 4f **Feature importance.** Source: `chapter_feature_importance_draft` (CURRENT).

5. **Limitations.** Source: `chapter_limitations_draft` (fixed) + `chapter_partner_limitation_draft` (CURRENT).
   - Backbone is the binding limitation; placement reduces to the same cause.
   - Honest where it loses: Polyphest on easy, GRAMPA-iter at extreme fractionation.

6. **Future work.** [WRITE from the discussion.]
   - Train on fractionated data (extend the fractionation win) — the strongest one.
   - Iterative MUL-tree decoding for nested events (a MUL-tree-trained model), motivated
     by GRAMPA-iter's edge at extreme fractionation.
   - Learned gene-tree encoding vs hand-crafted features (the archived Legacy model).
   - Backbone rebuild via phasing (from the diagnostic).
   - Calibrated / val-tuned decode threshold.

7. **Conclusion (short). [WRITE NEW]** A fast learned method, competitive accuracy at a
   fraction of the runtime, robust where copy-number methods break, with an honest map of
   which method to prefer under which conditions.

## What is DONE vs TODO

- Current and ready: setup, diagnostic, retention, feature-importance, partner-limitation.
- Fixed, ready: results-benchmark (table + runtime), limitations.
- Drafted: decode two-mode subsection (`chapter_decode_draft.md`).
- To write new: introduction, future-work, conclusion.
- To fill with numbers: 4a and 4e finalized 5-replicate kernel-ploidy numbers; the 4b
  runtime table (sacct wall times + --profile).

## Suggested writing order

1. Method + decode spectrum (the core novelty; you understand it cold).
2. Results 4a/4b/4c (benchmark + runtime + completion) — slot finalized numbers as they land.
3. Diagnostic 4d and fractionation 4e are already drafted — light edits only.
4. Limitations + future work (mostly assembly).
5. Introduction and conclusion last, once the body is fixed.
