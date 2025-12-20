# Meeting Preparation Guide

## Quick Summary of What You Have

You've run a comprehensive analysis comparing **6 phylogenetic network reconstruction methods** across **21 real networks** under **3 ILS conditions** (low, medium, high).

### Methods Compared:
- **GRAMPA**
- **Polyphest** (3 variants: p50, p70, p90)
- **MPSUGAR**
- **PADRE**

### Key Metrics Measured (Lower is better for all):
1. **edit_distance** - Overall network similarity (0=perfect, 1=terrible)
2. **num_rets_diff** - How many reticulations off from truth
3. **ploidy_diff.dist** - Polyploid detection accuracy
4. **ret_leaf_jaccard.dist** - Did they find the right species in reticulations?
5. **ret_sisters_jaccard.dist** - Did they get the topology right?

---

## Understanding Your Results

### Original Summary Files (simulations/analysis/summary/)
These are comprehensive but overwhelming:
- **aggregated_metrics.csv** - 2,000+ rows of detailed statistics
- **level1_detailed_per_network/** - 15 CSV files with per-network breakdowns
- **level2_method_rankings.csv** - Method rankings across networks
- **level3_network_correlations.csv** - What makes networks harder?

### Simplified Meeting Summaries (simulations/analysis/meeting_summaries/)
**START HERE!** Much clearer for presentation:
- **00_KEY_FINDINGS.txt** - One-page summary of everything important
- **02_best_methods.csv** - Which method wins for each metric
- **01_method_rankings.csv** - Detailed rankings
- **03_per_network_edit_distance.csv** - Which networks are easy/hard

---

## Key Findings to Present (ILS Low)

### 1. Overall Winner
**Polyphest (p70 variant) is the best overall method**
- Average rank: 2.00 (best across all metrics)
- Wins in 3 out of 5 metrics
- Very consistent performance

### 2. Method Performance Ranking (Best to Worst)
1. **polyphest_p70** (avg rank: 2.00) - BEST
2. **polyphest_p90** (avg rank: 2.40)
3. **polyphest_p50** (avg rank: 2.60)
4. **padre** (avg rank: 3.00)
5. **grampa** (avg rank: 5.20)
6. **mpsugar** (avg rank: 5.80) - WORST

### 3. Best Method per Task
- **Overall structure** (edit_distance): PADRE (0.619)
- **Counting reticulations** (num_rets_diff): Polyphest p70 (perfect!)
- **Finding polyploids** (ploidy_diff): Polyphest p50 (perfect!)
- **Identifying reticulated species**: Polyphest p70 (0.0015)
- **Getting topology right**: Polyphest p70 (0.0234)

### 4. Easiest vs Hardest Networks
**Easiest** (all methods work well):
- Wisecaver_2023 (0.387)
- Ding_2023 (0.443)
- Popp_2005 (0.590)

**Hardest** (all methods struggle):
- Soza_2014 (0.969)
- Diaz-Perez_2018 (0.966)
- Sessa_2012b (0.943)

### 5. What Makes Networks Hard?
**GRAMPA and MPSUGAR struggle with:**
- More reticulations (r=0.999, p<0.0001)
- More polyploids (r=0.698, p<0.001)
- More species (r=0.626, p<0.01)

**Polyphest is more robust** - weaker correlations with network properties

---

## What About Level 3 (Network Correlations)?

### What It Shows:
Level 3 answers: **"What network properties make reconstruction harder?"**

For example:
- `grampa | num_rets_diff | num_reticulations -> r=0.999`
  - Translation: "GRAMPA gets worse at counting reticulations as networks have more reticulations"
  - This is a STRONG correlation (r=0.999) and highly significant (p<0.0001)

- `polyphest_p70 | edit_distance | num_species -> r=0.598`
  - Translation: "Polyphest p70 has moderate difficulty with larger networks"
  - This is a MODERATE correlation (r=0.598) and significant (p=0.019)

### Key Insight:
**GRAMPA and MPSUGAR** are sensitive to network complexity (high correlations)
**Polyphest variants** are more robust (lower correlations)

---

## What About Secondary Metrics (FP, FN, TP)?

These are just **technical breakdowns** of the primary metrics:
- **FP (False Positive)** - What they found but shouldn't have
- **FN (False Negative)** - What they missed
- **TP (True Positive)** - What they correctly found

**For your meeting: Focus on the primary metrics (.dist), not the FP/FN/TP breakdowns.**

The `.dist` metrics already incorporate all this information in an easy-to-understand format.

---

## Is the Current Summary Approach the Best?

### What's Good:
1. **Comprehensive** - Covers all important comparisons
2. **Multi-level** - Can drill down from high-level to details
3. **Statistically sound** - Mean ± std across 5 replicates
4. **Well-documented** - Clear metrics and interpretations

### What Could Be Improved:
1. **Too many files** - Original output is overwhelming
2. **Secondary metrics clutter** - FP/FN/TP add confusion
3. **No cross-configuration comparison** - Have to look at 3 separate summaries

### What I Created for You:
**Simplified meeting summaries** that focus on:
- One-page key findings
- Clear "winner" identification
- Easy-to-read tables
- Only the metrics that matter

---

## Recommended Presentation Flow

### Slide 1: Study Overview
- 6 methods tested
- 21 real networks
- 3 ILS conditions
- 5 metrics evaluated

### Slide 2: Overall Winner (ILS Low)
```
Method Performance Ranking:
1. Polyphest p70 (avg rank: 2.00) ← BEST
2. Polyphest p90 (avg rank: 2.40)
3. Polyphest p50 (avg rank: 2.60)
4. PADRE (avg rank: 3.00)
5. GRAMPA (avg rank: 5.20)
6. MPSUGAR (avg rank: 5.80)
```

### Slide 3: Best Method per Task
Show the 02_best_methods.csv as a table

### Slide 4: Network Difficulty
- Show easiest vs hardest networks
- Mention what makes networks hard (correlations)

### Slide 5: Effect of ILS Level
Compare results across low/medium/high ILS
(Look at all three 00_KEY_FINDINGS.txt files)

### Slide 6: Key Insights
- Polyphest is most reliable
- GRAMPA/MPSUGAR struggle with complexity
- Some networks are universally hard

---

## Next Steps Before Your Meeting

1. **Read all three 00_KEY_FINDINGS.txt files**
   - conf_ils_low_10M
   - conf_ils_medium_10M
   - conf_ils_high_10M

2. **Note how rankings change with ILS level**
   - Does Polyphest stay best?
   - Does GRAMPA get worse with more ILS?

3. **Prepare to discuss**
   - Why is Polyphest better?
   - What do the failed comparisons mean? (check 04_comparison_report.txt)
   - Which method would you recommend for real data?

4. **Potential supervisor questions**
   - "Why did some methods fail on certain networks?"
   - "How does ILS affect method performance?"
   - "Is the sample size (5 replicates) sufficient?"
   - "Should we test on more networks?"

---

## Files to Review (Priority Order)

### MUST READ (10 minutes):
1. `meeting_summaries/conf_ils_low_10M/00_KEY_FINDINGS.txt`
2. `meeting_summaries/conf_ils_medium_10M/00_KEY_FINDINGS.txt`
3. `meeting_summaries/conf_ils_high_10M/00_KEY_FINDINGS.txt`

### GOOD TO KNOW (20 minutes):
4. `meeting_summaries/conf_ils_low_10M/02_best_methods.csv`
5. `meeting_summaries/conf_ils_low_10M/04_comparison_report.txt` (check failure rate)
6. Original `SUMMARY_PIPELINE_GUIDE.md` (understanding methodology)

### IF TIME PERMITS:
7. `meeting_summaries/*/03_per_network_edit_distance.csv` (which networks are hard)
8. `summary/conf_ils_low_10M/level3_network_correlations.csv` (correlation details)

---

## Common Questions and Answers

### Q: Why are there so many Polyphest variants?
A: They differ in the percentile threshold (p50, p70, p90) used to filter branches.
   p70 seems to be the sweet spot for your data.

### Q: Why did some comparisons fail?
A: Check 04_comparison_report.txt. Common reasons:
   - Polyphest returned empty networks for some cases (division by zero)
   - Some methods didn't complete for all replicates

### Q: Should I use all these metrics?
A: Focus on **edit_distance** for overall performance and **num_rets_diff** for
   reticulation accuracy. The others provide additional insights.

### Q: What's the difference between Level 2 and Level 3?
A: - Level 2: "Which method is best overall?"
   - Level 3: "What network properties make methods fail?"

### Q: Is the summary approach good enough?
A: Yes! It's comprehensive and well-structured. The simplified meeting summaries
   make it easier to present. Consider creating cross-configuration comparisons
   for future work.

---

## Suggested Improvements for Future

1. **Create a comparison across ILS levels**
   - Single table showing how each method's rank changes with ILS

2. **Add visualization scripts**
   - Heatmaps of method performance per network
   - Bar charts of average rankings

3. **Statistical testing**
   - Wilcoxon signed-rank test to compare methods
   - Bootstrap confidence intervals

4. **Simplify output files**
   - Skip FP/FN/TP secondary metrics
   - Create a single "executive summary" PDF

---

## Final Advice

**For your meeting:**
1. Lead with the KEY_FINDINGS.txt - it tells the whole story
2. Have the detailed CSVs ready for follow-up questions
3. Be prepared to explain what each metric means
4. Discuss implications: "Which method should we use for real data?"

**Remember:**
- Your analysis is solid and comprehensive
- The simplified summaries make it presentable
- Focus on the story: "Polyphest works best, here's why"

Good luck with your meeting!
