# Alignment with `a plan.txt`

The manuscript follows the plan's six-part structure, with the substantive
results organized around inferential, effect-size and precision stability.
The original plan is not rewritten to conceal implementation differences.

## Coverage

1. **Introduction:** preprocessing flexibility, its relevance to reading
   tasks, and the distinction between empirical sensitivity and known-truth
   simulation evidence.
2. **Preprocessing:** all 13 requested exclusion rules and all six nominal
   scales; explicit grouping, boundary conventions, scaled MAD and
   median-relative bounds. Natural and base-10 logarithms are the same
   inferential specification, leaving 65 distinct workflows.
3. **Dataset inventory:** citations/hosts, task, factorial design,
   participants, eligible observations, and paper-derived original
   exclusion/transformation policies. Shared participants and publications
   are not treated as independent experiment replications.
4. **Analysis:** supported maximal participant/item mixed models,
   explicit sum contrasts, correlation removal followed by variance-guided
   component removal until convergence and nonsingularity. Estimates,
   SEs, t statistics, degrees of freedom, p values, confidence intervals,
   retained observations and random-effect variances are exported.
5. **Results:** significance crossings, direction changes, standardized
   effect ranges, SE and interval-width ratios; main-effect/interaction
   summaries; task, region and design comparisons; associations with
   skewness, participant/item counts and reference-effect magnitude.
   The revised results additionally isolate exclusion from transformation
   by comparing matched workflows differing in one decision dimension.
6. **Discussion:** specific interpretation of the most consequential
   decisions, relevant methodological literature, scope/limitations,
   and practical recommendations as the last discussion subsection.

## Explicit departures and qualifications

- **Coverage is 79 primary analytical inputs, not 80 completed independent
  studies.** One unusually expensive input was excluded for runtime; ten
  previously prepared inputs are inventoried but not fitted as primary data.
  There are 111 source-selected region jobs. This restriction is disclosed.
- **Degrees of freedom are Satterthwaite approximations** from lmerTest,
  rather than an unspecified t-reference approximation.
- **Significant-interaction follow-ups are conditional contrasts of the
  accepted model, not newly fitted nested-contrast models.** They retain
  its covariance and degrees-of-freedom calculation. This differs literally
  from the plan's request for an additional fitted nested-contrast model.
  Follow-ups are separately exported, not counted as unconditional primary
  multiverse tests.
- **Stability uses a declared common reference, not every paper's original
  significance classification.** Original papers used heterogeneous
  eligibility rules and models. The reference and its limits are explicit;
  the new matched comparisons also avoid dependence on one anchor.
- **Random-effect estimates are fully exported**, but are not pooled into
  one common-scale population-variance effect. Covariance structures and
  outcome scales change across workflows; treating those estimates as one
  invariant variance would be misleading.
- **Dataset-level effect ranges remain reproducible in numerical exports.**
  The main paper replaces a crowded, selected-one-coefficient-per-input
  plot with decision-level sensitivity. All coefficient–region ranges
  remain summarized, without selecting just one effect from each input.
- The paper does not establish true/false positives or an optimal
  preprocessing procedure. The archive contains observations, not labels
  identifying a known cognitive generating process.
