# Reading-time corpus: source and input audit

This ledger reconstructs the collection directly from the 44 study folders,
including archived tables and original-RT columns alongside processed outcomes.
It does not rely on an unavailable task-agent list.

## Current verified extraction

- 80 primary **analytical inputs**, from 38 publications.
- 10 previously processed analytical inputs are listed separately, **not** included in the primary RDS/model run.
- 112 selected dataset-by-region model jobs; 8,736 nominal pipelines (13 policies × 6 scales).
- Compact primary RDS inputs: 1.39 MB total.

## Completed manuscript run

The extraction totals above describe the available corpus, not the fitted
sample. The final paper uses **79 primary inputs from 37 publications** and
**111 dataset–region jobs**. One primary input (`gastmann2026_spr`) was
excluded for exceptional runtime, as recorded in
`../analysis/excluded_datasets.json`. All 111 included jobs completed;
the ten prepared inputs remain inventoried but are not included in the
paper's numerical analyses.

An analytical input is a modelable design/paradigm, **not necessarily an independent
participant cohort**. Different paradigms or item subsets from one publication may
share participants. Do not describe the input count as that many independent
experiments, papers, or samples. The initial bibliography's approximately 80
entries and the previously mentioned 49 datasets are not target counts.

## Files

- `eligible_analysis_inputs.csv`: exact source, source hash, status, design and
  sample/point counts, plus the article-reported RT outlier and transformation
  policies and the exact multiverse policies; primary rows point to reusable RDS files.
- `paper_rt_policies.csv`: dataset-by-dataset source-paper policy audit, including
  task scope, evidence, source document, and explicit notes where the article
  does not report a rule. A source-paper method is not assumed to have been applied
  to the current multiverse input.
- `folder_coverage.csv`: a decision for **every** local study folder.
- `region_selection.csv`: selected regions and the source-based selection evidence.
- Primary RDS location: `generated-artifacts/brm_multiverse/primary_inputs/`,
  relative to the workspace root.

Each RDS contains `data` and `metadata`. Data retain participant/material grouping,
factor codes, region and positive original-scale RTs. There is no policy-specific
outlier filtering or transformation in these inputs. Invalid/missing observations
and native task eligibility are handled before the common input is written.
The RDS round trip and each region's data/design fingerprint are checked.

The published-paper policy columns describe what each source study reports;
the multiverse columns separately enumerate the screens and scales used here.
The 13 screens are applied within condition × region, with all eligible
observations (including the focal participant) contributing to each reference.
Secondary author-prepared inputs are listed but are not fitted in the primary run.

## Region and processing safeguards

Only critical sites are selected by default. Postcritical sites require evidence
of an experimental main effect in the original report; an interaction alone,
examining a region, or a nuisance age effect is insufficient. Some source papers
define several distinct critical sites or later spillover positions. These are
kept separately, not silently averaged.

Important corrections: ellipsis uses the combined *did too* region rather than
*too* alone; Tucker replications use the authors' participant-cohort definitions
and retain subject gender/number as experimental factors; Kim critical positions
are 7/6/8 across its experiments; Karimi word positions come from the actual word
label rather than its sentence-length field. Changed inputs cannot reuse the
old model cache.

Model results are checkpointed after each policy, with input-specific fingerprints.
A partial checkpoint is never reported as a completed region. Completed model
counts and failed pipelines are read from the run outputs, not inferred from
this extraction ledger. Release verification checks source hashes, compact-input
hashes, region fingerprints, model-output arithmetic, summary totals and figure
points. The publication class and bibliography style are copied without changes.
