# Reading-time preprocessing multiverse — OSF release (moved to GitHub after OSF removed Projects)

## Contents

- `data/`: **79 primary analytical inputs**, original-scale eligible RTs,
  factor mappings and input/checkpoint fingerprints; no secondary inputs.
- `checkpoints/`: **111 completed region jobs**. These small saved results
  avoid repeating the expensive model fits.
- `code/analysis/`: model, summary, figure and verification code.
- `code/source_adapters/`: source-specific extraction code and region plan.
- `provenance/`: verified inventory, factor levels and exact contrast coding,
  full-rank fixed designs and supported random structures, condition-cell counts,
  source-family dependency hashes, region evidence, timings and audit reports.

There are 37 source publications and 513,693 eligible region-level observations.
Analytical inputs and regions are not independent participant cohorts.
There are 78 nominal workflows per region but only 65 distinct specifications:
natural and base-10 logarithms give identical inferential results.

## Verify the download

    python3 reproduce.py --check

`MANIFEST.csv` contains SHA-256 and byte size for every release file except
the manifest itself. No network access is needed.

## Note
The original files included a script to generate the paper, but they are removed for the duration of review.

## Runtime warning and a complete fresh fit

The completed jobs recorded **107.1 cumulative
job-hours**; the longest region job recorded
**12.3 hours**.
These are summed worker elapsed times, **not CPU time or wall-clock duration**.
Jobs ran in parallel; recorded resumptions are included, abandoned/excluded
runs and unrecorded setup are not. A fresh fit can take multiple days,
depending on hardware and parallelism. This is an empirical multiverse
reanalysis, not a Monte Carlo simulation.

To deliberately discard the shortcut and fit into a new empty build:

    python3 reproduce.py --refit --workers 1

Use a sensible worker count for your computer. The full-fit runner uses
fork-based parallelism and is intended for macOS/Linux. The default command
uses the supplied checkpoints. Reproducing tables/figures does NOT require
refitting the models.

## Scope of verification

Fresh source-adapter reconstruction was compared observation by observation
with all 79 compact inputs. Participant/item counts, complete factorial
cells, factor mappings, selected regions and fitted-input identities were
checked. `provenance/provenance_verification.json` states the scope and limits.
`provenance/coding.csv` records the actual contrast values by level, not
guessed semantic labels for numeric source codes. Two-level factors use
-0.5/+0.5 in the declared order; larger factors retain K-1 sum contrasts.
`provenance/structures.csv` maps supported random-slope columns to fixed
predictors separately for participants and items in every region.
Source-family hashes deliberately include companion inputs and annotations,
sometimes files used by other cohorts in the same source family.

This verifies the implemented mappings; it is not an independent rereading
of every original article or proof that all modelling assumptions are correct.
The original 90-entry audit also inventories 10 prepared/secondary inputs
and one primary input excluded for runtime. These are not part of the
79-input results or data here. Historical `original_removal` and
`original_transform` fields in extraction metadata are not authoritative:
use the paper-derived policies in `provenance/source_audit/paper_rt_policies.csv`.
In particular, two old Fujita metadata labels conflict with the published
reciprocal transformation; the inventory uses the paper evidence.

## Data provenance, privacy and reuse

Compact inputs retain the source-supplied grouping identifiers needed to
reproduce the models. An automated screen found no email-like or
human-name-like subject identifiers; this heuristic is not a privacy
certification. Participant counts must not be summed across shared cohorts.

Original source archives and publisher PDFs are not duplicated. The
manuscript cites source studies and hosts; source paths and hashes are in
`provenance/`. Repeating the original extraction additionally requires that
source collection and the source-paper text used for region evidence.
The saved audit proves the local reconstruction checks, not that source
files absent from this release were re-read on your computer.

Original data, publisher files and dependencies retain their original rights
and terms. This release does not assert new redistribution permissions or
apply a blanket licence to third-party material. Confirm the relevant source
terms before public upload. Author-specific funding and conflict declarations
in the manuscript also remain for the author to confirm.
