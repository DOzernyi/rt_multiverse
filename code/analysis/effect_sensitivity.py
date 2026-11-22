"""Summarize all 78 workflows against no exclusion / untransformed RTs.

Run from the package or workspace root. No model is fitted or changed.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path("generated-artifacts/brm_multiverse")
OUT = ROOT / "effect_sensitivity"
POLICIES = [
    "none", "sd1.5", "sd2", "sd2.5", "sd3", "fixed100_2000",
    "fixed100_3000", "fixed200_2000", "fixed200_3000", "median",
    "mad2", "mad2.5", "mad3",
]
SCALES = ["identity", "ln", "log10", "inverse", "sqrt", "boxcox"]
KEYS = ["dataset", "region", "term"]
JOBS = ["dataset", "region"]
WORKFLOW = ["policy", "transformation"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(ROOT / "coefficients.csv.gz")
    c = c[c.term != "(Intercept)"].copy()
    assert len(c) == 380 * 78
    assert not c.duplicated(KEYS + WORKFLOW).any()
    assert c[["SE", "standardized_SE"]].gt(0).all().all()
    # Native SE / standardized SE recovers the exact pipeline outcome SD.
    c["standardized_CI_width"] = (
        (c.CI_high - c.CI_low) * c.standardized_SE / c.SE
    )
    ref = c[(c.policy == "none") & (c.transformation == "identity")]
    ref = ref[KEYS + ["standardized_estimate", "standardized_SE",
                      "standardized_CI_width"]]
    c = c.merge(ref, on=KEYS, suffixes=("", "_reference"), validate="many_to_one")
    c["effect_shift"] = abs(c.standardized_estimate - c.standardized_estimate_reference)
    c["SE_change"] = 100 * (c.standardized_SE / c.standardized_SE_reference - 1)
    c["CI_width_change"] = 100 * (
        c.standardized_CI_width / c.standardized_CI_width_reference - 1
    )
    c["SE_signed_change"] = c.SE_change
    c["CI_width_signed_change"] = c.CI_width_change
    c["SE_change"] = abs(c.SE_change)
    c["CI_width_change"] = abs(c.CI_width_change)
    c["interaction_shift"] = abs(
        abs(c.standardized_estimate) - abs(c.standardized_estimate_reference)
    )

    # Independent-slope models get arbitrary group suffixes in VarCorr.
    # Match by original grouping factor and named design term(s), never row order.
    v = pd.read_csv(ROOT / "random_variances.csv.gz")
    v = v[v.grp != "Residual"].copy()
    v["group"] = v.grp.str.replace(r"\.\d+$", "", regex=True)
    assert set(v.group) == {"subject", "item"}
    v["component"] = [
        g + "|" + "|".join(sorted([a, a if pd.isna(b) else b]))
        for g, a, b in zip(v.group, v.var1, v.var2)
    ]
    assert not v.duplicated(JOBS + WORKFLOW + ["component"]).any()
    # log10 is an exact positive rescaling of ln. VarCorr is already in
    # pipeline outcome-variance units, so its standardized entries are identical.
    assert "log10" not in set(v.transformation)
    v = pd.concat([v, v[v.transformation == "ln"].assign(transformation="log10")],
                  ignore_index=True)
    covariance = {}
    for key, z in v.groupby(JOBS + WORKFLOW):
        covariance[key] = z.set_index("component").vcov.to_dict()
    job_sources = c[JOBS + ["citation"]].drop_duplicates().set_index(JOBS).citation
    random_rows = []
    support_rows = []
    for (dataset, region, policy, scale), entries in covariance.items():
        reference = covariance[(dataset, region, "none", "identity")]
        common = set(entries) & set(reference)
        union = set(entries) | set(reference)
        assert common and union
        citation = job_sources.loc[(dataset, region)]
        common_shift = np.mean([abs(entries[k] - reference[k]) for k in common])
        row = dict(dataset=dataset, region=region, citation=citation,
                   policy=policy, transformation=scale,
                   random_covariance_shift=common_shift,
                   random_structure_change=100 * (1 - len(common) / len(union)))
        random_rows.append(row)
        support_rows.append(dict(
            **{k: row[k] for k in JOBS + ["citation"] + WORKFLOW},
            shared_entries=len(common), reference_entries=len(reference),
            target_entries=len(entries), added_entries=len(set(entries)-set(reference)),
            removed_entries=len(set(reference)-set(entries)),
        ))
    random_values = pd.DataFrame(random_rows)
    support = pd.DataFrame(support_rows)
    assert len(random_values) == 111 * 78
    support.to_csv(OUT / "random_component_support.csv", index=False)

    sources = []
    counts = []
    for metric in ["effect_shift", "SE_change", "CI_width_change",
                   "random_covariance_shift", "interaction_shift",
                   "random_structure_change", "SE_signed_change",
                   "CI_width_signed_change"]:
        d = random_values if metric.startswith("random_") else c
        if metric == "interaction_shift":
            d = d[d.term.str.contains(":", regex=False)]
        # Equal coefficient/region weight within input, equal input weight within
        # publication, and finally equal publication weight.
        by_input = d.groupby(WORKFLOW + ["dataset", "citation"])[metric].mean()
        by_source = by_input.groupby(WORKFLOW + ["citation"]).mean().rename("value")
        sources.append(by_source.reset_index().assign(metric=metric))
        for (policy, scale), z in d.groupby(WORKFLOW):
            counts.append(dict(
                policy=policy, transformation=scale, metric=metric,
                elementary_values=len(z), datasets=z.dataset.nunique(),
                region_jobs=len(z[JOBS].drop_duplicates()), sources=z.citation.nunique(),
            ))
    sources = pd.concat(sources, ignore_index=True)
    sources.to_csv(OUT / "publication_values.csv", index=False)
    counts = pd.DataFrame(counts)
    rng = np.random.default_rng(20261120)
    resamples = {}
    rows = []
    for (metric, policy, scale), z in sources.groupby(["metric"] + WORKFLOW):
        z = z.sort_values("citation")
        signature = tuple(z.citation)
        if signature not in resamples:
            resamples[signature] = rng.integers(0, len(z), size=(1000, len(z)))
        values = z.value.to_numpy()
        means = values[resamples[signature]].mean(axis=1)
        rows.append(dict(metric=metric, policy=policy, transformation=scale,
                         mean=float(values.mean()),
                         lower=float(np.quantile(means, .025)),
                         upper=float(np.quantile(means, .975))))
    summary = pd.DataFrame(rows).merge(
        counts, on=["metric"] + WORKFLOW, validate="one_to_one"
    )
    assert len(summary) == 8 * 13 * 6
    assert np.isfinite(summary[["mean", "lower", "upper"]]).all().all()
    assert np.allclose(summary.query("policy=='none' and transformation=='identity'")["mean"], 0)
    a = summary[summary.transformation == "ln"].set_index(["metric", "policy"])
    b = summary[summary.transformation == "log10"].set_index(["metric", "policy"])
    assert np.allclose(a["mean"], b["mean"], rtol=1e-12, atol=1e-12)
    directions = summary.metric.isin(["SE_signed_change", "CI_width_signed_change"])
    summary[directions].to_csv(OUT / "direction_values.csv", index=False)
    summary[~directions].to_csv(OUT / "figure_values.csv", index=False)
    caption = """# Preprocessing effects across all 78 requested specifications

## Read the figure

The reference is **no additional exclusion and untransformed RTs**, after
the common eligibility rules used in the manuscript. Each panel has 13
exclusion rows and six transformation columns. The outlined first row shows
transformations alone; the outlined first column shows exclusion alone.
All other cells combine exclusion and transformation. Natural and base-10
logs appear separately as requested, but have identical standardized outputs.
They are not independent inferential choices.
SD rows are symmetric mean ± k SD rules. The median-relative rule retains
0.5M–2.8M, where M is the condition–region median; MAD is scaled by 1.4826.
Fixed bounds are in milliseconds. These definitions are unchanged from the
completed manuscript analyses.

- **A: fixed-effect change.** Mean absolute change in the signed standardized
  coefficient, across all 380 non-intercept coefficient–region combinations.
  Units are the SD of each pipeline's transformed, retained outcome.
- **B: SE change.** Mean absolute percentage change in standardized coefficient
  SE, relative to the corresponding untrimmed, untransformed SE. Absolute
  changes prevent increases and decreases from cancelling across effects.
  Signed changes are provided separately in `direction_values.csv`.
  Near-zero effect estimates are never used as percentage denominators.
- **C: interval-width change.** Mean absolute percentage change in standardized
  95% confidence-interval width, including the fitted t degrees of freedom.
  This summarizes width, not movement of the interval's endpoints.
- **D: random-effect estimate change.** Mean absolute difference in shared
  VarCorr entries (variances AND covariances), excluding residual variance.
  Entries are in each pipeline's outcome-variance units (SD squared).
  Grouping suffixes are canonicalized by participant/item and named design
  term(s); rows are never matched by position. A removed or newly added entry
  has no comparable estimate and is NOT imputed as zero.
- **E: interaction-magnitude change.** Mean absolute difference between the
  magnitudes of corresponding standardized interaction coefficients.
  This uses the 150 interaction coefficient–region combinations, not the
  main effects. Multilevel contrasts are basis-dependent, not omnibus effects.
- **F: random-effect structure change.** Percentage of entries not shared
  between the two fitted random-effects covariance structures:
  100 × (1 − intersection size / union size). Both added and removed entries
  count; this is not the percentage of models with any structural change.

## Aggregation and uncertainty

Values are averaged within analytical input, then within publication,
then across publications with equal publication weights. A–C use 37
publications; E uses only publications contributing interaction coefficients.
Random-effect summaries first average comparable entries within a region job.
All five genuine scales and the algebraic log10 duplicate are shown, without
recounting that duplicate as independent evidence.

The heatmaps show mean magnitudes, not signed directions. `figure_values.csv`
also provides 95% percentile
publication-cluster bootstrap intervals (1,000 resamples, fixed seed),
source counts and the number of contributing coefficient/region values.
`publication_values.csv` gives the exact source-level quantities;
`random_component_support.csv` gives shared, added and removed entry counts.
Reference-cell changes are exactly zero. No model was refitted.

## Interpretation

Standardization removes arbitrary unit conversions, not the scientific
differences between nonlinear-scale estimands. These are descriptive shifts
relative to a stated reference, not estimates of bias or statistical accuracy.
Panel D is conditional on components estimable in both fitted structures;
panel F must be read alongside it. Small D with large F is not evidence of
unchanged random effects.

## Reproduce

From the bundle root, run:

    python3 "BRM impact/analysis/effect_sensitivity.py"
    Rscript "BRM impact/analysis/plot_effect_sensitivity.R"
    Rscript "BRM impact/analysis/verify_effect_sensitivity.R"

Python requires numpy and pandas. R uses base packages only. The two completed
output tables needed as inputs are included; the full model/data reproduction
remains in the manuscript reproducibility package. Generated outputs are in
`generated-artifacts/brm_multiverse/effect_sensitivity/`.
"""
    (OUT / "FIGURE_README.md").write_text(caption)
    print(summary.groupby("metric")["mean"].agg(["min", "max"]).to_string())
    print("Wrote", sum(~directions), "figure cells and", sum(directions),
          "signed-direction summaries to", OUT)


if __name__ == "__main__":
    main()
