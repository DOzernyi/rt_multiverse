"""Independently recompute manuscript totals and figure points with pandas."""
from pathlib import Path
import json
import itertools
import numpy as np
import pandas as pd

OUT = Path("generated-artifacts/brm_multiverse")
s = json.loads((OUT/"summary.json").read_text())
status = pd.read_csv(OUT/"pipeline_status.csv.gz")
c = pd.read_csv(OUT/"coefficients.csv.gz")
c = c[(c.term != "(Intercept)") & (c.transformation != "log10")].copy()
keys = ["dataset", "region", "term"]
base = c[(c.policy == "none") & (c.transformation == "ln")][keys+["p","formula"]]
base = base.rename(columns={"p":"reference_p", "formula":"reference_formula"})
c = c.merge(base, on=keys, validate="many_to_one")
c["disagreement"] = (c.p < .05) != (c.reference_p < .05)
pubmean = lambda d, v: d.groupby(["dataset","citation"])[v].mean().groupby("citation").mean().mean()
g = c.groupby(keys)
unstable = g.p.agg(lambda x: (x < .05).nunique() > 1).sum()
direction = g.estimate.agg(lambda x: np.sign(x).nunique() > 1).sum()
assert len(status) == s["pipelines"]["nominal"]
assert status.success.sum() == s["pipelines"]["accepted_nominal"]
assert len(g) == s["pipelines"]["coefficient_regions"]
assert int(unstable) == s["pipelines"]["unstable"]
assert int(direction) == s["pipelines"]["direction_changes"]
assert abs(pubmean(c,"disagreement")-s["pipelines"]["disagreement"]["mean"]) < 1e-12
same = c[c.formula == c.reference_formula]
assert abs(pubmean(same,"disagreement")-s["pipelines"]["same_structure_disagreement"]["mean"]) < 1e-12
ranges = g.standardized_estimate.max()-g.standardized_estimate.min()
se_ratio = g.standardized_SE.max()/g.standardized_SE.min()
for values, key in [(ranges,"effect_range"),(se_ratio,"precision_ratio")]:
    for q, label in [(.25,"q25"),(.5,"median"),(.75,"q75")]:
        assert abs(values.quantile(q)-s[key][label]) < 1e-12
inference = pd.read_csv(OUT/"inference_summary.csv")
for _, row in inference.iterrows():
    d = c[(c.policy == row.policy) & (c.transformation == row.transformation)]
    assert abs(pubmean(d,"disagreement")-row["mean"]) < 1e-12
retention = pd.read_csv(OUT/"retention_summary.csv")
for _, row in retention.iterrows():
    d = status[(status.policy == row.policy)&(status.transformation == "ln")]
    assert abs(pubmean(d,"proportion_excluded")-row["mean"]) < 1e-12

# Independent matched-pair decomposition: reconstruct directly from coefficients.
policies = ["none","sd1.5","sd2","sd2.5","sd3","fixed100_2000",
            "fixed100_3000","fixed200_2000","fixed200_3000","median",
            "mad2","mad2.5","mad3"]
scales = ["identity","ln","inverse","sqrt","boxcox"]
idx = keys + ["citation"]
p = c.pivot(index=idx, columns=["policy","transformation"], values="p")
t = c.pivot(index=idx, columns=["policy","transformation"], values="t")
f = c.pivot(index=idx, columns=["policy","transformation"], values="formula")
assert p.shape == (380,65) and p.notna().all().all()
source_frames = []
for dimension, levels, held in [
    ("exclusion",policies,scales),("transformation",scales,policies)
]:
    for a,b in itertools.combinations(levels,2):
        frames = []
        for h in held:
            ka,kb = ((a,h),(b,h)) if dimension=="exclusion" else ((h,a),(h,b))
            sa,sb = p[ka]<.05,p[kb]<.05
            same = f[ka]==f[kb]
            v = pd.DataFrame({
                "disagreement":(sa!=sb).astype(float),
                "significant_gain":((~sa)&sb).astype(float),
                "significant_loss":(sa&(~sb)).astype(float),
                "absolute_t_shift":abs(t[ka]-t[kb]),
                "same_formula_fraction":same.astype(float),
                "same_formula_disagreement":(sa!=sb).astype(float).where(same)
            })
            frames.append(v)
        v = pd.concat(frames)
        v = v.groupby(level=["dataset","citation"]).mean().groupby("citation").mean()
        v = v.reset_index().assign(dimension=dimension,level_a=a,level_b=b)
        source_frames.append(v)
source_values = pd.concat(source_frames,ignore_index=True)
saved_sources = pd.read_csv(OUT/"decision_pair_sources.csv")
join_keys = ["dimension","level_a","level_b","citation"]
checked = source_values.merge(saved_sources,on=join_keys,
                              suffixes=("_python","_r"),validate="one_to_one")
assert len(checked)==len(saved_sources)==len(source_values)
metrics = ["disagreement","significant_gain","significant_loss","absolute_t_shift",
           "same_formula_fraction","same_formula_disagreement"]
for m in metrics:
    assert np.allclose(checked[m+"_python"],checked[m+"_r"],
                       atol=1e-12,rtol=1e-12,equal_nan=True),m
pair_summary = pd.read_csv(OUT/"decision_pair_summary.csv")
for _,row in pair_summary.iterrows():
    v = source_values[(source_values.dimension==row.dimension)&
                      (source_values.level_a==row.level_a)&
                      (source_values.level_b==row.level_b)]
    assert abs(v.disagreement.mean()-row["mean"])<1e-12
    assert row.lower<=row["mean"]<=row.upper
marginals = pd.read_csv(OUT/"decision_marginals.csv")
for _,row in marginals.iterrows():
    v = source_values[source_values.dimension==row.dimension]
    if row.level!="ALL":
        v = v[(v.level_a==row.level)|(v.level_b==row.level)]
    expected = v.groupby("citation")[metrics].mean().mean()
    assert abs(expected.disagreement-row["mean"])<1e-12
    assert abs(expected.absolute_t_shift-row.absolute_t_shift)<1e-12
    assert row.lower<=row["mean"]<=row.upper
for dimension in ["exclusion","transformation"]:
    row = marginals[(marginals.dimension==dimension)&(marginals.level=="ALL")].iloc[0]
    assert abs(row["mean"]-s["decision_sensitivity"][dimension]["mean"])<1e-12
assert s["decision_sensitivity"]["comparisons"] == 197600
delta = s["decision_sensitivity"]["exclusion"]["mean"]-s["decision_sensitivity"]["transformation"]["mean"]
assert abs(delta-s["decision_sensitivity"]["exclusion_minus_transformation"]["mean"])<1e-12
width = (c.CI_high-c.CI_low)*c.standardized_SE/c.SE
widths = c.assign(width=width).groupby(keys).width.agg(lambda a:a.max()/a.min())
saved_effects = pd.read_csv(OUT/"effect_stability.csv").set_index(keys)
assert np.allclose(widths.sort_index(),
                   saved_effects.CI_width_ratio.sort_index(),atol=1e-11,rtol=1e-12)
kind_summary = pd.read_csv(OUT/"effect_kind_summary.csv")
for _,row in kind_summary.iterrows():
    z = saved_effects[saved_effects.effect_kind==row.effect_kind]
    assert len(z)==row.coefficient_regions
    assert abs(z.standardized_range.median()-row.range_median)<1e-12
    assert abs(z.SE_ratio.median()-row.SE_ratio_median)<1e-12
    assert abs(z.CI_width_ratio.median()-row.CI_width_ratio_median)<1e-12
comparisons = pd.read_csv(OUT/"coefficient_comparisons.csv.gz")
assert len(comparisons) == len(c)
text = (OUT/"compiled/generated/results_text.tex").read_text()
assert f"{int(unstable):,}" in text
body = Path("BRM impact/analysis/manuscript_body.tex").read_text()
assert r"\label{fig:estimate-sensitivity}" in body
assert r"\ref{fig:estimate-sensitivity}" in text
new_figure = OUT/"effect_sensitivity"
assert len(pd.read_csv(new_figure/"figure_values.csv")) == 468
assert (new_figure/"verification.txt").read_text().startswith("PASSED:")
for suffix in [".eps",".pdf"]:
    assert (OUT/"compiled/figures"/("preprocessing_effects"+suffix)).read_bytes() == (
        new_figure/("preprocessing_effects"+suffix)).read_bytes()
runtime = json.loads((OUT/"release_audit/runtime_summary.json").read_text())
job_times = pd.read_csv(OUT/"release_audit/jobs.csv")
assert abs(job_times.elapsed_seconds.sum()-runtime["recorded_job_elapsed_seconds"])<1e-6
assert abs(runtime["recorded_job_elapsed_seconds"]-s["computational_seconds"])<1e-6
footnote = (OUT/"compiled/generated/runtime_footnote.tex").read_text()
assert f"{runtime['recorded_job_elapsed_hours']:.1f}" in footnote
assert "not CPU time or end-to-end wall time" in footnote
discussion = body.split(r"\section{Discussion}",1)[1].split(r"\section{Reproducibility",1)[0]
assert discussion.rsplit(r"\subsection{",1)[1].startswith("Practical advice and recommendations}")
report = {
    "passed": True,
    "independent_engine": "Python/pandas, separate from R summarization",
    "coefficient_regions": len(g),
    "unstable": int(unstable),
    "direction_changes": int(direction),
    "publication_weighted_disagreement": float(pubmean(c,"disagreement")),
    "inference_figure_points_checked": len(inference),
    "retention_figure_points_checked": len(retention),
    "matched_decision_pairs_checked": len(pair_summary),
    "matched_decision_source_rows_checked": len(saved_sources),
    "matched_decision_marginals_checked": len(marginals),
    "matched_pair_comparisons": s["decision_sensitivity"]["comparisons"],
    "checks": ["pipeline totals", "significance and direction changes",
               "publication weighting", "same-formula restriction",
               "effect-range and precision quartiles", "figure point values",
               "all matched decision comparisons and publication weights",
               "main-effect and interaction confidence-interval width summaries",
               "six-panel estimate-sensitivity figure included unchanged",
               "recommendations are final discussion subsection"]
}
(OUT/"independent_summary_verification.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report,indent=2))
