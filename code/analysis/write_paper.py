"""Build the empirical manuscript only from verified, completed primary outputs."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import zipfile
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
ANALYSIS = ROOT / "analysis"
TEMPLATE = ROOT / "_BRN_SN__multiverse_preprocessing"
OUT = Path(os.environ.get("BRM_MULTIVERSE_OUT", WORKSPACE / "generated-artifacts/brm_multiverse"))
TITLE = "Preprocessing decisions and statistical inference in psycholinguistic reading tasks: a multiverse analysis"


def tex(value):
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
                    "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in str(value))


def pct(value):
    return f"{100*value:.1f}" + r"\%"


def interval(x):
    return (pct(x["mean"]) + r" (publication-cluster bootstrap 95\% interval "
            + pct(x["lower"]) + "--" + pct(x["upper"]) + ")")


def extract_entries(bib, wanted):
    # Registry identifier spelling differs from the existing bibliography key;
    # this is the same verified article, not an invented reference.
    bib = bib.replace("@article{fujita_reassessing_2026,",
                      "@article{fujita_reassessment_2026,")
    records = {}
    for match in re.finditer(r"(?m)^@\w+\s*\{\s*([^,\s]+)\s*,", bib):
        key = match.group(1)
        if key not in wanted:
            continue
        begin = bib.index("{", match.start())
        depth = 0
        for end in range(begin, len(bib)):
            if bib[end] == "{" and bib[end-1] != "\\":
                depth += 1
            elif bib[end] == "}" and bib[end-1] != "\\":
                depth -= 1
                if depth == 0:
                    entry = bib[match.start():end+1]
                    # DOI is sufficient for these entries; omit duplicate
                    # publisher URLs that the supplied APA style prints unbroken.
                    if re.search(r"(?im)^\s*doi\s*=", entry):
                        entry = re.sub(r"(?im)^\s*(?:url|urldate)\s*=\s*\{[^}]*\},?\s*\n", "", entry)
                    records[key] = entry
                    break
    if wanted - records.keys():
        raise ValueError("Missing references: " + repr(wanted - records.keys()))
    return "\n\n".join(records[k] for k in sorted(records)) + "\n"


def host(m):
    identifier = m["source"].split("/")[1].split("-")[0]
    if re.fullmatch(r"[a-z0-9]{5}", identifier):
        return rf"\href{{https://osf.io/{identifier}/}}{{{identifier}}}"
    if identifier == "github":
        return r"\href{https://vboyce.github.io/Maze/}{GitHub}"
    if identifier == "figshare":
        return r"\href{https://figshare.com/projects/18823}{Figshare}"
    return tex(identifier)


def table(registry):
    lines = [r"\scriptsize",
        r"Each row is an analytical input, not necessarily an independent cohort. "
        r"Observation counts sum only source-selected analysis regions, before multiverse exclusions. "
        r"The original-policy columns describe the source paper; discrepancies with archived "
        r"code and exact excerpts are preserved in \texttt{paper\_rt\_policies.csv}.",
        r"\begin{longtable}{@{}p{3.05cm}p{.85cm}p{.65cm}p{.85cm}rrp{5.1cm}p{4.4cm}@{}}",
        r"\caption{Primary inputs, source hosts, samples, and source-reported preprocessing.}\\",
        r"\toprule Reference / input & Host & Task & Design & $N$ & Obs. & Original exclusion policy & Original transformation/model scale\\\midrule",
        r"\endfirsthead",
        r"\toprule Reference / input & Host & Task & Design & $N$ & Obs. & Original exclusion policy & Original transformation/model scale\\\midrule",
        r"\endhead\bottomrule\endfoot"]
    for m in registry:
        lines.append(rf"{m['display_id']}: \citep{{{m['citation']}}}, {tex(m['label'])} & "
                     rf"{host(m)} & {m['task']} & {m['design'].replace('x',r'$\times$')} & "
                     rf"{m['participants']} & {m['observations']:,} & "
                     rf"{tex(m['paper_rt_outlier_policy'])} & "
                     rf"{tex(m['paper_rt_transform_or_model_scale'])}\\")
    lines += [r"\end{longtable}\normalsize",
              r"\subsection*{Design and provenance mappings}",
              r"\scriptsize\begin{longtable}{@{}p{5cm}p{16cm}@{}}\toprule Input & Factor levels and common eligibility\\\midrule\endhead"]
    for m in registry:
        mapping = "; ".join(f"{f} = {m['factors'][f]}: " + ", ".join(v)
                            for f, v in m["levels"].items())
        lines.append(m["display_id"] + ": " + tex(m["dataset"]) + " & " + tex(mapping + ". " + m["eligibility"]) + r"\\")
    lines.append(r"\bottomrule\end{longtable}\normalsize")
    return "\n".join(lines) + "\n"


def result_text(s):
    p = s["pipelines"]
    r = pd.read_csv(OUT / "retention_summary.csv")
    subgroup = pd.read_csv(OUT / "subgroup_summary.csv")
    q = s["precision_ratio"]
    e = s["effect_range"]
    decisions = s["decision_sensitivity"]
    pairs = pd.read_csv(OUT / "decision_pair_summary.csv")
    marginals = pd.read_csv(OUT / "decision_marginals.csv")
    kinds = pd.read_csv(OUT / "effect_kind_summary.csv")
    largest_cell = pd.read_csv(OUT / "inference_summary.csv").sort_values(
        "mean", ascending=False).iloc[0]
    names = {
        "none": "no additional exclusion", "sd1.5": r"mean \(\pm1.5\) SD",
        "sd2": r"mean \(\pm2\) SD", "sd2.5": r"mean \(\pm2.5\) SD",
        "sd3": r"mean \(\pm3\) SD", "fixed100_2000": "100--2,000 ms",
        "fixed100_3000": "100--3,000 ms", "fixed200_2000": "200--2,000 ms",
        "fixed200_3000": "200--3,000 ms", "median": "median-relative bounds",
        "mad2": r"median \(\pm2\) MAD", "mad2.5": r"median \(\pm2.5\) MAD",
        "mad3": r"median \(\pm3\) MAD", "identity": "untransformed milliseconds",
        "ln": "natural logarithms", "inverse": "the inverse transformation",
        "sqrt": "square roots", "boxcox": "Box--Cox"}
    def pair(dimension, a, b):
        return pairs[(pairs.dimension == dimension) &
                     (((pairs.level_a == a) & (pairs.level_b == b)) |
                      ((pairs.level_a == b) & (pairs.level_b == a)))].iloc[0]
    strongest = {
        dim: marginals[(marginals.dimension == dim) &
                       (marginals.level != "ALL")].sort_values("mean", ascending=False).iloc[0]
        for dim in ("exclusion", "transformation")}
    largest_pairs = {
        dim: pairs[pairs.dimension == dim].sort_values("mean", ascending=False).iloc[0]
        for dim in ("exclusion", "transformation")}
    lines = [r"\subsection{Coverage and fit availability}",
        rf"The completed primary multiverse comprises {s['corpus']['region_jobs']} region-level jobs "
        rf"and {p['nominal']:,} nominal outputs. Of these, {p['accepted_nominal']:,} "
        rf"yielded accepted inferential results. After collapsing the equivalent logarithms, "
        rf"{p['accepted_distinct']:,} of {p['distinct']:,} distinct workflows were accepted "
        rf"({pct(p['accepted_distinct']/p['distinct'])}); "
        rf"{p['failed_distinct']:,} were unsuccessful. There were {p['fit_attempts']:,} "
        r"recorded mixed-model fitting attempts, excluding the analytically derived base-10 "
        r"outputs. Failed workflows are unavailable inference, not evidence of no effect.",
        r"\subsection{Significance stability}",
        rf"Across successful distinct pipelines, {p['unstable']:,} of "
        rf"{p['coefficient_regions']:,} coefficient--region combinations "
        rf"({pct(p['unstable']/p['coefficient_regions'])}) included both "
        rf"significant and nonsignificant results. The publication-weighted rate of "
        rf"disagreement with the no-additional-exclusion natural-log reference was "
        rf"{interval(p['disagreement'])}. This is an average over the declared workflows, "
        r"not a probability that a source finding is incorrect.",
        rf"The accepted reference covered {p['reference_coefficients']:,} "
        rf"coefficient--region combinations. Restricting comparisons to an identical "
        rf"final random-effects formula gave {interval(p['same_structure_disagreement'])} "
        r"disagreement. This restriction is conditional on fit availability; it does not "
        r"identify the causal contribution of random-effects reduction.",
        r"Figure~\ref{fig:inference} displays the complete exclusion-by-scale surface. "
        r"Each cell averages coefficients within datasets and then gives publications "
        r"equal weight. A low disagreement rate means similarity to the chosen reference, "
        r"not superior recovery of a true effect.",
        rf"The largest departure from that reference occurred for "
        rf"{names[largest_cell['policy']]} combined with "
        rf"{names[largest_cell['transformation']]} "
        rf"({pct(largest_cell['mean'])}). This combined departure must not be "
        r"mistaken for the isolated influence of either decision.",
        r"\subsection{Which preprocessing decisions change inference most?}",
        rf"Changing only the exclusion rule, with the scale held fixed, produced "
        rf"{interval(decisions['exclusion'])} disagreement in \(p<.05\) classification. "
        rf"Changing only the transformation, with the retained sample held fixed, "
        rf"produced {interval(decisions['transformation'])}. The paired "
        rf"publication-weighted difference was "
        rf"{100*decisions['exclusion_minus_transformation']['mean']:.1f} percentage points "
        rf"(95\% cluster interval "
        rf"{100*decisions['exclusion_minus_transformation']['lower']:.1f}--"
        rf"{100*decisions['exclusion_minus_transformation']['upper']:.1f}). "
        r"Thus exclusion choices were more consequential on this metric, averaged "
        r"over this declared decision space; transformation was not negligible. "
        r"Figure~\ref{fig:decisions} summarizes matched comparisons rather than "
        r"confounding both decisions in departures from one reference.",
        rf"Among individual rules, {names[strongest['exclusion']['level']]} "
        rf"had the highest average disagreement with the other 12 rules "
        rf"({interval(strongest['exclusion'])}, averaging over five scales). "
        rf"Among scales, {names[strongest['transformation']['level']]} had the highest "
        rf"average disagreement with the other four scales "
        rf"({interval(strongest['transformation'])}, averaging over 13 rules). "
        r"These rankings describe how often a decision changes the conclusion, "
        r"not how well it recovers a true effect.",
        rf"The largest exclusion-rule contrast was "
        rf"{names[largest_pairs['exclusion']['level_a']]} versus "
        rf"{names[largest_pairs['exclusion']['level_b']]} "
        rf"({pct(largest_pairs['exclusion']['mean'])}); the largest scale contrast was "
        rf"{names[largest_pairs['transformation']['level_a']]} versus "
        rf"{names[largest_pairs['transformation']['level_b']]} "
        rf"({pct(largest_pairs['transformation']['mean'])}). "
        rf"For familiar anchored comparisons, changing from no exclusion to "
        rf"mean \(\pm1.5\) SD changed {pct(pair('exclusion','none','sd1.5')['mean'])} "
        rf"of classifications, whereas adopting median-relative bounds changed "
        rf"{pct(pair('exclusion','none','median')['mean'])}, holding scale fixed. "
        rf"Switching from natural logarithms to the inverse changed "
        rf"{pct(pair('transformation','ln','inverse')['mean'])}, whereas switching "
        rf"to Box--Cox changed {pct(pair('transformation','ln','boxcox')['mean'])}, "
        r"holding exclusion fixed. Averaging masks variation across contrasts and "
        r"does not imply that these workflows agree for every linguistic effect.",
        rf"The directions of change differed. Relative to no exclusion, "
        rf"mean \(\pm1.5\) SD changed "
        rf"{pct(pair('exclusion','none','sd1.5')['significant_gain'])} "
        rf"of all matched tests from nonsignificant to significant and "
        rf"{pct(pair('exclusion','none','sd1.5')['significant_loss'])} "
        r"in the reverse direction. Relative to natural logarithms, the inverse "
        rf"changed {pct(pair('transformation','ln','inverse')['significant_loss'])} "
        rf"from significant to nonsignificant and "
        rf"{pct(pair('transformation','ln','inverse')['significant_gain'])} "
        r"in the reverse direction. These are percentages of all matched tests, "
        r"not probabilities conditional on initially being significant. More "
        r"significant tests do not establish greater power, and fewer do not "
        r"establish better false-positive control.",
        rf"The mean absolute change in the coefficient's signed \(t\) statistic "
        rf"was {decisions['exclusion']['absolute_t_shift']:.2f} for an exclusion-rule "
        rf"change and {decisions['transformation']['absolute_t_shift']:.2f} for a "
        r"scale change. This supplements the binary threshold metric; it is "
        r"not an estimate of Type I or Type II error.",
        rf"Exact final formulas agreed in "
        rf"{pct(decisions['exclusion']['same_formula_fraction'])} of exclusion "
        rf"comparisons and {pct(decisions['transformation']['same_formula_fraction'])} "
        rf"of scale comparisons. Within the selected same-formula comparisons, "
        rf"the corresponding disagreement summaries were "
        rf"{pct(decisions['exclusion']['same_formula_disagreement'])} and "
        rf"{pct(decisions['transformation']['same_formula_disagreement'])}. "
        r"These differently selected subsets cannot isolate a causal covariance-model "
        r"contribution, but show that disagreement is not confined to pairs with "
        r"different formula strings.",
        r"\subsection{Removal, effect size, and precision}"]
    lo = r.loc[r["policy"] != "none"].sort_values("mean").iloc[0]
    hi = r.sort_values("mean").iloc[-1]
    med = r.loc[r["policy"] == "median"].iloc[0]
    lines.append(
        rf"Publication-weighted additional exclusion ranged from {pct(lo['mean'])} "
        rf"under \texttt{{{tex(lo['policy'])}}} to {pct(hi['mean'])} "
        rf"under \texttt{{{tex(hi['policy'])}}}, excluding the no-removal reference. "
        rf"The median-relative rule removed {pct(med['mean'])} "
        rf"(95\% cluster interval {pct(med['lower'])}--{pct(med['upper'])}; "
        r"Figure~\ref{fig:retention}). These percentages concern the available eligible "
        r"observations, not all trials originally collected.")
    lines.append(
        rf"The within-coefficient range of standardized estimates had median "
        rf"{e['median']:.3f} (interquartile range {e['q25']:.3f}--{e['q75']:.3f}). "
        rf"Estimated direction changed in {p['direction_changes']:,} of "
        rf"{p['coefficient_regions']:,} coefficient--region combinations "
        rf"({pct(p['direction_changes']/p['coefficient_regions'])}). "
        r"Figure~\ref{fig:precision} summarizes all coefficient--region ranges, "
        r"rather than selecting one coefficient per dataset. All individual estimates "
        r"and confidence limits remain available in the numerical supplement. "
        r"Changing the outcome scale changes the estimand and its SD denominator, "
        r"so these ranges are descriptive rather than estimates of bias.")
    lines.append(
        rf"The maximum-to-minimum standardized SE ratio within a coefficient--region "
        rf"combination had median {q['median']:.2f} "
        rf"(interquartile range {q['q25']:.2f}--{q['q75']:.2f}; "
        r"Figure~\ref{fig:precision}). A ratio of two means that the largest "
        r"standardized SE is twice the smallest. This measure describes sensitivity "
        r"of uncertainty to the complete workflow, including model reduction; "
        r"it is not a direct comparison of raw and logarithmic SE units.")
    lines += [
        r"\begin{table}[t]\centering\small",
        r"\caption{Effect and precision stability by coefficient type. Entries are "
        r"medians over coefficient--region combinations, not pooled effects. "
        r"Range is the maximum minus minimum standardized estimate; SE and interval "
        r"ratios compare the largest with the smallest value over 65 workflows.}",
        r"\begin{tabular}{lrrrr}\toprule Type & \(n\) & Range & SE ratio & CI-width ratio\\\midrule"]
    for _, x in kinds.iterrows():
        lines.append(rf"{tex(x['effect_kind'])} & {int(x['coefficient_regions'])} & "
                     rf"{x['range_median']:.3f} & {x['SE_ratio_median']:.2f} & "
                     rf"{x['CI_width_ratio_median']:.2f}\\")
    lines += [r"\bottomrule\end{tabular}\label{tab:effect-kinds}\end{table}",
              r"Table~\ref{tab:effect-kinds} separates main-effect and interaction "
              r"coefficients. Interval widths are expressed in each pipeline's "
              r"outcome-SD units, using the fitted degrees of freedom. Full native-scale "
              r"intervals and random-effect variances accompany the numerical outputs. "
              r"Because covariance structures are adapted separately, changes in "
              r"random-effect estimates combine retention, scale, and component removal; "
              r"we do not interpret them as change in one invariant population variance."]
    lines.append(
        r"Figure~\ref{fig:estimate-sensitivity} complements these range summaries "
        r"with a complete exclusion-by-transformation display of effect, SE, "
        r"confidence-interval-width, random-effect and interaction-magnitude changes. "
        r"Its no-removal, untransformed reference differs deliberately from the "
        r"natural-log reference used for significance classifications: the first "
        r"row then isolates transformation alone, the first column exclusion alone, "
        r"and the remaining cells their combinations. The paired covariance panels "
        r"separate shifts in shared estimates from changes in which entries are "
        r"estimated. The two logarithmic bases have identical standardized values "
        r"and are not counted as independent choices.")
    lines.append(
        rf"Estimated Box--Cox parameters ranged from {s['lambda']['min']:.3f} "
        rf"to {s['lambda']['max']:.3f}; {s['lambda']['boundaries']} were within "
        r".01 of a search endpoint (Figure~\ref{fig:lambda}). The distributions "
        r"vary across the retained samples. The profile is not a test that "
        r"mixed-model residuals are normal, nor a reason to select transformations "
        r"according to which produces the preferred inferential outcome.")
    lines.append(r"\subsection{Task, contrast, region, and design comparisons}")
    lines += [r"\begin{table*}[t]\centering\small",
              r"\caption{Descriptive subgroup stability. Sources may contribute to more than one row. Intervals resample publications, not pipelines.}",
              r"\begin{tabular}{llrrrr}\toprule Group & Level & Coefficient--regions & Sources & Disagreement (\%) & 95\% interval\\\midrule"]
    for _, x in subgroup.iterrows():
        lines.append(rf"{tex(x['group'])} & {tex(x['level'])} & "
                     rf"{int(x['coefficient_regions'])} & {int(x['sources'])} & "
                     rf"{100*x['mean']:.1f} & {100*x['lower']:.1f}--{100*x['upper']:.1f}\\")
    lines += [r"\bottomrule\end{tabular}\label{tab:subgroups}\end{table*}",
              r"Table~\ref{tab:subgroups} distinguishes "
              r"task, main-effect versus interaction coefficients, analysis region, and "
              r"factorial design. These groups overlap in source composition, and "
              r"designs and tasks are not randomly assigned across publications. "
              r"Consequently, differences are not adjusted task effects. Coefficients "
              r"for multilevel factors are basis-dependent contrasts, not omnibus tests.",
              r"\subsection{Exploratory dataset-level moderators}"]
    for variable, label in [("skewness", "Initial RT skewness"),
                            ("participants", "Participant count"),
                            ("items", "Item count"),
                            ("reference_magnitude", "Absolute reference-effect magnitude")]:
        x = s["moderators"][variable]
        lines.append(rf"{label} correlated with dataset-level disagreement at "
                     rf"\(\rho={x['rho']:.2f}\) (Spearman; publication-cluster "
                     rf"bootstrap 95\% interval {x['lower']:.2f} to {x['upper']:.2f}; "
                     rf"{x['datasets']} inputs).")
    lines.append(
        r"These associations are exploratory and not mutually adjusted. The "
        r"reference magnitude is itself estimated from these data; its association "
        r"with threshold crossing is partly coupled to the definition of disagreement. "
        r"No moderator establishes a mechanism or validates a cutoff.")
    return "\n\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archived-inputs", action="store_true",
                        help="Use shipped compact inputs and source-hash ledger without the original source archive.")
    args = parser.parse_args()
    audit = json.loads((OUT / "verification_results.json").read_text())
    if not audit["passed"]:
        raise RuntimeError("Results verification failed")
    registry = json.loads((OUT / "datasets.json").read_text())
    exclusions = json.loads((ANALYSIS / "excluded_datasets.json").read_text())
    excluded = {m["dataset"] for m in exclusions}
    registry = [m for m in registry if m["tier"] == "primary" and m["dataset"] not in excluded]
    for i, m in enumerate(registry, 1):
        m["display_id"] = f"D{i:02d}"
    policies = pd.read_csv(ROOT / "corpus_audit/paper_rt_policies.csv").fillna("").set_index("dataset")
    ledger = pd.read_csv(ROOT / "corpus_audit/eligible_analysis_inputs.csv").set_index("dataset")
    source_checks = []
    archived_checks = {}
    if args.archived_inputs:
        archived_checks = {x["dataset"]: x for x in
                           json.loads((OUT / "source_verification.json").read_text())}
    for m in registry:
        # Require current source audit, not the old filename-based tier assignment.
        assert ledger.loc[m["dataset"], "tier"] == "primary"
        assert ledger.loc[m["dataset"], "sha256"] == m["sha256"]
        source = ROOT / m["source"]
        if args.archived_inputs:
            assert archived_checks[m["dataset"]]["sha256"] == m["sha256"]
        else:
            assert hashlib.sha256(source.read_bytes()).hexdigest() == m["sha256"], str(source)
        source_checks.append({"dataset": m["dataset"], "source": m["source"], "sha256": m["sha256"]})
        m.update(policies.loc[m["dataset"]].to_dict())
    if not args.archived_inputs:
        (OUT / "source_verification.json").write_text(json.dumps(source_checks, indent=2))
    s = json.loads((OUT / "summary.json").read_text())
    runtime = json.loads((OUT / "release_audit/runtime_summary.json").read_text())
    assert abs(runtime["recorded_job_elapsed_seconds"] - s["computational_seconds"]) < .01
    assert len(registry) == s["corpus"]["datasets"] == audit["datasets"]
    assert audit["region_jobs"] == s["corpus"]["region_jobs"]
    directory = OUT / "compiled"
    directory.mkdir(parents=True, exist_ok=True)
    for name in ["generated", "figures"]:
        (directory / name).mkdir(exist_ok=True)
    runtime_footnote = (
        rf"\footnote{{Refitting the complete multiverse is computationally expensive. "
        rf"The {runtime['complete_region_jobs']} completed region jobs record "
        rf"{runtime['recorded_job_elapsed_hours']:.1f} hours of cumulative job elapsed time; "
        rf"the longest job took {runtime['longest_job_seconds']/3600:.1f} hours. "
        r"These are summed worker durations, not CPU time or end-to-end wall time: "
        r"jobs overlapped in parallel. Recorded resumptions are included, but abandoned "
        r"or excluded runs and unrecorded setup are not. A new full fit can require "
        r"multiple days, depending on hardware and parallelism. Supplied checkpoints "
        r"allow tables and figures to be regenerated without repeating those fits.}"
    )
    (directory / "generated/runtime_footnote.tex").write_text(runtime_footnote + "\n")
    original = (TEMPLATE / "paper_new.tex").read_text()
    class_line = next(x for x in original.splitlines() if x.startswith(r"\documentclass"))
    preamble = class_line + "\n" + r"""
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{mathtools,amssymb,booktabs,longtable,graphicx,pdflscape}
\usepackage{natbib}
\usepackage{xurl}
\usepackage{tikz}
\usetikzlibrary{arrows.meta,positioning,fit}
\newcommand{\rtprocess}{\texttt{rt\_process()}}
\hypersetup{colorlinks=true,citecolor=blue,linkcolor=blue,urlcolor=blue}
\raggedbottom
"""
    frontmatter = original.split(r"\begin{document}", 1)[1].split(r"\section{Introduction}", 1)[0]
    frontmatter = frontmatter.replace(r"\title[title]{title}",
                                    rf"\title[Preprocessing multiverse]{{{TITLE}}}")
    b = s["pipelines"]["disagreement"]
    abstract = (
        f"Reading-time preprocessing changes both the observations entering an "
        f"analysis and the scale of the estimated effects. We crossed 13 exclusion "
        f"rules with six transformation specifications in {len(registry)} source-audited "
        f"primary analytical inputs from {s['corpus']['publications']} publications "
        f"using self-paced reading or Maze tasks. Critical regions and "
        f"source-justified post-critical regions yielded {s['corpus']['region_jobs']} "
        f"analysis jobs. Equivalent logarithms were collapsed for stability summaries, "
        f"leaving 65 distinct workflows per job. Design-supported mixed models "
        f"followed an explicit convergence and singularity reduction sequence. "
        f"{s['pipelines']['unstable']} of {s['pipelines']['coefficient_regions']} "
        f"coefficient--region combinations crossed the unadjusted significance "
        f"threshold across successful workflows. Publication-weighted disagreement "
        f"with a no-additional-exclusion natural-log reference was "
        f"{100*b['mean']:.1f}\\% (95\\% cluster-bootstrap interval "
        f"{100*b['lower']:.1f}--{100*b['upper']:.1f}\\%). "
        f"In matched comparisons, changing only exclusion altered "
        f"{100*s['decision_sensitivity']['exclusion']['mean']:.1f}\\% of "
        f"classifications, compared with "
        f"{100*s['decision_sensitivity']['transformation']['mean']:.1f}\\% "
        f"when changing only scale. The 1.5-SD rule and inverse transformation "
        f"were the most disruptive choices within their respective dimensions. "
        f"We separately quantify direction, standardized effect-size and precision "
        f"stability and provide all pipeline failures. These results characterize "
        f"dependence on the declared analytical space, not the superiority of "
        f"one preprocessing rule. We recommend explicit bounds and grouping, "
        f"scale-specific interpretation, and transparent sensitivity reporting.")
    frontmatter = re.sub(r"\\abstract\{\s*\}", lambda _: r"\abstract{" + abstract + "}", frontmatter)
    frontmatter = frontmatter.replace(r"\keywords{}", r"\keywords{reading times, preprocessing, multiverse, mixed-effects models, self-paced reading, Maze}")
    corpus = (
        rf"The primary reanalysis contains {len(registry)} analytical inputs from "
        rf"{s['corpus']['publications']} publications: {s['corpus']['SPR']} SPR "
        rf"and {s['corpus']['Maze']} Maze inputs. These yield "
        rf"{s['corpus']['region_jobs']} dataset--region jobs "
        rf"({s['corpus']['critical_jobs']} critical and "
        rf"{s['corpus']['postcritical_jobs']} source-justified post-critical jobs), "
        rf"totalling {s['corpus']['observations']:,} eligible region-level observations. "
        r"The folder audit identified 80 primary analytical inputs and 10 "
        r"previously prepared inputs; one primary input "
        r"(\texttt{gastmann2026\_spr}) was excluded after exceptionally long "
        r"computation, before the final aggregation. This computational exclusion "
        r"is not evidence of a null result. The 10 prepared inputs are inventoried "
        r"separately but are not included in the numerical results reported here.")
    (directory / "generated/corpus_text.tex").write_text(corpus + "\n")
    (directory / "generated/results_text.tex").write_text(result_text(s))
    p = s["pipelines"]
    kind = pd.read_csv(OUT / "effect_kind_summary.csv").set_index("effect_kind")
    discussion = (
        r"\subsection{Which decisions matter most?}" "\n\n"
        rf"The clearest practical result is that exclusion choices changed "
        rf"inferential classifications more often than transformation choices "
        rf"in this multiverse. In matched comparisons, exclusion changes "
        rf"produced {pct(s['decision_sensitivity']['exclusion']['mean'])} disagreement "
        rf"and scale changes produced "
        rf"{pct(s['decision_sensitivity']['transformation']['mean'])}; "
        rf"the difference was "
        rf"{100*s['decision_sensitivity']['exclusion_minus_transformation']['mean']:.1f} "
        r"percentage points. The most disruptive rule was mean \(\pm1.5\) SD; "
        r"the most disruptive transformation was the inverse. These are not rankings "
        r"of validity. They identify the decisions for which a change of conclusion "
        r"most needs to be shown and explained. The comparison holds one preprocessing "
        r"dimension fixed, but retains the declared model-reduction procedure in both "
        r"workflows. It therefore compares operational analysis choices, not a pure "
        r"causal effect of deletion or scale with covariance parameters held constant."
        "\n\n"
        r"The two leading choices act differently. The 1.5-SD rule removes "
        r"observations using a relatively narrow condition-specific interval, "
        r"changing the retained tails and the information available to estimate "
        r"effects and variation. The inverse, \(-1000/x\), strongly compresses long "
        r"responses and is more locally sensitive to differences among shorter "
        r"responses. The empirical rankings are consistent with those mechanisms, "
        r"but do not establish that either set of affected observations was invalid. "
        r"Nor is removal percentage itself an adequate explanation: changing scale "
        r"can alter inference without removing any further observation, and "
        r"the placement of excluded observations matters as well as their number."
        "\n\n"
        r"The direction of the anchored changes is especially informative. "
        r"The 1.5-SD rule more often created than removed a significant result "
        r"relative to no exclusion, whereas the inverse more often removed "
        r"than created significance relative to natural logarithms. These "
        r"opposing patterns preclude treating all preprocessing as a generic "
        r"increase or decrease in evidence. They also make outcome-driven "
        r"selection consequential: an analyst seeking more threshold crossings "
        r"could prefer one branch, whereas an analyst seeking fewer could prefer "
        r"another. Neither preference is validated by these empirical data, "
        r"because the true status of the effects is unknown."
        "\n\n"
        r"Natural-log and Box--Cox workflows were comparatively similar, whereas "
        r"raw and inverse workflows were the most discrepant scale pair. Thus, "
        r"the useful scale sensitivity is not the nominal contrast between two "
        r"logarithmic bases, which give identical tests. It is the contrast between "
        r"scales with materially different weighting of the response distribution "
        r"and interpretations of factorial effects. This directly supports "
        r"scale-specific interpretation rather than a search for the transformation "
        r"that supplies the desired significance decision "
        r"\citep{lo_transform_2015}."
        "\n\n"
        rf"A separate result is the difference between how often a contrast "
        rf"can change and how often an individual workflow differs from a reference. "
        rf"{pct(p['unstable']/p['coefficient_regions'])} of coefficient--region "
        rf"combinations crossed \(p=.05\) somewhere in the declared space, whereas "
        rf"the publication-weighted disagreement rate across workflows was "
        rf"{pct(p['disagreement']['mean'])}. The former is an ``at least once'' "
        r"criterion over 65 workflows and is necessarily sensitive to the breadth "
        r"of that space. It must not be read as the percentage of source conclusions "
        r"that failed to replicate. The latter also depends on the chosen reference: "
        r"it measures departure from a common analysis, not error against a known truth."
        "\n\n"
        rf"Effect and uncertainty summaries reinforce this distinction. Direction "
        rf"changed for {pct(p['direction_changes']/p['coefficient_regions'])} of "
        rf"coefficient--region combinations, but a sign reversal near zero is not "
        r"evidence for a substantial opposite effect. Estimates and uncertainty "
        r"should therefore be inspected rather than interpreting every "
        r"sign change as a theoretical contradiction. Likewise, the median "
        rf"standardized SE ratio of {s['precision_ratio']['median']:.2f} shows "
        r"that precision can change even when a binary conclusion does not. "
        r"Because these are nonlinear scales with pipeline-specific denominators, "
        r"the analysis does not identify a single invariant millisecond effect "
        r"whose bias can be estimated without a generating model."
        "\n\n"
        rf"Interaction coefficients had a larger median standardized range "
        rf"({kind.loc['Interaction','range_median']:.3f}) than main-effect "
        rf"coefficients ({kind.loc['Main effect','range_median']:.3f}), whereas "
        rf"their median maximum-to-minimum SE ratios were almost identical "
        rf"({kind.loc['Interaction','SE_ratio_median']:.2f} and "
        rf"{kind.loc['Main effect','SE_ratio_median']:.2f}). "
        r"This distinction matters for claims about selective interference or "
        r"modulation: the estimated pattern of a factorial contrast may change "
        r"even when the relative change in its uncertainty resembles that of a "
        r"main effect. The comparison remains descriptive, because the two "
        r"coefficient sets differ in designs and source composition and because "
        r"standardization does not make nonlinear-scale interactions identical "
        r"scientific quantities."
        "\n\n"
        r"The task comparison provides no basis for declaring either SPR or Maze "
        r"inherently immune to preprocessing. Only seven publications contribute "
        r"Maze inputs, compared with 33 contributing SPR inputs, with overlap "
        r"between source sets. The uncertainty for Maze is consequently wide. "
        r"Post-critical regions are additionally selected on source-reported effects, "
        r"so their descriptive instability cannot be compared with critical regions "
        r"as if both were sampled by the same rule. Source composition, region "
        r"selection, and experimental design remain alternative explanations."
        "\n\n"
        r"Initial skewness was positively associated with disagreement, consistent "
        r"with longer tails making retention and scale decisions consequential. "
        r"This is not evidence that skewness alone determines validity: participant "
        r"counts, item counts, constructions, and task properties covary, and the "
        r"moderator analyses are not adjusted causal comparisons. In particular, "
        r"the absence of a precise item-count association is not evidence that "
        r"items do not matter for mixed-model uncertainty."
        "\n\n"
        rf"All {p['distinct']:,} distinct workflows ultimately yielded accepted "
        rf"inference, but only {p['accepted_first']:,} were accepted at the first "
        r"maximal attempt. Successful completion is therefore not equivalent to "
        r"support for an unchanged maximal covariance model. The reduction trace "
        r"is essential: dropping small variance components, including an intercept "
        r"when selected by the declared rule, can change uncertainty and can "
        r"remove a source of clustering. The same-formula comparison is useful "
        r"as a sensitivity description, but does not replace a preregistered "
        r"fixed-structure analysis or prove nominal Type I error control."
        "\n\n")
    (directory / "generated/discussion_text.tex").write_text(discussion)
    (directory / "generated/dataset_tables.tex").write_text(table(registry))
    shutil.copyfile(ANALYSIS / "workflow.tex", directory / "generated/workflow.tex")
    gs = shutil.which("gs") or str(next(Path("/nix/store").glob("*ghostscript*/bin/gs")))
    # Remove only obsolete generated figure copies, never source figures.
    for stem in ("effect_ranges", "subgroups"):
        for suffix in (".eps", ".pdf"):
            (directory / "figures" / (stem + suffix)).unlink(missing_ok=True)
    names = ["retention", "inference_stability", "decision_sensitivity", "lambda", "precision"]
    for name in names:
        src = OUT / "figures" / (name + ".eps")
        (directory / "figures" / src.name).write_bytes(src.read_bytes())
        subprocess.run([gs, "-q", "-dBATCH", "-dNOPAUSE", "-dEPSCrop", "-sDEVICE=pdfwrite",
                        "-sOutputFile=" + str(directory / "figures" / (name + ".pdf")), str(src)], check=True)
    extra = OUT / "effect_sensitivity"
    panel_values = pd.read_csv(extra / "figure_values.csv")
    if len(panel_values) != 468 or not (extra / "verification.txt").read_text().startswith("PASSED:"):
        raise RuntimeError("The preprocessing-effects figure has not passed verification")
    for suffix in (".eps", ".pdf"):
        (directory / "figures" / ("preprocessing_effects" + suffix)).write_bytes(
            (extra / ("preprocessing_effects" + suffix)).read_bytes())
    content = preamble + r"\begin{document}" + frontmatter + (ANALYSIS / "manuscript_body.tex").read_text() + "\n" + r"\end{document}" + "\n"
    (directory / "paper.tex").write_text(content)
    wanted = {m["citation"] for m in registry}
    for match in re.finditer(r"\\cite[pt]?\{([^}]+)\}", content):
        wanted.update(k.strip() for k in match.group(1).split(","))
    for name in ["sn-jnl.cls", "sn-apacite.bst"]:
        shutil.copyfile(TEMPLATE / name, directory / name)
    bib = ((TEMPLATE / "references.bib").read_text() + "\n"
           + (ANALYSIS / "methodology_refs.bib").read_text())
    (directory / "multiverse_refs.bib").write_text(extract_entries(bib, wanted))
    texbin = Path("/nix/store/wg7ib51s6xjplja05bxfhnd2qgfa34p1-texlive-combined-2024/bin")
    commands = [("pdflatex", "-interaction=nonstopmode", "-halt-on-error", "paper.tex"),
                ("bibtex", "paper"),
                ("pdflatex", "-interaction=nonstopmode", "-halt-on-error", "paper.tex"),
                ("pdflatex", "-interaction=nonstopmode", "-halt-on-error", "paper.tex")]
    logs = []
    for exe, *args in commands:
        result = subprocess.run([shutil.which(exe) or str(texbin / exe), *args], cwd=directory,
                                capture_output=True, text=True, errors="replace")
        logs.append(result.stdout + result.stderr)
        (directory / "compile.log").write_text("\n".join(logs))
        if result.returncode:
            raise RuntimeError("Compilation failed:\n" + logs[-1][-5000:])
    log = (directory / "paper.log").read_text(errors="replace")
    if re.search(r"Citation .* undefined|There were undefined references", log):
        raise RuntimeError("Unresolved references")
    for name in ["sn-jnl.cls", "sn-apacite.bst"]:
        assert (directory/name).read_bytes() == (TEMPLATE/name).read_bytes()
    package = OUT / "multiverse_overleaf.zip"
    with zipfile.ZipFile(package, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in ["paper.tex", "multiverse_refs.bib", "sn-jnl.cls", "sn-apacite.bst"]:
            archive.write(directory/name, name)
        for folder in ["generated", "figures"]:
            for path in sorted((directory/folder).rglob("*")):
                if path.is_file():
                    archive.write(path, str(path.relative_to(directory)))
    (OUT/"delivery.json").write_text(json.dumps({
        "paper_pdf":str((directory/"paper.pdf").relative_to(WORKSPACE)),
        "overleaf_zip":str(package.relative_to(WORKSPACE)),
        "datasets":len(registry),"region_jobs":audit["region_jobs"],
        "publisher_files_unchanged":True},indent=2))
    print("Compiled", directory/"paper.pdf")


if __name__ == "__main__":
    main()
