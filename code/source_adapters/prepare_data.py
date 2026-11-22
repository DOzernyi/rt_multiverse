"""Explicit, source-specific adapters. All generated data stay outside sources."""
from pathlib import Path
import hashlib
import json
import re
import csv
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "preprocessing_robustness"
OUT = ROOT.parent / "generated-artifacts/brm_multiverse"


def source(prefix, filename):
    folders = [p for p in DATA.iterdir() if p.is_dir() and p.name.startswith(prefix)]
    if len(folders) != 1:
        raise ValueError((prefix, folders))
    return folders[0] / filename


def read(path, **kwargs):
    if "sep" not in kwargs:
        with path.open(errors="replace") as stream:
            header = stream.readline()
        try:
            kwargs["sep"] = csv.Sniffer().sniff(header, delimiters=",;\t").delimiter
        except csv.Error:
            kwargs["sep"] = ","
    return pd.read_csv(path, encoding_errors="replace", **kwargs)


def save_dataset(frame, *, dataset, citation, label, task, tier, path,
                 subject, item, rt, factors, region, regions, position=None,
                 original_removal="", original_transform="", eligibility="",
                 source_script="", roi_evidence=None):
    """One observation per participant/item/ROI; no data-dependent ROI choice."""
    frame = frame.copy()
    frame[".position"] = pd.to_numeric(frame[position], errors="coerce") if position else 1
    def region_key(value):
        value = str(value).strip()
        return re.sub(r"\.0+$", "", value) if re.fullmatch(r"[+-]?\d+\.0+", value) else value
    regions = {region_key(k):v for k,v in regions.items()}
    frame[".region"] = frame[region].map(region_key)
    frame = frame[frame[".region"].isin(regions)].copy()
    frame["RT"] = pd.to_numeric(frame[rt], errors="coerce")
    names = []
    for n, column in enumerate(factors, 1):
        name = f"F{n}"
        frame[name] = frame[column].astype("string").str.strip()
        names.append(name)
    frame["subject"] = frame[subject].astype("string")
    frame["item"] = frame[item].astype("string")
    valid = np.isfinite(frame.RT) & (frame.RT > 0)
    for key in names + ["subject", "item"]:
        valid &= frame[key].notna() & ~frame[key].isin(["", "nan", "NA"])
    invalid_count = int((~valid).sum())
    frame = frame.loc[valid].copy()
    # A label such as "spillover" sometimes covers multiple words. Select the
    # first physical position within each source item/condition/ROI, not a
    # fastest/slowest response and not an average of differently timed words.
    if position:
        first = frame.groupby(["item"] + names + [".region"])[".position"].transform("min")
        frame = frame[frame[".position"].eq(first)].copy()
    frame["region"] = frame[".region"].map(regions)
    # A fully repeated design can legitimately present the same material in
    # multiple conditions to one person. Never discard these observations.
    repeated_cells = int(frame.duplicated(["subject", "item", "region"] + names).sum())
    levels = {n: sorted(frame[n].unique().tolist()) for n in names}
    if any(len(v) < 2 for v in levels.values()):
        raise ValueError(f"{dataset}: missing factor levels {levels}")
    expected = int(np.prod([len(v) for v in levels.values()]))
    if len(frame[names].drop_duplicates()) != expected:
        raise ValueError(f"{dataset}: incomplete factorial design")
    target = OUT / "canonical" / (dataset + ".csv")
    target.parent.mkdir(parents=True, exist_ok=True)
    frame[["subject", "item", "region", "RT"] + names].to_csv(target, index=False)
    # Keep only pseudonymous author-supplied IDs in canonical data.
    return dict(dataset=dataset, citation=citation, label=label, task=task,
                tier=tier, data=str(target.relative_to(OUT)),
                source=str(path.relative_to(ROOT)),
                sha256=hashlib.file_digest(path.open("rb"), "sha256").hexdigest(),
                factors=dict(zip(names, factors)), levels=levels,
                design="x".join(str(len(v)) for v in levels.values()),
                participants=int(frame.subject.nunique()), items=int(frame.item.nunique()),
                observations=len(frame), invalid_or_missing=invalid_count,
                repeated_subject_item_condition_cells=repeated_cells,
                region_counts=frame.groupby("region").size().to_dict(),
                roi_evidence=roi_evidence or {},
                skewness=float(frame.RT.skew()),
                original_removal=original_removal, original_transform=original_transform,
                eligibility=eligibility, source_script=source_script)


def main(write=True):
    OUT.mkdir(parents=True, exist_ok=True)
    records, excluded = [], []

    def add(frame, **kwargs):
        try:
            records.append(save_dataset(frame, **kwargs))
        except Exception as exc:
            excluded.append(dict(dataset=kwargs["dataset"], reason=str(exc),
                                 source=str(kwargs["path"].relative_to(ROOT))))

    # The author scripts specify exactly 80 participants in fixed row blocks.
    for exp, length, filename, critical in [(1, 192, "ex1_raw_data_rt.csv", 6),
                                           (2, 218, "ex2_raw_data.csv", 4)]:
        p = source("aw8t3", filename)
        d = read(p)
        if len(d) != 80 * length:
            raise ValueError("Raw Fujita row blocks disagree with supplied author script")
        d["subject_reconstructed"] = np.repeat(np.arange(1, 81), length)
        d["ambiguity"] = np.where(d.d6.isin(["gp1", "gp2"]), "ambiguous", "unambiguous")
        d["gender"] = np.where(d.d6.isin(["gp1", "gp3"]), "match", "mismatch")
        add(d, dataset=f"fujita2021_e{exp}", citation="fujita_parsing_2021",
            label=f"E{exp}", task="SPR", tier="primary", path=p,
            subject="subject_reconstructed", item="d7", rt="d10",
            factors=["ambiguity", "gender"], region="d8",
            regions={str(critical): "critical", str(critical+1): "postcritical"},
            original_removal="100--10000 ms", original_transform="natural log",
            eligibility="Raw row blocks reconstructed exactly as in author script; no RT trimming",
            source_script=f"ex{exp}_code_for_raw_data.R")

    # Tier decisions are overridden by the source-column audit in the complete
    # corpus builder; filenames alone cannot establish previous RT exclusions.
    families = [
        ("78vmz", "fujita_onlinex_2024", 2, "SPR", ["ambiguity", "gender"],
         "region_g", {"disambiguation":"critical", "post_disambiguation":"postcritical"}, "phrase", "rt"),
        ("nhc2p", "fujita_memory_2024", 2, "Maze", ["grammaticality", "distractor"],
         "region_g", {"critical":"critical", "spillover":"postcritical"}, "phrase", "rt"),
        ("pykua", "fujita_interference_2023", 4, "SPR", ["target", "distractor"],
         "region_g", {"critical":"critical", "spillover":"postcritical"}, "phrase", "rt"),
        ("qfnsg", "fujita_active_2026", 3, "SPR", None, "region",
         {"Disambiguating region":"critical", "Post-disambiguating region":"postcritical"}, "segment", "RT"),
        ("t42b6", "fujita_real-time_2026", 3, None, ["interpretation", "gender"],
         "region", {"Too region":"critical", "Spillover region":"postcritical"}, "segment", "RT")
    ]
    for prefix, cite, count, task, factors, reg, rois, pos, rt in families:
        for exp in range(1, count+1):
            filename = f"ex{exp}_df.csv" if prefix == "pykua" else f"Ex{exp}CleanData.csv"
            p = source(prefix, filename)
            fs = factors or ["ambiguity", ["control", "wh", "dependency"][exp-1]]
            add(read(p), dataset=f"{cite}_e{exp}", citation=cite, label=f"E{exp}",
                task=task or ("SPR" if exp == 3 else "Maze"), tier="secondary", path=p,
                subject="subject", item="item", rt=rt, factors=fs, region=reg,
                regions=rois, position=pos,
                original_removal="Author-prepared input; see source-script audit",
                original_transform="natural log",
                eligibility="Previously prepared input; irreversible source exclusions may remain",
                source_script="script_for_data_analysis.R" if prefix == "pykua" else "AnalysisCode.R")
    for filename, task in [("LmazeCleanData.csv", "Maze"), ("SPRCleanData.csv", "SPR")]:
        p = source("kftg7", filename)
        add(read(p), dataset=f"fujita2025_{task.lower()}", citation="fujita_l-maze_2025",
            label=task, task=task, tier="secondary", path=p,
            subject="subject", item="item", rt="RT", factors=["ambiguity", "gender"],
            region="region_g", regions={"Disambiguating region":"critical",
                                       "Post-disambiguating region":"postcritical"}, position="segment",
            original_removal="100--10000 ms after supplied cleaning", original_transform="natural log",
            eligibility="Author-cleaned input", source_script="AnalysisCode.R")
    p = source("9dp6e", "CleanData.csv")
    add(read(p), dataset="fujita2025_active", citation="fujita_active_2025", label="E1",
        task="SPR", tier="secondary", path=p, subject="subject", item="item", rt="RT",
        factors=["dependency","number"], region="region",
        regions={"Critical region":"critical","Spillover region":"postcritical"}, position="segment",
        original_removal="Author-prepared input; see source-script audit",
        original_transform="natural log", eligibility="Author-cleaned input", source_script="AnalysisCode.R")

    # The authors' archived analysis scripts define A/B by participant
    # cohorts, NOT by the subject-feature experimental manipulation.
    # In particular, both genders/numbers occur within each replication.
    for exp in range(1,6):
        p = source("figshare", f"experiment{exp}.csv")
        d = read(p)
        cohorts=[None] if exp in [1,3] else ["a","b"]
        for cohort in cohorts:
            v=d
            if cohort is not None:
                cutoff={2:128,4:112,5:200}[exp]
                in_a=d.SubjID.isin([f"S{i}" for i in range(1,cutoff+1)])
                v=d[in_a if cohort=="a" else ~in_a]
            v=v[v.ExpID.eq({1:"gender",2:"gender2",3:"ambiguity",
                            4:"fempire",5:"gennum"}[exp])]
            # E3 additionally manipulates attractor ambiguity; retain it.
            fs = ["Grammaticality","Match"] + (["Ambiguity"] if exp == 3 else [])
            if exp in [2,4,5]:
                fs+=[{2:"SubjGen",4:"SubjNum",5:"SubjPhi"}[exp]]
            add(v, dataset=f"tucker2021_e{exp}{cohort or ''}",
                citation="tucker_attraction_2021",
                label=f"E{exp}{cohort or ''}",
                task="SPR", tier="primary", path=p, subject="SubjID", item="ItemNum",
                rt="RT", factors=fs, region="Region", regions={"Verb":"critical","Verb+1":"postcritical"},
                original_removal="By-region 1% winsorization (article)", original_transform="untransformed descriptives",
                eligibility="Original RT column, before winsorization; experimental regions only",
                source_script=f"AnalysisScripts.zip: 0{exp}-Experiment{exp}-Analysis-Frontiers.R")

    p = source("czq8x", "Formatted_Data_Maze.csv")
    d = read(p)
    d = d[d.cond.isin(["BothTheEver","MatrixTheEver","MatrixNoEver",
                      "BothTheStill","MatrixTheStill","MatrixNoStill"])].copy()
    # Correct target selections only; missing post-error regions are structural.
    d = d[d.correct.astype(str).str.lower().isin(["yes","true","1"])]
    d["polarity"] = np.where(d.cond.str.endswith("Ever"), "NPI-ever", "PPI-still")
    d["negation"] = d.cond.str.replace(r"(Ever|Still)$", "", regex=True)
    add(d, dataset="orth2025_maze", citation="orth_positive_2025", label="E3", task="Maze",
        tier="primary", path=p, subject="subj", item="item", rt="rt",
        factors=["polarity","negation"], region="wordnumber",
        regions={"10":"critical","11":"postcritical"},
        original_removal="100--10000 ms; three unprepared trials removed (script)",
        original_transform="natural log", eligibility="Original RT before script's RT bounds; correct choices",
        source_script="Stats_Maze.R")

    # Saul source retains raw RT alongside modRT created only in its analysis.
    p = source("hdzsr", "DO_SPR.csv")
    d = read(p)
    d = d[d.trialType.eq("experimental-trial") & ~d.subject.astype(str).isin(["24","32"])]
    add(d, dataset="saul2025_do", citation="saul_interference_2025", label="DO", task="SPR",
        tier="primary", path=p, subject="subject", item="item", rt="RT",
        factors=["Distractor","Resumption"], region="wordN",
        regions={"8":"critical","9":"postcritical"},
        original_removal=">100 ms; participant mean+2.5 SD winsorization",
        original_transform="natural log",
        eligibility="Original trial RT; author exclusions of participants 24/32 for RT-based reasons held fixed",
        source_script="DO_SPR.qmd")

    for exp, filename, reg, rois in [
        (1, "Experiment 1/SPR-attraction.csv", "tag2", {"refl":"critical","spill":"postcritical"}),
        (2, "Experiment 2/SPR-grammaticality.csv", "tag", {"refl":"critical","spill":"postcritical"})]:
        p = source("qz93x", filename)
        d = read(p)
        if exp == 1:
            d = d[d.condition.isin([a+"-"+b for a in ["match","mismatch"] for b in ["gram","ungram","none"]])].copy()
            d["distractor"] = d.condition.str.split("-").str[0]
            d["verb_cue"] = d.condition.str.split("-").str[1]
            fs = ["verb_cue","distractor"]
        else:
            d = d[d.condition.str.startswith("SV",na=False)].copy()
            d["SV_agreement"] = d.condition.str.split("-").str[0]
            d["VR_agreement"] = d.condition.str.split("-").str[1]
            fs = ["SV_agreement","VR_agreement"]
        add(d, dataset=f"keshev2024_e{exp}", citation="keshev_representation_2024",
            label=f"E{exp}", task="SPR",tier="primary",path=p,
            subject="subject",item="item",rt="RT",factors=fs,region=reg,regions=rois,position="wordN",
            original_removal=("Adjusted trim column supplied; E1 bounds not specified in loaded script"
                              if exp == 1 else "120 ms minimum; participant mean+2.5 SD winsorization"),
            original_transform="natural log",eligibility="Original RT retained separately from adjusted trim",
            source_script=f"Experiment {exp}/E{exp}_analysis.Rmd")

    for edition in ["spr","l_maze","g_maze","one_b_maze","gulo_maze"]:
        export = OUT / "boyce" / (edition + ".csv")
        if not export.exists():
            continue
        p = source("github-boyce",f"experiment/Data/Processed/{edition}.rds")
        add(read(export), dataset=f"boyce2020_{edition}",citation="boyce_maze_2020",
            label={"spr":"SPR","l_maze":"L-maze","g_maze":"G-maze",
                   "one_b_maze":"Jozefowicz A-maze","gulo_maze":"Gulordava A-maze"}[edition],
            task="SPR" if edition=="spr" else "Maze",tier="primary",path=p,
            subject="subject",item="group",rt="rt",factors=["type"],region="mod_word_number",
            regions={"0":"critical","1":"postcritical"},
            original_removal="Positive, recorded RT; correct Maze selections; SPR accuracy >.8",
            original_transform="natural log",
            eligibility="Original positive RT, author-native eligibility, relative-clause source-label join",
            source_script="experiment/Code/results.Rmd")

    # Save a machine-verifiable registry rather than silently skipping adapters.
    if write:
        (OUT / "datasets.json").write_text(json.dumps(records, indent=2))
        (OUT / "adapter_exclusions.json").write_text(json.dumps(excluded, indent=2))
        pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v, dict)} for r in records]).to_csv(
            OUT / "dataset_registry.csv", index=False)
    print(f"Prepared {len(records)} datasets: " +
          str(pd.Series([r['tier'] for r in records]).value_counts().to_dict()))
    for r in records:
        print(r["dataset"], r["tier"], r["design"], r["participants"], r["observations"])
    for r in excluded:
        print("REQUIRES RESOLUTION", r)
    return records, excluded


if __name__ == "__main__":
    main()
