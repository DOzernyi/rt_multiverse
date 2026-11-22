"""Source-specific expansion of the reading-time corpus; no numerical fits.

Each adapter extracts original-scale RT, not the author's adjusted outcome.
Unresolved inputs are recorded explicitly rather than treated as null results.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import re
import zipfile
import numpy as np
import pandas as pd
import prepare_data as base
from audit_raw_sources import openpyxl_module

ROOT, DATA, OUT = base.ROOT, base.DATA, base.OUT
records, issues = [], []
save_original = base.save_dataset
REGION_PLAN=json.loads((ROOT/"analysis/source_region_plan.json").read_text())


def region_proof(dataset):
    plan=REGION_PLAN.get(dataset)
    if plan is None:
        return {}
    papers=sorted(OUT.glob(plan["paper_prefix"]+"*.txt"),
                  key=lambda p:p.stat().st_size,reverse=True)
    if not papers:
        raise ValueError(f"{dataset}: no original-paper text for region evidence")
    path=papers[0]
    lines=path.read_text(errors="replace").splitlines()
    ranges=plan["lines"]
    if plan.get("anchor_regex"):
        found=False
        for paper in papers:
            candidate=paper.read_text(errors="replace").splitlines()
            hits=[i for i,line in enumerate(candidate) if re.search(plan["anchor_regex"],line)]
            if hits:
                path=paper;lines=candidate
                ranges=[[max(1,hits[0]-12),hits[0]+5]]
                found=True;break
        if not found:
            raise ValueError(f"{dataset}: the source-paper statistic was not located")
    passages=[]
    for a,b in ranges:
        passages.append("\n".join(f"{i+1}: {lines[i]}" for i in range(a-1,min(b,len(lines)))))
    return {"source":path.name,"main_effect":plan["reason"],
            "original_numbered_passages":passages}

# These files demonstrably retain durations outside the cuts subsequently
# applied by the provided author scripts. This replaces filename-based tiers.
RAW_FAMILIES = {"fujita_memory_2024", "fujita_active_2026",
                "fujita_real-time_2026", "fujita_l-maze_2025", "fujita_active_2025"}


def workbook(path, sheet=0):
    openpyxl_module()
    return pd.read_excel(path, sheet_name=sheet)


def archived(path, member, **kwargs):
    with zipfile.ZipFile(path) as archive:
        return pd.read_csv(archive.open(member),encoding_errors="replace",**kwargs)


def citation(prefix):
    folder = base.source(prefix, "").name
    return "keshev_representation_2026" if prefix=="gdjne" else folder.split("-",1)[1]


def evidence(prefix, text):
    return {"critical": {"source":base.source(prefix,"").name,"definition":text}}


def add(d, *, prefix, filename, dataset, subject, item, rt, factors,
        region=None, regions=None, task="SPR", label=None, note="",
        tier="primary", roi_evidence=None, position=None, source_script=""):
    p=base.source(prefix,filename)
    plan=REGION_PLAN.get(dataset,{})
    if "wide_extra" in plan:
        parts=[]
        for column,roi in {rt:"critical",**plan["wide_extra"]}.items():
            part=d.copy();part[".raw_region_RT"]=part[column];part[".roi"]=roi;parts.append(part)
        d=pd.concat(parts,ignore_index=True)
        rt=".raw_region_RT";region=".roi"
        regions={r:r for r in ["critical",*plan["wide_extra"].values()]}
    if region is None:
        d=d.copy()
        d[".roi"]="critical"
        region=".roi"; regions={"critical":"critical"}
    if plan.get("extra_regions"):
        regions={**regions,**plan["extra_regions"]}
    if plan:
        roi_evidence=dict(roi_evidence or evidence(prefix,note))
        for roi in [*plan.get("extra_regions",{}).values(),*plan.get("wide_extra",{}).values()]:
            roi_evidence[roi]=region_proof(dataset)
    try:
        r=save_original(d,dataset=dataset,citation=citation(prefix),
            label=label or dataset,task=task,tier=tier,path=p,
            subject=subject,item=item,rt=rt,factors=factors,region=region,
            regions=regions,position=position,source_script=source_script,
            eligibility="Recorded original-scale RT; positive finite observations; "+note,
            original_removal="See original script and folder audit",
            original_transform="See original script and folder audit",
            roi_evidence=roi_evidence or evidence(prefix,note))
        expected=set(plan.get("extra_regions",{}).values())|set(plan.get("wide_extra",{}).values())
        if expected-set(r["region_counts"]):
            raise ValueError(f"Planned postcritical selectors had no observations: {expected-set(r['region_counts'])}")
        records.append(r)
        print("EXTRACTED",dataset,r["observations"],r["design"],flush=True)
    except Exception as exc:
        issues.append({"dataset":dataset,"folder":p.parent.name,"reason":str(exc)})
        print("UNRESOLVED",dataset,str(exc),flush=True)


def legacy_save(d, **kwargs):
    """Retain supported adapters, correcting source tiers and ROI errors."""
    cite=kwargs["citation"]
    if cite in RAW_FAMILIES:
        kwargs["tier"]="primary"
        kwargs["eligibility"]="Original RT column retained before author-script RT exclusions"
        if cite=="fujita_memory_2024":
            kwargs["source_script"]="AnalysisCodeBF.R; AnalysisCodeEffectEstimate.R"
    # Critical-only is the default. Postcritical inclusion needs a separate,
    # paper-supported selection, not an automatic second job for every source.
    plan=REGION_PLAN.get(kwargs["dataset"],{})
    allowed=set(plan.get("extra_regions",{}).values())
    kwargs["regions"]={k:v for k,v in kwargs["regions"].items()
                       if v=="critical" or v in allowed}
    kwargs["roi_evidence"]={"critical":{
        "source":kwargs.get("source_script",""),
        "definition":"Source-defined critical label/position: "+str(kwargs["regions"])}}
    for roi in plan.get("extra_regions",{}).values():
        kwargs["roi_evidence"][roi]=region_proof(kwargs["dataset"])
    if cite=="fujita_real-time_2026":
        d=d.copy()
        exp=int(kwargs["dataset"][-1])
        positions=[11,12] if exp==3 else [10,11]
        eligible=d[d.segment.isin(positions)].copy()
        eligible[".valid"]=np.isfinite(eligible.RT)&eligible.RT.gt(0)
        grouped=eligible.groupby(["subject","item"],dropna=False)
        total=grouped.RT.sum(min_count=2)
        complete=grouped[".valid"].sum().eq(2) & grouped.size().eq(2)
        total=total.where(complete)
        d=d[d.region.isin(["Did region",*plan.get("extra_regions",{})])].copy()
        d.loc[d.region.eq("Did region"),"RT"]=pd.MultiIndex.from_frame(
            d.loc[d.region.eq("Did region"),["subject","item"]]).map(total)
        kwargs["regions"]={"Did region":"critical",**plan.get("extra_regions",{})}
        kwargs["roi_evidence"]={"critical":{
            "source":"AnalysisCode.R: RTc construction and ellipsis-region models",
            "definition":f"Sum of original RT at segments {positions}; both components required.",
            "correction":"Earlier extraction incorrectly used Too region alone."}}
        for roi in plan.get("extra_regions",{}).values():
            kwargs["roi_evidence"][roi]=region_proof(kwargs["dataset"])
    result=save_original(d,**kwargs)
    # Existing adapters' critical observations/designs are unchanged apart
    # from the explicitly corrected ellipsis sum. This one-time migration
    # additionally checks source hash, factors, levels and n in the R runner.
    result["legacy_critical_unchanged"]=(cite not in [
        "fujita_real-time_2026","tucker_attraction_2021"])
    return result


def existing():
    base.save_dataset=legacy_save
    previous, failures=base.main(write=False)
    records.extend(previous)
    issues.extend(failures)
    base.save_dataset=save_original


def straightforward():
    p=base.source("396hw","Data.xlsx")
    add(workbook(p),prefix="396hw",filename="Data.xlsx",dataset="dlugosz2025_gender",
        subject="Participant",item="Item",rt="RT",factors=["Transparency","Congruency"],
        note="Target-noun RT supplied directly; proficiency not included in the lighter factorial model.",
        source_script="Analysis.R")
    for exp in (1,2):
        name=f"spr_dat_{'existing' if exp==1 else 'novel'}_compounds.csv"
        add(base.read(base.source("5r93v",name)),prefix="5r93v",filename=name,
            dataset=f"benjamin2023_{'existing' if exp==1 else 'novel'}",
            subject="subject_id",item="compound_id",rt="rt",factors=["context"],
            note="Target compound RT; focused supportive versus nonsupportive context model. "
                 "Continuous entropy/frequency covariates are not included in this alternative model.")
    name="Tan_Foltz_reading_times_region4.csv"
    add(base.read(base.source("4szpx",name)),prefix="4szpx",filename=name,
        dataset="tan2020_attachment",subject="sub",item="sent_id",rt="rt_raw",
        factors=["sent_type","quest_type","language"],
        note="Paper: disambiguating region 4; rt_raw, not rt_fin.")
    name="4 - main - data/readingtimes_88pp.csv"
    d=base.read(base.source("bdg7u",name))
    add(d,prefix="bdg7u",filename=name,dataset="slaats2026_agreement",
        subject="UserId",item="set_number",rt="RT",factors=["agreement","surprisal"],
        region="word_index",regions={"5.0":"critical","6.0":"postcritical_1",
                                     "7.0":"postcritical_2","8.0":"postcritical_3"},
        note="Paper Table 1: target is sixth word (zero-based word_index=5).",
        roi_evidence={
            "critical":{"source":"Paper, target-word method and Table 1","definition":"word_index=5"},
            "postcritical_1":{"source":"Paper, results","definition":"word_index=6",
                "main_effect":"Effect of agreement in all analysed regions, target and three spillovers."},
            "postcritical_2":{"source":"Paper, results","definition":"word_index=7",
                "main_effect":"Effect of agreement in all analysed regions."},
            "postcritical_3":{"source":"Paper, results","definition":"word_index=8",
                "main_effect":"Effect of agreement in all analysed regions."}})
    name="Data/SPR_data.csv"
    d=base.read(base.source("czvb9",name))
    parts=[]
    for col,roi in [("noun_RT_orig","critical"),("noun_spill1_RT_orig","postcritical_1"),
                    ("noun_spill2_RT_orig","postcritical_2")]:
        part=d.copy();part["raw_RT"]=part[col];part[".roi"]=roi;parts.append(part)
    add(pd.concat(parts,ignore_index=True),prefix="czvb9",filename=name,
        dataset="haeuser2024_age",subject="exp_subject_id",item="itemnr2",rt="raw_RT",
        factors=["condition","group"],region=".roi",
        regions={r:r for r in ["critical","postcritical_1","postcritical_2"]},
        note="Use *_orig RT backups, not winsorized noun_RT columns.",
        roi_evidence={
            "critical":{"source":"Paper, statistical analysis","definition":"Predictable/unpredictable noun."},
            "postcritical_1":{"source":"Paper, Noun+1 results","definition":"noun_spill1_RT_orig",
                "main_effect":"Main effect of predictability, p=.002."},
            "postcritical_2":{"source":"Paper, Noun+2 results","definition":"noun_spill2_RT_orig",
                "main_effect":"Significant effect of predictability, p<.001."}})
    for kind in ["Def","Indef"]:
        name=f"Data{kind}.xlsx"
        d=workbook(base.source("fp86z",name))
        add(d,prefix="fp86z",filename=name,dataset=f"dlugosz2025_{kind.lower()}",
            subject="Participant",item="Item",rt="TriggerNoun",
            factors=["Grammaticality","Group"],
            note="Critical trigger-noun column, not the word column Subject. "
                 "9999 is retained pending documented timeout interpretation; no ad-hoc RT cut applied.")
    name="Data and scripts for analysis/SPR_data.csv"
    d=base.read(base.source("j86m9",name));d=d[d.Item_type.eq("exp")]
    add(d,prefix="j86m9",filename=name,dataset="gastmann2026_spr",
        subject="Participant",item="Item",rt="RT_phrase3",
        factors=["Setting","Sent_type","Plausibility","Cog_stat"],
        note="Author SPR_analysis.Rmd models RT_phrase3 for experimental SRC/ORC items.",
        source_script="Data and scripts for analysis/SPR_analysis.Rmd")
    name="data_spr.csv";d=base.read(base.source("kcqt9",name))
    for kind,selector,fs,crit,post in [
        ("number",d.CataphorNumber.notna(),["CataphorNumber","NumberMatch"],6,7),
        ("gender",d.GenderMatch.notna(),["GenderMatch"],7,8)]:
        add(d[selector],prefix="kcqt9",filename=name,dataset=f"giskes2023_{kind}",
            subject="subj",item="item",rt="RT",factors=fs,region="region",
            regions={str(crit):"critical",str(post):"postcritical"},
            source_script="AbstractPrediction_plotsandmodels.R",
            note="Number and independent gender item sets; no filler pooling.",
            roi_evidence={
                "critical":{"source":"Paper and AbstractPrediction_plotsandmodels.R",
                            "definition":f"Region {crit}: "+("V2 verb" if kind=="number" else "name")},
                "postcritical":{"source":"Paper, respective results section",
                    "definition":f"Region {post}",
                    "main_effect":"Match effect Pr=.999 (number); gender mismatch effect Pr>.999 (gender)."}})
    name="AllData_OSF.xlsx"
    add(workbook(base.source("na6vy",name),"READING"),prefix="na6vy",filename=name,
        dataset="haigh2022_context",subject="participant",item="TargetWord",rt="RT",
        factors=["SentenceContext","TaskDemands","AgeGroup"],
        note="READING sheet: target-word RT only; repeated conditions retained.")
    name="data_SPR.csv";d=base.read(base.source("rnbjm",name))
    d["Type"]=d.Type.str.strip()
    for kind in ["Ability","Epistemic","Permission"]:
        add(d[d.Type.eq(kind)],prefix="rnbjm",filename=name,
            dataset=f"mifka2024_{kind.lower()}",subject="Participant",item="Sentence",
            rt="RT",factors=["Congruency","Group"],region="Segment",regions={"1":"critical"},
            note="Modal is Segment 1, per paper's segment definitions. "
                 "Modality-specific subdesigns share participants; not independent cohorts.")
    name="CleanData.csv";d=base.read(base.source("rh4xz",name));parts=[]
    for field,value,roi in [("regionA","reflexive","critical_reflexive"),
                            ("regionB","disambiguation","critical_disambiguation")]:
        part=d[d[field].eq(value)].copy();part[".roi"]=roi;parts.append(part)
    add(pd.concat(parts),prefix="rh4xz",filename=name,dataset="fujita2023_predictive",
        subject="subject",item="item",rt="RT",factors=["Wh","Gender"],task="Maze",
        region=".roi",regions={"critical_reflexive":"critical_reflexive",
                               "critical_disambiguation":"critical_disambiguation"},
        note="Two independently defined critical sites in regionA and regionB; no post-region auto-addition.")
    name="CleanDataAdditionalExperiment.csv"
    add(base.read(base.source("9dp6e",name)),prefix="9dp6e",filename=name,
        dataset="fujita2025_active_additional",subject="subject",item="item",rt="RT",
        factors=["dependency"],region="region",regions={"Critical region":"critical"},
        note="Additional experiment's dependency contrast, not a duplicate of the main experiment.")
    for variant in ["AVCSPR","NAVCSPR","NCSPR","PCSPR"]:
        for paradigm,field,reg in [("GP","ambiguity","Disambiguating region"),
                                   ("SV","grammaticality","Verb region")]:
            name=f"ReadingTime/{variant}_{paradigm}_RT.csv"
            add(base.read(base.source("p2rck",name)),prefix="p2rck",filename=name,
                dataset=f"fujita2026_reassessment_{variant.lower()}_{paradigm.lower()}",
                subject="subject",item="item",rt="RT",factors=[field],region="region",
                regions={reg:"critical"},
                note="Separate GP/SV paradigms and task implementations; cohorts may overlap.")


def archives_and_positioned():
    name="datasets/rawData/Archive.zip";p=base.source("rsdp5",name)
    for exp,critical in [(1,7),(2,6),(3,8)]:
        d=archived(p,f"exp{exp}_raw.csv")
        d=d[d.Clause.isin(["appos","restr"])].copy()
        # Original functions.R::assign_region specifies different critical
        # positions in E1/E2/E3. Applying E1's position to all three is wrong.
        add(d,prefix="rsdp5",filename=name,dataset=f"kim2024_e{exp}",
            subject="Subj",item="ItemNo",rt="RT",
            factors=["Clause","Distractor","Grammaticality"],region="WordNo",
            regions={str(critical):"critical"},source_script="functions.R::assign_region",
            note=f"Archive member exp{exp}_raw.csv; critical WordNo={critical}.")
    for exp in (1,2,3):
        name=f"01-data/tidy/exp{exp}_data_tidy.csv.zip"
        d=archived(base.source("n38tk",name),f"exp{exp}_data_tidy.csv")
        d=d[d.type.eq("critical")].copy()
        add(d,prefix="n38tk",filename=name,dataset=f"palleschi2026_e{exp}",
            subject="px",item="item",rt="raw_rt",factors=["lifetime","tense"],
            region="region",regions={"verb":"critical"},
            note="Raw word duration from individual experiment archive; not transformed or pooled export.")
    prior_name="A1_E1: Semantic Richness, Engagement, & Aging/Karimi et al_2023_RT Data.csv"
    prior=base.read(base.source("tcz5n",prior_name))
    add(prior,prefix="tcz5n",filename=prior_name,
        dataset="karimi2025_prior_semantic_data",label="Archived Karimi et al. (2023) data",
        subject="sub",item="item",rt="rt",factors=["condition"],
        region="word_ID_final",regions={"0":"critical"},
        source_script="A1_E1: Semantic Richness, Engagement, & Aging/Re-analysis of Karimi et al_2023_R1_OSF.Rmd",
        note="Separate archived 2023 experiment reanalysed in the 2025 article, "
             "not the newly collected age-group dataset; original milliseconds, critical noun at relative position zero.")
    for exp,subdir,agecol in [
        (1,"A1_E1: Semantic Richness, Engagement, & Aging/spr_data_OSF.csv","Age_Group"),
        (3,"A2_E1_Modifier Content Retrieval Interference and Aging/Statistical Analyses/_R15_Aim2_Exp1_spr_data_OSF_Final.csv","AgeGroup")]:
        d=base.read(base.source("tcz5n",subdir))
        # WordNum is total sentence length, NOT the current word position.
        # The source script extracts the zero-based suffix of w_0, w_1, ...
        word_position=pd.to_numeric(d.word.str.extract(r"(\d+)$",expand=False),errors="coerce")+1
        d=d[word_position.eq(pd.to_numeric(d.TargetPos,errors="coerce"))].copy()
        d["RT_ms"]=pd.to_numeric(d.RT,errors="coerce")*1000
        add(d,prefix="tcz5n",filename=subdir,dataset=f"karimi2025_e{exp}",
            subject="participant",item="item",rt="RT_ms",factors=["Condition",agecol],
            note="TargetPos locates the critical target noun; PsychoPy seconds explicitly converted to ms.")
    name="SPR_data.csv";d=base.read(base.source("qudny",name))
    d=d[d.Item_type.astype(str).str.lower().isin(["exp","experimental"])].copy()
    d["RT_target"]=[row.get(f"RT_phrase{int(row.Max_phrase)}",np.nan)
        if pd.notna(row.Max_phrase) else np.nan for _,row in d.iterrows()]
    add(d,prefix="qudny",filename=name,dataset="brockmeyer2026_prediction",
        subject="Participant",item="item",rt="RT_target",factors=["Condition"],
        source_script="SPR_analysis.Rmd: target selected by Max_phrase",
        note="Variable sentence-final target selected using Max_phrase, not a fixed phrase column.")


def remaining_tables():
    d=base.read(base.source("t9nf3","spr_data_long.csv"))
    add(d[d.word.eq("noun")],prefix="t9nf3",filename="spr_data_long.csv",
        dataset="haeuser2022_prediction",subject="exp_subject_id",item="itemnr2",
        rt="RT",factors=["condition"],tier="secondary",
        note="Prepared long table: source script sets RTs <100 or >3000 ms to NA, "
             "without a preserved original-RT column. Critical noun recovered; "
             "not eligible for primary raw-RT multiverse or primary RDS export.")
    # Raw RT joined to categorical stimulus annotations only: no prepared
    # outcome or RT-based row selection is imported from the lookup.
    lookup=base.read(base.source("4fpru","pandora_spr_774_crit.csv"))
    lookup=lookup[["condition","syn","sem"]].drop_duplicates()
    if lookup.condition.duplicated().any():
        raise ValueError("Schoknecht condition codes do not map uniquely")
    for exp in (1,2):
        name=f"pandora_repS{exp}_data_cleaned_short.csv"
        d=base.read(base.source("4fpru",name))
        d=d[d.label.eq("experimental")].merge(lookup,on="condition",how="left",validate="many_to_one")
        add(d,prefix="4fpru",filename=name,dataset=f"schoknecht2025_rep_s{exp}",
            subject="fullname",item="item",rt="rt",factors=["syn","sem"],region="region",
            regions={"critical":"critical"},
            note="Replication source RT; pooled export used ONLY to decode a/b/c/d stimulus factors.")
    name="data/loco_SPR_rawdata.csv"
    d=base.read(base.source("m7dqc",name));d=d[d.cond.isin(["loco","noco"])]
    add(d,prefix="m7dqc",filename=name,dataset="schoknecht2026_local_coherence",
        subject="subject",item="item",rt="rt",factors=["cond"],region="region",
        regions={"verb":"critical"},source_script="prep_dfs_for_models.R",
        note="Experimental local/nonlocal-coherence items; verb is the disambiguation site.")
    for edition,name in [("a","saluki-spr-4afc.csv"),("b","saluki-spr-yn.csv")]:
        d=base.read(base.source("gdjne",name))
        d=d[d.exp.eq("experimental-trial")].copy()
        pos=pd.to_numeric(d.wordN,errors="coerce")
        item=pd.to_numeric(d.item,errors="coerce")
        adjust=(item.isin([6,10])&pos.gt(2))|(item.eq(23)&pos.gt(7))|(item.eq(32)&pos.gt(5))
        d["aligned_region"]=pos-adjust.astype(int)
        fs=["cue","agreement"] if edition=="a" else ["Cue","Agreement"]
        add(d,prefix="gdjne",filename=name,dataset=f"keshev2026_e3{edition}",
            subject="subject",item="item",rt="RT",factors=fs,region="aligned_region",
            regions={"9":"critical"} if edition=="a" else {"9.0":"critical"},
            source_script=f"saluki-spr-{edition}.Rmd: region tagging",
            note="Auxiliary at aligned region i=9. Item-specific source-script offsets for 6,10,23,32.")
    name="DataForOSF.zip";parts=[]
    for member,age in [("DO.csv","older"),("DY.csv","younger")]:
        d=archived(base.source("cwyk6",name),member)
        d=d[d["Zone.Type"].eq("reading") & d["Screen.Name"].eq("Screen 2")
            & d.Experiment.eq("EitherOr")].copy()
        d["AgeGroup"]=age
        d["participant"]=age+":"+d["Participant.Private.ID"].astype(str)
        parts.append(d)
    add(pd.concat(parts,ignore_index=True),prefix="cwyk6",filename=name,
        dataset="cutter2023_either",subject="participant",item="Item",
        rt="Reaction.Time",factors=["Condition","AgeGroup"],region="Attempt",
        regions={"3.0":"critical"},source_script="EitherAnalysisScriptOSF.R",
        note="Author script: reading-zone Screen 2, EitherOr experiment, Attempt 3 target. No author RT cuts.")
    name="Statistical analysis/E2_dataset_master.xlsx"
    d=workbook(base.source("he6cf",name),"SPR.long")
    add(d,prefix="he6cf",filename=name,dataset="ishihara2025_e2",
        subject="PID",item="Item",rt="RT",factors=["Group","Grammaticality"],
        note="Trial-level SPR.long outcome; School/Home administration groups. "
             "CR is correctness, NOT a region index; SOE is not the original RT.")
    for pronoun,name in [("null","SPR_Null_Critical_df.csv"),("overt","SPR_Overt_Critical_df.csv")]:
        d=base.read(base.source("gxpwd",name))
        add(d,prefix="gxpwd",filename=name,dataset=f"cairncross2026_{pronoun}",
            subject="Participant",item="Item",rt="ReadingTime",factors=["Antecedent"],
            tier="secondary",note="Available critical-window RT; source paper documents prior RT exclusions. "
            "This prepared source is not exported as an untrimmed primary input.")
    name="DATA/SPR_data.csv"
    add(base.read(base.source("e56k4",name)),prefix="e56k4",filename=name,
        dataset="haeuser2025_reduced",subject="exp_subject_id",item="itemnr2",
        rt="noun_RT",factors=["condition","group"],tier="secondary",
        note="Available noun RT; no *_orig backup for restoring potential earlier critical-region screening.")
    # Orth's acquisition CSVs have numbered comment headers and no CSV header.
    cols=["time","subj","controller","number","element","cond","item","wordnumber",
          "option1","option2","side","correct","rt","fulltext","withredo"]
    for edition in ["L","A"]:
        name=f"{edition}_Maze/{edition}_Maze_Result.csv"
        rows=[]
        with base.source("naj9k",name).open(errors="replace") as f:
            for row in csv.reader(line for line in f if not line.startswith("#")):
                if len(row)>=15:
                    rows.append(row[:15])
        d=pd.DataFrame(rows,columns=cols)
        d=d[d.cond.isin(["Expt-A","Expt-B","Expt-C","Expt-D"])
            & d.correct.eq("yes")].copy()
        d["RC_type"]=np.where(d.cond.isin(["Expt-A","Expt-B"]),"SRC","ORC")
        d["locality"]=np.where(d.cond.isin(["Expt-A","Expt-C"]),"local","nonlocal")
        add(d,prefix="naj9k",filename=name,dataset=f"orth2026_{edition.lower()}maze",
            subject="subj",item="item",rt="rt",factors=["RC_type","locality"],task="Maze",
            region="wordnumber",regions={"2":"critical"},
            source_script=f"{edition}_Maze/Stats_{edition}_Maze_Updated.R",
            note="Relative-pronoun site wordnumber=2, correct first choice only. "
                 "The separate Appendix archive is not silently added as a new cohort.")


def aligned_and_serialized():
    # Reproduce the authors' explicit nil-condition exceptions, without any
    # RT-dependent selection or importing their already trimmed region tables.
    for exp in [2,4]:
        name=f"experiment {exp}/results/results_all.csv"
        d=base.read(base.source("2xcgu",name))
        d=d[d.controller.eq("DashedSentence")&d.item_type.eq("target")].copy()
        d["pos"]=pd.to_numeric(d.quest_or_wordnum,errors="coerce")
        anchors=d[d.ans_or_word.isin(["too","tooc"])].groupby(["item_no","dist"]).pos.first()
        d["anchor"]=pd.MultiIndex.from_frame(d[["item_no","dist"]]).map(anchors)
        script=(f"experiment {exp}/src/"+("03-regions.R" if exp==2 else "02-regions.R"))
        text=base.source("2xcgu",script).read_text()
        for a,b,pos in re.findall(
                r'dist\s*==\s*"([^"]+)"\s*&\s*item\s*==\s*(\d+)\s*~\s*(\d+)',text):
            dist=pd.to_numeric(d.dist,errors="coerce") if exp==2 else d.dist
            val=float(a) if exp==2 else a
            d.loc[d.trigger.eq("nil")&d.item_no.eq(int(b))&dist.eq(val),"anchor"]=int(pos)
        for item,dist,pos in re.findall(
                r'item\s*==\s*(\d+)\s*&\s*dist\s*==\s*"([^"]+)"\s*~\s*(\d+)',text):
            d.loc[d.trigger.eq("nil")&d.item_no.eq(int(item))&d.dist.eq(dist),"anchor"]=int(pos)
        d["aligned_region"]=d.pos-d.anchor
        add(d,prefix="2xcgu",filename=name,dataset=f"winkowski2024_e{exp}",
            subject="subj",item="item_no",rt="rt_or_correct",factors=["trigger","dist"],
            region="aligned_region",regions={"0.0":"critical","1.0":"postcritical"},
            source_script=script,note="Too versus corresponding nil position; exact author alignment exceptions.",
            roi_evidence={
                "critical":{"source":script,"definition":"Author-relative region 0."},
                "postcritical":{"source":"Paper, E2/E4 reading-time results",
                    "definition":"Author-relative region 1.",
                    "main_effect":"Trigger effect negative credible interval at postcritical region."}})
    name="experiment-raw-data/emoji-homophones-data.csv"
    d=base.read(base.source("d34y5",name))
    groups=[
        ['nachtisch','beleuchtung','stadt','waldtier','alkohol','feld'],
        ['grad','prinzessin','blatt','nagetier','bank','tränen','ausscheidung'],
        ['tür','computer','akrobatik','sitzmöbel'],
        ['stall','kuh','längenumrechnung','meeressäuger'],
        ['klempner','obstsalat','demo','aktien','politik','schach'],
        ['papier','wettrennen','schmuggel']]
    codes={"A":["match","match","word","homophone","word","homophone"],
           "B":["homophone","word","match","word","homophone","match"],
           "C":["word","homophone","homophone","match","match","word"]}
    d["condition"]=pd.NA
    for group,labels in codes.items():
        for words,label in zip(groups,labels):
            which=d.participant_group.eq(group)&d.itemname.str.contains("|".join(words),na=False)
            d.loc[which,"condition"]=label
    d=d[d.condition.notna()].copy()
    def critical(value):
        try:
            values=[float(v) for v in str(value).split("|")]
            return values[-4] if len(values)>=4 else np.nan
        except ValueError:
            return np.nan
    d["RT_critical"]=d.reaction_times.map(critical)
    add(d,prefix="d34y5",filename=name,dataset="scheffler2022_emoji",
        subject="participant_id",item="itemname",rt="RT_critical",factors=["condition"],
        note="Original notebook get_critical returns the fourth-from-last RT; "
             "condition mapping copied from its six item lists and A/B/C assignment.",
        source_script="analysis-code/emoji-homophones-lmm.ipynb")


def r_objects():
    name="Experiment 1 (SPR)/lc_spr_trialinfo2.Rda"
    export=OUT/"source_exports/paape_spr2.csv"
    if export.exists():
        d=base.read(export)
        d=d[pd.to_numeric(d.no,errors="coerce").eq(pd.to_numeric(d.critno,errors="coerce"))]
        add(d,prefix="f8qwh",filename=name,dataset="paape2025_spr_semantic",
            subject="subj",item="item",rt="rt",factors=["coh","red"],
            note="Second source object only, untrimmed rt. no==critno selects the verb. "
                 "First source object is combined in author script but has prior RT bounds.",
            source_script="Experiment 1 (SPR)/analyses-spr-new.Rmd")
    else:
        issues.append({"dataset":"paape2025_spr_semantic",
                       "reason":"Run export_source_objects.R before rebuilding the corpus."})


def raw_logs_and_annotations():
    name="Data/data(200)_final.txt.zip"
    d=archived(base.source("3tjwv",name),"data(200)_final.txt",sep="\t")
    d=d[d.Type.eq("RCIS")&d["Controller.name"].eq("DashedSentence")].copy()
    lookup=archived(base.source("3tjwv","Data/data(200)res_final.txt.zip"),
                    "data(200)res_final.txt",sep="\t")
    factors=lookup[["Item.number","item","clause.order","information.order"]].drop_duplicates()
    if factors["Item.number"].duplicated().any():
        raise ValueError("Wilke stimulus-factor lookup is not unique")
    cl=archived(base.source("3tjwv","Data/clause.data(200)res_final.txt.zip"),
                "clause.data(200)res_final.txt",sep="\t")
    for frame in (d,cl):
        frame["Word.number"]=pd.to_numeric(frame["Word.number"],errors="coerce")
    loc=cl[["Item.number","Word.number","crit"]].drop_duplicates()
    if loc.duplicated(["Item.number","Word.number"]).any():
        raise ValueError("Wilke clause-region lookup is not unique")
    d=d.merge(factors,on="Item.number",how="left",validate="many_to_one").merge(
        loc,on=["Item.number","Word.number"],how="left",validate="many_to_one")
    add(d,prefix="3tjwv",filename=name,dataset="wilke2025_appositives",
        subject="subj",item="item",rt="rt",factors=["clause.order","information.order"],
        region="crit",regions={"ARC":"critical_arc","matrix":"critical_matrix"},
        note="Untrimmed RT in author-retained 200-participant cohort. Prepared files "
             "supply only categorical item/word-position annotations, never RT or row eligibility.",
        source_script="Data/ISclauses_res.rtData.R")
    parts=[]
    folder=base.source("7zvhg","")
    def collect(stream):
        collected=[]
        for row in csv.reader(line for line in stream if not line.startswith("#")):
            if len(row)!=33 or row[7]!="Controller-DashedSentence":
                continue
            if row[20] not in ["High-cloze","Low-cloze"]:
                continue
            try:
                pos=float(row[9]); target=float(row[23]); rt=float(row[29])
            except ValueError:
                continue
            if pos==target:
                collected.append({"subject":row[1],"item":row[19],
                    "Condition":row[20],"Task":row[18],"RT":rt,
                    "recorded_word":row[10],"target_word":row[22]})
        return pd.DataFrame(collected)
    for path in sorted(folder.glob("results_list*.csv")):
        with path.open(errors="replace") as stream:
            parts.append(collect(stream))
    for path in sorted(folder.glob("results_list*.zip")):
        with zipfile.ZipFile(path) as z:
            for member in z.namelist():
                if member.endswith(".csv") and not member.startswith("__"):
                    with z.open(member) as stream:
                        parts.append(collect(io.TextIOWrapper(stream,encoding="utf-8",errors="replace")))
    if parts:
        d=pd.concat(parts,ignore_index=True)
        fs=["Condition"]+(["Task"] if d.Task.nunique()>1 else [])
        # Stimulus identity check is independent of RT: a wrong index convention
        # must be detected before including any model.
        normalise=lambda s:s.astype(str).str.casefold().str.replace(r"\W","",regex=True)
        matches=normalise(d.recorded_word).eq(normalise(d.target_word))
        if len(d)==0 or matches.mean()<.99:
            issues.append({"dataset":"ito2026_predictability",
                "reason":f"Target_position identity verification failed: {matches.mean()} matches."})
        else:
            add(d,prefix="7zvhg",filename="results_listA1.csv",dataset="ito2026_predictability",
                subject="subject",item="item",rt="RT",factors=fs,
                note="All four list logs, PennController-DashedSentence rows; "
                     "column 30 reading RT at Target_position; verified against actual target word. "
                     "The source hash records A1; companion list hashes are in the input manifest.")


def write_registry():
    ids=[r["dataset"] for r in records]
    if len(ids)!=len(set(ids)):
        raise ValueError("Duplicate dataset identifiers")
    (OUT/"datasets.json").write_text(json.dumps(records,indent=2))
    (OUT/"adapter_exclusions.json").write_text(json.dumps(issues,indent=2))
    pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,dict)}
                  for r in records]).to_csv(OUT/"dataset_registry.csv",index=False)
    print("TOTAL",len(records),"PRIMARY",sum(r["tier"]=="primary" for r in records),
          "ISSUES",len(issues),flush=True)


if __name__=="__main__":
    existing()
    straightforward()
    archives_and_positioned()
    remaining_tables()
    aligned_and_serialized()
    r_objects()
    raw_logs_and_annotations()
    write_registry()
