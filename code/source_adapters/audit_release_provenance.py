"""Reconstruct source adapters in an isolated output tree and verify all inputs.

The source data and fitted inputs are read-only. Generated CSVs are temporary.
"""
from pathlib import Path
import contextlib
from decimal import Decimal
import hashlib
import importlib
import json
import re
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd

WORKSPACE = Path(__file__).resolve().parents[2]
PROJECT = WORKSPACE / "BRM impact"
OUT = WORKSPACE / "generated-artifacts/brm_multiverse"
AUDIT = OUT / "release_audit"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def numeric_label(s):
    if re.fullmatch(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?", s):
        n = Decimal(s)
        return format(n.normalize(), "f") if n else "0"
    return s


def main():
    started = time.perf_counter()
    AUDIT.mkdir(parents=True, exist_ok=True)
    stage = AUDIT / "source_reconstruction"
    local_out = stage / "generated-artifacts/brm_multiverse"
    local_out.mkdir(parents=True, exist_ok=True)
    stage_project = stage / "BRM impact"
    stage_project.mkdir(exist_ok=True)
    link = stage_project / "preprocessing_robustness"
    if not link.exists():
        link.symlink_to(PROJECT / "preprocessing_robustness", target_is_directory=True)
    # Paper text is needed only to assert the source-backed post-critical rule.
    for path in OUT.glob("*.txt"):
        if "-" in path.name:
            (local_out/path.name).write_bytes(path.read_bytes())
    with (AUDIT / "source_reconstruction.log").open("w") as log:
        for name in ["prepare_boyce.R", "export_source_objects.R"]:
            # These two readers need base R only, not a new package-closure lookup.
            subprocess.run(["Rscript", str(PROJECT / "analysis" / name)],
                           cwd=stage, stdout=log, stderr=log, check=True)
    import prepare_data as base
    base.OUT = local_out
    build = importlib.import_module("rebuild_corpus")
    opened = set()
    data_root = (PROJECT / "preprocessing_robustness").resolve()
    def observe(event, args):
        if event == "open" and args and isinstance(args[0], (str, bytes)):
            path = Path(args[0].decode() if isinstance(args[0],bytes) else args[0]).resolve()
            if path.is_file() and path.is_relative_to(data_root):
                opened.add(path)
    sys.addaudithook(observe)
    with (AUDIT / "source_reconstruction.log").open("a") as log:
        with contextlib.redirect_stdout(log):
            build.existing()
            build.straightforward()
            build.archives_and_positioned()
            build.remaining_tables()
            build.aligned_and_serialized()
            build.r_objects()
            build.raw_logs_and_annotations()
            build.write_registry()
    # R binary-data readers run in a child process, outside Python's audit hook.
    opened.add(data_root / "f8qwh-paape_local_2025/Experiment 1 (SPR)/lc_spr_trialinfo2.Rda")
    boyce = data_root / "github-boyce_maze_2020/experiment"
    opened.add(boyce / "Labelled/rel_clause_labelled.rds")
    for name in ["spr","l_maze","g_maze","one_b_maze","gulo_maze"]:
        opened.add(boyce / "Data/Processed" / (name + ".rds"))

    registry = json.loads((OUT / "datasets.json").read_text())
    excluded = {x["dataset"] for x in json.loads((PROJECT / "analysis/excluded_datasets.json").read_text())}
    included = [m for m in registry if m["tier"] == "primary" and m["dataset"] not in excluded]
    reconstructed = {m["dataset"]: m for m in json.loads((local_out / "datasets.json").read_text())}
    inventory = pd.read_csv(PROJECT / "corpus_audit/eligible_analysis_inputs.csv").set_index("dataset")
    policies = pd.read_csv(PROJECT / "corpus_audit/paper_rt_policies.csv").set_index("dataset")
    regions = pd.read_csv(PROJECT / "corpus_audit/region_selection.csv")
    results = []
    region_rows = []
    dependencies = []
    file_hashes = {p: digest(p) for p in sorted(opened)}
    for m in included:
        r = reconstructed[m["dataset"]]
        for key in ["citation","task","tier","design","participants","items","observations",
                    "factors","region_counts","source","sha256"]:
            assert r[key] == m[key], (m["dataset"], key, r[key], m[key])
        for f in m["factors"]:
            assert {numeric_label(x) for x in r["levels"][f]} == {
                numeric_label(x) for x in m["levels"][f]}, (m["dataset"],f)
        cols = ["subject","item","region",*m["factors"]]
        a = pd.read_csv(local_out/r["data"],dtype=str,keep_default_na=False)
        b = pd.read_csv(AUDIT/"rds_exports"/(m["dataset"]+".csv"),dtype=str,keep_default_na=False)
        for frame in (a,b):
            for col in cols:
                frame[col] = frame[col].map(numeric_label)
            frame["RT"] = pd.to_numeric(frame.RT,errors="raise")
        a = a.sort_values(cols+["RT"]).reset_index(drop=True)
        b = b.sort_values(cols+["RT"]).reset_index(drop=True)
        assert a[cols].equals(b[cols]), (m["dataset"],"row/group/factor correspondence")
        assert np.allclose(a.RT,b.RT,rtol=5e-13,atol=1e-9), (m["dataset"],"RTs")
        max_diff = float(np.max(abs(a.RT-b.RT)))
        for ledger in (inventory,policies):
            assert ledger.loc[m["dataset"],"citation"] == m["citation"]
        for key in ["design","task","source","sha256","participants","items","observations","input_md5"]:
            assert str(inventory.loc[m["dataset"],key]) == str(m[key]), (m["dataset"],key)
        selected = regions[regions.dataset == m["dataset"]]
        assert set(selected.region) == set(m["region_counts"])
        for _, row in selected.iterrows():
            assert int(row.observations) == m["region_counts"][row.region]
            assert isinstance(row.selection_evidence,str) and row.selection_evidence.strip()
            region_rows.append(row.to_dict())
        source = PROJECT/m["source"]
        source_hash = file_hashes.get(source.resolve()) or digest(source)
        assert source_hash == m["sha256"]
        document = str(policies.loc[m["dataset"],"source_document"])
        documents = [p for p in (WORKSPACE/Path(document), PROJECT/Path(document)) if p.is_file()]
        # Existing evidence may name several PDF locations in a single field.
        document_exists = bool(documents) or any(
            pdf.name in document for parent in source.parents
            if parent.is_dir() and parent.parent == data_root
            for pdf in parent.rglob("*.pdf"))
        folder = source.relative_to(data_root).parts[0]
        family_files = [p for p in file_hashes if p.relative_to(data_root).parts[0] == folder]
        assert source.resolve() in file_hashes or m["dataset"].startswith(("boyce","paape"))
        for p in family_files:
            dependencies.append(dict(dataset=m["dataset"],source_path=str(p.relative_to(PROJECT)),
                sha256=file_hashes[p],bytes=p.stat().st_size,
                scope="Files read by this source-family's reconstruction; may include companion cohorts.",
                role="lead input" if p.resolve()==source.resolve() else "source-family dependency"))
        results.append(dict(dataset=m["dataset"],citation=m["citation"],design=m["design"],
            participants=m["participants"],items=m["items"],observations=len(a),
            reconstructed=True,factor_and_row_mapping_verified=True,
            maximum_RT_roundtrip_difference_ms=max_diff,lead_source_hash_verified=True,
            paper_document_present=document_exists,
            source_family_files_hashed=len(family_files)))
    pd.DataFrame(results).to_csv(AUDIT/"source_reconstruction_checks.csv",index=False)
    pd.DataFrame(dependencies).to_csv(AUDIT/"source_dependencies.csv",index=False)
    pd.DataFrame(region_rows).to_csv(AUDIT/"verified_region_selection.csv",index=False)
    report = dict(passed=True,inputs=len(included),observations=sum(m["observations"] for m in included),
        source_files_hashed=len(file_hashes),source_family_dependency_rows=len(dependencies),
        elapsed_seconds=time.perf_counter()-started,
        maximum_RT_roundtrip_difference_ms=max(x["maximum_RT_roundtrip_difference_ms"] for x in results),
        scope="Fresh reconstruction with the supplied source adapters, independent comparison to compact fitted inputs and audit ledgers.",
        limitations=[
            "Adapter reconstruction verifies implemented mappings, not an independent rereading of every source paper.",
            "Source-family dependencies are deliberately overinclusive, not a claim every listed file contributes rows to every input.",
            "Equivalent numeric identifier spellings are normalized only for comparison; fitted inputs are not changed.",
            "Source document availability does not establish redistribution permission."])
    (AUDIT/"provenance_verification.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
