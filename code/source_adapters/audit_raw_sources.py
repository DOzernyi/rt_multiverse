"""Read-only, archive-aware audit of every supplied study folder.

Reports schemas and aggregate RT statistics only: never publishes response rows
or participant identifiers and never extracts archives over original inputs.
"""
from pathlib import Path
import csv
import argparse
import glob
import hashlib
import io
import json
import re
import sys
import warnings
import zipfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "preprocessing_robustness"
OUT = ROOT.parent / "generated-artifacts/brm_multiverse/raw_source_audit"
RT_PATTERN = re.compile(
    r"(^rt$|(^|[._ ])rt([._ 0-9]|$)|^rt[0-9]+$|reading.?time|reaction.?time|"
    r"^d10$|^(noun|spill|verb).*rt$|^rt_or_correct$)", re.I)
CODE_PATTERN = re.compile(
    r"read[._](csv|table|excel|delim)|readxl|readRDS|load\(|"
    r"trim|winsori|outlier|raw.?data|log.?rt|log\([^)]*[Rr][Tt]|"
    r"filter\(|subset\(|[Rr][Tt]\s*[<>]|[Rr][Tt]\s*/\s*1000")
CODE_EXT = {".r", ".rmd", ".qmd", ".py", ".md", ".ipynb"}
TABLE_EXT = {".csv", ".tsv", ".txt"}


def rt_columns(columns, locator):
    names = [str(c) for c in columns]
    if "fp86z-" in locator and any(n in locator for n in ("DataDef.xlsx","DataIndef.xlsx")):
        return [c for c in names if c in ("Subject","Verb","Article","TriggerNoun","Preposition","Noun")]
    if "ye6dj-" in locator and locator.endswith("RT_data.tsv"):
        return [c for c in names if c in ("FFdur","RPdur","FPRT","RRT","TRT")]
    if "reliability_spr" in locator:
        return [c for c in names if re.match(r"^[AEP]_", c)]
    return [c for c in names if RT_PATTERN.search(c)]


def moments(values):
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return {"finite": 0}
    return {"finite": int(len(v)), "positive": int((v > 0).sum()),
            "zero": int((v == 0).sum()), "negative": int((v < 0).sum()),
            "minimum": float(v.min()), "maximum": float(v.max()),
            "below_100": int((v < 100).sum()),
            "above_10000": int((v > 10000).sum()),
            "sum": float(v.sum())}


def combine(a, b):
    if not b.get("finite"):
        return a
    if not a.get("finite"):
        return b
    return {k: (min(a[k], b[k]) if k == "minimum" else
                max(a[k], b[k]) if k == "maximum" else a[k]+b[k]) for k in a}


def numeric_values(series):
    direct = pd.to_numeric(series, errors="coerce")
    if direct.notna().sum() == series.notna().sum():
        return direct
    # Serialized trial-wise RT arrays are retained as arrays, not means.
    result = direct[direct.notna()].tolist()
    for value in series[direct.isna()].dropna():
        text = str(value).strip()
        if text.startswith(("[", "(")) and text.endswith(("]", ")")):
            result.extend(re.findall(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", text))
        elif "|" in text and re.fullmatch(r"[-+\d.eE|\s]+",text):
            result.extend(text.split("|"))
    return result


def profile_ibex(opener, locator, maze=False):
    stats, controllers, rows = {}, set(), 0
    reading_index = None
    with opener() as stream:
        def lines():
            nonlocal reading_index
            for line in io.TextIOWrapper(stream,encoding="utf-8-sig",errors="replace"):
                if line.startswith("#"):
                    match=re.match(r"#\s*(\d+)\.\s*Reading time[. ]",line,re.I)
                    if match:
                        reading_index=int(match.group(1))-1
                else:
                    yield line
        for row in csv.reader(lines()):
            if len(row)<3:
                continue
            rows += 1
            controller=row[2]
            controllers.add(controller)
            index=12 if maze else reading_index
            if index is not None and len(row)>index and (maze or "Dashed" in controller):
                m=moments([row[index]])
                if m.get("finite"):
                    stats["rt_column_"+str(index+1)]=combine(stats.get("rt_column_"+str(index+1),{}),m)
    return {"locator":locator,"format":"Ibex acquisition log","rows":rows,
            "controllers":sorted(controllers),"rt_columns":list(stats),
            "rt_statistics":stats,"note":("RT column 13 comes from the supplied cols definition."
              if maze else "Reading-time field identified from the embedded numbered comment headers.")}


def profile_csv(opener, locator):
    with opener() as stream:
        sample = stream.read(65536)
    text = sample.decode("utf-8-sig", errors="replace")
    if not text.strip():
        return {"locator": locator, "empty": True}
    try:
        separator = csv.Sniffer().sniff(text, delimiters=",;\t").delimiter
    except csv.Error:
        separator = "\t" if "\t" in text.splitlines()[0] else ","
    options = dict(sep=separator, encoding="utf-8-sig", encoding_errors="replace",
                   on_bad_lines="error", low_memory=False)
    try:
        with opener() as stream:
            head = pd.read_csv(stream, nrows=8, **options)
        fields = rt_columns(head.columns, locator)
        record = {"locator": locator, "format": "delimited", "separator": separator,
                  "columns": list(map(str, head.columns)), "rt_columns": fields}
        if not fields:
            return record
        stats = {}
        rows = 0
        selected = fields + [c for c in head.columns if re.fullmatch(
            r"condition|cond|type|region|tag|tag2|exp|experiment|controller|task|expt",
            str(c), flags=re.I) and c not in fields]
        levels = {c: set() for c in selected if c not in fields}
        with opener() as stream:
            for chunk in pd.read_csv(stream, usecols=selected, chunksize=25000, **options):
                rows += len(chunk)
                for field in fields:
                    stats[field] = combine(stats.get(field, {}), moments(numeric_values(chunk[field])))
                for field in levels:
                    if len(levels[field]) < 40:
                        levels[field].update(map(str, chunk[field].dropna().unique()[:40]))
        record.update(rows=rows, rt_statistics=stats,
                      code_levels={k: sorted(v)[:40] for k, v in levels.items()})
        return record
    except Exception as error:
        return {"locator": locator, "format": "delimited",
                "error": f"{type(error).__name__}: {error}"[:500],
                "rt_header_hint": bool(RT_PATTERN.search(text.splitlines()[0]))}


def openpyxl_module():
    try:
        import openpyxl
    except ImportError:
        for pattern in ["*openpyxl*/lib/python*/site-packages",
                        "*et*xmlfile*/lib/python*/site-packages"]:
            sys.path.extend(glob.glob("/nix/store/" + pattern))
        import openpyxl
    return openpyxl


def profile_workbook(path):
    path = path.resolve()
    openpyxl = openpyxl_module()
    result = []
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet in workbook:
        rows = sheet.iter_rows(values_only=True)
        initial = []
        for _ in range(12):
            try:
                initial.append(next(rows))
            except StopIteration:
                break
        if not initial:
            continue
        # Metadata can precede the real table header.
        header_index = max(range(len(initial)), key=lambda i:
                           len(rt_columns(initial[i], str(path))) * 1000 +
                           sum(isinstance(v, str) and bool(v.strip()) for v in initial[i]))
        header = [str(x) if x is not None else "" for x in initial[header_index]]
        fields = rt_columns(header, str(path))
        record = {"locator": str(path.relative_to(DATA))+"::sheet="+sheet.title,
                  "format": "xlsx", "header_row": header_index+1,
                  "columns": header, "rt_columns": fields}
        if fields:
            indices = [(i, field) for i, field in enumerate(header) if field in fields]
            stats, count, batch = {}, 0, []
            import itertools
            for row in itertools.chain(initial[header_index+1:], rows):
                if not any(x is not None for x in row):
                    continue
                count += 1
                batch.append(row)
                if len(batch) == 10000:
                    for i, field in indices:
                        values = [r[i] if i < len(r) else None for r in batch]
                        stats[field] = combine(stats.get(field, {}), moments(values))
                    batch = []
            for i, field in indices:
                values = [r[i] if i < len(r) else None for r in batch]
                stats[field] = combine(stats.get(field, {}), moments(values))
            record.update(rows=count, rt_statistics=stats)
        result.append(record)
    workbook.close()
    return result


def script_evidence(text, locator):
    return [{"locator": locator, "line": i, "text": line.strip()[:500]}
            for i, line in enumerate(text.splitlines(), 1)
            if CODE_PATTERN.search(line)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--retry-workbooks", action="store_true")
    parser.add_argument("--special-cases", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.special_cases:
        profiles=json.loads((OUT/"table_profiles.json").read_text())
        updates=[]
        for path in DATA.glob("7zvhg-*/*.csv"):
            updates.append(profile_ibex(lambda p=path:p.open("rb"),str(path.relative_to(DATA))))
        for path in DATA.glob("7zvhg-*/*.zip"):
            with zipfile.ZipFile(path) as z:
                for member in z.namelist():
                    updates.append(profile_ibex(lambda z=z,m=member:z.open(m),
                        str(path.relative_to(DATA))+"::"+member))
        for path in DATA.glob("naj9k-*/*/*Result.csv"):
            updates.append(profile_ibex(lambda p=path:p.open("rb"),str(path.relative_to(DATA)),maze=True))
        for path in DATA.glob("d34y5-*/experiment-raw-data/*.csv"):
            updates.append(profile_csv(lambda p=path:p.open("rb"),str(path.relative_to(DATA))))
        for path in DATA.glob("ye6dj-*/RT_data.tsv"):
            updates.append(profile_csv(lambda p=path:p.open("rb"),str(path.relative_to(DATA))))
        for record in updates:
            record["study_folder"]=record["locator"].split("/")[0]
        names={r["locator"] for r in updates}
        profiles=[r for r in profiles if r["locator"] not in names]+updates
        (OUT/"table_profiles.json").write_text(json.dumps(profiles,indent=2,ensure_ascii=False))
        print("DONE: inspected headerless logs, pipe-separated word times and eye-tracking fields")
        return
    if args.retry_workbooks:
        profiles = json.loads((OUT/"table_profiles.json").read_text())
        profiles = [x for x in profiles if x.get("format") != "xlsx"]
        for path in sorted(DATA.rglob("*.xlsx")):
            records = profile_workbook(path)
            for record in records:
                record["study_folder"] = path.relative_to(DATA).parts[0]
            profiles.extend(records)
        (OUT/"table_profiles.json").write_text(json.dumps(profiles,indent=2,ensure_ascii=False))
        print("DONE: inspected all workbooks")
        return
    profiles, scripts, archives, folders = [], [], [], []
    warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")
    for folder in sorted(p for p in DATA.iterdir() if p.is_dir()):
        print("AUDITING", folder.name, flush=True)
        folder_files = sorted(p for p in folder.rglob("*") if p.is_file())
        for path in folder_files:
            locator = str(path.relative_to(DATA))
            suffix = path.suffix.lower()
            if suffix in CODE_EXT:
                scripts.extend(script_evidence(path.read_text(errors="replace"), locator))
            elif suffix in TABLE_EXT:
                record = profile_csv(lambda p=path: p.open("rb"), locator)
                record["study_folder"] = folder.name
                profiles.append(record)
            elif suffix == ".xlsx":
                try:
                    records = profile_workbook(path)
                    for record in records:
                        record["study_folder"] = folder.name
                    profiles.extend(records)
                except Exception as error:
                    profiles.append({"study_folder": folder.name, "locator": locator,
                                     "format": "xlsx", "error": str(error)[:400]})
            elif suffix == ".zip":
                with zipfile.ZipFile(path) as archive:
                    members = [m for m in archive.infolist() if not m.is_dir()]
                    archives.append({"locator": locator, "study_folder": folder.name,
                                     "members": [{"name": m.filename, "bytes": m.file_size}
                                                 for m in members]})
                    for member in members:
                        if "__MACOSX" in member.filename or "/.git/" in member.filename:
                            continue
                        ext = Path(member.filename).suffix.lower()
                        inner = locator+"::"+member.filename
                        if ext in CODE_EXT and member.file_size < 2_000_000:
                            scripts.extend(script_evidence(archive.read(member).decode(
                                "utf-8", errors="replace"), inner))
                        elif ext in TABLE_EXT and member.file_size:
                            record = profile_csv(lambda n=member.filename, z=archive: z.open(n), inner)
                            record["study_folder"] = folder.name
                            profiles.append(record)
            elif not suffix and "Data/Raw/" in locator:
                profiles.append({"study_folder": folder.name, "locator": locator,
                                 "format": "extensionless acquisition log",
                                 "bytes": path.stat().st_size,
                                 "note": "Use the supplied original parser; not a generic CSV."})
        folders.append({"study_folder": folder.name, "files_checked": len(folder_files),
                        "profiled_tables": sum(p["study_folder"] == folder.name for p in profiles),
                        "script_evidence_lines": sum(s["locator"].startswith(folder.name+"/") for s in scripts)})
        # Incremental audit checkpoint; source files and model inputs are untouched.
        for name, data in [("table_profiles.json", profiles), ("script_evidence.json", scripts),
                           ("archive_inventory.json", archives), ("folder_coverage.json", folders)]:
            (OUT/name).write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print("DONE", len(folders), "folders", len(profiles), "table/sheet profiles",
          len(archives), "archives", len(scripts), "script evidence lines", flush=True)


if __name__ == "__main__":
    main()
