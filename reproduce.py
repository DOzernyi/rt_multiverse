"""Portable entry point for the clean OSF release (Python standard library)."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time


def copy_file(source,target):
    source=Path(source);target=Path(target)
    target.write_bytes(source.read_bytes())
    return str(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify every release-file hash; no R or TeX required.")
    parser.add_argument("--refit", action="store_true", help="Fit all models afresh instead of using checkpoints; potentially multiple days.")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers for --refit; default 1.")
    args = parser.parse_args()
    release = Path(__file__).resolve().parent
    with (release/"MANIFEST.csv").open(newline="") as stream:
        manifest = list(csv.DictReader(stream))
    for row in manifest:
        path = release/row["path"]
        if not path.is_file():
            raise SystemExit("Missing release file: "+row["path"])
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream,"sha256").hexdigest()
        if digest != row["sha256"]:
            raise SystemExit("Changed release file: "+row["path"])
    print(f"Verified {len(manifest)} release-file hashes.",flush=True)
    if args.check:
        return
    if args.workers<1:
        raise SystemExit("--workers must be positive")
    if not shutil.which("Rscript"):
        raise SystemExit("Rscript is required. See README.md for R packages and TeX/Ghostscript dependencies.")
    build_root = release/"_build"
    build_root.mkdir(exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="reproduce-",dir=build_root))
    project = stage/"BRM impact"
    analysis = project/"analysis"
    shutil.copytree(release/"code/analysis",analysis,copy_function=copy_file)
    shutil.copytree(release/"code/source_adapters",analysis,dirs_exist_ok=True,copy_function=copy_file)
    shutil.copytree(release/"provenance/source_audit",project/"corpus_audit",copy_function=copy_file)
    template = project/"_BRN_SN__multiverse_preprocessing"
    template.mkdir()
    copy_file(release/"code/template/paper_new.tex",template/"paper_new.tex")
    for name in ["sn-jnl.cls","sn-apacite.bst"]:
        copy_file(release/"manuscript"/name,template/name)
    copy_file(release/"manuscript/multiverse_refs.bib",template/"references.bib")
    out = stage/"generated-artifacts/brm_multiverse"
    out.mkdir(parents=True)
    shutil.copytree(release/"data/primary_inputs",out/"primary_inputs",copy_function=copy_file)
    copy_file(release/"data/datasets.json",out/"datasets.json")
    copy_file(release/"provenance/source_verification.json",out/"source_verification.json")
    if not args.refit:
        shutil.copytree(release/"checkpoints",out/"fits",copy_function=copy_file)
    print("Reproducing in:",stage,flush=True)
    if args.refit:
        print("WARNING: full refitting is enabled; the original run recorded 107.1 cumulative job-hours, not wall time.",flush=True)
    commands = []
    if args.refit:
        commands.append(["Rscript",str(analysis/"multiverse.R"),str(args.workers),"primary"])
    for name in ["summarize.R","verify_results.R","verify_inventory.R"]:
        commands.append(["Rscript",str(analysis/name)])
    commands += [
        [sys.executable,str(analysis/"effect_sensitivity.py")],
        ["Rscript",str(analysis/"plot_effect_sensitivity.R")],
        ["Rscript",str(analysis/"verify_effect_sensitivity.R")],
        [sys.executable,str(analysis/"write_paper.py"),"--archived-inputs"],
        [sys.executable,str(analysis/"verify_summary.py")],
    ]
    start = time.perf_counter()
    with (stage/"reproduction.log").open("w") as log:
        for command in commands:
            print("Running",Path(command[1]).name,flush=True)
            result = subprocess.run(command,cwd=stage,stdout=log,stderr=log)
            if result.returncode:
                raise SystemExit(f"Reproduction stopped with exit {result.returncode}. See {stage/'reproduction.log'}")
    duration = time.perf_counter()-start
    (stage/"run_summary.json").write_text(json.dumps(
        dict(elapsed_seconds=duration,refitted=args.refit,workers=args.workers,
             paper=str(out/"compiled/paper.pdf")),indent=2))
    print(f"Finished in {duration:.1f} seconds. PDF: {out/'compiled/paper.pdf'}",flush=True)


if __name__=="__main__":
    main()
