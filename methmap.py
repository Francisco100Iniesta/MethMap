#!/usr/bin/env python3
"""Portable entry point for the MethMap tissue reference pipeline (stdlib only)."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

VERSION = "0.1.0"
ROOT = Path(__file__).resolve().parent
CANONICAL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
METRICS = [f"p{i}" for i in range(5, 96, 5)] + ["hits", "coverage_pct", "score", "iqr"]
SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def fail(message):
    raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def run(command, log=None, env=None):
    command = [str(x) for x in command]
    print("+ " + shlex.join(command), flush=True)
    if log:
        with Path(log).open("a") as handle:
            handle.write("+ " + shlex.join(command) + "\n")
            handle.flush()
            subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT, env=env)
    else:
        subprocess.run(command, check=True, env=env)


def load_samples(manifest, tissue, input_format):
    samples = []
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["sample_id", "tissue", "path"]:
            fail("Sample TSV must have exactly these columns: sample_id, tissue, path")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                fail(f"Malformed manifest row {reader.line_num}")
            if row["tissue"] != tissue:
                continue
            sid = row["sample_id"]
            if not SAFE_NAME.fullmatch(sid):
                fail(f"Unsafe sample ID: {sid!r}")
            path = Path(row["path"]).expanduser()
            if not path.is_absolute():
                path = manifest.parent / path
            path = path.resolve()
            if not path.is_file():
                fail(f"Input not found for {sid}: {path}")
            samples.append(dict(sample_id=sid, tissue=tissue, path=str(path),
                                input_format=input_format, bytes=path.stat().st_size))
    if not samples:
        fail(f"No samples for tissue {tissue!r}")
    for key in ["sample_id", "path"]:
        if len({x[key] for x in samples}) != len(samples):
            fail(f"Duplicate {key} in selected tissue")
    return samples


def read_fai(path, chromosomes):
    sizes = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                fail("Expected a FASTA .fai index with at least five columns")
            name, length = fields[0], int(fields[1])
            if name in sizes or length < 1:
                fail(f"Invalid or duplicate FAI entry: {name}")
            sizes[name] = length
    selected = chromosomes.split(",") if chromosomes else [c for c in CANONICAL if c in sizes]
    if not selected or len(selected) != len(set(selected)):
        fail("No canonical chromosomes found, or repeated --chromosomes values")
    for name in selected:
        if not SAFE_NAME.fullmatch(name) or name not in sizes:
            fail(f"Chromosome missing from reference or unsafe name: {name!r}")
    return {c: sizes[c] for c in selected}


def check_tools(input_format, bigwig):
    wanted = ["Rscript"] + (["modkit", "samtools"] if input_format == "bam" else [])
    if bigwig:
        wanted.append("bedGraphToBigWig")
    missing = [name for name in wanted if shutil.which(name) is None]
    if missing:
        fail("Missing executables: " + ", ".join(missing) + ". Activate the environment.yml environment.")
    subprocess.run(["Rscript", "-e", "stopifnot(requireNamespace('data.table', quietly=TRUE))"], check=True)


def tool_versions():
    versions = {"methmap": VERSION, "python": sys.version.split()[0]}
    for name, command in {
        "R": ["Rscript", "--version"], "samtools": ["samtools", "--version"],
        "modkit": ["modkit", "--version"], "bedGraphToBigWig": ["bedGraphToBigWig"],
        "data.table": ["Rscript", "-e", "cat(as.character(packageVersion('data.table')))"]
    }.items():
        if shutil.which(command[0]):
            result = subprocess.run(command, capture_output=True, text=True, errors="replace")
            lines = (result.stdout + result.stderr).strip().splitlines()
            versions[name] = lines[0] if lines else "unknown"
    return versions


def check_bam(path, sizes):
    if not any(p.is_file() for p in [Path(str(path) + ".bai"), path.with_suffix(".bai"),
                                    Path(str(path) + ".csi"), path.with_suffix(".csi")]):
        fail(f"BAM index missing: {path}; use samtools index")
    run(["samtools", "quickcheck", "-v", path])
    header = subprocess.check_output(["samtools", "view", "-H", str(path)], text=True)
    bam_sizes, sorted_bam = {}, False
    for line in header.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            sorted_bam = "SO:coordinate" in fields
        elif fields[0] == "@SQ":
            values = dict(field.split(":", 1) for field in fields[1:])
            bam_sizes[values["SN"]] = int(values["LN"])
    if not sorted_bam:
        fail(f"BAM must declare coordinate sorting: {path}")
    for chrom, size in sizes.items():
        if bam_sizes.get(chrom) != size:
            fail(f"BAM/reference sequence name or length mismatch: {chrom} in {path}")


def filter_bedmethyl(path, destination, sizes, minimum):
    """Stream unphased Modkit bedMethyl into one normalized table per chromosome."""
    destination.mkdir(parents=True)
    handles = {c: (destination / f"{c}.tsv").open("w") for c in sizes}
    qc = {"input_rows": 0, "selected_5mc_rows": 0, "below_coverage": 0, "retained_loci": 0}
    try:
        for handle in handles.values():
            handle.write("chr\tstart\tend\tpercent\tcoverage\n")
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as source:
            for line_number, line in enumerate(source, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.split()
                if len(fields) != 18:
                    fail(f"{path}:{line_number}: expected 18-column Modkit bedMethyl")
                qc["input_rows"] += 1
                if fields[0] not in sizes or fields[3] != "m":
                    continue
                qc["selected_5mc_rows"] += 1
                try:
                    start, end, coverage = int(fields[1]), int(fields[2]), int(fields[9])
                    counts = [int(x) for x in fields[11:18]]
                except ValueError:
                    fail(f"{path}:{line_number}: invalid numeric value")
                if not (0 <= start < end <= sizes[fields[0]]) or end != start + 1:
                    fail(f"{path}:{line_number}: expected a one-base interval within reference bounds")
                if fields[5] not in ("+", "-"):
                    fail(f"{path}:{line_number}: strand must be + or -; combined-strand input is unsupported")
                if coverage < 0 or min(counts) < 0 or sum(counts[:3]) != coverage:
                    fail(f"{path}:{line_number}: invalid valid-coverage/count fields")
                if coverage < minimum:
                    qc["below_coverage"] += 1
                    continue
                try:
                    percent = float(fields[10])
                except ValueError:
                    fail(f"{path}:{line_number}: invalid methylation percentage")
                if not math.isfinite(percent) or not 0 <= percent <= 100:
                    fail(f"{path}:{line_number}: methylation must be a finite percentage in [0, 100]")
                if coverage and abs(percent - 100 * counts[0] / coverage) > 0.11:
                    fail(f"{path}:{line_number}: percentage/count mismatch; column 11 must use 0–100")
                handles[fields[0]].write(f"{fields[0]}\t{start}\t{end}\t{percent}\t{coverage}\n")
                qc["retained_loci"] += 1
    finally:
        for handle in handles.values():
            handle.close()
    return qc


def join_tracks(staging, sizes, bigwig):
    tracks = staging / "tracks"
    tracks.mkdir()
    empty = []
    for metric in METRICS:
        target = tracks / f"{metric}.bedGraph"
        with target.open("wb") as output:
            for chrom in sorted(sizes):
                with (staging / "per_chromosome" / metric / f"{chrom}.bedGraph").open("rb") as source:
                    shutil.copyfileobj(source, output)
        if target.stat().st_size and bigwig:
            run(["bedGraphToBigWig", target, staging / "chrom.sizes", tracks / f"{metric}.bw"],
                staging / "logs" / "bigwig.log")
        elif not target.stat().st_size:
            empty.append(metric)
    return empty


def execute(args):
    reference = args.reference.expanduser().resolve()
    fai = Path(str(reference) + ".fai")
    if not reference.is_file() or not fai.is_file():
        fail("Supply an existing uncompressed FASTA and its .fai index (samtools faidx reference.fa)")
    if reference.suffix in (".gz", ".bgz"):
        fail("This release requires an uncompressed reference FASTA")
    if not SAFE_NAME.fullmatch(args.tissue):
        fail("Tissue must contain letters, numbers, underscore, dash or dot")
    sizes = read_fai(fai, args.chromosomes)
    samples = load_samples(args.samples.expanduser().resolve(), args.tissue, args.input_format)
    out = args.out.expanduser().resolve()
    if out.exists():
        fail(f"Output already exists; use a new directory: {out}")
    base = None
    base_path = args.base.expanduser().resolve() if args.command == "add" else None
    if base_path:
        base = json.loads((base_path / "run.json").read_text())
        if base.get("schema_version") != 1 or base.get("status") != "complete":
            fail("Base must be a completed reference from this MethMap package")
        for name, value in [("tissue", args.tissue), ("assembly", args.assembly),
                            ("min_coverage", args.min_coverage), ("chromosomes", sizes)]:
            if base.get(name) != value:
                fail(f"Base/new run mismatch: {name}; use the same settings")
        if base.get("modkit_policy") != {"seed": args.seed, "filter_percentile": 0.1}:
            fail("Base/new run Modkit settings mismatch")
        for key in ["sample_id", "path"]:
            previous = {x[key] for x in base["samples"]}
            if previous.intersection(x[key] for x in samples):
                fail(f"A new sample reuses an existing {key}")
    total_samples = len(samples) + (len(base["samples"]) if base else 0)
    if args.min_samples > total_samples:
        fail("--min-samples cannot exceed cohort size")
    plan = dict(command=args.command, tissue=args.tissue, assembly=args.assembly,
                reference=str(reference), chromosomes=sizes, new_samples=samples,
                total_samples=total_samples, output=str(out), min_coverage=args.min_coverage,
                min_samples=args.min_samples, bigwig=not args.no_bigwig,
                threads=args.threads, jobs=args.jobs)
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return
    check_tools(args.input_format, not args.no_bigwig)
    if args.input_format == "bam":
        for sample in samples:
            check_bam(Path(sample["path"]), sizes)
    print("Computing reference SHA-256...", flush=True)
    reference_hash = digest(reference)
    if base and reference_hash != base["reference"]["sha256"]:
        fail("Reference FASTA differs from the base reference")
    out.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=out.name + ".partial-", dir=out.parent))
    print(f"Working directory: {staging}", flush=True)
    try:
        for name in ["logs", "reference", "manifests", "samples", "qc"]:
            (staging / name).mkdir()
        (staging / "chrom.sizes").write_text("".join(f"{c}\t{sizes[c]}\n" for c in sorted(sizes)))
        r_env = {**os.environ, "METHMAP_R_THREADS": str(max(1, args.threads // args.jobs))}
        sample_qc = list(base.get("sample_qc", [])) if base else []
        for sample in samples:
            sid = sample["sample_id"]
            source = Path(sample["path"])
            if args.input_format == "bam":
                source = staging / "samples" / f"{sid}.bedmethyl"
                run(["modkit", "pileup", sample["path"], source, "--modified-bases", "5mC", "5hmC",
                     "--reference", reference, "--threads", args.threads,
                     "--seed", args.seed, "--filter-percentile", "0.1",
                     "--log-filepath", staging / "logs" / f"{sid}.modkit.log"],
                    staging / "logs" / f"{sid}.pileup.log")
            qc = filter_bedmethyl(source, staging / "samples" / sid, sizes, args.min_coverage)
            sample_qc.append({"sample_id": sid, **qc})
            print(f"{sid}: {qc['retained_loci']} retained 5mC loci", flush=True)
            run(["Rscript", ROOT / "scripts" / "procesador.R", staging / "samples" / sid],
                staging / "logs" / f"{sid}.normalize.log", env=r_env)
            if args.input_format == "bam" and not args.keep_intermediates:
                source.unlink()

        def chromosome_job(chrom):
            manifest = staging / "manifests" / f"{chrom}.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_id", "path"])
                for sample in samples:
                    writer.writerow([sample["sample_id"], staging / "samples" / sample["sample_id"] / f"{chrom}.rds"])
            target = staging / "reference" / f"{chrom}.rds"
            if base_path:
                command = ["Rscript", ROOT / "scripts" / "Add_info.R", base_path / "reference" / f"{chrom}.rds", manifest, target]
            else:
                command = ["Rscript", ROOT / "scripts" / "baseconstructor.R", chrom, manifest, target]
            log = staging / "logs" / f"{chrom}.reference.log"
            run(command, log, env=r_env)
            run(["Rscript", ROOT / "scripts" / "rds_to_bed.R", target, staging, args.min_samples], log, env=r_env)

        with ThreadPoolExecutor(max_workers=min(args.jobs, len(sizes))) as pool:
            list(pool.map(chromosome_job, sizes))
        qc_rows = []
        for chrom in sizes:
            with (staging / "qc" / f"{chrom}.tsv").open() as handle:
                qc_rows.extend(csv.DictReader(handle, delimiter="\t"))
        if not sum(int(x["observed_loci"]) for x in qc_rows):
            fail("No eligible 5mC sites; check chromosome names, code m and coverage threshold")
        with (staging / "qc" / "chromosomes.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(qc_rows)
        with (staging / "qc" / "samples.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(sample_qc[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(sample_qc)
        empty_tracks = join_tracks(staging, sizes, not args.no_bigwig)
        state = {**plan, "schema_version": 1, "status": "complete", "version": VERSION,
                 "created_utc": datetime.now(timezone.utc).isoformat(),
                 "samples": (base["samples"] if base else []) + samples,
                 "reference": {"path": str(reference), "sha256": reference_hash, "fai_sha256": digest(fai)},
                 "modkit_policy": {"seed": args.seed, "filter_percentile": 0.1},
                 "sample_qc": sample_qc, "empty_tracks": empty_tracks,
                 "tool_versions": tool_versions(),
                 "parent": {"path": str(base_path), "run_sha256": digest(base_path / "run.json")} if base else None}
        (staging / "run.json").write_text(json.dumps(state, indent=2) + "\n")
        shutil.rmtree(staging / "manifests")
        if not args.keep_intermediates:
            shutil.rmtree(staging / "samples")
            shutil.rmtree(staging / "per_chromosome")
        if out.exists():
            fail(f"Output was created by another process during this run: {out}")
        staging.rename(out)
        print(f"MethMap completed: {out}", flush=True)
    except Exception:
        print(f"Run failed. Diagnostic files retained at {staging}; final output was not created.", file=sys.stderr)
        raise


def parser():
    p = argparse.ArgumentParser(description="MethMap: tissue-specific 5mC reference maps, focused on hg38.")
    p.add_argument("--version", action="version", version=f"MethMap {VERSION}")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor", help="Show installed executables and versions")
    for mode in ["build", "add"]:
        q = sub.add_parser(mode, help="Build a reference" if mode == "build" else "Extend a reference into a new directory")
        q.add_argument("--samples", required=True, type=Path, help="TSV: sample_id, tissue, path; paths relative to TSV")
        q.add_argument("--tissue", required=True)
        q.add_argument("--reference", required=True, type=Path, help="Uncompressed FASTA; adjacent .fai required")
        q.add_argument("--out", required=True, type=Path, help="New output directory (must not exist)")
        if mode == "add":
            q.add_argument("--base", required=True, type=Path, help="Completed MethMap output directory")
        q.add_argument("--input-format", choices=["bam", "bedmethyl"], default="bam")
        q.add_argument("--assembly", default="hg38", help="Assembly label recorded in metadata; not a sequence identity check")
        q.add_argument("--chromosomes", help="Comma-separated contigs; default chr1..chr22,chrX,chrY present in FAI")
        q.add_argument("--min-coverage", type=positive, default=10)
        q.add_argument("--min-samples", type=positive, default=1, help="Minimum observed samples for exported tracks")
        q.add_argument("--threads", type=positive, default=4)
        q.add_argument("--jobs", type=positive, default=1, help="Concurrent chromosomes; memory multiplies with this value")
        q.add_argument("--seed", type=int, default=42, help="Seed for Modkit threshold sampling")
        q.add_argument("--no-bigwig", action="store_true", help="Export bedGraph and statistics only")
        q.add_argument("--keep-intermediates", action="store_true")
        q.add_argument("--dry-run", action="store_true", help="Validate manifest and FAI, print plan; do not run tools")
    return p


def main():
    args = parser().parse_args()
    try:
        if args.command == "doctor":
            print(json.dumps({"executables": {x: shutil.which(x) for x in
                  ["Rscript", "samtools", "modkit", "bedGraphToBigWig"]}, "versions": tool_versions()}, indent=2))
        else:
            execute(args)
    except (ValueError, OSError, subprocess.CalledProcessError, KeyError) as error:
        print(f"MethMap error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
