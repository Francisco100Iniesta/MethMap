"""Numerical and failure-path checks using only synthetic records."""
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("methmap", ROOT / "methmap.py")
mm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mm)


def bedrow(start=1, percent=50, coverage=10, code="m", chrom="chr1"):
    count = round(percent * coverage / 100)
    return "\t".join(map(str, [chrom, start, start+1, code, coverage, "+", start,
        start+1, "255,0,0", coverage, percent, count, coverage-count, 0, 0, 0, 0, 0])) + "\n"


class FilterTests(unittest.TestCase):
    def test_valid_coverage_threshold_and_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p / "input.bed").write_text(bedrow(percent=0) + bedrow(start=5, coverage=8) +
                bedrow(start=9, code="h") + bedrow(start=13, chrom="chrUn"))
            qc = mm.filter_bedmethyl(p / "input.bed", p / "out", {"chr1": 100}, 10)
            self.assertEqual(qc, dict(input_rows=4, selected_5mc_rows=2, below_coverage=1, retained_loci=1))
            self.assertIn("chr1\t1\t2\t0.0\t10", (p / "out/chr1.tsv").read_text())

    def test_invalid_coordinates_and_percentage_fail(self):
        for row in [bedrow(start=100), bedrow().replace("\t50\t5\t", "\t0.5\t5\t"),
                    bedrow().replace("\t50\t5\t", "\tnan\t5\t")]:
            with self.subTest(row=row), tempfile.TemporaryDirectory() as temp:
                p = Path(temp)
                (p / "input.bed").write_text(row)
                with self.assertRaises(ValueError):
                    mm.filter_bedmethyl(p / "input.bed", p / "out", {"chr1": 100}, 10)

    def test_gzip_and_unordered_input(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            with gzip.open(p / "in.bed.gz", "wt") as h:
                h.write(bedrow(start=9) + bedrow(start=1))
            qc = mm.filter_bedmethyl(p / "in.bed.gz", p / "out", {"chr1": 100}, 10)
            self.assertEqual(qc["retained_loci"], 2)

    def test_no_valid_calls_are_filtered_before_percent_validation(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p / "input.bed").write_text(bedrow(coverage=0).replace("\t50\t0\t", "\tnan\t0\t"))
            qc = mm.filter_bedmethyl(p / "input.bed", p / "out", {"chr1": 100}, 10)
            self.assertEqual(qc["below_coverage"], 1)
            self.assertEqual(qc["retained_loci"], 0)


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is required for integration checks")
class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="methmap test ")
        cls.work = Path(cls.temp.name)
        relocated = cls.work / "source code"
        shutil.copytree(ROOT / "scripts", relocated / "scripts")
        shutil.copyfile(ROOT / "methmap.py", relocated / "methmap.py")
        cls.base = cls.work / "base"
        cls.updated = cls.work / "updated"
        cls.common = ["--tissue", "DEMO", "--assembly", "synthetic", "--reference",
                      str(ROOT / "examples/demo/synthetic.fa"), "--input-format", "bedmethyl", "--no-bigwig", "--jobs", "2"]
        cls.call("build", "--samples", ROOT / "examples/demo/samples.tsv", "--out", cls.base, *cls.common,
                 entry=relocated / "methmap.py")
        cls.old_hash = hashlib.sha256((cls.base / "reference/chr1.rds").read_bytes()).hexdigest()
        cls.call("add", "--base", cls.base, "--samples", ROOT / "examples/demo/samples_to_add.tsv",
                 "--out", cls.updated, "--min-samples", "2", *cls.common, entry=relocated / "methmap.py")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @staticmethod
    def call(*args, ok=True, entry=None):
        result = subprocess.run([sys.executable, str(entry or ROOT / "methmap.py"), *map(str, args)],
                                capture_output=True, text=True)
        if (result.returncode == 0) != ok:
            raise AssertionError(result.stdout + result.stderr)
        return result

    def stats(self, out):
        with gzip.open(out / "statistics/chr1.tsv.gz", "rt") as h:
            return {int(row["start"]): row for row in csv.DictReader(h, delimiter="\t")}

    def test_known_quantiles_and_zero_values(self):
        row = self.stats(self.base)[1]
        for field, value in {"p5": 5, "p25": 25, "p50": 50, "p75": 75, "p95": 95,
                             "hits": 6, "coverage_pct": 100, "score": 100, "iqr": 50}.items():
            self.assertAlmostEqual(float(row[field]), value)

    def test_missingness_new_loci_and_base_immutability(self):
        rows = self.stats(self.updated)
        self.assertEqual(float(rows[1]["p50"]), 60)
        self.assertEqual(float(rows[1001]["hits"]), 1)
        self.assertAlmostEqual(float(rows[1001]["coverage_pct"]), 100/7)
        self.assertEqual(float(rows[1001]["p50"]), 80)
        self.assertNotIn("\t1001\t", (self.updated / "tracks/p50.bedGraph").read_text())
        self.assertEqual(hashlib.sha256((self.base / "reference/chr1.rds").read_bytes()).hexdigest(), self.old_hash)
        self.assertEqual(len(json.loads((self.updated / "run.json").read_text())["samples"]), 7)
        with gzip.open(self.updated / "statistics/chr2.tsv.gz", "rt") as h:
            new_chromosome = list(csv.DictReader(h, delimiter="\t"))
        self.assertEqual(len(new_chromosome), 1)
        self.assertEqual(float(new_chromosome[0]["p50"]), 60)
        self.assertEqual(int(new_chromosome[0]["hits"]), 1)

    def test_empty_chromosome_and_sorted_tracks(self):
        with gzip.open(self.base / "statistics/chr2.tsv.gz", "rt") as h:
            self.assertEqual(len(h.read().splitlines()), 1)
        positions = [int(x.split("\t")[1]) for x in (self.base / "tracks/p50.bedGraph").read_text().splitlines()]
        self.assertEqual(positions, sorted(set(positions)))

    def test_duplicate_cohort_sample_rejected(self):
        result = self.call("add", "--base", self.base, "--samples", ROOT / "examples/demo/samples.tsv",
                           "--out", self.work / "duplicate", *self.common, ok=False)
        self.assertIn("reuses an existing sample_id", result.stderr)
        self.assertFalse((self.work / "duplicate").exists())

    def test_existing_output_protected(self):
        result = self.call("build", "--samples", ROOT / "examples/demo/samples.tsv",
                           "--out", self.base, *self.common, ok=False)
        self.assertIn("Output already exists", result.stderr)

    def test_duplicate_locus_fails_without_final_output(self):
        source = self.work / "duplicated.bedmethyl"
        source.write_text(bedrow() * 2)
        manifest = self.work / "duplicates.tsv"
        manifest.write_text("sample_id\ttissue\tpath\nduplicate\tDEMO\tduplicated.bedmethyl\n")
        output = self.work / "bad"
        self.call("build", "--samples", manifest, "--out", output, *self.common, ok=False)
        self.assertFalse(output.exists())
        partial = next(self.work.glob("bad.partial-*"))
        self.assertIn("Duplicate coordinates", (partial / "logs/duplicate.normalize.log").read_text())

    def test_fasta_change_rejected_even_at_same_lengths(self):
        fasta = self.work / "altered.fa"
        fasta.write_text((ROOT / "examples/demo/synthetic.fa").read_text().replace("ACGT", "TCGT", 1))
        shutil.copyfile(ROOT / "examples/demo/synthetic.fa.fai", str(fasta) + ".fai")
        result = self.call("add", "--base", self.base, "--samples", ROOT / "examples/demo/samples_to_add.tsv",
                           "--out", self.work / "mismatch", *self.common, "--reference", fasta, ok=False)
        self.assertIn("FASTA differs", result.stderr)


@unittest.skipUnless(all(shutil.which(t) for t in ["Rscript", "modkit", "samtools", "bedGraphToBigWig"]),
                     "Full BAM integration requires the environment.yml tools")
class BamIntegrationTests(unittest.TestCase):
    def test_modified_bam_to_bigwig(self):
        with tempfile.TemporaryDirectory(prefix="methmap-bam-") as temp:
            p = Path(temp)
            ref = ROOT / "examples/demo/synthetic.fa"
            for sid, methylated_reads in [("A", 3), ("B", 9)]:
                sam = p / f"{sid}.sam"
                lines = ["@HD\tVN:1.6\tSO:coordinate\n", "@SQ\tSN:chr1\tLN:2000\n", "@SQ\tSN:chr2\tLN:500\n"]
                for i in range(12):
                    ml = ",".join(["255,0" if i < methylated_reads else "0,0"] * 25)
                    mm_tag = "C+mh," + ",".join(["0"] * 25) + ";"
                    fields = [f"read{i}", "0", "chr1", "1", "60", "100M", "*", "0", "0", "ACGT"*25, "I"*100,
                              "MM:Z:" + mm_tag, "ML:B:C," + ml]
                    lines.append("\t".join(fields) + "\n")
                sam.write_text("".join(lines))
                subprocess.run(["samtools", "view", "-b", "-o", str(p / f"{sid}.bam"), str(sam)], check=True)
                subprocess.run(["samtools", "index", str(p / f"{sid}.bam")], check=True)
            manifest = p / "samples.tsv"
            manifest.write_text("sample_id\ttissue\tpath\nA\tDEMO\tA.bam\nB\tDEMO\tB.bam\n")
            ReferenceTests.call("build", "--samples", manifest, "--tissue", "DEMO", "--reference", ref,
                                "--assembly", "synthetic", "--out", p / "out", "--threads", "2")
            with gzip.open(p / "out/statistics/chr1.tsv.gz", "rt") as h:
                rows = list(csv.DictReader(h, delimiter="\t"))
            self.assertEqual(len(rows), 25)
            self.assertAlmostEqual(float(rows[0]["p50"]), 50)
            self.assertAlmostEqual(float(rows[0]["p5"]), 27.5)
            self.assertEqual(len(list((p / "out/tracks").glob("*.bw"))), 23)
            for path in (p / "out/tracks").glob("*.bw"):
                self.assertEqual(path.read_bytes()[:4], bytes.fromhex("26fc8f88"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
