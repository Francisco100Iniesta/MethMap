#!/usr/bin/env python3
"""Create deterministic synthetic bedMethyl records. No biological sample data."""
from pathlib import Path
import math

DEST = Path(__file__).resolve().parent / "demo"


def row(chrom, start, percent, coverage=100, code="m"):
    modified = round(percent * coverage / 100)
    actual = 100 * modified / coverage
    return "\t".join(map(str, [chrom, start, start + 1, code, coverage, "+", start,
                    start + 1, "255,0,0", coverage, round(actual, 2), modified,
                    coverage - modified, 0, 0, 0, 0, 0])) + "\n"


def main():
    DEST.mkdir(exist_ok=True)
    indexes = []
    with (DEST / "synthetic.fa").open("wb") as handle:
        for chrom, sequence in [("chr1", "ACGT" * 500), ("chr2", "ACGT" * 125)]:
            handle.write(f">{chrom}\n".encode())
            indexes.append(f"{chrom}\t{len(sequence)}\t{handle.tell()}\t50\t51\n")
            for start in range(0, len(sequence), 50):
                handle.write((sequence[start:start+50] + "\n").encode())
    (DEST / "synthetic.fa.fai").write_text("".join(indexes))
    for sample in range(1, 8):
        rows = [row("chr1", 1, min(100, (sample - 1) * 20))]
        for i in range(150):
            if (i + 3 * sample) % 23 == 0:
                continue  # Missing observation, not a methylation value of zero.
            level = 48 + 29 * math.sin(i / 22) + (sample - 3.5) * (2 + 9 * math.exp(-((i-75)/18)**2))
            level += 4 * math.sin(i * .8 + sample)
            rows.append(row("chr1", 101 + 4 * i, round(max(0, min(100, level)))))
        rows += [row("chr1", 901, 50, 8), row("chr1", 905, 50, code="h"), row("chrUn_demo", 1, 50)]
        if sample == 7:
            rows.append(row("chr1", 1001, 80))
            rows.append(row("chr2", 1, 60))
        (DEST / f"sample{sample:02}.bedmethyl").write_text("".join(reversed(rows)))
    (DEST / "samples.tsv").write_text("sample_id\ttissue\tpath\n" + "".join(
        f"sample{i:02}\tDEMO\tsample{i:02}.bedmethyl\n" for i in range(1, 7)))
    (DEST / "samples_to_add.tsv").write_text("sample_id\ttissue\tpath\nsample07\tDEMO\tsample07.bedmethyl\n")


if __name__ == "__main__":
    main()
