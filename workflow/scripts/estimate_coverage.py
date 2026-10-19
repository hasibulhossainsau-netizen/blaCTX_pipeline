#!/usr/bin/env python3
"""
estimate_coverage.py
---------------------
Coverage = total trimmed bases (R1 + R2) / assembly length (bp).
Writes a single-line text file: "<coverage_x>\t<total_bases>\t<assembly_length>"

Usage:
  python estimate_coverage.py --r1 R1.fastq.gz --r2 R2.fastq.gz \
      --assembly assembly.fasta --output coverage.txt
"""

import argparse
import gzip
import sys
from pathlib import Path


def count_bases_fastq_gz(path: Path) -> int:
    total = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:  # sequence line
                total += len(line.strip())
    return total


def count_assembly_length(fasta_path: Path) -> int:
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--r1", required=True, type=Path, help="Trimmed R1 FASTQ (.gz)")
    ap.add_argument("--r2", default=None, type=Path, help="Trimmed R2 FASTQ (.gz); omit for single-end")
    ap.add_argument("--assembly", required=True, type=Path, help="Assembly FASTA")
    ap.add_argument("--output", required=True, type=Path, help="Output coverage text file")
    args = ap.parse_args()

    inputs = [args.r1, args.assembly] + ([args.r2] if args.r2 else [])
    for path in inputs:
        if not path.is_file():
            sys.exit(f"ERROR: input file not found: {path}")

    try:
        total_bases = count_bases_fastq_gz(args.r1)
        if args.r2:
            total_bases += count_bases_fastq_gz(args.r2)
        assembly_len = count_assembly_length(args.assembly)
    except (OSError, EOFError, UnicodeDecodeError) as exc:
        sys.exit(f"ERROR: could not read input: {exc}")

    if assembly_len == 0:
        sys.exit(f"ERROR: assembly has zero length: {args.assembly}")
    coverage = total_bases / assembly_len

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(f"{coverage:.2f}\t{total_bases}\t{assembly_len}\n")
    print(f"Coverage: {coverage:.2f}x  (bases={total_bases}, assembly_len={assembly_len})")


if __name__ == "__main__":
    main()
