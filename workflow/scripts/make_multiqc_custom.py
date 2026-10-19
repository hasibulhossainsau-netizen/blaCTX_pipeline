#!/usr/bin/env python3
"""
make_multiqc_custom.py
----------------------
Build MultiQC custom-content tables from the pipeline's own outputs:

  assembly_qc_mqc.tsv   per-sample BUSCO %, coverage, contig count, assembly length, QC pass
  amr_target_mqc.tsv    per-sample RGI / ResFinder hit counts for the target gene

Inputs are looked up by name in --input-dir:
  <sample>.filter_status.json, <sample>.rgi.txt, <sample>.ResFinder_results.txt

Usage:
  python make_multiqc_custom.py --input-dir . --target-gene blaCTX-M-15 --outdir .
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_screening import parse_resfinder_target_hits, parse_rgi_target_hits  # noqa: E402

ASSEMBLY_HEADER = """\
# id: 'assembly_qc'
# section_name: 'Assembly QC'
# description: 'BUSCO completeness, read depth and contiguity used for assembly filtering.'
# plot_type: 'table'
"""

AMR_HEADER = """\
# id: 'amr_target'
# section_name: 'Target gene screening ({gene})'
# description: 'Number of {gene} hits per tool (assemblies that failed QC are not screened).'
# plot_type: 'table'
"""


def write_table(path: Path, header: str, columns: list, rows: list) -> None:
    lines = [header.rstrip("\n"), "Sample\t" + "\t".join(columns)]
    for sample, values in rows:
        lines.append(sample + "\t" + "\t".join("" if v is None else str(v) for v in values))
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--input-dir", required=True, type=Path, help="Directory holding per-sample result files")
    ap.add_argument("--target-gene", default="blaCTX-M-15", help="Gene to count AMR hits for")
    ap.add_argument("--outdir", required=True, type=Path, help="Where to write the *_mqc.tsv files")
    args = ap.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"ERROR: input directory not found: {args.input_dir}")
    status_files = sorted(args.input_dir.glob("*.filter_status.json"))
    if not status_files:
        sys.exit(f"ERROR: no *.filter_status.json files in {args.input_dir}")

    assembly_rows, amr_rows = [], []
    for status_file in status_files:
        sample = status_file.name[: -len(".filter_status.json")]
        try:
            status = json.loads(status_file.read_text())
        except json.JSONDecodeError as exc:
            sys.exit(f"ERROR: malformed JSON in {status_file}: {exc}")
        assembly_rows.append((sample, [
            status.get("busco_completeness"),
            status.get("coverage_x"),
            status.get("n_contigs"),
            status.get("total_assembly_length"),
            status.get("pass"),
        ]))

        rgi = args.input_dir / f"{sample}.rgi.txt"
        rf = args.input_dir / f"{sample}.ResFinder_results.txt"
        amr_rows.append((sample, [
            len(parse_rgi_target_hits(rgi, args.target_gene)) if rgi.exists() else None,
            len(parse_resfinder_target_hits(rf, args.target_gene)) if rf.exists() else None,
        ]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_table(args.outdir / "assembly_qc_mqc.tsv", ASSEMBLY_HEADER,
                ["BUSCO complete (%)", "Coverage (x)", "Contigs", "Assembly length (bp)", "QC pass"],
                assembly_rows)
    write_table(args.outdir / "amr_target_mqc.tsv", AMR_HEADER.format(gene=args.target_gene),
                ["RGI hits", "ResFinder hits"], amr_rows)
    print(f"Wrote MultiQC tables for {len(status_files)} samples to {args.outdir}")


if __name__ == "__main__":
    main()
