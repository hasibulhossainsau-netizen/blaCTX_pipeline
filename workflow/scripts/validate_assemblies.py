#!/usr/bin/env python3
"""
validate_assemblies.py
------------------------
1. `filter-one`     - validate a single assembly against BUSCO/coverage/contig
                      thresholds; write a JSON status file consumed downstream.
2. `promote`        - copy a passing assembly into the final pass/ directory
                      (or write an empty placeholder for a failed one).
3. `build-summary` - aggregate filter status + RGI + ResFinder results for the
                      whole cohort into a single pandas DataFrame, exported as
                      TSV + XLSX, with target gene presence/absence flagged.

Usage:
  python validate_assemblies.py filter-one --sample-id KPN001 --fasta assembly.fasta \
      --busco-summary short_summary.specific.txt --coverage-file coverage.txt \
      --busco-min 95 --coverage-min 50 --max-contigs 2 --min-contig-len 500 \
      --outdir results/04_filtered_assemblies --status-json KPN001.filter_status.json

  python validate_assemblies.py promote --status-json KPN001.filter_status.json \
      --fasta assembly.fasta --output KPN001.pass.fasta

  python validate_assemblies.py build-summary --samples-tsv config/samples.tsv \
      --filter-dir results/04_filtered_assemblies --rgi-dir results/05_amr/rgi \
      --resfinder-dir results/05_amr/resfinder --target-gene blaCTX-M-15 \
      --out-tsv results/06_metadata/cohort_metadata_summary.tsv \
      --out-xlsx results/06_metadata/cohort_metadata_summary.xlsx
"""

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

# amr_screening.py lives alongside this script
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amr_screening import (  # noqa: E402
    parse_rgi_target_hits,
    parse_resfinder_target_hits,
)


# =============================================================================
# Assembly parsing helpers
# =============================================================================

def parse_fasta_contigs(fasta_path: Path, min_contig_len: int) -> list[int]:
    """Return list of contig lengths >= min_contig_len."""
    lengths = []
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_len > 0:
        lengths.append(current_len)
    return [l for l in lengths if l >= min_contig_len]


def parse_busco_completeness(busco_summary_path: Path) -> float | None:
    """
    Extract 'C:<pct>%[S:<pct>%,D:<pct>%],F:...,M:...,n:...' line from a BUSCO
    short_summary*.txt file and return the total Complete percentage.
    """
    if not busco_summary_path.exists():
        return None
    text = busco_summary_path.read_text()
    match = re.search(r"C:([\d.]+)%", text)
    return float(match.group(1)) if match else None


def parse_coverage_file(coverage_path: Path) -> float | None:
    if not coverage_path.exists():
        return None
    line = coverage_path.read_text().strip().split("\t")
    try:
        return float(line[0])
    except (ValueError, IndexError):
        return None


# =============================================================================
# filter-one
# =============================================================================

def cmd_filter_one(args) -> None:
    contig_lengths = parse_fasta_contigs(args.fasta, args.min_contig_len)
    n_contigs = len(contig_lengths)
    total_length = sum(contig_lengths)

    busco_pct = parse_busco_completeness(args.busco_summary)
    coverage_x = parse_coverage_file(args.coverage_file)

    reasons = []
    busco_pass = busco_pct is not None and busco_pct >= args.busco_min
    if not busco_pass:
        reasons.append(f"BUSCO completeness {busco_pct} < {args.busco_min}")

    # Coverage is exempted (treated as pass) if genuinely unavailable (WGS_ASM samples)
    if coverage_x is None:
        coverage_pass = True
        reasons.append("Coverage NA (no read data; exempted, verify manually)")
    else:
        coverage_pass = coverage_x >= args.coverage_min
        if not coverage_pass:
            reasons.append(f"Coverage {coverage_x}x < {args.coverage_min}x")

    contig_pass = n_contigs <= args.max_contigs
    if not contig_pass:
        reasons.append(f"Contig count {n_contigs} > {args.max_contigs}")

    overall_pass = busco_pass and coverage_pass and contig_pass

    status = {
        "sample_id": args.sample_id,
        "busco_completeness": busco_pct,
        "coverage_x": coverage_x,
        "n_contigs": n_contigs,
        "total_assembly_length": total_length,
        "pass": overall_pass,
        "fail_reasons": reasons,
    }

    args.status_json.parent.mkdir(parents=True, exist_ok=True)
    args.status_json.write_text(json.dumps(status, indent=2))
    print(json.dumps(status, indent=2))


# =============================================================================
# promote
# =============================================================================

def cmd_promote(args) -> None:
    status = json.loads(args.status_json.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if status["pass"]:
        args.output.write_text(args.fasta.read_text())
        print(f"PASS: promoted {args.fasta} -> {args.output}")
    else:
        # Write empty placeholder so the DAG completes; failure is recorded
        # in filter_status.json and surfaced in the metadata summary.
        args.output.write_text("")
        print(f"FAIL: {status['sample_id']} did not pass filters "
              f"({'; '.join(status['fail_reasons'])}); wrote empty placeholder.")


# =============================================================================
# build-summary
# =============================================================================

def cmd_build_summary(args) -> None:
    samples_df = pd.read_csv(args.samples_tsv, sep="\t", comment="#")
    records = []

    for _, row in samples_df.iterrows():
        sample_id = str(row["sample_id"])

        status_path = args.filter_dir / f"{sample_id}.filter_status.json"
        status = json.loads(status_path.read_text()) if status_path.exists() else {}

        rgi_txt = args.rgi_dir / f"{sample_id}.rgi.txt"
        rgi_hits = parse_rgi_target_hits(rgi_txt, args.target_gene) if rgi_txt.exists() else []

        resfinder_txt = args.resfinder_dir / sample_id / "ResFinder_results_tab.txt"
        if not resfinder_txt.exists():
            resfinder_txt = args.resfinder_dir / sample_id / "ResFinder_results.txt"
        
        resfinder_hits = parse_resfinder_target_hits(resfinder_txt, args.target_gene) if resfinder_txt.exists() else []

        rgi_detected = len(rgi_hits) > 0
        resfinder_detected = len(resfinder_hits) > 0
        # Concordant only if BOTH tools detected the gene
        concordant_call = bool(rgi_detected and resfinder_detected)

        records.append({
            "sample_id": sample_id,
            "biosample": row.get("biosample", "NA"),
            "country": row.get("country", "NA"),
            "collection_year": row.get("collection_year", "NA"),
            "expected_st": row.get("expected_st", "NA"),
            "source_type": row.get("source_type", "NA"),
            "busco_completeness": status.get("busco_completeness"),
            "coverage_x": status.get("coverage_x"),
            "n_contigs": status.get("n_contigs"),
            "total_assembly_length": status.get("total_assembly_length"),
            "qc_pass": status.get("pass", False),
            "qc_fail_reasons": "; ".join(status.get("fail_reasons", [])) if status and status.get("fail_reasons") else "no_status_file",
            f"{args.target_gene}_rgi_detected": rgi_detected,
            f"{args.target_gene}_rgi_identity": rgi_hits[0]["identity"] if rgi_hits else None,
            f"{args.target_gene}_resfinder_detected": resfinder_detected,
            f"{args.target_gene}_resfinder_identity": resfinder_hits[0]["identity"] if resfinder_hits else None,
            f"{args.target_gene}_concordant_call": concordant_call,
        })

    summary_df = pd.DataFrame.from_records(records)

    # Cohort-level QC stats for a quick sanity check in logs
    n_total = len(summary_df)
    n_pass = int(summary_df["qc_pass"].sum())
    n_ctx_m15 = int(summary_df[f"{args.target_gene}_concordant_call"].sum())
    print(f"Cohort summary: {n_pass}/{n_total} samples passed QC; "
          f"{n_ctx_m15}/{n_total} concordant {args.target_gene} calls (RGI + ResFinder).")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.out_tsv, sep="\t", index=False)
    summary_df.to_excel(args.out_xlsx, index=False, engine="openpyxl")
    print(f"Wrote {args.out_tsv} and {args.out_xlsx}")


# =============================================================================
# CLI
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("filter-one")
    p1.add_argument("--sample-id", required=True)
    p1.add_argument("--fasta", required=True, type=Path)
    p1.add_argument("--busco-summary", required=True, type=Path)
    p1.add_argument("--coverage-file", required=True, type=Path)
    p1.add_argument("--busco-min", required=True, type=float)
    p1.add_argument("--coverage-min", required=True, type=float)
    p1.add_argument("--max-contigs", required=True, type=int)
    p1.add_argument("--min-contig-len", required=True, type=int)
    p1.add_argument("--outdir", required=True, type=Path)
    p1.add_argument("--status-json", required=True, type=Path)
    p1.set_defaults(func=cmd_filter_one)

    p2 = sub.add_parser("promote")
    p2.add_argument("--status-json", required=True, type=Path)
    p2.add_argument("--fasta", required=True, type=Path)
    p2.add_argument("--output", required=True, type=Path)
    p2.set_defaults(func=cmd_promote)

    p3 = sub.add_parser("build-summary")
    p3.add_argument("--samples-tsv", required=True, type=Path)
    p3.add_argument("--filter-dir", required=True, type=Path)
    p3.add_argument("--rgi-dir", required=True, type=Path)
    p3.add_argument("--resfinder-dir", required=True, type=Path)
    p3.add_argument("--target-gene", default="blaCTX-M-15")
    p3.add_argument("--out-tsv", required=True, type=Path)
    p3.add_argument("--out-xlsx", required=True, type=Path)
    p3.set_defaults(func=cmd_build_summary)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()