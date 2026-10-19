#!/usr/bin/env python3
"""
amr_screening.py
----------------

Subprocess wrappers around RGI (CARD) and ResFinder for AMR gene detection in
assembled K. pneumoniae contigs, with specific attention to blaCTX-M-15.

Important:
    CARD/RGI may use different naming conventions for CTX-M enzymes.
    Therefore, this script distinguishes:

        1. exact blaCTX-M-15 detection
        2. CTX-M family detection
        3. other CTX-M alleles

    A CTX-M family hit is NOT automatically considered an exact
    blaCTX-M-15 hit.

Subcommands:
    rgi
        Run RGI against assembled contigs.

    resfinder
        Run ResFinder against assembled contigs.

    parse-hits
        Parse RGI/ResFinder outputs and report CTX-M family / exact
        blaCTX-M-15 matches.

Usage:
    python amr_screening.py rgi \
        --fasta pass.fasta \
        --card-json resources/card_data/card.json \
        --aligner BLAST \
        --out-prefix results/05_amr/rgi/KPN001.rgi \
        --threads 4

    python amr_screening.py resfinder \
        --fasta pass.fasta \
        --db-path resources/resfinder_db \
        --species Klebsiella \
        --outdir results/05_amr/resfinder/KPN001

    python amr_screening.py parse-hits \
        --rgi-txt results/05_amr/rgi/KPN001.rgi.txt \
        --resfinder-txt results/05_amr/resfinder/KPN001/ResFinder_results.txt
"""

import argparse
import csv
import json
import logging
import re
import subprocess
from pathlib import Path


# =============================================================================
# Logging
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

log = logging.getLogger("amr_screening")


# =============================================================================
# Constants
# =============================================================================

TARGET_GENE_DEFAULT = "blaCTX-M-15"

# Canonical representation used internally
TARGET_CANONICAL = "ctxm15"


# =============================================================================
# General utilities
# =============================================================================

def run_cmd(cmd, **kwargs):
    """
    Execute external command and raise an error on failure.
    """

    log.info("RUN: %s", " ".join(map(str, cmd)))

    result = subprocess.run(
        cmd,
        text=True,
        capture_output=True,
        **kwargs
    )

    if result.returncode != 0:
        log.error(
            "Command failed (rc=%d)\nSTDOUT:\n%s\nSTDERR:\n%s",
            result.returncode,
            result.stdout[-2000:],
            result.stderr[-2000:]
        )

        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            result.stdout,
            result.stderr
        )

    return result


def normalize_gene_name(name: str) -> str:
    """
    Normalize a gene/enzyme name for comparison.

    Examples:
        blaCTX-M-15  -> blactxm15
        CTX-M-15     -> ctxm15
        blaCTX_M_15  -> blactxm15
    """

    if not name:
        return ""

    return re.sub(
        r"[^a-z0-9]",
        "",
        name.lower()
    )


def extract_ctxm_allele(name: str) -> str | None:
    """
    Extract CTX-M allele number from a gene/hit name.

    Examples:
        blaCTX-M-15 -> 15
        CTX-M-15 -> 15
        CTX-M-55 -> 55

    Returns:
        allele number as string, or None.
    """

    if not name:
        return None

    cleaned = name.lower()

    match = re.search(
        r"ctx[\s_-]*m[\s_-]*(\d+)",
        cleaned
    )

    return match.group(1) if match else None


def classify_ctxm_hit(name: str, target_gene: str = TARGET_GENE_DEFAULT) -> str:
    """
    Classify a hit.

    Returns one of:

        exact_target
        other_ctxm
        ctxm_family
        not_ctxm

    This function deliberately avoids treating every CTX-M hit as
    blaCTX-M-15.
    """

    if not name:
        return "not_ctxm"

    target_norm = normalize_gene_name(target_gene)
    hit_norm = normalize_gene_name(name)

    # -------------------------------------------------------------------------
    # Exact normalized match
    # -------------------------------------------------------------------------

    if target_norm and target_norm in hit_norm:
        return "exact_target"

    # Target allele
    target_allele = extract_ctxm_allele(target_gene)

    # -------------------------------------------------------------------------
    # CTX-M family
    # -------------------------------------------------------------------------

    hit_allele = extract_ctxm_allele(name)

    if hit_allele is not None:

        if target_allele is not None and hit_allele == target_allele:
            return "exact_target"

        return "other_ctxm"

    # -------------------------------------------------------------------------
    # Generic CTX-M family hit where allele cannot be extracted
    # -------------------------------------------------------------------------

    if re.search(r"ctx[\s_-]*m", name.lower()):
        return "ctxm_family"

    return "not_ctxm"


# =============================================================================
# RGI
# =============================================================================

def run_rgi(
    fasta: Path,
    card_json: Path,
    aligner: str,
    out_prefix: Path,
    threads: int
) -> None:
    """
    Run RGI.

    Note:
        card_json is retained as a CLI parameter for compatibility with the
        existing pipeline. RGI itself uses its configured CARD database.
    """

    out_prefix.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    txt_path = Path(str(out_prefix) + ".txt")
    json_path = Path(str(out_prefix) + ".json")

    # -------------------------------------------------------------------------
    # Empty FASTA guard
    # -------------------------------------------------------------------------

    if not fasta.exists() or fasta.stat().st_size == 0:

        log.warning(
            "Empty or missing FASTA provided for RGI (%s). "
            "Creating empty output files.",
            fasta
        )

        txt_path.touch()
        json_path.touch()

        return

    # -------------------------------------------------------------------------
    # RGI command
    # -------------------------------------------------------------------------

    cmd = [
        "rgi",
        "main",

        "--input_sequence",
        str(fasta),

        "--output_file",
        str(out_prefix),

        "--input_type",
        "contig",

        "--alignment_tool",
        aligner,

        "-n",
        str(threads),

        "--clean",
    ]

    run_cmd(cmd)

    log.info(
        "RGI complete: %s",
        txt_path
    )


def parse_rgi_hits(
    rgi_txt: Path,
    target_gene: str = TARGET_GENE_DEFAULT
) -> list[dict]:
    """
    Parse RGI tabular output.

    Returns ALL CTX-M-related hits, while explicitly classifying whether each
    hit is:

        exact_target
        other_ctxm
        ctxm_family
        not_ctxm

    Only CTX-M-related hits are returned.
    """

    hits = []

    if not rgi_txt.exists() or rgi_txt.stat().st_size == 0:

        log.warning(
            "RGI output missing or empty: %s",
            rgi_txt
        )

        return hits

    with open(
        rgi_txt,
        newline="",
        encoding="utf-8"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            best_hit = (
                row.get("Best_Hit_ARO", "")
                or row.get("ARO", "")
                or row.get("Gene", "")
                or ""
            )

            match_type = classify_ctxm_hit(
                best_hit,
                target_gene
            )

            # Ignore non-CTX-M hits
            if match_type == "not_ctxm":
                continue

            allele = extract_ctxm_allele(best_hit)

            hits.append({
                "gene": best_hit,

                "target_match_type": match_type,

                "ctxm_allele": allele,

                "is_exact_target": (
                    match_type == "exact_target"
                ),

                "is_ctxm_family": True,

                "contig": row.get(
                    "Contig",
                    ""
                ),

                "identity": row.get(
                    "Best_Identities",
                    ""
                ),

                "coverage": row.get(
                    "Coverage",
                    ""
                ),

                "cutoff": row.get(
                    "Cut_Off",
                    ""
                ),

                "model_type": row.get(
                    "Model_type",
                    ""
                ),

                "drug_class": row.get(
                    "Drug Class",
                    row.get(
                        "Drug_Class",
                        ""
                    )
                ),
            })

    return hits


def parse_rgi_target_hits(
    rgi_txt: Path,
    target_gene: str = TARGET_GENE_DEFAULT
) -> list[dict]:
    """
    Return only exact target hits.

    This is the function downstream code should use when calculating
    blaCTX-M-15 prevalence.
    """

    hits = parse_rgi_hits(
        rgi_txt,
        target_gene
    )

    return [
        hit
        for hit in hits
        if hit["target_match_type"] == "exact_target"
    ]


def parse_rgi_ctxm_family_hits(
    rgi_txt: Path,
    target_gene: str = TARGET_GENE_DEFAULT
) -> list[dict]:
    """
    Return all CTX-M family hits.
    """

    return parse_rgi_hits(
        rgi_txt,
        target_gene
    )


# =============================================================================
# ResFinder
# =============================================================================

def run_resfinder(
    fasta: Path,
    db_path: Path,
    species: str,
    outdir: Path
) -> None:
    """
    Run ResFinder.
    """

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    res_txt = outdir / "ResFinder_results.txt"

    # -------------------------------------------------------------------------
    # Empty FASTA guard
    # -------------------------------------------------------------------------

    if not fasta.exists() or fasta.stat().st_size == 0:

        log.warning(
            "Empty or missing FASTA provided for ResFinder (%s). "
            "Creating empty output file.",
            fasta
        )

        res_txt.touch()

        return

    cmd = [
        "python3",
        "-m",
        "resfinder",

        "-ifa",
        str(fasta),

        "-o",
        str(outdir),

        "-db_res",
        str(db_path),

        "-s",
        species,

        "--acquired",

        "-l",
        "0.60",

        "-t",
        "0.90",
    ]

    run_cmd(cmd)

    log.info(
        "ResFinder complete: %s",
        outdir
    )


def parse_resfinder_hits(
    resfinder_txt: Path,
    target_gene: str = TARGET_GENE_DEFAULT
) -> list[dict]:
    """
    Parse ResFinder output.

    Only CTX-M-related hits are returned.

    Exact blaCTX-M-15 is classified as:
        target_match_type = exact_target

    Other CTX-M alleles are classified as:
        target_match_type = other_ctxm
    """

    hits = []

    if (
        not resfinder_txt.exists()
        or resfinder_txt.stat().st_size == 0
    ):

        log.warning(
            "ResFinder output missing or empty: %s",
            resfinder_txt
        )

        return hits

    with open(
        resfinder_txt,
        newline="",
        encoding="utf-8"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            gene_name = (
                row.get(
                    "Resistance gene",
                    ""
                )
                or row.get(
                    "Gene",
                    ""
                )
                or row.get(
                    "Name",
                    ""
                )
                or ""
            )

            match_type = classify_ctxm_hit(
                gene_name,
                target_gene
            )

            if match_type == "not_ctxm":
                continue

            allele = extract_ctxm_allele(
                gene_name
            )

            hits.append({
                "gene": gene_name,

                "target_match_type": match_type,

                "ctxm_allele": allele,

                "is_exact_target": (
                    match_type == "exact_target"
                ),

                "is_ctxm_family": True,

                "identity": row.get(
                    "Identity",
                    ""
                ),

                "coverage": row.get(
                    "Coverage",
                    ""
                ),

                "contig": row.get(
                    "Contig",
                    row.get(
                        "Position in contig",
                        ""
                    )
                ),

                "phenotype": row.get(
                    "Phenotype",
                    ""
                ),
            })

    return hits


def parse_resfinder_target_hits(
    resfinder_txt: Path,
    target_gene: str = TARGET_GENE_DEFAULT
) -> list[dict]:
    """
    Return only exact target hits.
    """

    hits = parse_resfinder_hits(
        resfinder_txt,
        target_gene
    )

    return [
        hit
        for hit in hits
        if hit["target_match_type"] == "exact_target"
    ]


# =============================================================================
# Summary helper
# =============================================================================

def summarize_target_detection(
    rgi_txt: Path | None = None,
    resfinder_txt: Path | None = None,
    target_gene: str = TARGET_GENE_DEFAULT
) -> dict:
    """
    Generate a compact detection summary.

    IMPORTANT:
        exact_target_detected means exact allele-level detection.

        ctxm_family_detected means a CTX-M family hit was found, regardless
        of exact allele.
    """

    result = {
        "target_gene": target_gene,

        "rgi_exact_target_detected": False,
        "rgi_ctxm_family_detected": False,

        "resfinder_exact_target_detected": False,
        "resfinder_ctxm_family_detected": False,

        "rgi_exact_hits": [],
        "rgi_ctxm_hits": [],

        "resfinder_exact_hits": [],
        "resfinder_ctxm_hits": [],

        "exact_target_concordant": None,
    }

    # -------------------------------------------------------------------------
    # RGI
    # -------------------------------------------------------------------------

    if rgi_txt is not None:

        rgi_hits = parse_rgi_hits(
            rgi_txt,
            target_gene
        )

        result["rgi_ctxm_hits"] = rgi_hits

        result["rgi_ctxm_family_detected"] = (
            len(rgi_hits) > 0
        )

        exact_rgi = [
            hit
            for hit in rgi_hits
            if hit["target_match_type"] == "exact_target"
        ]

        result["rgi_exact_hits"] = exact_rgi

        result["rgi_exact_target_detected"] = (
            len(exact_rgi) > 0
        )

    # -------------------------------------------------------------------------
    # ResFinder
    # -------------------------------------------------------------------------

    if resfinder_txt is not None:

        resfinder_hits = parse_resfinder_hits(
            resfinder_txt,
            target_gene
        )

        result["resfinder_ctxm_hits"] = (
            resfinder_hits
        )

        result["resfinder_ctxm_family_detected"] = (
            len(resfinder_hits) > 0
        )

        exact_rf = [
            hit
            for hit in resfinder_hits
            if hit["target_match_type"] == "exact_target"
        ]

        result["resfinder_exact_hits"] = exact_rf

        result["resfinder_exact_target_detected"] = (
            len(exact_rf) > 0
        )

    # -------------------------------------------------------------------------
    # Concordance
    # -------------------------------------------------------------------------

    if (
        rgi_txt is not None
        and resfinder_txt is not None
    ):

        result["exact_target_concordant"] = (
            result["rgi_exact_target_detected"]
            ==
            result["resfinder_exact_target_detected"]
        )

    return result


# =============================================================================
# CLI
# =============================================================================

def main():

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    sub = ap.add_subparsers(
        dest="cmd",
        required=True
    )

    # -------------------------------------------------------------------------
    # RGI
    # -------------------------------------------------------------------------

    p_rgi = sub.add_parser(
        "rgi",
        help="Run RGI against assembled contigs."
    )

    p_rgi.add_argument(
        "--fasta",
        required=True,
        type=Path
    )

    p_rgi.add_argument(
        "--card-json",
        required=True,
        type=Path
    )

    p_rgi.add_argument(
        "--aligner",
        default="BLAST"
    )

    p_rgi.add_argument(
        "--out-prefix",
        required=True,
        type=Path
    )

    p_rgi.add_argument(
        "--threads",
        type=int,
        default=4
    )

    # -------------------------------------------------------------------------
    # ResFinder
    # -------------------------------------------------------------------------

    p_rf = sub.add_parser(
        "resfinder",
        help="Run ResFinder."
    )

    p_rf.add_argument(
        "--fasta",
        required=True,
        type=Path
    )

    p_rf.add_argument(
        "--db-path",
        required=True,
        type=Path
    )

    p_rf.add_argument(
        "--species",
        default="Klebsiella"
    )

    p_rf.add_argument(
        "--outdir",
        required=True,
        type=Path
    )

    # -------------------------------------------------------------------------
    # Parse
    # -------------------------------------------------------------------------

    p_parse = sub.add_parser(
        "parse-hits",
        help="Parse RGI/ResFinder CTX-M hits."
    )

    p_parse.add_argument(
        "--rgi-txt",
        type=Path
    )

    p_parse.add_argument(
        "--resfinder-txt",
        type=Path
    )

    p_parse.add_argument(
        "--target-gene",
        default=TARGET_GENE_DEFAULT
    )

    args = ap.parse_args()

    # -------------------------------------------------------------------------
    # Commands
    # -------------------------------------------------------------------------

    if args.cmd == "rgi":

        run_rgi(
            args.fasta,
            args.card_json,
            args.aligner,
            args.out_prefix,
            args.threads
        )

    elif args.cmd == "resfinder":

        run_resfinder(
            args.fasta,
            args.db_path,
            args.species,
            args.outdir
        )

    elif args.cmd == "parse-hits":

        result = summarize_target_detection(
            rgi_txt=args.rgi_txt,
            resfinder_txt=args.resfinder_txt,
            target_gene=args.target_gene
        )

        print(
            json.dumps(
                result,
                indent=2
            )
        )


if __name__ == "__main__":
    main()