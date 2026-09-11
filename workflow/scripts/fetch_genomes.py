#!/usr/bin/env python3
"""
fetch_genomes.py
-----------------
Fetch SRA read sets (PE/SE) or WGS assemblies for the K. pneumoniae cohort.
Prioritizes ENA Direct Download -> SRA Toolkit (prefetch+fasterq-dump) -> NCBI Entrez.
"""

import argparse
import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger("fetch_genomes")

MAX_RETRIES = 3
RETRY_BACKOFF_SEC = 15


def run_cmd(cmd: list[str], **kwargs) -> subprocess.CompletedProcess:
    log.info("RUN: %s", " ".join(cmd))
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kwargs)


def which_or_none(tool: str) -> str | None:
    return shutil.which(tool)


def is_valid_file(file_path: Path, min_bytes: int = 1000) -> bool:
    return file_path.exists() and file_path.stat().st_size > min_bytes


def with_retries(fn, *args, **kwargs):
    last_exc = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return fn(*args, **kwargs)
        except subprocess.CalledProcessError as e:
            last_exc = e
            log.warning(
                "Attempt %d/%d failed (rc=%d): %s",
                attempt, MAX_RETRIES, e.returncode, e.stderr.strip()[:500] if e.stderr else "",
            )
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF_SEC * attempt)
        except Exception as e:
            last_exc = e
            log.warning("Attempt %d/%d failed with exception: %s", attempt, MAX_RETRIES, e)
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF_SEC * attempt)
    raise last_exc


def normalize_layout(layout: str) -> str:
    """Normalize layout string to 'PE' or 'SE'."""
    layout_upper = layout.upper()
    if "PE" in layout_upper:
        return "PE"
    elif "SE" in layout_upper:
        return "SE"
    return "PE"


# -----------------------------------------------------------------------------
# SRA read fetching strategies
# -----------------------------------------------------------------------------

def fetch_sra_ena_ftp(accession: str, sample_id: str, outdir: Path, layout: str) -> bool:
    """
    Priority 1: Direct HTTP/FTP download from ENA using Portal API links.
    Eliminates fasterq-dump conversion overhead and improves speed.
    """
    layout = normalize_layout(layout)
    downloader = "wget" if which_or_none("wget") else ("curl" if which_or_none("curl") else None)
    if not downloader:
        log.warning("Neither wget nor curl is installed; skipping ENA direct download.")
        return False

    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp&format=json"
    r1_out = outdir / f"{sample_id}_R1.fastq.gz"
    r2_out = outdir / f"{sample_id}_R2.fastq.gz"

    try:
        log.info("Querying ENA Portal API for %s...", accession)
        req = urllib.request.urlopen(api_url, timeout=15)
        data = json.loads(req.read().decode("utf-8"))

        if not data or "fastq_ftp" not in data[0] or not data[0]["fastq_ftp"]:
            log.warning("ENA API returned no FASTQ FTP links for %s", accession)
            return False

        ftp_urls = data[0]["fastq_ftp"].split(";")

        if layout == "PE" and len(ftp_urls) >= 2:
            url_1 = f"https://{ftp_urls[0]}"
            url_2 = f"https://{ftp_urls[1]}"

            def _download_pe():
                if downloader == "wget":
                    run_cmd(["wget", "-q", "-O", str(r1_out), url_1])
                    run_cmd(["wget", "-q", "-O", str(r2_out), url_2])
                else:
                    run_cmd(["curl", "--fail", "--silent", "--show-error", "-o", str(r1_out), url_1])
                    run_cmd(["curl", "--fail", "--silent", "--show-error", "-o", str(r2_out), url_2])

            with_retries(_download_pe)
            if is_valid_file(r1_out) and is_valid_file(r2_out):
                log.info("ENA PE download successful for %s", accession)
                return True

        elif layout == "SE" and len(ftp_urls) >= 1:
            url_se = f"https://{ftp_urls[0]}"

            def _download_se():
                if downloader == "wget":
                    run_cmd(["wget", "-q", "-O", str(r1_out), url_se])
                else:
                    run_cmd(["curl", "--fail", "--silent", "--show-error", "-o", str(r1_out), url_se])

            with_retries(_download_se)
            if is_valid_file(r1_out):
                log.info("ENA SE download successful for %s", accession)
                return True

    except Exception as e:
        log.warning("ENA Direct download failed for %s: %s", accession, e)

    # Clean up broken files on failure
    r1_out.unlink(missing_ok=True)
    if layout == "PE":
        r2_out.unlink(missing_ok=True)
    return False


def fetch_sra_prefetch(accession: str, sample_id: str, outdir: Path, threads: int, layout: str) -> bool:
    """
    Priority 2: Fallback to NCBI sra-tools (prefetch + fasterq-dump).
    """
    layout = normalize_layout(layout)
    if which_or_none("prefetch") is None or which_or_none("fasterq-dump") is None:
        log.warning("sra-tools (prefetch/fasterq-dump) not found; skipping fallback.")
        return False

    cache_dir = outdir / "_sra_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    def _prefetch():
        return run_cmd(["prefetch", accession, "-O", str(cache_dir)])

    try:
        with_retries(_prefetch)
        sra_target = cache_dir / accession
        if not sra_target.exists():
            sra_target = accession

        def _fasterq():
            cmd = ["fasterq-dump", str(sra_target), "-O", str(outdir), "-e", str(threads), "--skip-technical"]
            if layout == "PE":
                cmd.append("--split-files")
            return run_cmd(cmd)

        with_retries(_fasterq)
    except Exception as e:
        log.warning("sra-tools failed for %s: %s", accession, e)
        shutil.rmtree(cache_dir, ignore_errors=True)
        return False

    r1_dest = outdir / f"{sample_id}_R1.fastq.gz"

    if layout == "PE":
        r1_src = outdir / f"{accession}_1.fastq"
        r2_src = outdir / f"{accession}_2.fastq"
        r2_dest = outdir / f"{sample_id}_R2.fastq.gz"

        if not (is_valid_file(r1_src) and is_valid_file(r2_src)):
            shutil.rmtree(cache_dir, ignore_errors=True)
            return False

        _gzip_and_rename(r1_src, r1_dest)
        _gzip_and_rename(r2_src, r2_dest)
        success = is_valid_file(r1_dest) and is_valid_file(r2_dest)

    else:  # SE Layout
        r1_src = outdir / f"{accession}.fastq"
        if not is_valid_file(r1_src):
            r1_src = outdir / f"{accession}_1.fastq"
            if not is_valid_file(r1_src):
                shutil.rmtree(cache_dir, ignore_errors=True)
                return False

        _gzip_and_rename(r1_src, r1_dest)
        success = is_valid_file(r1_dest)

    shutil.rmtree(cache_dir, ignore_errors=True)
    return success


def fetch_sra_entrez_direct(accession: str, sample_id: str, outdir: Path, email: str, api_key: str, layout: str) -> bool:
    """
    Priority 3: Final fallback to NCBI Entrez Direct (efetch).
    """
    layout = normalize_layout(layout)
    if which_or_none("efetch") is None:
        log.warning("efetch not found; skipping NCBI Entrez fallback.")
        return False

    out_fastq = outdir / f"{sample_id}_raw.fastq"
    cmd = ["efetch", "-db", "sra", "-id", accession, "-format", "fastq"]

    try:
        with open(out_fastq, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, text=True, env=_entrez_env(email, api_key))
    except subprocess.CalledProcessError:
        out_fastq.unlink(missing_ok=True)
        return False

    if not is_valid_file(out_fastq):
        out_fastq.unlink(missing_ok=True)
        return False

    r1_dest = outdir / f"{sample_id}_R1.fastq.gz"

    if layout == "PE":
        r2_dest = outdir / f"{sample_id}_R2.fastq.gz"
        _split_interleaved_fastq(out_fastq, r1_dest, r2_dest)
        out_fastq.unlink(missing_ok=True)
        return is_valid_file(r1_dest) and is_valid_file(r2_dest)
    else:
        _gzip_and_rename(out_fastq, r1_dest)
        return is_valid_file(r1_dest)


def _entrez_env(email: str, api_key: str) -> dict:
    env = os.environ.copy()
    if email:
        env["NCBI_EMAIL"] = email
    if api_key:
        env["NCBI_API_KEY"] = api_key
    return env


def _gzip_and_rename(src: Path, dest_gz: Path) -> None:
    with open(src, "rb") as f_in, gzip.open(dest_gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    src.unlink(missing_ok=True)


def _split_interleaved_fastq(interleaved: Path, r1_out: Path, r2_out: Path) -> None:
    with open(interleaved) as fh, gzip.open(r1_out, "wt") as r1, gzip.open(r2_out, "wt") as r2:
        record = []
        toggle = True
        for line in fh:
            record.append(line)
            if len(record) == 4:
                target = r1 if toggle else r2
                target.writelines(record)
                record = []
                toggle = not toggle


def fetch_sra(accession: str, sample_id: str, outdir: Path, threads: int, email: str, api_key: str, layout: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    norm_layout = normalize_layout(layout)

    r1_check = outdir / f"{sample_id}_R1.fastq.gz"

    if norm_layout == "PE":
        r2_check = outdir / f"{sample_id}_R2.fastq.gz"
        if is_valid_file(r1_check) and is_valid_file(r2_check):
            log.info("Valid PE files already exist for %s, skipping download.", sample_id)
            return
    else:
        if is_valid_file(r1_check):
            log.info("Valid SE file already exists for %s, skipping download.", sample_id)
            return

    # ENA Direct Download -> SRA Toolkit -> NCBI Entrez
    methods = [
        ("ENA Direct Download", lambda: fetch_sra_ena_ftp(accession, sample_id, outdir, norm_layout)),
        ("SRA Toolkit (prefetch+fasterq-dump)", lambda: fetch_sra_prefetch(accession, sample_id, outdir, threads, norm_layout)),
        ("NCBI Entrez (efetch)", lambda: fetch_sra_entrez_direct(accession, sample_id, outdir, email, api_key, norm_layout)),
    ]

    for name, method in methods:
        log.info("Attempting fetch strategy: %s for accession %s", name, accession)
        try:
            if method():
                log.info("Successfully fetched %s (%s) using %s", accession, sample_id, name)
                return
        except Exception:
            log.exception("Strategy '%s' failed for %s, attempting next fallback...", name, accession)

    raise RuntimeError(f"All fetch methods exhausted for accession {accession} (sample {sample_id})")


# -----------------------------------------------------------------------------
# WGS assembly fetching
# -----------------------------------------------------------------------------

def fetch_wgs_acc_download(accession: str, outdir: Path) -> Path | None:
    if which_or_none("ncbi-acc-download") is None:
        return None
    def _download():
        return run_cmd(["ncbi-acc-download", "--format", "fasta", "-o", str(outdir / f"{accession}.fasta"), accession])
    with_retries(_download)
    fasta_path = outdir / f"{accession}.fasta"
    return fasta_path if is_valid_file(fasta_path) else None


def fetch_wgs_datasets(accession: str, outdir: Path) -> Path | None:
    if which_or_none("datasets") is None:
        return None
    def _download():
        return run_cmd(["datasets", "download", "genome", "accession", accession, "--include", "genome", "--filename", str(outdir / f"{accession}.zip")])
    with_retries(_download)
    zip_path = outdir / f"{accession}.zip"
    if not zip_path.exists():
        return None
    run_cmd(["unzip", "-o", str(zip_path), "-d", str(outdir)])
    candidates = list(outdir.rglob("*.fna"))
    return candidates[0] if candidates else None


def fetch_wgs(accession: str, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = fetch_wgs_acc_download(accession, outdir) or fetch_wgs_datasets(accession, outdir)
    if fasta is None:
        raise RuntimeError(f"Could not fetch WGS assembly {accession}")
    dest = outdir / "assembly.fasta"
    shutil.move(str(fasta), str(dest))
    log.info("Fetched WGS assembly %s -> %s", accession, dest)


# -----------------------------------------------------------------------------
# CLI Entry Point
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["sra", "wgs"], required=True)
    ap.add_argument("--layout", choices=['PE', 'SE', 'SRA_PE', 'SRA_SE'], default="PE", help="Read layout (PE or SE)")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    if args.accession in ("NA", "", None):
        log.error("No accession provided for sample %s", args.sample_id)
        sys.exit(1)

    try:
        if args.mode == "sra":
            fetch_sra(args.accession, args.sample_id, args.outdir, args.threads, args.email, args.api_key, args.layout)
        else:
            fetch_wgs(args.accession, args.outdir)
    except Exception as e:
        log.error("FATAL: failed to fetch %s (%s): %s", args.accession, args.sample_id, e)
        sys.exit(1)


if __name__ == "__main__":
    main()