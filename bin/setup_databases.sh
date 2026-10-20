#!/usr/bin/env bash
# Download the reference databases the pipeline needs into resources/ (git-ignored).
#   CARD (RGI)      -> resources/card_data/card.json
#   ResFinder       -> resources/resfinder/db_resfinder   (+ db_pointfinder, db_disinfinder)
#   BUSCO lineage   -> resources/busco_db/enterobacterales_odb10
#   Bakta (light)   -> resources/db-light
# Usage: bin/setup_databases.sh [-r resources_dir] [-d "card resfinder busco bakta"] [-f]
# Needs: curl, tar, git; bakta needs `bakta_db` (conda: bakta env; or --bakta-db handled via wget fallback).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="${ROOT}/resources"
WHICH="card resfinder busco bakta"
FORCE=0
BUSCO_LINEAGE="enterobacterales_odb10"
BAKTA_ZENODO="https://zenodo.org/records/14916843/files/db-light.tar.xz"  # Bakta DB v6.0 light

while getopts "r:d:fh" o; do
  case "$o" in
    r) RES="$OPTARG" ;;
    d) WHICH="$OPTARG" ;;
    f) FORCE=1 ;;
    *) sed -n '2,8p' "$0"; exit 0 ;;
  esac
done

log() { printf '[setup_databases] %s\n' "$*"; }
need() { command -v "$1" >/dev/null || { log "missing required tool: $1"; exit 1; }; }
skip() {
  if [[ -e "$1" && $FORCE -eq 0 ]]; then
    log "$2 present ($1); use -f to refresh"
    return 0
  fi
  return 1
}
mkdir -p "$RES"
need curl; need tar

setup_card() {
  skip "$RES/card_data/card.json" CARD && return
  mkdir -p "$RES/card_data"
  local t; t="$(mktemp -d)"
  log "downloading CARD (latest)"
  curl -fL "https://card.mcmaster.ca/latest/data" -o "$t/card.tar.bz2"
  tar -xjf "$t/card.tar.bz2" -C "$RES/card_data" ./card.json
  rm -rf "$t"
}

setup_resfinder() {
  skip "$RES/resfinder/db_resfinder" ResFinder && return
  need git
  mkdir -p "$RES/resfinder"; cd "$RES/resfinder"
  for d in resfinder pointfinder disinfinder; do
    rm -rf "db_$d"
    log "cloning db_$d"
    git clone --depth 1 "https://bitbucket.org/genomicepidemiology/${d}_db.git" "db_$d"
  done
  # KMA index (needed by ResFinder when run in read mode; harmless for assemblies)
  if command -v kma_index >/dev/null; then
    (cd db_resfinder && python3 INSTALL.py non_interactive 2>/dev/null || true)
  else
    log "kma_index not found; skipping KMA indexing (only required for read-based runs)"
  fi
}

setup_busco() {
  skip "$RES/busco_db/$BUSCO_LINEAGE" BUSCO && return
  mkdir -p "$RES/busco_db"
  local t; t="$(mktemp -d)"
  log "downloading BUSCO lineage $BUSCO_LINEAGE"
  curl -fL "https://busco-data.ezlab.org/v5/data/lineages/${BUSCO_LINEAGE}.2024-01-08.tar.gz" -o "$t/l.tar.gz"
  tar -xzf "$t/l.tar.gz" -C "$RES/busco_db"
  rm -rf "$t"
}

setup_bakta() {
  skip "$RES/db-light/bakta.db" Bakta && return
  need tar
  local t; t="$(mktemp -d)"
  log "downloading Bakta light DB (~2 GB)"
  if command -v bakta_db >/dev/null; then
    bakta_db download --output "$t" --type light
    rm -rf "$RES/db-light"; mv "$t"/db-light "$RES/db-light"
  else
    curl -fL "$BAKTA_ZENODO" -o "$t/db.tar.xz"
    tar -xJf "$t/db.tar.xz" -C "$RES"
  fi
  if command -v amrfinder_update >/dev/null; then
    amrfinder_update --force_update --database "$RES/db-light/amrfinderplus-db" || true
  fi
  rm -rf "$t"
}

for d in $WHICH; do
  case "$d" in
    card) setup_card ;; resfinder) setup_resfinder ;; busco) setup_busco ;; bakta) setup_bakta ;;
    *) log "unknown database '$d' (card|resfinder|busco|bakta)"; exit 1 ;;
  esac
done
log "done. Default paths match nextflow.config; override with --card_json/--resfinder_db/--busco_db/--bakta_db."
