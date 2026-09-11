#!/usr/bin/env python3
"""
parse_selection_results.py

Combines, per locus, the PAML codeml summary (site models + BEB codons) and
the HyPhy SLAC/MEME summary (per-codon omega + p-values) into one flat,
codon-resolved table — the format expected by downstream structural mapping
onto blaCTX-M-15 models (OmegaFold/FoldX).

A codon is called under positive selection if EITHER:
  - PAML BEB posterior >= threshold (from M2a/M8), or
  - HyPhy SLAC flags p_positive <= threshold, or
  - HyPhy MEME flags episodic positive selection (p_value <= threshold)

Usage:
    python parse_selection_results.py \
        --results-dir results/selection \
        --loci blaCTX-M-15 rpoB fusA ompK35 ompK36 \
        --out-csv results/selection/combined_selection_summary.csv \
        --out-json results/selection/combined_selection_summary.json
"""

import argparse
import csv
import json
import os


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def build_locus_table(locus: str, results_dir: str) -> list:
    codeml_path = os.path.join(results_dir, "selection", locus, "codeml", "summary.json")
    hyphy_path = os.path.join(results_dir, "selection", locus, "hyphy", "slac_meme_summary.json")

    codeml_data = load_json(codeml_path) or {}
    hyphy_data = load_json(hyphy_path) or {}

    beb_by_pos = {}
    for site in codeml_data.get("positively_selected_codons", []):
        beb_by_pos.setdefault(site["codon_position"], []).append(site)

    hyphy_by_pos = {s["codon_position"]: s for s in hyphy_data.get("per_codon", [])}

    all_positions = sorted(set(beb_by_pos) | set(hyphy_by_pos))
    rows = []
    for pos in all_positions:
        beb_hits = beb_by_pos.get(pos, [])
        hy = hyphy_by_pos.get(pos, {})

        paml_prob = max((h["posterior_probability"] for h in beb_hits), default=None)
        paml_model = ",".join(sorted({h["source_model"] for h in beb_hits})) or None
        paml_omega = next((h["omega_estimate"] for h in beb_hits if h.get("omega_estimate")), None)

        paml_positive = bool(beb_hits)
        slac_positive = hy.get("slac_classification") == "positive_selection"
        meme_positive = bool(hy.get("meme_episodic_positive_selection"))
        consensus = paml_positive or slac_positive or meme_positive

        n_methods_agreeing = sum([paml_positive, slac_positive, meme_positive])

        rows.append({
            "locus": locus,
            "codon_position": pos,
            "paml_beb_posterior": paml_prob,
            "paml_source_model": paml_model,
            "paml_omega_estimate": paml_omega,
            "paml_positive_selection": paml_positive,
            "slac_omega": hy.get("slac_omega"),
            "slac_p_positive": hy.get("slac_p_positive"),
            "slac_p_purifying": hy.get("slac_p_purifying"),
            "slac_classification": hy.get("slac_classification"),
            "meme_omega_episodic": hy.get("meme_omega_episodic"),
            "meme_p_value": hy.get("meme_p_value"),
            "meme_positive_selection": meme_positive,
            "consensus_positive_selection": consensus,
            "n_methods_agreeing": n_methods_agreeing,
            "selection_class": "positive" if consensus else "purifying_or_neutral",
        })

    # locus-level LRT results, attached once per locus for context
    lrt = codeml_data.get("LRTs", {})
    return rows, lrt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--loci", nargs="+", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    all_rows = []
    locus_lrts = {}
    for locus in args.loci:
        rows, lrt = build_locus_table(locus, args.results_dir)
        all_rows.extend(rows)
        locus_lrts[locus] = lrt

    output_dir = os.path.dirname(args.out_csv)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fieldnames = [
        "locus", "codon_position",
        "paml_beb_posterior", "paml_source_model", "paml_omega_estimate", "paml_positive_selection",
        "slac_omega", "slac_p_positive", "slac_p_purifying", "slac_classification",
        "meme_omega_episodic", "meme_p_value", "meme_positive_selection",
        "consensus_positive_selection", "n_methods_agreeing", "selection_class",
    ]
    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    n_positive = sum(1 for r in all_rows if r["consensus_positive_selection"])
    out_json = {
        "loci_analyzed": args.loci,
        "n_codons_total": len(all_rows),
        "n_codons_under_positive_selection": n_positive,
        "locus_level_LRTs": locus_lrts,
        "codons": all_rows,
    }
    with open(args.out_json, "w") as fh:
        json.dump(out_json, fh, indent=2)

    print(f"[parse_selection_results] {len(all_rows)} codon rows across "
          f"{len(args.loci)} loci -> {args.out_csv} / {args.out_json}")
    print(f"[parse_selection_results] {n_positive} codons flagged under "
          f"consensus positive selection")


if __name__ == "__main__":
    main()
