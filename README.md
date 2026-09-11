# blaCTX Pipeline

Nextflow workflow for reproducible *Klebsiella pneumoniae* whole-genome
analysis, with emphasis on assembly quality control and detection of
extended-spectrum beta-lactamase genes, including **blaCTX-M-15**.

## Overview

The pipeline automates raw read processing, assembly, quality filtering, and dual-tool AMR screening:
1. **Fetch & Preprocess:** Retrieves SRA/ENA reads and performs QC/trimming (FastQC, Trimmomatic).
2. **Assemble & QC:** Assembles genomes (Shovill/SKESA) and evaluates quality (BUSCO, read-depth coverage).
3. **Filter:** Filters assemblies based on completeness, coverage, and contig thresholds.
4. **AMR Screening:** Screens passing assemblies using **RGI/CARD** and **ResFinder**.
5. **Report:** Generates a unified cohort-level metadata summary (`TSV` / `XLSX`).

## Prerequisites

- **OS:** Linux / WSL2
- **Workflow Manager:** Nextflow (>= 22.10)
- **Environment Manager:** Conda / Mamba

## Directory structure

```
.
├── main.nf                    # Nextflow DSL2 entry point
├── nextflow.config            # Workflow parameters & execution settings 
├── config/
│   ├── config.yaml            # Project thresholds & DB paths
│   └── samples.tsv            # Sample manifest
├── workflow/
│   ├── envs/                  # Conda environment YAMLs
│   ├── scripts/               # Utility Python scripts
├── resources/                 # Local reference DBs
```

## Setup & Credentials 

Before running, export your NCBI credentials for metadata fetching: 

```bash
export NCBI_EMAIL="you@example.org"
export NCBI_API_KEY="replace-with-a-secret"
```

## Quick Start 

Execute the complete workflow 

```bash
nextflow run main.nf \
  --samples config/samples.tsv \
  --results results \
  --threads 8 \
  --target_gene blaCTX-M-15 \
  -with-conda \
  -resume
```

Note: To inspect standalone Python script contracts, run python workflow/scripts/<script_name>.py --help.

## Outputs

| Directory | Contents |
|---|---|
| `results/00_raw_reads/` | Raw FASTQ files |
| `results/01_qc/fastqc_raw/` | Raw-read FastQC reports |
| `results/01_qc/trimmed/` | Paired trimmed reads |
| `results/02_assembly/` | Per-sample assemblies and coverage files |
| `results/03_busco/` | BUSCO short summaries |
| `results/04_filtered_assemblies/` | Filter status JSON and promoted assemblies |
| `results/05_amr/rgi/` | RGI tabular and JSON outputs |
| `results/05_amr/resfinder/` | ResFinder outputs |
| `results/06_metadata/` | Cohort TSV and XLSX summary |


## Citation & License 

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

Please cite Nextflow, Shovill, BUSCO, RGI/CARD, and ResFinder when publishing results.

