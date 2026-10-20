# Sample manifest (`samples.tsv`)

Tab-separated, one row per sequencing run. TSV cannot hold comments, so the format is documented here.

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes* | Unique name used for all output files. *Falls back to `Run` if absent. |
| `sra_accession` | yes* | SRA/ENA run accession (`SRR`/`ERR`/`DRR`). *Falls back to `Run` if absent. |
| `source_type` | yes | `SRA_PE`/`PE` (paired-end, Shovill+SKESA) or `SRA_SE`/`SE` (single-end, `spades.py --isolate`). Mixed layouts are fine. |
| anything else | no | Carried through unchanged into `results/06_metadata/cohort_metadata_summary.{tsv,xlsx}`. |

The bundled example has 8 *K. pneumoniae* Illumina runs (5 PE, 3 SE), all checked against the ENA Portal API
(species, layout, 280-520 Mb each; the SE runs are ~80-95x). With 8 samples the phylogeny stage (needs >= 4 passing
genomes) also runs. A full SRA RunInfo export also works as a manifest, as long as the columns above are present.
