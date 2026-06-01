---
layout: default
---

# PEARS

**P**ipeline for g**e**ne-fusion se**a**rching in **R**na **S**ingle-cell sequences

PEARS is a Nextflow DSL2 pipeline that detects gene fusions at single-cell resolution from 10x scRNA-seq and Visium HD spatial transcriptomics data. It combines three complementary fusion-calling approaches — [FUSCIA](https://github.com/ding-lab/fuscia), [Flexiplex](https://github.com/DavidsonGroup/flexiplex), and [Arriba](https://github.com/suhrig/arriba) — and assigns cell barcodes to each detected fusion event, producing per-cell fusion calls.

## Pipeline overview

1. **Reference preparation** — Downloads genome FASTA and GTF annotation (or uses pre-built references).
2. **Fusion target generation** — Builds search targets from a known fusions list using the reference annotation.
3. **Alignment** — Aligns reads with STARsolo (chimeric-aware) and produces a BAM and single-cell count matrix.
4. **Fusion detection** — Calls fusions using FUSCIA, Flexiplex, and Arriba in parallel.
5. **Formatting** — Consolidates results into per-tool CSVs and a combined `combined_fusions.csv`. For Visium HD data, spatial bin barcodes are written to `combined_fusions_spatial.csv`.

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.10)
- [Conda](https://docs.conda.io/) (environments are built automatically)

## Installation

```bash
git clone https://github.com/DavidsonGroup/pears.git
```

## Usage

Running locally:

```bash
nextflow run /path/to/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-v3" \
  --genome_version "GRCh38+GENCODE44" \
  --out_dir "pears_output" \
  -profile "local" \
  -resume
```

Running on a SLURM cluster (recommended for large datasets):

```bash
nextflow run /path/to/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-v3" \
  --genome_version "GRCh38+GENCODE44" \
  --out_dir "pears_output" \
  -profile "slurm" \
  -resume
```

Running on Visium HD spatial transcriptomics data:

```bash
nextflow run /path/to/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-visiumHD" \
  --genome_version "GRCh38+GENCODE44" \
  --out_dir "pears_output" \
  -profile "slurm" \
  -resume
```

The `-resume` flag allows you to continue from the last successful step if the pipeline is interrupted.

> **Tip:** Nextflow also supports running directly from GitHub without cloning first: `nextflow run DavidsonGroup/pears [arguments]`.

## Protocol presets

`--protocol` sets the barcode whitelist and UMI length for the given 10x chemistry.

| Preset | Chemistry | UMI length |
|---|---|---|
| `10x-3prime-v2` | 3' Gene Expression v2 | 10 bp |
| `10x-3prime-v3` | 3' Gene Expression v3/v3.1 | 12 bp |
| `10x-3prime-v4` | 3' Gene Expression v4 | 12 bp |
| `10x-5prime-v2` | 5' Gene Expression v1/v2 | 10 bp |
| `10x-5prime-v3` | 5' Gene Expression v3 | 12 bp |
| `10x-3prime-visiumHD` | Visium HD spatial transcriptomics | 9 bp |

## Key arguments

| Argument | Default | Description |
|---|---|---|
| `--fastq_r1` | — | Path to Read 1 FASTQ file(s) (gzipped). |
| `--fastq_r2` | — | Path to Read 2 FASTQ file(s) (gzipped). |
| `--known_fusions_list` | — | CSV of known/candidate fusions to search for. |
| `--protocol` | — | 10x chemistry preset (see table above). |
| `--genome_version` | `GRCh38+GENCODE44` | Reference genome to download. |
| `--out_dir` | `pears_output` | Output directory. |
| `--discover_fusions` | `false` | Discover novel fusions via Arriba in addition to known fusions. |
| `--min_arriba_support` | `20000` | Minimum reads for a novel Arriba fusion to be included. |
| `--arriba_exclusion_file` | — | Path to a gzipped Arriba blacklist (`.tsv.gz`) to filter likely false positives. See [Arriba releases](https://github.com/suhrig/arriba/releases) for pre-built blacklists. |
| `--visium_bin_size` | `8` | *(Visium HD only)* Bin size in microns (`2`, `8`, or `16`) for spatial barcode conversion. |
| `--cpus` | `16` | CPUs per process. |
| `--memory` | `128 GB` | Memory per process. |
| `--time` | `48h` | Wall-time limit per process. |
| `-profile` | — | `local` or `slurm`. |

For the full argument reference see the [README](https://github.com/DavidsonGroup/pears/blob/main/README.md).

## Known fusions list format

The `--known_fusions_list` input is a CSV with the following columns:

| Column | Description |
|---|---|
| `fusion genes` | Fusion pair separated by `--` (e.g. `BCAS4--BCAS3`). |
| `chrom1` | Chromosome of gene 1. |
| `base1` | Breakpoint position of gene 1. |
| `strand1` | Strand of gene 1 (`+` or `-`). |
| `chrom2` | Chromosome of gene 2. |
| `base2` | Breakpoint position of gene 2. |
| `strand2` | Strand of gene 2 (`+` or `-`). |

This format is compatible with [JAFFA](https://github.com/Oshlack/JAFFA) output. Additional columns are ignored.

```
fusion genes,chrom1,base1,strand1,chrom2,base2,strand2,classification
BCAS4--BCAS3,chr20,50795173,+,chr17,61368327,+,HighConfidence
RPS6KB1--VMP1,chr17,59914703,+,chr17,59839768,+,HighConfidence
```

## Output

Results are written to `--out_dir`:

| File | Description |
|---|---|
| `combined_fusions.csv` | All tools merged: UMI counts per tool and total per (fusion, cell) pair. |
| `combined_fusions_spatial.csv` | *(Visium HD only)* Combined fusions with SpaceRanger-format spatial barcodes (e.g. `s_008um_00241_00258-1`). |
| `fuscia_fusion_calls.csv` | Per-cell fusion calls from FUSCIA. |
| `flexiplex_fusion_calls.csv` | Per-cell fusion calls from Flexiplex. |
| `arriba_fusion_calls.csv` | Per-cell fusion calls from Arriba. |
| `STARsolo/` | BAM alignment and single-cell count matrix. |
| `arriba_out/` | Arriba fusion table and per-fusion barcode files. |
| `fusion_targets.csv` | Generated fusion target sequences. |
| `nextflow_report.html` | Nextflow execution report. |

## Credits

Adapted from [FUSCIA](https://github.com/ding-lab/fuscia) (Steven Foltz, 2019) and [Flexiplex](https://github.com/DavidsonGroup/flexiplex) (Davidson et al., 2022).
