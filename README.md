# PEARS

**P**ipeline for g**e**ne-fusion se**a**rching in **R**na **S**ingle-cell sequences

PEARS is a Nextflow DSL2 pipeline that detects gene fusions at single-cell resolution from 10X scRNA-seq data. It combines three complementary fusion-calling approaches — [FUSCIA](https://github.com/sfoltz/fuscia), [Flexiplex](https://github.com/DavidsonGroup/flexiplex), and [Arriba](https://github.com/suhrig/arriba) — and assigns cell barcodes to each detected fusion event, producing per-cell fusion calls.

## Pipeline overview

1. **Reference preparation** — Downloads genome FASTA and GTF annotation (or uses pre-built references).
2. **Fusion target generation** — Builds search targets from a known fusions list using the reference annotation.
3. **Alignment** — Aligns reads with STARsolo (chimeric-aware) and produces a BAM and single-cell count matrix.
4. **Fusion detection** - Calls fusions using three tools in parallel:
   - **FUSCIA**
   - **Flexiplex**
   - **Arriba**
5. **Formatting** — Consolidates results into three CSV files of per-cell fusion calls (`fuscia_fusion_calls.csv`, `flexiplex_fusion_calls.csv`, `arriba_fusion_calls.csv`).

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.10)
- [Conda](https://docs.conda.io/) (environments are built automatically from `env/pears_env.yml`)

## Usage

Running locally:

```bash
nextflow run DavidsonLab/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-v3" \
  --genome_version "GRCh38+GENCODE44" \
  -profile "local"
```

Running on SLURM cluster:

```bash
nextflow run DavidsonLab/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-v3" \
  --genome_version "GRCh38+GENCODE44" \
  -profile "slurm"
```

## Arguments

### Input files (required)

| Argument | Default | Description |
|---|---|---|
| `--fastq_r1` | — | Glob pattern or path to Read 1 FASTQ files (gzipped). |
| `--fastq_r2` | — | Glob pattern or path to Read 2 FASTQ files (gzipped). |
| `--known_fusions_list` | — | CSV file of known/candidate fusions to search for (see [Known fusions list format](#known-fusions-list-format)). |

### Protocol and read structure

Either `--protocol` **or** both `--barcode_include_list` and `--umi_len` must be provided.

`--protocol` is the recommended way to configure the read structure — it automatically selects the correct barcode whitelist and UMI length for the given 10x chemistry. If `--barcode_include_list` or `--umi_len` are also provided, they **override** the corresponding value set by `--protocol`.

| Argument | Default | Description |
|---|---|---|
| `--protocol` | — | 10x Chromium chemistry preset (see table below). |
| `--barcode_include_list` | *set by `--protocol`* | Path to a custom barcode whitelist file (can be gzipped). **Overrides** the whitelist set by `--protocol`. |
| `--umi_len` | *set by `--protocol`* | UMI length in bases. **Overrides** the value set by `--protocol`. |

### Reference genome

Provide **all three** of `--ref_fasta`, `--ref_gtf`, and `--star_genome_index` to use pre-built references — this skips genome download and STAR index building. If any of the three is missing, the pipeline downloads the genome specified by `--genome_version` and builds the index automatically. In that case `--genome_version` is the only relevant argument; the other three are ignored.

| Argument | Default | Description |
|---|---|---|
| `--genome_version` | `GRCh38+GENCODE44` | Genome build to download. Available versions: `GRCh38+GENCODE40` through `GRCh38+GENCODE49`. **Ignored** when all three pre-built reference arguments are provided. |
| `--ref_fasta` | — | Path to a pre-built genome FASTA. Must be provided together with `--ref_gtf` and `--star_genome_index`. |
| `--ref_gtf` | — | Path to a pre-built GTF annotation. Must be provided together with `--ref_fasta` and `--star_genome_index`. |
| `--star_genome_index` | — | Path to a pre-built STAR genome index directory. Must be provided together with `--ref_fasta` and `--ref_gtf`. |

### Tool parameters

| Argument | Default | Description |
|---|---|---|
| `--flexiplex_searchlen` | `20` | Length of fusion junction sequence to search for (2x actual overlap). |
| `--flexiplex_demultiplex_options` | *auto-generated* | Flexiplex demultiplexing options string. When not set, auto-generated as `-b "?{barcode_len}" -u "?{umi_len}" -e 1 -f 0` where barcode length is read from the whitelist file and UMI length comes from `--protocol` or `--umi_len`. **Overrides** the auto-generated value when set explicitly. |
| `--fuscia_mapqual` | `30` | Minimum mapping quality for FUSCIA read extraction. |
| `--fuscia_up` | `1000` | Upstream search distance (bp) when no gene annotation is available. |
| `--fuscia_down` | `1000` | Downstream search distance (bp) when no gene annotation is available. |

### Output

| Argument | Default | Description |
|---|---|---|
| `--out_dir` | `pears_output` | Directory for all pipeline outputs. |

### Profiles

Use `-profile` to select the execution environment:

| Profile | Description |
|---|---|
| `local` | Run on the local machine. |
| `slurm` | Submit jobs to a SLURM cluster. |

## Known fusions list format

The `--known_fusions_list` input is a CSV with the following required columns:

| Column | Description |
|---|---|
| `fusion genes` | Fusion gene pair separated by `--` (e.g. `BCAS4--BCAS3`). |
| `chrom1` | Chromosome of gene 1 (e.g. `chr20`). |
| `base1` | Breakpoint position of gene 1. |
| `strand1` | Strand of gene 1 (`+` or `-`). |
| `chrom2` | Chromosome of gene 2 (e.g. `chr17`). |
| `base2` | Breakpoint position of gene 2. |
| `strand2` | Strand of gene 2 (`+` or `-`). |

Additional columns (e.g. `classification`) are ignored. This format is compatible with [JAFFA](https://github.com/Oshlack/JAFFA) output. Fusions involving mitochondrial genes (`MT-`) are automatically filtered out.

Example:

```
fusion genes,chrom1,base1,strand1,chrom2,base2,strand2,classification
BCAS4--BCAS3,chr20,50795173,+,chr17,61368327,+,HighConfidence
RPS6KB1--VMP1,chr17,59914703,+,chr17,59839768,+,HighConfidence
SLC25A24--NBPF6,chr1,108161182,-,chr1,108470597,+,HighConfidence
```

## Output

Results are written to `--out_dir` (default `pears_output/`):

| File/Directory | Description |
|---|---|
| `fuscia_fusion_calls.csv` | Per-cell fusion calls (cell barcode, UMI, fusion) from FUSCIA. |
| `flexiplex_fusion_calls.csv` | Per-cell fusion calls (cell barcode, UMI, fusion) from Flexiplex. |
| `arriba_fusion_calls.csv` | Per-cell fusion calls (cell barcode, UMI, fusion) from Arriba. |
| `STARsolo/` | BAM alignment, index, and single-cell count matrix. |
| `fuscia_out/` | Per-fusion FUSCIA discordant read files. |
| `flexiplex_out/` | Per-fusion Flexiplex barcode files. |
| `arriba_out/` | Arriba fusion table and per-fusion barcode files. |
| `fusion_targets.csv` | Generated fusion target coordinates and sequences. |
| `nextflow_report.html` | Nextflow execution report. |
| `nextflow_trace.txt` | Nextflow process trace log. |

## Credits

Adapted from [FUSCIA](https://github.com/sfoltz/fuscia) (Steven Foltz, 2019) and [Flexiplex](https://github.com/DavidsonGroup/flexiplex) (Davidson et al., 2022).
