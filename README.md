# PEARS

**P**ipeline for g**e**ne-fusion se**a**rching in **R**na **S**ingle-cell sequences

PEARS is a Nextflow DSL2 pipeline that detects gene fusions at single-cell resolution from 10X scRNA-seq data. It combines three complementary fusion-calling approaches — [FUSCIA](https://github.com/ding-lab/fuscia), [Flexiplex](https://github.com/DavidsonGroup/flexiplex), and [Arriba](https://github.com/suhrig/arriba) — and assigns cell barcodes to each detected fusion event, producing per-cell fusion calls.

## Pipeline overview

1. **Reference preparation** — Downloads genome FASTA and GTF annotation (or uses pre-built references).
2. **Fusion target generation** — Builds search targets from a known fusions list using the reference annotation.
3. **Alignment** — Aligns reads with STARsolo (chimeric-aware) and produces a BAM and single-cell count matrix.
4. **Fusion detection** - Calls fusions using [FUSCIA](https://github.com/ding-lab/fuscia), [Flexiplex](https://github.com/DavidsonGroup/flexiplex), and [Arriba](https://github.com/suhrig/arriba) in parallel.
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
  -profile "local" \
  -resume
```

Running on SLURM cluster:

```bash
nextflow run DavidsonLab/pears \
  --fastq_r1 "/path/to/Reads_R1.fastq.gz" \
  --fastq_r2 "/path/to/Reads_R2.fastq.gz" \
  --known_fusions_list "known_fusions.csv" \
  --protocol "10x-3prime-v3" \
  --genome_version "GRCh38+GENCODE44" \
  -profile "slurm" \
  -resume
```

Running on SLURM will be much faster due to parallelization of the fusion detection step. The `-resume` flag allows you to continue from the last successful step if the pipeline is interrupted.

## Reusing references

> [!TIP]
> After running the pipeline for the first time, you can extract the references and index from the `references/` and `STAR_index/` folders in the output directory. Note that the STAR index is built with a specific `--sjdbOverhang` parameter according to the [read length of the input FASTQ](https://www.biostars.org/p/93883/). If you want to reuse the same index, make sure that the read lengths of the new FASTQ files are the same as the original ones. If they differ, you should rebuild the STAR index with the appropriate `--sjdbOverhang` (read length - 1) for optimal performance.

See [Pre-built or reusing reference overrides](#pre-built-or-reusing-reference-overrides) for the relevant arguments.

## Arguments

### Basic

| Argument | Default | Description |
|---|---|---|
| `--fastq_r1` | — | Glob pattern or path to Read 1 FASTQ files (gzipped). |
| `--fastq_r2` | — | Glob pattern or path to Read 2 FASTQ files (gzipped). |
| `--known_fusions_list` | — | CSV file of known/candidate fusions to search for (see [Known fusions list format](#known-fusions-list-format)). |
| `--protocol` | — | 10x Chromium chemistry preset (see [Protocol presets](#protocol-presets)). Sets the barcode whitelist and UMI length automatically. |
| `--genome_version` | `GRCh38+GENCODE44` | Genome build to download. Available versions: `GRCh38+GENCODE40` through `GRCh38+GENCODE49`. |
| `--out_dir` | `pears_output` | Directory for all pipeline outputs. |
| `-profile` | — | Execution environment: `local` or `slurm`. |

### Protocol presets

`--protocol` sets the barcode whitelist and UMI length for the given 10x chemistry. These values can be individually overridden with `--barcode_include_list` and `--umi_len` (see [Read structure overrides](#read-structure-overrides)).

| Preset | Chemistry | UMI length | Barcode whitelist |
|---|---|---|---|
| `10x-3prime-v2` | 3' Gene Expression v2 | 10 bp | 737K-august-2016 |
| `10x-3prime-v3` | 3' Gene Expression v3/v3.1 | 12 bp | 3M-february-2018 |
| `10x-3prime-v4` | 3' Gene Expression v4 | 12 bp | 3M-3pgex-may-2023 |
| `10x-5prime-v2` | 5' Gene Expression v1/v2 | 10 bp | 737K-august-2016 |
| `10x-5prime-v3` | 5' Gene Expression v3 | 12 bp | 3M-5pgex-jan-2023 |

### Read structure overrides

`--barcode_include_list` and `--umi_len` override the corresponding values set by `--protocol`. If `--protocol` is omitted entirely, **both** must be provided.

| Argument | Default | Overrides |
|---|---|---|
| `--barcode_include_list` | *set by `--protocol`* | Barcode whitelist. Path to a custom whitelist file (can be gzipped). |
| `--umi_len` | *set by `--protocol`* | UMI length in bases. |

### Pre-built or reusing reference overrides

By default, the pipeline downloads the genome specified by `--genome_version` and builds the STAR index automatically. To skip this, provide **all three** arguments below — `--genome_version` is then ignored.

| Argument | Default | Overrides |
|---|---|---|
| `--ref_fasta` | *downloaded* | Genome FASTA. Must be provided together with `--ref_gtf` and `--star_genome_index`. |
| `--ref_gtf` | *downloaded* | GTF annotation. Must be provided together with `--ref_fasta` and `--star_genome_index`. |
| `--star_genome_index` | *built from download* | STAR genome index directory. Must be provided together with `--ref_fasta` and `--ref_gtf`. |

### Tool parameters

| Argument | Default | Description |
|---|---|---|
| `--flexiplex_searchlen` | `20` | Length of fusion junction sequence to search for (2x actual overlap). |
| `--flexiplex_demultiplex_options` | *auto-generated* | Flexiplex demultiplexing options string. When not set, auto-generated as `-b "?{barcode_len}" -u "?{umi_len}" -e 1 -f 0` where barcode length is read from the whitelist file and UMI length comes from `--protocol` or `--umi_len`. Setting this explicitly overrides the auto-generated value. |
| `--fuscia_mapqual` | `30` | Minimum mapping quality for FUSCIA read extraction. |
| `--fuscia_up` | `1000` | Upstream search distance (bp) when no gene annotation is available. |
| `--fuscia_down` | `1000` | Downstream search distance (bp) when no gene annotation is available. |

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
| `fuscia_fusion_calls.csv` | Per-cell fusion calls from FUSCIA. |
| `flexiplex_fusion_calls.csv` | Per-cell fusion calls from Flexiplex. |
| `arriba_fusion_calls.csv` | Per-cell fusion calls from Arriba. |
| `STARsolo/` | BAM alignment, index, and single-cell count matrix. |
| `fuscia_out/` | Per-fusion FUSCIA discordant read files. |
| `flexiplex_out/` | Per-fusion Flexiplex barcode files. |
| `arriba_out/` | Arriba fusion table and per-fusion barcode files. |
| `fusion_targets.csv` | Generated fusion target coordinates and sequences. |
| `nextflow_report.html` | Nextflow execution report. |
| `nextflow_trace.txt` | Nextflow process trace log. |

### Fusion calls CSV format

The three fusion calls CSVs (`fuscia_fusion_calls.csv`, `flexiplex_fusion_calls.csv`, `arriba_fusion_calls.csv`) share the same format:

| Column | Description |
|---|---|
| `cell_barcode` | 10x cell barcode (CB tag) identifying the single cell. The trailing `-1` suffix is stripped. |
| `molecular_barcode` | Unique Molecular Identifier (UMI / UB tag) distinguishing distinct RNA molecules from PCR duplicates within the same cell. |
| `fusion` | Detected gene fusion, formatted as `GENE1--GENE2` (names taken from the `--known_fusions_list` input). |

Example:

```
cell_barcode,molecular_barcode,fusion
CCACAAAAGGTTCTTG,CAGGGATCAGTA,JPH1--NCOA2
CAGGGCTCACTTGGGC,TGATAGGAATCG,JPH1--NCOA2
GTGTGGCGTGGCCCAT,GGTAATCAGCAA,KIAA1429--RP11-586K2.1
```

Each row represents one unique observation of a fusion transcript in a specific cell. Rows are deduplicated — each (cell_barcode, molecular_barcode, fusion) triple appears at most once. Multiple rows with different cell barcodes for the same fusion indicate independent cells harbouring that fusion. Multiple rows with the same cell barcode but different UMIs for the same fusion indicate multiple distinct fusion transcript molecules captured in that cell, providing stronger evidence. Fusions detected by more than one tool in the same cell are higher confidence and can be identified by cross-referencing the three CSVs.

## Credits

Adapted from [FUSCIA](https://github.com/ding-lab/fuscia) (Steven Foltz, 2019) and [Flexiplex](https://github.com/DavidsonGroup/flexiplex) (Davidson et al., 2022).
