#!/usr/bin/env python3
"""
Transfer cell barcodes and UMIs from the pipeline-wide flexiplex demultiplexing
table onto a STAR BAM file.

Reads are demultiplexed from R1 before alignment (modules/demultiplex.nf), STAR
aligns R2 only, and this script writes the barcode back onto the alignment
records as the CB and UB tags that fuscia expects. Every alignment record of a
read is tagged, chimeric and supplementary records included, which is what
fusion calling needs and what aligner-side demultiplexing did not give us
(see modules/demultiplex.nf).

By default only reads overlapping the fusion target regions are tagged, which
keeps the read-ID lookup table small. Use --all to tag every read in the BAM
(the lookup table is then held in memory for the whole run, which for a full
10x library needs tens of GB of RAM).
"""

import argparse
import sys

import pandas as pd
import pysam


def target_regions(targets_csv, pad):
    """Yield (contig, start, end) for both breakpoints of every fusion target."""
    df = pd.read_csv(targets_csv)
    regions = []
    for _, row in df.iterrows():
        # gene1/gene2 hold the far edge of the gene body (or base +/- the
        # fuscia_up / fuscia_down fallback) and base1/base2 the breakpoint, so
        # each pair spans the same interval fuscia fetches - see runFuscia.
        for chrom, a, b in ((row["chrom1"], row["gene1"], row["base1"]),
                            (row["chrom2"], row["gene2"], row["base2"])):
            try:
                start = min(int(a), int(b)) - pad
                end = max(int(a), int(b)) + pad
            except (TypeError, ValueError):
                print(f"WARNING: could not read a region from "
                      f"{row.get('fusion_genes', '?')} ({chrom}, {a}, {b}) - "
                      "reads over this breakpoint will not be tagged",
                      file=sys.stderr)
                continue
            regions.append((str(chrom), max(0, start), end))
    return regions


def resolve_contig(contig, bam_contigs):
    """Match a contig name against the BAM header, allowing for a chr prefix."""
    if contig in bam_contigs:
        return contig
    alt = contig[3:] if contig.startswith("chr") else "chr" + contig
    if alt in bam_contigs:
        return alt
    return None


def names_in_regions(bam_path, regions):
    """Collect the query names of all reads overlapping the given regions."""
    names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_contigs = set(bam.references)
        for contig, start, end in regions:
            resolved = resolve_contig(contig, bam_contigs)
            if resolved is None:
                print(f"WARNING: contig {contig} not in BAM header, skipping",
                      file=sys.stderr)
                continue
            end = min(end, bam.get_reference_length(resolved))
            for read in bam.fetch(resolved, max(0, start), end):
                names.add(read.query_name)
    return names


def load_barcodes(table_path, keep_names=None):
    """Read the flexiplex barcode table into {read_id: (barcode, umi)}."""
    barcodes = {}
    with open(table_path) as fh:
        for line in fh:
            if line.startswith("Read\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            read_id = fields[0]
            if keep_names is not None and read_id not in keep_names:
                continue
            barcodes[read_id] = (fields[1], fields[4])
    return barcodes


def main():
    parser = argparse.ArgumentParser(
        description="Write CB/UB tags onto a BAM from a flexiplex barcode table."
    )
    parser.add_argument("--bam", required=True, help="Input (indexed) BAM file.")
    parser.add_argument("--barcodes", required=True,
                        help="Pipeline-wide flexiplex reads_barcodes.txt table.")
    parser.add_argument("--output", required=True, help="Output BAM file.")
    parser.add_argument("--targets",
                        help="fusion_targets.csv - only reads overlapping these "
                             "regions are tagged.")
    parser.add_argument("--pad", type=int, default=1000,
                        help="Padding (bp) added either side of each target region.")
    parser.add_argument("--all", action="store_true",
                        help="Tag every read in the BAM (memory hungry).")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for BAM compression/decompression.")
    args = parser.parse_args()

    if not args.all and not args.targets:
        parser.error("either --targets or --all is required")

    keep_names = None
    if not args.all:
        regions = target_regions(args.targets, args.pad)
        keep_names = names_in_regions(args.bam, regions)
        print(f"{len(keep_names)} reads in {len(regions)} target regions",
              file=sys.stderr)

    barcodes = load_barcodes(args.barcodes, keep_names)
    print(f"{len(barcodes)} of these reads have a barcode assignment",
          file=sys.stderr)

    n_records = n_tagged = 0
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam_in:
        with pysam.AlignmentFile(args.output, "wb", template=bam_in,
                                 threads=args.threads) as bam_out:
            for read in bam_in:
                n_records += 1
                hit = barcodes.get(read.query_name)
                if hit is not None:
                    read.set_tag("CB", hit[0], value_type="Z")
                    read.set_tag("UB", hit[1], value_type="Z")
                    n_tagged += 1
                bam_out.write(read)

    print(f"Tagged {n_tagged} of {n_records} alignment records", file=sys.stderr)


if __name__ == "__main__":
    main()
