#!/usr/bin/env python3
"""
Look up cell barcodes / UMIs for sets of read IDs in the pipeline-wide
flexiplex demultiplexing table (see modules/demultiplex.nf).

This replaces the per-fusion flexiplex demultiplexing that used to be run on
the fusion-supporting reads: every read is now demultiplexed once, before
alignment, and downstream steps only join on the read ID.

All fusions are handled in a single pass over the table, which is large: one
scan for the whole run rather than one scan (and one scheduler job) per fusion.

Each input file is a list of read IDs named <fusion>_read_ids.txt; the matching
rows are written to barcodes_<fusion>_reads_barcodes.txt. The output has the
same columns as a flexiplex reads_barcodes.txt table (Read, CellBarcode,
FlankEditDist, BarcodeEditDist, UMI) so that format_barcodes.py can consume it
unchanged.
"""

import argparse
import os
import sys

HEADER = "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"
IN_SUFFIX = "_read_ids.txt"
OUT_SUFFIX = "_reads_barcodes.txt"
OUT_PREFIX = "barcodes_"


def output_name(read_ids_path):
    """barcodes_<fusion>_reads_barcodes.txt for a <fusion>_read_ids.txt input."""
    base = os.path.basename(read_ids_path)
    if base.endswith(IN_SUFFIX):
        base = base[:-len(IN_SUFFIX)]
    else:
        base = os.path.splitext(base)[0]
    return OUT_PREFIX + base + OUT_SUFFIX


def load_read_ids(path):
    """Read IDs, one per line. Leading '@' and any trailing comment are stripped."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            read_id = line.strip()
            if not read_id:
                continue
            if read_id.startswith("@"):
                read_id = read_id[1:]
            ids.add(read_id.split()[0])
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Extract barcode table rows for one or more lists of read IDs."
    )
    parser.add_argument("--barcodes", required=True,
                        help="Pipeline-wide flexiplex reads_barcodes.txt table.")
    parser.add_argument("read_ids", nargs="+",
                        help="Files of read IDs, one per line, one file per fusion.")
    args = parser.parse_args()

    # read ID -> the fusions wanting it, so the table is only scanned once
    wanted = {}
    for i, path in enumerate(args.read_ids):
        for read_id in load_read_ids(path):
            wanted.setdefault(read_id, []).append(i)
    print(f"Looking up {len(wanted)} read IDs from {len(args.read_ids)} "
          f"fusions in {args.barcodes}", file=sys.stderr)

    out_paths = [output_name(p) for p in args.read_ids]
    counts = [0] * len(args.read_ids)
    handles = [open(p, "w") for p in out_paths]
    try:
        for fh in handles:
            fh.write(HEADER + "\n")

        with open(args.barcodes) as table:
            header = table.readline()
            if header and not header.startswith("Read\t"):
                # No header line in the table - treat the first line as data.
                table.seek(0)
            for line in table:
                read_id = line.split("\t", 1)[0]
                for i in wanted.get(read_id, ()):
                    handles[i].write(line)
                    counts[i] += 1
    finally:
        for fh in handles:
            fh.close()

    for path, count in zip(out_paths, counts):
        print(f"  {count}\t{path}", file=sys.stderr)
    print(f"Wrote {sum(counts)} barcode assignments", file=sys.stderr)


if __name__ == "__main__":
    main()
