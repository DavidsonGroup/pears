#!/usr/bin/env python3
"""
Turn per-fusion lists of read IDs into fusion calls, using the run-wide
barcode table from the demultiplexing step (see modules/demultiplex.nf).

Every read is demultiplexed once, before alignment, so the fusion callers do
not demultiplex anything themselves: they report which reads support which
fusion, and the cell barcode and UMI are looked up here. All fusions are done
in a single pass over the barcode table - one scan for the whole run rather
than one scan, and one scheduler job, per fusion.

Input files are the read ID lists written by the fusion callers, named
<fusion>_<chrom1>_<base1>_<chrom2>_<base2>_<tool>_read_ids.txt. The output is
the cell_barcode, molecular_barcode, fusion table that combine_fusions.py
consumes, so no separate formatting step is needed.
"""

import argparse
import csv
import os
import re
import sys

# Flexiplex rewrites read IDs: it appends the strand it matched on (_+ / _-),
# may prepend "<barcode>_<umi>#", and may add the FLAMES-style _1of2 counter.
STRAND_SUFFIX = re.compile(r"_[+-]([0-9]+of[0-9]+)?(_C)?$")


def normalise_read_id(read_id):
    """Recover the original read ID from a flexiplex-rewritten one."""
    read_id = read_id.strip()
    if not read_id:
        return ""
    if read_id.startswith("@"):
        read_id = read_id[1:]
    read_id = read_id.split()[0]
    if "#" in read_id:
        read_id = read_id.rsplit("#", 1)[1]
    return STRAND_SUFFIX.sub("", read_id)


def fusion_name(path):
    """The fusion is the first field of the read ID file name."""
    return os.path.basename(path).split("_")[0]


def load_read_ids(path):
    ids = set()
    with open(path) as fh:
        for line in fh:
            read_id = normalise_read_id(line)
            if read_id:
                ids.add(read_id)
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Build fusion calls from read IDs and the barcode table."
    )
    parser.add_argument("--barcodes", required=True,
                        help="Run-wide flexiplex reads_barcodes.txt table.")
    parser.add_argument("--output", required=True,
                        help="Output CSV: cell_barcode, molecular_barcode, fusion.")
    parser.add_argument("read_ids", nargs="+",
                        help="Read ID lists, one file per fusion.")
    args = parser.parse_args()

    # read ID -> the fusions it supports, so the table is only scanned once
    wanted = {}
    fusions = []
    for path in args.read_ids:
        fusion = fusion_name(path)
        fusions.append(fusion)
        for read_id in load_read_ids(path):
            wanted.setdefault(read_id, []).append(fusion)
    print(f"Looking up {len(wanted)} read IDs from {len(set(fusions))} fusions "
          f"in {args.barcodes}", file=sys.stderr)

    seen = set()
    rows = []
    n_found = 0
    with open(args.barcodes) as table:
        header = table.readline()
        if header and not header.startswith("Read\t"):
            table.seek(0)
        for line in table:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            hits = wanted.get(fields[0])
            if not hits:
                continue
            n_found += 1
            for fusion in hits:
                row = (fields[1], fields[4], fusion)
                if row not in seen:
                    seen.add(row)
                    rows.append(row)

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["cell_barcode", "molecular_barcode", "fusion"])
        writer.writerows(rows)

    print(f"{n_found} of those reads had a barcode; wrote {len(rows)} unique "
          f"cell/UMI/fusion rows to {args.output}", file=sys.stderr)
    if wanted and not n_found:
        print("WARNING: no read ID matched the barcode table - check that the "
              "read IDs have not been rewritten by the fusion caller",
              file=sys.stderr)


if __name__ == "__main__":
    main()
