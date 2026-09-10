#!/usr/bin/env python3
"""
Pull a named subset of reads out of one or more FASTQ files.

Used to cut R1 down to the reads that actually need a barcode before running
flexiplex: the reads over the fusion breakpoints plus the reads arriba and
flexiplex called as fusion supporting. Demultiplexing the whole library is
dominated by flexiplex's fallback for reads that do not match the barcode list
exactly, which costs an edit distance against every barcode in the list, so
restricting the read set is worth far more than adding threads.

Reads are matched on the read ID, i.e. everything up to the first whitespace
in the header, with any leading '@' ignored.
"""

import argparse
import gzip
import sys


def load_read_ids(path):
    ids = set()
    with open(path, "rb") as fh:
        for line in fh:
            read_id = line.strip()
            if not read_id:
                continue
            if read_id.startswith(b"@"):
                read_id = read_id[1:]
            ids.add(read_id.split(None, 1)[0])
    return ids


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def main():
    parser = argparse.ArgumentParser(
        description="Write out the FASTQ records whose read ID is in a list."
    )
    parser.add_argument("--reads", required=True,
                        help="File of read IDs, one per line.")
    parser.add_argument("--output", help="Output FASTQ (default: stdout).")
    parser.add_argument("fastq", nargs="+", help="Input FASTQ file(s), .gz allowed.")
    args = parser.parse_args()

    ids = load_read_ids(args.reads)
    print(f"Looking for {len(ids)} reads in {len(args.fastq)} FASTQ file(s)",
          file=sys.stderr)

    out = open(args.output, "wb") if args.output else sys.stdout.buffer
    n_reads = n_written = 0
    try:
        for path in args.fastq:
            with open_fastq(path) as fh:
                for header in fh:
                    try:
                        seq, plus, qual = next(fh), next(fh), next(fh)
                    except StopIteration:
                        sys.exit(f"ERROR: truncated FASTQ record in {path}")
                    n_reads += 1
                    if header[1:].split(None, 1)[0] in ids:
                        out.write(header + seq + plus + qual)
                        n_written += 1
    finally:
        if args.output:
            out.close()

    print(f"Wrote {n_written} of {n_reads} reads", file=sys.stderr)
    if n_written == 0:
        print("WARNING: no reads matched - check that the read IDs in the BAM "
              "and the FASTQ headers agree", file=sys.stderr)


if __name__ == "__main__":
    main()
