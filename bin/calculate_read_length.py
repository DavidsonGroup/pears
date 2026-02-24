#!/usr/bin/env python3
"""
Calculate the mode of R2 read lengths from FASTQ files.
Emits a warning if >10% of reads differ from the mode.

Usage:
    calculate_read_length.py <fastq_file> [<fastq_file> ...]

Output:
    Prints the mode read length to stdout.
    Warnings are printed to stderr.
"""

import gzip
import sys
from collections import Counter
from pathlib import Path


def open_fastq(filepath):
    """Open a FASTQ file, handling gzip if needed."""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def get_read_lengths(fastq_files, max_reads=10000):
    """Extract read lengths from the first max_reads of FASTQ files."""
    lengths = []
    reads_counted = 0

    for fq in sorted(fastq_files):
        if reads_counted >= max_reads:
            break
        with open_fastq(fq) as f:
            line_num = 0
            for line in f:
                line_num += 1
                # Sequence is on line 2 of each 4-line FASTQ record
                if line_num % 4 == 2:
                    lengths.append(len(line.strip()))
                    reads_counted += 1
                    if reads_counted >= max_reads:
                        break
    return lengths


def main():
    if len(sys.argv) < 2:
        print("Usage: calculate_read_length.py <fastq_file> [<fastq_file> ...]", file=sys.stderr)
        sys.exit(1)

    fastq_files = [Path(f) for f in sys.argv[1:]]

    # Validate files exist
    for fq in fastq_files:
        if not fq.exists():
            print(f"ERROR: File not found: {fq}", file=sys.stderr)
            sys.exit(1)

    # Get read lengths
    lengths = get_read_lengths(fastq_files, max_reads=10000)

    if not lengths:
        print("ERROR: Could not read any sequences from FASTQ files", file=sys.stderr)
        sys.exit(1)

    # Calculate mode
    length_counts = Counter(lengths)
    mode_length, mode_count = length_counts.most_common(1)[0]
    total_reads = len(lengths)

    # Calculate percentage of reads matching mode
    mode_percentage = (mode_count / total_reads) * 100
    non_mode_percentage = 100 - mode_percentage

    # Warn if >10% of reads differ from mode
    if non_mode_percentage > 10:
        print(f"WARNING: {non_mode_percentage:.1f}% of reads have length different from mode ({mode_length}bp)", file=sys.stderr)
        print(f"         Length distribution: {dict(length_counts.most_common(5))}", file=sys.stderr)

    print(f"R2 read length mode: {mode_length}bp ({mode_percentage:.1f}% of {total_reads} reads sampled)", file=sys.stderr)

    # Output the mode length to stdout
    print(mode_length)


if __name__ == "__main__":
    main()
