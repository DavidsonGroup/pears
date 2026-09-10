#!/usr/bin/env python3
"""
Assign a cell barcode and UMI to each read, writing a flexiplex-format
reads_barcodes.txt table.

The barcode and UMI sit at a fixed offset in R1, so they are read positionally
rather than searched for. Correction is done by enumerating the sequences one
error away from an observed barcode and looking those up in the whitelist,
rather than comparing each read against every whitelist barcode: the
neighbourhood of a 16bp barcode has at most 176 members whatever the size of
the whitelist, so the whitelist is touched twice by a streaming membership
test instead of once per read.

Correction prefers barcodes already seen in these reads. A barcode observed in
the data is a much more likely source than a whitelist barcode that appears
nowhere in it, so the whitelist is only consulted for observed sequences that
no other observed barcode explains. Where more than one candidate remains, the
read is assigned by combining how abundant each candidate is with how likely
the implied sequencing error is given the base quality, and is left unassigned
if that does not favour one candidate clearly enough.

With --barcode-list-out the barcode list is written as well (or instead), for
feeding flexiplex or for QC.
"""

import argparse
import gzip
import math
import sys
from collections import Counter

BASES = b"ACGT"
HEADER = "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"
# Posterior at which an ambiguous barcode is called rather than dropped
MIN_POSTERIOR = 0.9
# Probability assigned to an indel, which has no single quality score to read
INDEL_PRIOR = 1e-4


def open_maybe_gz(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def iter_reads(fastqs):
    """Yield (read_id, sequence, quality) for every record."""
    for path in fastqs:
        with open_maybe_gz(path) as fh:
            for header in fh:
                try:
                    seq, _, qual = next(fh), next(fh), next(fh)
                except StopIteration:
                    sys.exit(f"ERROR: truncated FASTQ record in {path}")
                yield header[1:].split(None, 1)[0], seq.rstrip(), qual.rstrip()


def neighbours(barcode):
    """Every barcode one error away from this observed sequence.

    The barcode is read at a fixed offset, so an indel shifts everything after
    it: a deletion pulls one base in from downstream and an insertion pushes
    the last base out. Covering all three error types means substituting each
    base, inserting a base (dropping the last observed one), and deleting a
    base (appending each possible base).
    """
    length = len(barcode)
    for i in range(length):
        for base in BASES:
            if barcode[i] != base:
                yield barcode[:i] + bytes([base]) + barcode[i + 1:]
            yield barcode[:i] + bytes([base]) + barcode[i:length - 1]
            yield barcode[:i] + barcode[i + 1:] + bytes([base])


def whitelist_members(whitelist_path, wanted):
    """The whitelist entries that are in `wanted`, in one streaming pass."""
    hits = set()
    if not wanted:
        return hits
    with open_maybe_gz(whitelist_path) as fh:
        for line in fh:
            barcode = line.split()[0] if line.strip() else b""
            if barcode in wanted:
                hits.add(barcode)
    return hits


def explain(observed, candidate, next_base):
    """How `candidate` could have produced this read, or None if it could not.

    Returns (substituted position or None, offset of the UMI). An indel inside
    the barcode shifts everything after it, so the UMI is not at a fixed
    offset: a base lost from the barcode pulls the UMI one base earlier, and a
    base gained pushes it one base later.
    """
    length = len(observed)
    if len(candidate) == length:
        diffs = [i for i, (a, b) in enumerate(zip(observed, candidate)) if a != b]
        if len(diffs) == 1:
            return diffs[0], length
    # a base was lost from the barcode: the candidate less one base is what we saw
    head = observed[:length - 1]
    for i in range(len(candidate)):
        if candidate[:i] + candidate[i + 1:] == head:
            return None, length - 1
    # a base was gained: drop one of the bases we saw and the next base follows
    if next_base:
        for i in range(length):
            if observed[:i] + observed[i + 1:] + next_base == candidate:
                return None, length + 1
    return None, None


def error_probability(substituted, quality):
    """How likely the implied sequencing error is."""
    if substituted is None:
        return INDEL_PRIOR
    phred = quality[substituted] - 33
    # a specific substitution is one of the three possible wrong bases
    return 10 ** (-phred / 10.0) / 3.0


def resolve(observed, candidates, quality, next_base, counts):
    """Pick between candidates, returning (barcode, umi_offset) or (None, None).

    Candidates come from a neighbourhood generated without the read's context,
    so each one is first checked against this read; the rest are weighted by
    how abundant they are and how likely the error they imply is.
    """
    scored = []
    for candidate in candidates:
        substituted, offset = explain(observed, candidate, next_base)
        if offset is None:
            continue
        scored.append(((counts.get(candidate, 0) + 1)
                       * error_probability(substituted, quality),
                       candidate, offset))
    if not scored:
        return None, None
    if len(scored) == 1:
        return scored[0][1], scored[0][2]
    total = sum(s for s, _, _ in scored)
    best = max(scored)
    if total <= 0 or best[0] / total < MIN_POSTERIOR:
        return None, None
    return best[1], best[2]


def main():
    parser = argparse.ArgumentParser(
        description="Demultiplex reads against a barcode whitelist."
    )
    parser.add_argument("--whitelist", required=True,
                        help="Barcode whitelist, one barcode per line.")
    parser.add_argument("--barcode-length", type=int, required=True,
                        help="Barcode length in bases, from the start of the read.")
    parser.add_argument("--umi-length", type=int,
                        help="UMI length in bases, immediately after the barcode. "
                             "Required with --table-out.")
    parser.add_argument("--table-out",
                        help="Write a flexiplex-format reads_barcodes.txt here.")
    parser.add_argument("--barcode-list-out",
                        help="Write the barcode list here.")
    parser.add_argument("--edit-distance", type=int, default=1, choices=[0, 1],
                        help="1 (default) corrects barcodes one error from a "
                             "known barcode; 0 keeps only exact matches.")
    parser.add_argument("fastq", nargs="+", help="Read file(s), .gz allowed.")
    args = parser.parse_args()

    if not args.table_out and not args.barcode_list_out:
        parser.error("at least one of --table-out or --barcode-list-out is required")
    if args.table_out and args.umi_length is None:
        parser.error("--table-out requires --umi-length")

    length = args.barcode_length

    # Pass 1: what barcodes are present, and how often
    counts = Counter()
    n_reads = 0
    for _, seq, _ in iter_reads(args.fastq):
        n_reads += 1
        if len(seq) >= length:
            counts[seq[:length]] += 1
    observed = set(counts)
    print(f"{n_reads} reads, {len(observed)} distinct barcode sequences",
          file=sys.stderr)

    exact = whitelist_members(args.whitelist, observed)
    print(f"{len(exact)} are in the whitelist exactly", file=sys.stderr)

    # Candidates for the sequences that are not, preferring barcodes we have
    # already seen over whitelist barcodes that appear nowhere in the data
    from_observed = {}
    from_whitelist = {}
    near = set()
    if args.edit_distance == 1:
        variants = set()
        unexplained = {}
        for barcode in observed - exact:
            candidates = set(neighbours(barcode))
            seen_here = sorted(candidates & exact)
            if seen_here:
                from_observed[barcode] = seen_here
            else:
                unexplained[barcode] = candidates
                variants |= candidates
        near = whitelist_members(args.whitelist, variants)
        for barcode, candidates in unexplained.items():
            hits = sorted(candidates & near)
            if hits:
                from_whitelist[barcode] = hits
        print(f"{len(observed) - len(exact)} are not: {len(from_observed)} are "
              f"one error from a barcode already seen, {len(from_whitelist)} "
              f"from {len(near)} whitelist barcodes, "
              f"{len(observed) - len(exact) - len(from_observed) - len(from_whitelist)} "
              f"have no candidate", file=sys.stderr)

    if args.barcode_list_out:
        final = sorted(exact | near)
        with open(args.barcode_list_out, "wb") as out:
            for barcode in final:
                out.write(barcode + b"\n")
        print(f"Wrote {len(final)} barcodes to {args.barcode_list_out}",
              file=sys.stderr)

    if not args.table_out:
        return

    # Pass 2: assign every read
    tally = Counter()
    with open(args.table_out, "w") as out:
        out.write(HEADER + "\n")
        for read_id, seq, qual in iter_reads(args.fastq):
            barcode = seq[:length]
            if barcode in exact:
                assigned, distance, umi_at = barcode, 0, length
                tally["exact"] += 1
            else:
                candidates = from_observed.get(barcode) or from_whitelist.get(barcode)
                if not candidates:
                    tally["no candidate"] += 1
                    continue
                assigned, umi_at = resolve(barcode, candidates, qual,
                                           seq[length:length + 1], counts)
                if assigned is None:
                    tally["ambiguous"] += 1
                    continue
                distance = 1
                tally["corrected"] += 1
            umi = seq[umi_at:umi_at + args.umi_length]
            out.write(f"{read_id.decode()}\t{assigned.decode()}\t0\t{distance}\t"
                      f"{umi.decode()}\n")

    assigned_total = tally["exact"] + tally["corrected"]
    print(f"Assigned {assigned_total} of {n_reads} reads "
          f"({tally['exact']} exact, {tally['corrected']} corrected); "
          f"{tally['ambiguous']} ambiguous, {tally['no candidate']} with no "
          f"candidate barcode", file=sys.stderr)
    if assigned_total == 0:
        sys.exit("ERROR: no read was assigned a barcode - check --protocol and "
                 "the read structure")


if __name__ == "__main__":
    main()
