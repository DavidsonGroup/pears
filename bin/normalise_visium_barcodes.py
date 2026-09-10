#!/usr/bin/env python3
"""
Normalise the two-stage VisiumHD flexiplex output into a single barcode table.

The VisiumHD spot barcode is split across two positions in R1, so flexiplex is
run twice. The first pass matches the second half of the barcode and writes it
into the read ID (flexiplex's default -l true behaviour), the second pass
matches the UMI and the first half of the barcode and writes the
reads_barcodes.txt table. The Read column of that table therefore looks like

    <BC2>_#<original read ID>_+1of1

This script splits that apart and writes a table keyed on the original read ID,
with the two barcode halves joined as BC2_BC1 - the same composite barcode that
format_barcodes.py --type flexiplex_hd used to build, and the format expected by
convert_barcodes_spatial.py.
"""

import argparse
import re
import sys

HEADER = "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"

# flexiplex appends e.g. _+1of1, _-2of2 or _+1of1_C to the read ID
SUFFIX_RE = re.compile(r"_[+-]\d+of\d+(_C)?$")


def split_read_field(field):
    """Return (barcode_from_read_id, original_read_id) for a flexiplex Read field."""
    if "#" not in field:
        return None, field
    barcode_part, read_id = field.split("#", 1)
    # barcode_part is "<barcode>_<umi>"; the first pass carries no UMI
    barcode = barcode_part.split("_", 1)[0]
    read_id = SUFFIX_RE.sub("", read_id)
    return barcode, read_id


def main():
    parser = argparse.ArgumentParser(
        description="Normalise two-stage VisiumHD flexiplex output."
    )
    parser.add_argument("--input", required=True,
                        help="reads_barcodes.txt from the second flexiplex pass.")
    parser.add_argument("--output", required=True,
                        help="Normalised barcode table.")
    args = parser.parse_args()

    n_in = n_out = 0
    with open(args.input) as fh, open(args.output, "w") as out:
        out.write(HEADER + "\n")
        for line in fh:
            if line.startswith("Read\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            n_in += 1
            bc2, read_id = split_read_field(fields[0])
            if bc2 is None:
                print(f"WARNING: no barcode in read ID '{fields[0]}', skipping",
                      file=sys.stderr)
                continue
            out.write("\t".join([read_id, f"{bc2}_{fields[1]}",
                                 fields[2], fields[3], fields[4]]) + "\n")
            n_out += 1

    print(f"Normalised {n_out} of {n_in} barcode assignments to {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
