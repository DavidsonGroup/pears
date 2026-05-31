#!/usr/bin/env python3
"""
Convert VisiumHD two-part barcodes (BC1_BC2) in combined_fusions.csv to
standard SpaceRanger spatial barcode format: s_008um_XXXXX_YYYYY-1.

The VisiumHD whitelist maps each barcode sequence (by line position, 1-indexed)
to a coordinate on the capture array (1-3350 on each axis). The spatial bin
barcode is derived by dividing each axis position by the bin divisor and
zero-padding to 5 digits.

Bin size -> divisor:
  2um  -> 1   (s_002um_...)
  8um  -> 4   (s_008um_...)
  16um -> 8   (s_016um_...)
"""

import argparse
import sys
import pandas as pd


BIN_CONFIG = {
    2:  (1,  "002um"),
    8:  (4,  "008um"),
    16: (8,  "016um"),
}


def load_whitelist(path):
    """Return dict {barcode_seq: 1-based_position} from whitelist (one seq per line)."""
    mapping = {}
    with open(path) as fh:
        for i, line in enumerate(fh, start=1):
            seq = line.strip()
            if seq:
                mapping[seq] = i
    return mapping


def convert_barcode(bc, whitelist, divisor, prefix):
    """Convert a BC1_BC2 string to s_{prefix}_XXXXX_YYYYY-1, or None if unmappable."""
    parts = bc.split("_", 1)
    if len(parts) != 2:
        return None
    bc1, bc2 = parts
    pos1 = whitelist.get(bc1)
    pos2 = whitelist.get(bc2)
    if pos1 is None or pos2 is None:
        return None
    x = round(pos1 / divisor)
    y = round(pos2 / divisor)
    return f"s_{prefix}_{x:05d}_{y:05d}-1"


def main():
    parser = argparse.ArgumentParser(
        description="Convert VisiumHD BC1_BC2 barcodes to spatial bin format."
    )
    parser.add_argument("--input",    required=True, help="combined_fusions.csv")
    parser.add_argument("--whitelist", required=True, help="Uncompressed barcode whitelist")
    parser.add_argument("--bin-size", required=True, type=int, choices=[2, 8, 16],
                        help="VisiumHD bin size in microns (2, 8, or 16)")
    parser.add_argument("--output",   required=True, help="Output CSV path")
    args = parser.parse_args()

    divisor, prefix = BIN_CONFIG[args.bin_size]

    print(f"Loading whitelist: {args.whitelist}", file=sys.stderr)
    whitelist = load_whitelist(args.whitelist)
    print(f"  {len(whitelist)} barcodes loaded", file=sys.stderr)

    df = pd.read_csv(args.input)
    if "cell_barcode" not in df.columns:
        sys.exit(f"ERROR: cell_barcode column not found in {args.input}")

    converted = df["cell_barcode"].apply(
        lambda bc: convert_barcode(bc, whitelist, divisor, prefix)
    )
    n_ok   = converted.notna().sum()
    n_fail = converted.isna().sum()
    print(f"Converted {n_ok} barcodes; {n_fail} could not be mapped", file=sys.stderr)

    df["cell_barcode"] = converted
    df = df.dropna(subset=["cell_barcode"])
    df.to_csv(args.output, index=False)
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
