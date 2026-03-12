#!/usr/bin/env python3
## Written by ChatGPT March 2026

import argparse
import csv
from collections import defaultdict


def read_method_file(path, method_name, per_pair):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)

        required = {"cell_barcode", "molecular_barcode", "fusion"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            cell = row["cell_barcode"].strip()
            umi = row["molecular_barcode"].strip()
            fusion = row["fusion"].strip()

            if not cell or not umi or not fusion:
                continue

            key = (fusion, cell)
            per_pair[key][method_name].add(umi)
            per_pair[key]["total"].add(umi)


def main():
    parser = argparse.ArgumentParser(
        description="Combine Arriba, Flexiplex, and Fuscia fusion calls per fusion/cell pair."
    )
    parser.add_argument("--arriba", required=True, help="Arriba CSV")
    parser.add_argument("--flexiplex", required=True, help="Flexiplex CSV")
    parser.add_argument("--fuscia", required=True, help="Fuscia CSV")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    args = parser.parse_args()

    per_pair = defaultdict(
        lambda: {
            "total": set(),
            "arriba": set(),
            "flexiplex": set(),
            "fuscia": set(),
        }
    )

    read_method_file(args.arriba, "arriba", per_pair)
    read_method_file(args.flexiplex, "flexiplex", per_pair)
    read_method_file(args.fuscia, "fuscia", per_pair)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "fusion",
                "cell_barcode",
                "total_umis",
                "arriba_umis",
                "flexiplex_umis",
                "fuscia_umis",
            ]
        )

        for (fusion, cell) in sorted(per_pair):
            writer.writerow(
                [
                    fusion,
                    cell,
                    len(per_pair[(fusion, cell)]["total"]),
                    len(per_pair[(fusion, cell)]["arriba"]),
                    len(per_pair[(fusion, cell)]["flexiplex"]),
                    len(per_pair[(fusion, cell)]["fuscia"]),
                ]
            )


if __name__ == "__main__":
    main()
