#!/usr/bin/env python3
"""
Consolidate barcode formatting for Fuscia, Flexiplex, and Arriba outputs.
"""

import argparse
import os
import pandas as pd


def format_fuscia(input_files):
    """Format Fuscia output files."""
    df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])

    for file in input_files:
        r = pd.read_table(file)
        if not r.empty:
            r['fusion'] = os.path.basename(file).split("_")[0]
            r = r[r['cell_barcode'] != '-']
            df = pd.concat([df, r], axis=0, ignore_index=True)

    df['cell_barcode'] = df['cell_barcode'].str.replace('-1', "")

    df = df.iloc[:, :-3]
    df = df.drop_duplicates()
    return df


def format_flexiplex_arriba(input_files):
    """Format Flexiplex or Arriba output files."""
    df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'fusion'])

    for file in input_files:
        if os.path.basename(file).startswith('barcodes'):
            r = pd.read_table(file)
            if not r.empty:
                df_temp = pd.DataFrame().assign(
                    cell_barcode=r['CellBarcode'],
                    molecular_barcode=r['UMI']
                )
                df_temp['fusion'] = os.path.basename(file).split("_")[1]
                df = pd.concat([df, df_temp], axis=0, ignore_index=True)

    df = df.drop_duplicates()
    return df


def main():
    parser = argparse.ArgumentParser(description='Format barcode files from fusion detection tools')
    parser.add_argument('--type', '-t', required=True,
                        choices=['fuscia', 'flexiplex', 'arriba'],
                        help='Type of input files to format')
    parser.add_argument('--output', '-o', required=True,
                        help='Output CSV file path')
    parser.add_argument('input_files', nargs='+',
                        help='Input files to process')

    args = parser.parse_args()

    if args.type == 'fuscia':
        df = format_fuscia(args.input_files)
    elif args.type in ['flexiplex', 'arriba']:
        df = format_flexiplex_arriba(args.input_files)
    else:
        raise ValueError("Unsupported type. Use 'fuscia', 'flexiplex', or 'arriba'.")

    df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
