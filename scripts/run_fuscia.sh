import pandas as pd
import sys
import os

def run_fuscia(masterdata, map_qual, fuscia, cur_dir)
    for index, row in r.iterrows():
        #if map_qual != NULL
        gene1 = f'{row["chrom1"]}:{min(row["gene1"], row["base1"])}-{max(row["gene1"], row["base1"])}'
        gene2 = f'{row["chrom2"]}:{min(row["gene2"], row["base2"])}-{max(row["gene2"], row["base2"])}'
        os.system(f'python {fuscia}/discover_chimeric_transcripts.py {cur_dir}/cellranger_output/outs/*.bam {gene1} {gene2} {cur_dir}/fuscia_out {row["fusion genes"]} {map_qual}')

