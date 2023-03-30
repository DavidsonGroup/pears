import pandas as pd
import sys
import os

def run_flexiplex(masterdata, reads, flexiplex, cur_dir):
    os.makedirs(f'{cur_dir}/flexiplex_output', exist_ok = True)
    r = pd.read_csv(masterdata)
    #os.system('pwd')
    #os.system(f'paste -d "&" ${reads}/*.fastq > jreads.fastq')
    path = f'{cur_dir}/flexiplex_output'
    for index, row in r.iterrows():
        fusion = row['fusion genes']
        if os.path.isfile(f'{path}/{fusion}_reads.fastq') == True:
            i = 1
            fusion_name = f'{fusion}_{i}'
            while os.path.isfile(f'{path}/{fusion}_{i}_reads.fastq') == True:
                i += 1
                fusion_name = f'{fusion}_{i}'
        else:
            fusion_name = fusion
        left = row['sequence1']
        right = row['sequence2']
        os.system(f'paste -d "&" {reads}/*.fastq | {flexiplex} -n {fusion_name} -l {left} -k {right} -r "" -f 1 -e 1 -u 0 -i false -s false > {path}/{fusion_name}_reads.fastq')
        os.system(f'{flexiplex} -l "" -r "&" -e 0 -f 0 -n {fusion_name} {path}/{fusion_name}_reads.fastq')
        os.system(f'{flexiplex} -l "" -k {fusion_name}_barcodes_counts.txt -b 16 -r "&" -e 0 -f 0 -n barcodes_{fusion_name} {path}/{fusion_name}_reads.fastq')
       

