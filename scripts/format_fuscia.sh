import pandas as pd
import os
import sys

def add_fusion_name(file): 	
    r = pd.read_table(file)
    if r.empty == False:
        r['fusion'] =  os.path.basename(file).split(".")[0]
        return r

def format_fuscia(fuscia_dir):
    df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'chrom', 'start', 'end','fusion'])
    for file in os.listdir(fuscia_dir):
        r = add_fusion_name(f'{path}{file}')
        df = pd.concat([df, r], axis = 0, ignore_index = True)
    df['cell_barcode'] = df['cell_barcode'].str.replace('-1', "") 
    return df
