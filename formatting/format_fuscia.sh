import pandas as pd
import os
import sys

def add_fusion_name(file): 	
    r = pd.read_table(file)
    if r.empty == False:
        r['fusion'] =  os.path.basename(file)[:-33]
        return r

path = sys.argv[1]
cur_dir = sys.argv[2]
df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'chrom', 'start', 'end','fusion'])
for file in os.listdir(path):
    r = add_fusion_name(f'{path}{file}')
    df = pd.concat([df, r], axis = 0, ignore_index = True)
df['cell_barcode'] = df['cell_barcode'].str.replace('-1', "") 
df.to_csv(f'{cur_dir}/master_fuscia.csv', index=False)
