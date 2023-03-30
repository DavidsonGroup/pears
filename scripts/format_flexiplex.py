import pandas as pd
import os
import sys

def add_fusion_name(file):
    r = pd.read_table(file)
    if r.empty == False:
        df_temp = pd.DataFrame().assign(cell_barcode = r['Read'], molecular_barcode = r['BarcodeEditDist'])
        df_temp['fusion'] =  os.path.basename(file).split("_")[1]
        return df_temp

def format_flexiplex(cur_dir):
    df = pd.DataFrame(columns = ['cell_barcode', 'molecular_barcode', 'fusion'])

    for file in os.listdir(path):
        if os.path.basename(file)[0:8] == 'barcodes':
            #print (os.path.basename(file)[10:-19])
            r = add_fusion_name(f'{path}/{file}')
            df = pd.concat([df, r], axis = 0, ignore_index = True)
        
    return df

