import pandas as pd
import sys
import os
import threading

def run_flexiplex(masterdata, reads, flexiplex, cur_dir, start, end):
    r = pd.read_csv(masterdata)
    #os.system('pwd')
    #os.system(f'paste -d "&" ${reads}/*.fastq > jreads.fastq')
    path = f'{cur_dir}/flexiplex_output'
    for index, row in r.iloc[start:end].iterrows():
        fusion_name = row['fusion genes']
        i = 1
        while os.path.isfile(f'{path}/{fusion_name}_reads_barcodes.txt') == True:
            fusion_name = f'{fusion_name}_{i}'
            i += 1
        left = row['sequence1']
        right = row['sequence2']
        os.system(f'paste -d "&" {reads}/*.fastq | {flexiplex} -n {fusion_name} -l {left} -k {right} -r "" -f 1 -e 1 -u 0 -s false > {path}/{fusion_name}_reads.fastq ')
        os.system(f'{flexiplex} -l "" -r "" -e 0 -f 0 -n {fusion_name} {path}/{fusion_name}_reads.fastq')
        os.system(f'{flexiplex} -l "" -k {fusion_name}_barcodes_counts.txt -b 16 -r "" -e 0 -f 0 -n barcodes_{fusion_name} {path}/{fusion_name}_reads.fastq')


def split_processing(items, num_splits):                                      
    split_size = len(items) // num_splits                                       
    threads = []                                                                
    for i in range(num_splits):                                                 
        # determine the indices of the list this thread will handle             
        start = i * split_size                                                  
        # special case on the last chunk to account for uneven splits           
        end = None if i+1 == num_splits else (i+1) * split_size                 
        # create the thread                                                     
        threads.append(                                                         
            threading.Thread(target=run_flexiplex, args=(masterdata, reads, flexiplex, cur_dir, start, end)))         
        threads[-1].start() # start the thread we just created                  

    # wait for all threads to finish                                            
    for t in threads:                                                           
        t.join()                                                                



reads = sys.argv[1]
flexiplex = sys.argv[2]
masterdata = sys.argv[3]
cur_dir = sys.argv[4]
os.makedirs(f'{cur_dir}/flexiplex_output', exist_ok = True)
split_processing(masterdata, 8)
