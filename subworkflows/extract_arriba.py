#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:08:53 2023

@author: wu.s
"""

import pandas as pd

arriba = pd.read_csv("/Users/wu.s/Documents/Research/Arriba/fusions.tsv", sep='\t')
print(arriba)
arriba_c = arriba[['#gene1','gene2','read_identifiers']]


fusions = {}
for index, row in arriba_c.iterrows():
    fusion = f'{row["#gene1"]}--{row["gene2"]}'
    reads = row["read_identifiers"].split(',')
    for item in reads:
        if fusion not in fusions:
            fusions[fusion] = [item]
        else:
            fusions[fusion].append(item)

print(fusions)
df = pd.DataFrame([(k, v) for k, lst in fusions.items() for v in lst], columns=['Key', 'Value'])
df.to_csv('/Users/wu.s/Desktop/9cl/Arribafusions_ID.csv', index=False)