#!/usr/bin/env python3
"""
Created on Mon Mar 20 20:11:22 2023

@author: wu.s
"""

#finding fuscia gene range
#input: gene_file(*/reference/gene/genes.gtf (refdata)), fusion_list (e.g. jjaffa output file)
#if region is bigger/smaller than breakpoint site = region
#output: {gene: [start, end]}
#issue: if name does not match (ENSG code might match) it is ignored.
import pandas as pd
import pybedtools
import os
from Bio.Seq import Seq
import sys

def get_gene_dict(shr_output):
    r = pd.read_csv(shr_output)
    f = r["fusion genes"].str.split('--', expand = True).to_numpy().flatten()

    gene_dict = {}

    for i in f:
        if i not in gene_dict:
            gene_dict[i] = []

    return gene_dict


def get_gene_range(gene_gtf, gene_dict):
    df = pd.read_csv(gene_gtf, delimiter="\t", usecols=[0, 2, 3, 4, 8] , skiprows=(5), header=None, names=['chrom', 'feature', 'start', 'end', 'attributes'])


    df['attributes'] = df['attributes'].str.rsplit('gene_name "').str.get(1)
    df['attributes'] = df['attributes'].str.rsplit('";').str.get(0)

    for index, row in df.iterrows():
        if row['feature'] == 'gene':
            if row['attributes'] in gene_dict and gene_dict[row['attributes']] == []:
                gene_dict[row['attributes']]= [row['chrom'], row['start'], row['end']]
            elif row['attributes'] in gene_dict and gene_dict[row['attributes']] != []:
                if row['start'] < gene_dict[row['attributes']][1]:
                    gene_dict[row['attributes']][1] = row['start']
                if row['end'] > gene_dict[row['attributes']][2]:
                    gene_dict[row['attributes']][1] = row['end']
    #gene_dict = dict( [(k,v) for (k,v) in gene_dict.items() if v] )

    return gene_dict

def gene_range(shr_output, gene_dict, up, down):

    r = pd.read_csv(shr_output)
    df = r[['fusion genes', 'chrom1', 'base1', 'strand1', 'chrom2', 'base2','strand2']].copy()
    df['fusion genes'] = df['fusion genes'].str.replace(':','--')
    df = df.assign(gene1=None, sequence1=None, gene2=None, sequence2=None, confidence=None)
    df = df[['fusion genes', 'chrom1', 'gene1', 'base1', 'sequence1', 'strand1', 'chrom2', 'gene2', 'base2', 'sequence2', 'strand2']]
    df = df[~df['fusion genes'].str.contains('MT-')]

    gene1 = []
    gene2 = []
    for index, row in df.iterrows():
        genes = row['fusion genes'].split("--")
        if gene_dict[genes[0]] != []:
            if row['strand1'] == '+' :
                gene1.append(gene_dict[genes[0]][1]) #start of gene 1
            else:
                gene1.append(gene_dict[genes[0]][2]) #end of gene 1
        else:
            if row['strand1'] == '+' :
                gene1.append(row['base1']-int(down))
            else:
                gene1.append(row['base1']+int(up))
        if gene_dict[genes[1]] != []:
            if row['strand1'] == '+' :
                gene2.append(gene_dict[genes[1]][2]) #end of gene 2
            else:
                gene2.append(gene_dict[genes[1]][1]) #start of gene 2
        else:
            if row['strand1'] == '+' :
                gene2.append(row['base2']+int(up))
            else:
                gene2.append(row['base2']-int(down))

    df['gene1'] = gene1
    df['gene2'] = gene2
    #df.set_index('fusion genes', inplace =True)

    return df

def get_sequence(gene, fasta):
    y = f'{gene[0]} {gene[1]} {gene[2]}'
    x = pybedtools.BedTool(y, from_string=True)
    a = x.sequence(fi=fasta)
    return open(a.seqfn).read().splitlines()[1]

def get_flexi_sequences(df, len_barcode, fasta):
    gene1_seq = []
    gene2_seq = []
    chromosomes = [f'chr{i}' for i in range(1, 24)] + ['chrX', 'chrY']
    for index, row in df.iterrows():
        if row['chrom1'] in chromosomes and row['chrom2'] in chromosomes:
            if row['strand1'] == '+':
                gene1seq = (get_sequence([row['chrom1'], row['base1'] - len_barcode, row['base1']], fasta))
            else:
                gene1seq = (get_sequence([row['chrom1'], row['base1']-1, row['base1'] + len_barcode - 1], fasta))
                gene1seq = Seq(gene1seq).reverse_complement()
            if row['strand2'] == '+':
                gene2seq = (get_sequence([row['chrom2'], row['base2']-1, row['base2'] + len_barcode - 1], fasta))
            else:
                gene2seq = (get_sequence([row['chrom2'], row['base2'] - len_barcode, row['base2']], fasta))
                gene2seq = Seq(gene2seq).reverse_complement()
        else:
            gene1seq = None
            gene2seq = None
        gene1_seq.append(str(gene1seq).upper())
        gene2_seq.append(str(gene2seq).upper())


    df['sequence1'] = gene1_seq
    df['sequence2'] = gene2_seq

    return df

def main():
    shr_output = sys.argv[1]
    gene = sys.argv[2]
    fasta = sys.argv[3]
    flexi_searchlen = sys.argv[4]
    out_dir = sys.argv[5]
    up = sys.argv[6]
    down = sys.argv[7]

    gene_dict = get_gene_dict(shr_output)
    add_gene_range = get_gene_range(gene, gene_dict)
    #with open(f'{out_dir}/generange.txt', 'w') as f:
        #for key in add_gene_range.keys():
            #f.write("%s, %s\n"%(key, add_gene_range[key]))
    fuscia_gene_range = gene_range(shr_output, add_gene_range, up, down)
    flexi_gene_range = get_flexi_sequences(fuscia_gene_range, int(flexi_searchlen), fasta)
    flexi_gene_range = flexi_gene_range.rename(columns={"fusion genes": "fusion_genes"})
    flexi_gene_range = flexi_gene_range.drop_duplicates() #remove duplicates in the fusion list
    flexi_gene_range.to_csv(f'{out_dir}/fusion_targets.csv', index = False)


if __name__ == "__main__":
    main()
