"""
Created on Mon Mar 20 20:11:22 2023

@author: wu.s
"""

#finding fuscia gene range
#input: gene_file(*/reference/gene/genes.gtf (refdata)), gene_list (jaffa output file)
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
    #r = r[r.classification == 'HighConfidence']
    f = r["fusion genes"].str.split(':', expand = True).to_numpy().flatten()

    gene_dict = {}

    for i in f:
        if i not in gene_dict:
            gene_dict[i] = []

    return gene_dict


def get_gene_range(gene_gtf, gene_dict):
    r = pd.read_csv(gene_gtf, delimiter="\t", skiprows=(5), header=None)
    df = pd.DataFrame(r, columns = [0, 3, 4, 8])

    df[8] = df[8].str.rsplit('gene_name "').str.get(1)
    df[8] = df[8].str.rsplit('";').str.get(0)

    for index, row in df.iterrows():
        if row[8] in gene_dict and gene_dict[row[8]] == []:
            gene_dict[row[8]]= [row[0], row[3], row[4]]
        elif row[8] in gene_dict and gene_dict[row[8]] != []:
            if row[3] < gene_dict[row[8]][1]:
                gene_dict[row[8]][0] = row[3]
            if row[4] > gene_dict[row[8]][2]:
                gene_dict[row[8]][1] = row[4]

    gene_dict = dict( [(k,v) for (k,v) in gene_dict.items() if v] )

    return gene_dict

def gene_range(shr_output, gene_dict):

    r = pd.read_csv(shr_output)
    df = r[['fusion genes', 'chrom1', 'base1', 'strand1', 'chrom2', 'base2','strand2', 'classification']].copy()
    df['fusion genes'] = df['fusion genes'].str.replace(':','--')
    df = df.assign(gene1=None, sequence1=None, gene2=None, sequence2=None, confidence=None)
    df = df[['fusion genes', 'chrom1', 'gene1', 'base1', 'sequence1', 'strand1', 'chrom2', 'gene2', 'base2', 'sequence2', 'strand2', 'confidence']]

    gene1 = []
    gene2 = []
    for index, row in df.iterrows():
        genes = row['fusion genes'].split("--")
        if genes[0] in gene_dict and genes[1] in gene_dict:
            gene1.append(gene_dict[genes[0]][1])
            gene2.append(gene_dict[genes[1]][2])
        else:
            df = df.drop(index)
    
    df['gene1'] = gene1
    df['gene2'] = gene2

    return df

def get_sequence(gene, fasta):
    y = f'{gene[0]} {gene[1]} {gene[2]}'
    x = pybedtools.BedTool(y, from_string=True)
    a = x.sequence(fi=fasta)
    return open(a.seqfn).read().splitlines()[1]

def get_flexi_sequences(df, len_barcode, fasta):

    gene1_seq = []
    gene2_seq = []
    for index, row in df.iterrows():
        gene1seq = (get_sequence([row['chrom1'], row['base1'] - len_barcode, row['base1']], fasta))
        gene2seq = (get_sequence([row['chrom2'], row['base2'], row['base2'] + len_barcode], fasta))
        if row['strand1'] == '-':
            gene1seq = Seq(gene1seq).reverse_complement()
        elif row['strand2'] == '-':
            gene2seq = Seq(gene2seq).reverse_complement()
        gene1_seq.append(gene1seq)
        gene2_seq.append(gene2seq)

    df['sequence1'] = gene1_seq
    df['sequence2'] = gene2_seq

    return df

def gen_masterdata(shr_output, reference, out_dir, flexi_seqlen):
    gene = reference + '/genes/genes.gtf'
    fasta = reference + '/fasta/genome.fa'
    gene_dict = get_gene_dict(shr_output)
    add_gene_range = get_gene_range(gene, gene_dict)
    fuscia_gene_range = gene_range(shr_output, add_gene_range)
    #20 is default, it can also be altered
    flexi_gene_range = get_flexi_sequences(fuscia_gene_range, flexi_seqlen, fasta)
    return flexi_gene_range
