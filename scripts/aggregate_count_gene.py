#!/usr/bin/env python3

import os

from scipy.io import mmread
import pandas as pd

prefix = '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/count/kallisto_0.49.0/human_CR_3.0.0/'
index = ['standard_1', 'standard_offlist_1']
postfix = 'default/10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo/20/run1/counts_unfiltered/'
cxg = 'cells_x_genes.mtx'

df = pd.DataFrame()

gen = [l.rstrip() for l in open(os.path.join(prefix, index[0], postfix, 'cells_x_genes.genes.txt'))]
bcs_consensus = set(l.rstrip() for l in open('/home/kristjan/kallisto_bf_analysis/benchmark/barcodes_consensus'))

df['gene'] = gen

for idx in index:
    bcs = (l.rstrip() for l in open(os.path.join(prefix, idx, postfix, 'cells_x_genes.barcodes.txt')))

    bcs = [i for i, b in enumerate(bcs) if b in bcs_consensus]

    a = mmread(os.path.join(prefix, idx, postfix, cxg)).todense()

    for i, gene in enumerate(gen):
        df.loc[df['gene'] == gene, [idx]] = a[bcs, i].sum()

df.to_csv('aggregate_count_gene.csv')
