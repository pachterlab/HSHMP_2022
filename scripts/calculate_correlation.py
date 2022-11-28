#!/usr/bin/env python3

import sys

import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

barcodes_path = '/home/kristjan/kallisto_bf_analysis/benchmark/barcodes_consensus'
res_prefix = '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/results_sim_vs_'
res_id = [
    'kallisto',
    'kallisto_offlist_filtered_barcodes',
    'salmon_splici_cr_like',
    'salmon_splici_cr_like_em',
    'salmon_standard_cr_like',
    'salmon_standard_cr_like_em',
    'salmon_standard_sketch_cr_like',
    'salmon_standard_sketch_cr_like_em',
    'star'
]

def parse_results(path: str, barcodes: set, pad: bool) -> list[float]:

    # print(f'parsing {path}')
    rhos = []
    with open(path) as fh:

        # Ignore first line
        fh.readline()
        # Get number of genes in intersection
        l = int(fh.readline().split()[1])

        in_list = False

        for idx, line in enumerate(fh):

            mod = idx % 5
            line = line.rstrip('\n')
            # Filter on barcodes
            if mod == 0:
                if line not in barcodes:
                    in_list = False
                else:
                    in_list = True
            elif in_list:
                if mod == 1:
                    vec1 = [int(i) for i in line.split(',')]
                    if pad:
                        vec1 = [*vec1, *[0] * (l - len(vec1))]
                elif mod == 2:
                    vec2 = [int(i) for i in line.split(',')]
                    if pad:
                        vec2 = [*vec2, *[0] * (l - len(vec2))]
                elif mod == 3:
                    rho, _ = spearmanr(vec1, vec2)
                    rhos.append(rho)
                    vec1 = []
                    vec2 = []

    return rhos

if __name__ == '__main__':

    bc = set(b.rstrip('\n') for b in open(barcodes_path))

    pad = len(sys.argv) > 1 and sys.argv[1] == 'pad'
    data = {}
    for dat in res_id:
        data[dat] = parse_results(f'{res_prefix}{dat}.txt', bc, pad)

    df = pd.DataFrame(data=data)
    if pad:
        df.to_csv('spearman_coef_full.csv')
    else:
        df.to_csv('spearman_coef.csv')
    # df = pd.melt(df, value_vars=res_id)
    # sns.displot(df, x='value', hue='variable')
    # plt.show()

    # print(df)

