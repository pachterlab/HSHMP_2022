#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

barcodes_path = '/home/kristjan/kallisto_bf_analysis/benchmark/barcodes_consensus'
res_prefix = '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/results_sim_vs_'
res_id = [
    # 'kallisto',
    'kallisto_offlist',
    'salmon_splici_cr_like',
    'salmon_splici_cr_like_em',
    'star'
]

def parse_results(path: str, barcodes: set) -> list[float]:

    print(f'parsing {path}')
    rhos = []
    with open(path) as fh:

        # Ignore first two lines
        fh.readline()
        fh.readline()
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
                elif mod == 2:
                    vec2 = [int(i) for i in line.split(',')]
                elif mod == 3:
                    rho, _ = spearmanr(vec1, vec2)
                    rhos.append(rho)
                    vec1 = []
                    vec2 = []

    return rhos

if __name__ == '__main__':

    bc = set(b.rstrip('\n') for b in open(barcodes_path))

    data = {}
    for dat in res_id:
        data[dat] = parse_results(f'{res_prefix}{dat}.txt', bc)

    for idx, l in data.items():
        print(idx, len(l))

    df = pd.DataFrame(data=data)
    df.to_csv('spearman_coef.csv')
    # df = pd.melt(df, value_vars=res_id)
    # sns.displot(df, x='value', hue='variable')
    # plt.show()

    # print(df)

