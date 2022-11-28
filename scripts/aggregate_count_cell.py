#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

barcodes_path = '/home/kristjan/kallisto_bf_analysis/benchmark/barcodes_consensus'
res_prefix = '/home/dsullivan/benchmarking/starsolo/STARsoloManuscript/results_sim_vs_'
res_id = [
    'kallisto',
    'kallisto_offlist_filtered_barcodes',
    # 'salmon_splici_cr_like',
    # 'salmon_splici_cr_like_em',
    # 'salmon_standard_cr_like',
    # 'salmon_standard_cr_like_em',
    # 'salmon_standard_sketch_cr_like',
    # 'salmon_standard_sketch_cr_like_em',
    # 'star'
]

def parse_results(path: str, col: str, bc: set, df: pd.DataFrame) -> pd.DataFrame:

    # print(f'parsing {path}')
    expr = []
    with open(path) as fh:

        # Ignore first line
        fh.readline()
        # Get number of genes in intersection
        l = int(fh.readline().split()[1])

        in_list = False
        barcode = ''

        for idx, line in enumerate(fh):

            mod = idx % 5
            line = line.rstrip('\n')
            # Filter on barcodes
            if mod == 0:
                if line not in bc:
                    in_list = False
                    barcode = ''
                else:
                    in_list = True
                    barcode = line
            elif in_list and mod == 2:
                df.loc[df['barcode'] == barcode, [col]] = sum(int(i) for i in line.split(','))

    return df

if __name__ == '__main__':

    bc = set(b.rstrip('\n') for b in open(barcodes_path))

    df = pd.DataFrame()
    df['barcode'] = list(bc)
    for dat in res_id:
        df[dat] = 0
        df = parse_results(f'{res_prefix}{dat}.txt', dat, bc, df)

    df.to_csv('aggregate_count_cell.csv')
