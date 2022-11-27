#!/usr/bin/env python3

import pandas as pd
from scipy.stats.mstats import mquantiles

spearman_path = '../spearman_coef.csv'

df = pd.read_csv(spearman_path)
print(mquantiles(df[df.columns[1:]], axis=0))

