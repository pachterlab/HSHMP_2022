#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

spearman_path = '../spearman_coef.csv'


df = pd.read_csv(spearman_path)
df = pd.melt(df, value_vars=df.columns[1:])

sns.set_theme()
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})
sns.set_style("whitegrid")
sns.displot(df[df['variable'].isin(['kallisto', 'kallisto_offlist_filtered_barcodes', 'salmon_splici_cr_like', 'salmon_splici_cr_like_em', 'star'])], x='value', hue='variable', kde=True)#, legend=False)
#  sns.displot(df, x='value', hue='variable', kde=True)#, legend=False)

plt.xlim(0.98, 1)
#  plt.legend(title='Method', loc='upper left', labels=['STARsolo',
#                                                       'Alevin-fry splici cr-like-em',
#                                                       'Alevin-fry splici cr-like',
#                                                       'kallisto (with off-list)',
#                                                       'kallisto'])
plt.xlabel('Per-cell Spearman coefficient')
plt.ylabel('Number of cells')
plt.show()
