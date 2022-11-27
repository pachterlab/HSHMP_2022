#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

spearman_path = '../spearman_coef.csv'


df = pd.read_csv(spearman_path)
df = pd.melt(df, value_vars=df.columns[1:])

sns.set_theme()
#  sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})
sns.set_style("whitegrid")
ax = sns.boxplot(df, x='value', y='variable', notch=True, showcaps=False, flierprops={"marker": "x"})#, legend=False)
#  ax.set_yticklabels(["kallisto", "kallisto (with off-list)", "Alevin-fry cr-like", "Alevin-fry cr-like-em", "STARsolo"], rotation=45)
ax.set(ylabel=None)
#  plt.xlim(0.65, 1)
#  plt.legend(title='Method', loc='upper left', labels=['STARsolo',
#                                                       'Alevin-fry splici cr-like-em',
#                                                       'Alevin-fry splici cr-like',
#                                                       'kallisto (with off-list)',
#                                                       'kallisto'])
plt.xlabel('Per-cell Spearman coefficient')
plt.show()
