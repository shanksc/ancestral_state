import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

usage = '<tsv> <out>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

#99.9 percentile

df = pd.read_csv(sys.argv[1], sep='\t', index_col=False)
print(df.head())

percentile = df['mismatches'].quantile(.999)
print(percentile)
#should already be sorted
#df = df[df['pbs'] > CUTOFF]
df.reset_index(inplace=True, drop=True)
df['i'] = df.index

plot = sns.relplot(data=df, x='i', y='mismatches', aspect=2.5, hue='chr', palette=sns.color_palette('colorblind'), legend=None, linewidth=0, s=11)
plt.margins(x=0.01, y=0.01)
chrom_df=df.groupby('chr')['i'].median()
plot.ax.set_xlabel('Chromsome')
plot.ax.set_xticks(chrom_df, labels=[x for x in range(1, len(chrom_df.index)+1)])
plot.ax.tick_params(axis='x', labelsize=9)
plt.ylabel('Mismatches in 5kb window')
plt.axhline(y=percentile, linestyle='--', color='black', linewidth=1.5, label='99.9th percentile')
plt.legend(loc='upper left')
plt.savefig(sys.argv[2])

