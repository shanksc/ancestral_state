import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.stats import ttest_rel
import numpy as np

usage = '<tsv> <out>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

df = pd.read_csv(sys.argv[1], delimiter='\t')


difference = df['conf_1'] - df['conf_2']

difference *= 100
print(len(difference))

#paired t test, conf_1 > conf_2 
print(ttest_rel(df['conf_1'], df['conf_2'], alternative='greater'))


#histogram
plt.figure(figsize = (12,6))

#sns.kdeplot(difference, bw = 0.5 , fill = True)
plt.hist(difference, bins=100, color='blue', alpha=0.7)
plt.xlabel('Difference in % high confidence bases (T2T - Ensembl)')
plt.ylabel('Frequency')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig(sys.argv[2]+'_hist_full.pdf', bbox_inches='tight')
plt.xlim(-.5,1)
plt.savefig(sys.argv[2]+'_hist.pdf',bbox_inches='tight')
plt.clf()

bins = np.linspace(0, 1, 100)
plt.hist(df['conf_1'], density=True, bins=bins, alpha=0.7, label='T2T')
plt.hist(df['conf_2'], density=True, bins=bins, alpha=0.7, label='EPO')
plt.legend(loc='upper left')
plt.xlabel('% High confidence bases')
plt.savefig(sys.argv[2]+'_hists.pdf', bbox_inches='tight')
