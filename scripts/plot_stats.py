import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

usage = '<stats> <epo stats> <out>'
if len(sys.argv) != 4:
    print(usage)
    sys.exit()


stats = pd.read_csv(sys.argv[1], sep='\t')
epo_stats = pd.read_csv(sys.argv[2], sep='\t')
#print(stats.head())
#print(epo_stats.head())

bases=['A','T','G','C']
total_bases = sum(np.sum(stats[base]) for base in stats)
total_bases_epo = sum(np.sum(epo_stats[base]) for base in stats)
print(f'total bases: {total_bases}')
print(f'total bases epo: {total_bases_epo}')

#sum bases for each chr
stats['sum'] = stats[bases].sum(axis=1)
epo_stats['sum'] = epo_stats[bases].sum(axis=1)
#print(stats.head())
#print(epo_stats.head())

stats.set_index('chr', inplace=True)
epo_stats.set_index('chr', inplace=True)
combined_stats = pd.concat([stats['sum'], epo_stats['sum']], axis=1, keys=['stats', 'epo_stats'])

fig, ax = plt.subplots(figsize=(12, 6))
bar_width = 0.3
bar_positions = np.arange(len(combined_stats.index))

# Bar for 'stats'
ax.bar(bar_positions - bar_width/2, combined_stats['stats'], bar_width, label='T2T Primate Alignment', color='#377eb8')

# Bar for 'epo_stats'
ax.bar(bar_positions + bar_width/2, combined_stats['epo_stats'], bar_width, label='EPO Alignment', color='#ff7f00')

# Set labels and title
ax.set_ylabel("Ancestral Bases Annotated")
ax.set_xlabel("Chromosome")
ax.set_xticks(bar_positions)
ax.set_xticklabels(combined_stats.index)
from matplotlib.ticker import FuncFormatter

def millions_formatter(x, pos):
    return f'{int(x/1e6)}M'


ax.yaxis.set_major_formatter(FuncFormatter(millions_formatter))

ax.legend()

plt.savefig(sys.argv[3], bbox_inches='tight', dpi=300)
