from Bio import SeqIO
import sys
import numpy as np
import matplotlib.pyplot as plt

usage = '<fasta 1> <fasta 2> <out>'
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

seq_1 = next(SeqIO.parse(sys.argv[1], 'fasta'))
seq_2 = next(SeqIO.parse(sys.argv[2], 'fasta'))


eqs = []
starts = []
#5kb windows for now, no step size
window_size = 5000
for i in range(0, len(seq_1.seq), window_size):

    window_seq_1 = seq_1.seq[i:i+window_size]
    window_seq_2 = seq_2.seq[i:i+window_size]
    eq=0
    tot=0
    for base1, base2 in zip(window_seq_1, window_seq_2):
        if base1 == '.' or base2 == '.' or base1 == '-' or base2 == '-' or base1 == 'N' or base2 == 'N':
            continue
    
        tot+=1
        if base2.upper() == base1.upper():
            eq += 1
    
    if tot - eq == 0:
        continue
    #probably just fine to look at mismatches alone
    eqs.append(tot - eq)
    #1 based start
    starts.append(i+1)

np.save(sys.argv[3]+'_mismatch', eqs)
np.save(sys.argv[3]+'_pos', starts)

idx = np.argmax(eqs)
#print(np.amax(eqs))
print(f'most mismatches in 5kb window {eqs[idx]}')
print(f'position of window {starts[idx]}')

#plot distribution
plt.hist(eqs, density=True, color='blue', bins=200, alpha=.7)
plt.xlim(left=0)
plt.xlabel('Mismatches in 5kb windows')
plt.ylabel('Density')
plt.savefig(sys.argv[3]+'_windows.pdf', bbox_inches='tight')

        


