from Bio import SeqIO
import sys
import numpy as np
import matplotlib.pyplot as plt

usage = '<fasta 1> <fasta 2>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

seq_1 = next(SeqIO.parse(sys.argv[1], 'fasta'))
seq_2 = next(SeqIO.parse(sys.argv[2], 'fasta'))

tot = 0 #total number of bases that have coverage
eq = 0 #bases that are equal
diffs = []
diff = 0
for base1, base2 in zip(seq_1.seq, seq_2.seq):
    if base1 == '.' or base2 == '.':
        continue
    
    #if tot > 0 and tot % 1000000 == 0:
        
        #print(eq/tot, flush=True)
    tot += 1
    if base2.upper() == base1.upper():
        eq += 1
        if diff > 0:
            diffs.append(diff)
        diff = 0
    else:
        diff += 1

print(np.mean(diffs))
print(np.std(diffs))
print(np.amax(diffs))
#np.save('diffs', diffs)
print(eq/tot)

#plt.hist(diffs, bins='rice', color='red')
#plt.xlabel('Length of differences')
#plt.ylabel('Count')
#plt.savefig('diffs.pdf')

