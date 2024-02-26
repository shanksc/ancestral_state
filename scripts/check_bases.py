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
positions = []
diff = 0
diff_str1 = ''
diff_str2 = ''
pos = 0 #1 based 
fst_diff = 0
#count number of bases that are uppercase where there is disagreement 
diff_upper_1 = 0
diff_upper_2 = 0
for base1, base2 in zip(seq_1.seq, seq_2.seq):
    pos += 1
    if base1 == '.' or base2 == '.' or base1 =='-' or base2 =='-' or base1 == 'N' or base2 == 'N':
        continue
    
    #if tot > 0 and tot % 1000000 == 0:
        
    tot += 1 #print(eq/tot, flush=True)
    if base2.upper() == base1.upper():
        eq += 1
        if diff > 0:
            #print(pos)
            #print(base2.upper())
            #print(base1.upper())
            #sys.exit()
            '''
            if diff >= 10:
                print(fst_diff)
                print(diff_str1)
                print(diff_str2)
                sys.exit()
            '''
            diffs.append(diff)
            positions.append(fst_diff)
        diff = 0
        '''
        diff_str1 = ''
        diff_str2 = ''
        '''
    else:
        #print(base1, base2)
        #print(pos)
            #sys.exit()
        #record start of new diff
        if diff == 0:
            fst_diff = pos
        diff += 1
        '''
        diff_str1 += base1
        diff_str2 += base2
        '''

tot_diffs = np.sum(diffs)
print(f'total diffs: {tot_diffs}')
print(f'u: {np.mean(diffs)}')
print(f'sd: {np.std(diffs)}')
print(f'max diff length: {np.amax(diffs)}')
print(f'position of max diff: {positions[np.argmax(diffs)]}')
#np.save('diffs', diffs)
print(eq/tot)
#print(f'base 1 uppercase {diff_upper_1/tot_diffs}')
#print(f'base 2 uppercase {diff_upper_2/tot_diffs}')

#plt.hist(diffs, bins='rice', color='red')
#plt.xlabel('Length of differences')
#plt.ylabel('Count')
#plt.savefig('diffs.pdf')

