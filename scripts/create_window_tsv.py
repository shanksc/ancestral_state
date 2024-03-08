import numpy as np
import pandas as pd
import sys

usage = '<dir> <out>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

arrays = []
path = sys.argv[1]

#write out tsv containing positions
with open(sys.argv[2], 'w') as f:
    f.write('chr\tstart\tmismatches\tconf_1\tconf_2\n')
    for n in range(1, 23):
        #10_window_mismatch.npy 10_window_pos.npy
        mismatches = np.load(f'{path}/{n}_window_mismatch.npy')
        starts = np.load(f'{path}/{n}_window_pos.npy')
        conf_1 = np.load(f'{path}/{n}_window_conf_1.npy')
        conf_2 = np.load(f'{path}/{n}_window_conf_2.npy')

        #write out mismatches and starts
        for i in range(starts.shape[0]):
            f.write(f'{n}\t{starts[i]}\t{mismatches[i]}\t{conf_1[i]}\t{conf_2[i]}\n')
