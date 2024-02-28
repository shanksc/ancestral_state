import numpy as np
import matplotlib.pyplot as plt
import sys

usage = '<dir> <out>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

arrays = []
path = sys.argv[1]

for file_name in [f'{path}/{n}.npy' for n in range(1,23)]:
    array = np.load(file_name)
    arrays.append(array)

result_array = np.concatenate(arrays, axis=0)

np.save(sys.argv[2]+'.npy', result_array)

plt.hist(result_array, density=True, bins=np.amax(result_array), color='red', alpha=0.7)
plt.xlabel('Length of continuous mismatches (bp)')
plt.ylabel('Density')
plt.xlim(1,7)
plt.savefig(sys.argv[2]+'.pdf', bbox_inches='tight')

