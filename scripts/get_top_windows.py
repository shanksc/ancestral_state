import pandas as pd
import sys

usage = '<tsv> <out>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

df = pd.read_csv(sys.argv[1], sep='\t')

cutoff = df['mismatches'].quantile(.99)
print(f'99th percentile: {cutoff}')

df = df[df['mismatches'] > cutoff]
df.to_csv(sys.argv[2], sep='\t', index=False)

