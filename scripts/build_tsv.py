import os
import sys
import re

usage = '<dir> <tsv>'
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

# Specify the directory containing the _stats.txt files
directory = sys.argv[1]

# Initialize an empty list to store the data
data = []

# Iterate through files in the directory
for filename in os.listdir(directory):
    if filename.endswith("_stats.txt") or filename.endswith(".fa.stats"):
        filepath = os.path.join(directory, filename)
        
        pattern = re.compile(r'(\d+)')
        if filename.endswith("_stats.txt"):
            pattern = re.compile(r'_chr(\d+)')

        match = pattern.search(filename)

        # Extract the matched number
        if match:
            chr_number = match.group(1)
            print(chr_number)
        else:
            continue 

        # Read the content of each file
        with open(filepath, 'r') as file:
            content = file.read()
            
            # Extract values from the content
            #values = [filename[:-10]]  # Extracting the prefix (assuming filenames are like 'something_stats.txt')
            
            pattern = re.compile(r'\b(\d+)\b')
            numbers = pattern.findall(content)
            data.append([chr_number] + [num for num in numbers])

# Write the data to a TSV file
tsv_filename = sys.argv[2]
with open(tsv_filename, 'w') as tsv_file:
    # Write header
    tsv_file.write('\t'.join(['chr', 'A', 'C', 'G', 'T', 'Periods']) + '\n')
    
    # Write data
    for values in sorted(data, key = lambda x: int(x[0])):
        tsv_file.write('\t'.join(values) + '\n')

print(f'TSV file "{tsv_filename}" has been created.')

