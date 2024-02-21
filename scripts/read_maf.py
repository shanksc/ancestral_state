from taffy.lib import AlignmentReader, TafIndex
import pathlib
import bisect
import sys
import re

usage = '<taf>'

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

def all_ancs_in_block(ancs, seqs):
    #max_len = len('Anc0.Anc0refChr22') 18
    valid_seqs = set([])
    count = 0
    for anc in ancs:
        for seq in seqs:
            # check that anc is substring
            if anc in seq:
                valid_seqs.add(seq)
                count += 1
                break

    return count == len(ancs), valid_seqs

#1000 genomes method. 
'''
1000 genomes:
we look at the human-chimp (a) ancestral sequence and compare it to the chimp
(b) and to the human-chimp-orang (c) ancestral sequences. We use the following
convention for the ancestral sequence:

In our case let (a) be the human-chimp-bonobo, (b) human-chimp-bonobo-gorilla, (c) human-chimp-bonobo-gorilla-orangutan

uppercase when all 3 –(a), (b) and (c)– sequences agree
lowercase when:
– there is no ancestral sequence for the ancestral sequence, i.e. there are
only two extant sequences in the alignment, but (a) and (b) agree.
– there is a gap in the sister sequence, but (a) and (c) agree.
– either (b) or (c) disagree with (a), but not both.
N when both (b) and (c) disagree with (a)
- (dash) when no there is no ancestral allele, this is a lineage-specific insertion
. (dot) when there is no alignment, i.e. no data.
'''
def get_anc_allele(col, ancs, valid_seq_to_idx):
    #ancs is in order of a, b, c 
    a, b, c = [col[valid_seq_to_idx[anc]] for anc in ancs]

    #unanimous 
    if a == b == c:
        return a.upper()
    #gap in sister/partial agreement
    elif a == c or a == b:
        return a.lower()
    #both b and c disagree with a
    elif b == c:
        return 'N'
    #lineage specific insertion
    return '-'
    
def remove_last_numbers(input_string):

    return re.sub(r'\d+$', '', input_string)
    
def process_block(block):
    
taf_file = sys.argv[1]
print(taf_file)
taf_index = TafIndex(taf_file + ".tai", is_maf=False)
#print(taf_index)
#presumably length would be max size of the contig and that should work?
target = 'hg38.chr20'
#make this is in order of a, b, c
ancs = ['Anc4.Anc4refChr', 'Anc3.Anc3refChr', 'Anc1.Anc1refChr']
with AlignmentReader(taf_file, taf_index=taf_index, sequence_name=target, start=0, length=int(20e6)) as mp:
    header = mp.get_header()
    sequence_names = set([]) 

    #we can write out the contig as we iterate through the blocks
    
    blocks = []
    #print(len(mp))
    for block in mp:
        #sequence_names.update(block.get_column_sequences()) 
        if block.row_number() > 3:
            for row in block:
                #for each start we can measure the previous gap and fill with '.'
                print(row.start())
                seqs = block.get_column_sequences()
                valid, valid_seqs = all_ancs_in_block(ancs, seqs)
                print(valid_seqs)
                valid_seq_to_idx = {}
                if valid:
                    indices = []
                    for i, seq in enumerate(seqs):
                        if seq in valid_seqs or seq in target:
                            valid_seq_to_idx[remove_last_numbers(seq)] = i

                    col = None
                    anc_alleles = []
                    for i in range(block.column_number()):
                        anc_alleles.append(get_anc_allele(block.get_column(i), ancs, valid_seq_to_idx))
                    
                    for allele in anc_alleles:
                        print(allele, end='')
                    #print(len(col))
                    #print(len(seqs))
                    sys.exit()
                    
                    
                    #print(f'Sequence names: {" ".join([block.g])}')
                    #blocks.append((block, row.start()))



            #print(f'Sequence names: {" ".join(block.get_column_sequences())}')
            #for i in range(block.column_number()):
                #we can iterate through this since thos corresponds to the order of sequences
                #col = block.get_column(i)
                

            #break
    #sort blocks by start 
    print(len(blocks))
    #sorted(blocks, key=lambda x: x[1])

    #we can check each block for containing the right ancestors 
    #for b in blocks[:10]:
    #    print(b[1])

    #print(blocks[0][1])
    #print(blocks[-1][1])
    #block = blocks[0][1]

    #block = blocks[0][0]
    #print(block)
    #for i in range(block.column_number()):
        #we can iterate through this since thos corresponds to the order of sequences
        #col = block.get_column(i)
        #print(col)










