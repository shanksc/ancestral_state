import argparse
import re
import time
import sys
from taffy.lib import AlignmentReader, TafIndex
import pandas as pd
#from Bio import SeqIO

BUFFER=100000

#may want to just use block obj to simplify later
class AncSeq:
    def __init__(self, seq, start, end, strand):
        self.seq = seq 
        self.start = start
        self.end = end 
        self.strand = strand
        #may want to add strand direction

        
def get_args():
    parser = argparse.ArgumentParser(description='check identity between alignments')
    #parser.add_argument('--ref', type=str, required=True, help='reference genome fasta (chm13, hg38)')
    parser.add_argument('--taf', type=str, required=True, help='taf (hg38/CHM13)')
    parser.add_argument('--ref', type=str, required=False, help='hg38/hs1', default='hg38')
    parser.add_argument('--tsv', type=str, required=True, help='tsv containing windows')
    parser.add_argument('--window-size', type=int, required=False, help='size of windows in bp', default=5000)
    parser.add_argument('--ancs', type=str, required=True, nargs='+', help='Anc4.Anc4.refChr, Anc3.Anc3.refChr etc.')
    #parser.add_argument('--chr', type=int, required=True, help='chr corresponding to taf')
    parser.add_argument('--contigs-info', type=str, required=True, help='tsv containing contig, size, chr corresponding to reference (hg38/CHM13)')
    parser.add_argument('--o', type=str, required=True, help='fasta to output for contig')

    return parser.parse_args()

            
#we need a simple heuristic for selecting anc sequences
def get_ancs_in_block(ancs, block, ref_row):
    ancestors_found = {}

    for row in block:
        for anc in ancs:
            #find all ancestor sequences of interest
            #need to understand how cactus comes up with different ancestral contigs
            if anc in row.sequence_name():
                #valid_seqs.append(row.sequence_name())
                if anc in ancestors_found:
                    ancestors_found[anc].append(row.sequence_name)
                else:
                    ancestors_found[anc] = [row.sequence_name]
                break

    return ancestors_found, len(ancestors_found) == len(ancs)

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
#return 1 if unanimous agreement 0 else
def get_anc_allele(col, ancs, anc_to_idx, row):

    #print(col)
    #print(ancs)
    #print(anc_to_idx)
    #look for any matching allele for each anc 
    a, b, c = ancs

    #print(a, b, c)
    a_alleles = [col[i].upper() for i in anc_to_idx[a]]
    b_alleles = [col[i].upper() for i in anc_to_idx[b]]
    c_alleles = [col[i].upper() for i in anc_to_idx[c]]


    #check for unanimous agreement
    for base in a_alleles:
        if base in b_alleles and base in c_alleles:
            #return base.upper()
            return 1
    return 0
    
#since we need an hg38 based fasta, we need only original aligned bases
def build_seq(anc_alleles, row):
    
    #get indices or original bases
    indices = []
    bases = row.bases()
    for i, base in enumerate(bases):
        if base != '-':
            indices.append(i)
    
    seq = ''.join([anc_alleles[idx] for idx in indices])
    assert len(seq) == row.length()

    return seq


#remove numbers from ancestor contig 
def remove_last_numbers(input_string):

    return re.sub(r'\d+$', '', input_string)


#iterate through blocks, getting the inferred ancestral sequence and the interval of the sequence
#we may want to refactor this to just use the block obj later.
def get_identity(taf_file, ancs, ref, windows, window_length, sizes):

    taf_index = TafIndex(taf_file + '.tai', is_maf=False)

    #windows here
    for n in ['19']:
        target = f'{ref}.chr{n}'
        #for window_start, _ in windows[n]:
            #print(target, window_start, window_length, sizes[n])
        with AlignmentReader(taf_file, taf_index, sequence_name=target, start=10000, length=sizes[n]) as mp:
            window_start, _ = windows[n].pop(0)
            equ = 0 #identical bases
            tot = 0 #total bases
                
            for block in mp:
                #print(block)
                col_names = block.get_column_sequences()
                row = block.first_row()
                #get next window
                if row.start() >= window_start + window_length:
                    identity = round(equ / tot, 5)
                    print(f'{n}\t{window_start}\t{identity}', flush=True)
                    equ = 0
                    tot = 0
                    if len(windows[n]) == 0:
                        break
                    window_start, _ = windows[n].pop(0)
                    #print(f'next window: {window_start}')
                    #print(row.start(), flush=True)

                #do we want percent identity between ancestors or between our primates?
                ancestor_names, valid = get_ancs_in_block(ancs, block, row)

                #print(row.length())
                if valid:  
                    #get identity between ancs
                    #print(ancestor_names)
                    anc_to_indices = {anc:[] for anc in ancs}
                    for i, name in enumerate(col_names):
                        prefix = remove_last_numbers(name)
                        if prefix in ancs:
                            anc_to_indices[prefix].append(i)
                    #print(anc_to_indices)

                    alleles = []
                    for i in range(block.column_number()):
                        #print(block.get_column(i))
                        #try:
                        col = block.get_column(i)
                        #except:
                            #continue
                        #print(col)
                        #this should always be true 
                        #if len(col) != row.length():
                            #print(col)
                            #print(row.length())
                            #sys.exit()
                        #only add it if it's in our original window
                        if row.start() + i >= window_start and row.start() + i < window_start + window_length:
                            #print(row.start())
                            #sys.exit()
                            alleles.append(get_anc_allele(block.get_column(i), ancs, anc_to_indices, row))
                        
                    equ += sum(alleles)
                    tot += len(alleles)
                    #print(equ/tot)
                    #sys.exit()
                    
            
            #identity = round(equ / tot, 5)
            #print(f'{n}\t{window_start}\t{identity}')
                    #maybe we treat this as identity == 0
                    #print(ancestor_names)

                    #raise ValueError("Block doesn't have complete ancestral coverage")




            
    '''
    with AlignmentReader(taf_file, taf_index=taf_index, sequence_name=target, start=start, length=window_length) as mp:
        print('iterating', flush=True)
        ancestral_seqs = []
        total_blocks = 0
        valid_blocks = 0
        for block in mp:
            
            #if total_blocks % 1000 == 0 and total_blocks > 0:
            #    print(f'blocks read: {total_blocks}')
            #The first row of each alignment block consists of a sequence from the reference genome.
            #The first row is on the forward (+) strand
            #this is in order of actual column sequence
            total_blocks += 1
            col_names = block.get_column_sequences()
            row = block.first_row()
            ancestor_names, valid = get_ancs_in_block(ancs, block, row)
            #valid, best_anc_seqs = get_best_ancs(valid_seqs, block
            #print(block)
            #sys.exit()
            if valid:
                valid_blocks+=1
                #print(valid_seqs, flush=True)
                
                valid_seq_to_idx = {}
                for i, seq in enumerate(seqs):
                    if seq in valid_seqs or seq in target:
                        valid_seq_to_idx[remove_last_numbers(seq)] = i

                #change this so that we look for across multiple ancestral sequences rather than those that have exact length
                anc_alleles = []
                for i in range(block.column_number()):
                    anc_alleles.append(get_anc_allele(block.get_column(i), ancs, valid_seq_to_idx))

                anc_seq = AncSeq(seq=build_seq(anc_alleles, row), start=row.start(), end=row.start()+row.length(), strand=row.strand())
                
                if abs((anc_seq.end - anc_seq.start) - len(anc_seq.seq)) > 1:
                    print(row)
                    print(row.length())
                    print(row.bases())
                    print(len(row.bases()))
                    print(anc_seq.start, anc_seq.end, len(anc_seq.seq), abs((anc_seq.end - anc_seq.start) - len(anc_seq.seq)))
                    print(block)
                    sys.exit()
                
                anc_to_indices = {anc:[] for anc in ancs}
                for i, name in enumerate(col_names):
                    prefix = remove_last_numbers(name)
                    if prefix in ancs:
                        anc_to_indices[prefix].append(i)

                anc_alleles = []
                for i in range(block.column_number()):
                    anc_alleles.append(get_anc_allele(block.get_column(i), ancs, anc_to_indices, row))
                
                anc_seq = AncSeq(seq=build_seq(anc_alleles, row), start=row.start(), end=row.start()+row.length(), strand=row.strand())
                ancestral_seqs.append(anc_seq)

    print(f'valid blocks: {valid_blocks}')
    print(f'total blocks: {total_blocks}')
    return ancestral_seqs
    '''

def get_windows(tsv, window_size):

    windows = {str(n):[] for n in range(1,23)}
    
    with open(tsv, 'r') as f:
        next(f) #skip header
        for line in f:
            chr, pos, _ = line.strip().split('\t')

            #make zero based
            start = int(pos)-1
            windows[chr].append((start, start + window_size))
    
    return windows
 

def get_chr_sizes(contigs_info):

    sizes = {str(n):0 for n in range(1,23)}
    with open(contigs_info, 'r') as f:
        for line in f:
            _, size, chr = line.strip().split('\t')
            if len(chr) == 5:
                chr = chr[-2:]
            else:
                chr = chr[-1]
            if chr in sizes:
                sizes[chr] = int(size)

    return sizes


def main():
    args = get_args()
    
    windows = get_windows(args.tsv, args.window_size)
    sizes = get_chr_sizes(args.contigs_info)
    get_identity(args.taf, args.ancs, args.ref, windows, args.window_size, sizes)


if __name__ == '__main__':
    main()
