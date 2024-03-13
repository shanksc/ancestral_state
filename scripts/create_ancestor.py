import argparse
import re
import time
import sys
from taffy.lib import AlignmentReader, TafIndex
#from Bio import SeqIO
#may want to just use block obj to simplify later
class AncSeq:
    def __init__(self, seq, start, end, strand):
        self.seq = seq 
        self.start = start
        self.end = end 
        self.strand = strand
        #may want to add strand direction

        
def get_args():
    parser = argparse.ArgumentParser(description='create fasta containing ancestral states')
    #parser.add_argument('--ref', type=str, required=True, help='reference genome fasta (chm13, hg38)')
    parser.add_argument('--taf', type=str, required=True, help='taf (hg38/CHM13)')
    parser.add_argument('--target', type=str, required=True, help='ex. hg38.chr20')
    parser.add_argument('--ancs', type=str, required=True, nargs='+', help='Anc4.Anc4.refChr, Anc3.Anc3.refChr etc.')
    #parser.add_argument('--chr', type=int, required=True, help='chr corresponding to taf')
    parser.add_argument('--contigs-info', type=str, required=True, help='tsv containing contig, size, chr corresponding to reference (hg38/CHM13)')
    parser.add_argument('--o', type=str, required=True, help='fasta to output for contig')

    return parser.parse_args()


def get_chr_size(target, contigs_info):

    target_chr = target.strip().split('.')[1]
    with open(contigs_info, 'r') as f:
        for line in f:
            _, size, chr = line.strip().split('\t')
            if target_chr in chr:
                print(chr, size)
                return int(size)

    return None

#for single-copy          
'''
def all_ancs_in_block(ancs, seqs):
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
'''
#we need a simple heuristic for selecting anc sequences
#for now we should be able to get 
def get_ancs_in_block(ancs, block, ref_row):

    #for now lets just select blocks that are the same size as the reference and on the postitive strand
    #this should exist for most blocks-check this
    valid_seqs = []
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
#don't really need ancestor names at all, just need to check that there's enought to be valid 
def get_anc_allele(col, ancs, anc_to_idx, row):

    #look for any matching allele for each anc 
    a, b, c = ancs

    a_alleles = [col[i].upper() for i in anc_to_idx[a]]
    b_alleles = [col[i].upper() for i in anc_to_idx[b]]
    c_alleles = [col[i].upper() for i in anc_to_idx[c]]

    canidates = []

    #may need to modify this to handle tie-breaking - we can use the reference/human anc 

    #check for unanimous agreement
    for base in a_alleles:
        if base in b_alleles and base in c_alleles:
            return base.upper()
            #canidates.append(base)
    
    #check for partial agreement with closest ancestor
    for base in a_alleles:
        if base in b_alleles:
            return base.lower()
    
    #check for partial agreement with farther ancestor
    for base in a_alleles:
        if base in c_alleles:
            return base.lower()
    
    #other ancestors agree but not with human-chimp-bonobo
    for base in b_alleles:
        if base in c_alleles:
            return 'N'
    
    #lineage specific
    return '-'

        
'''
#with this we have to select which anc seqs we want
def get_anc_allele(col, ancs, anc_to_indices):
    #ancs is in order of a, b, c 
    a, b, c = [col[anc_to_indices[anc]] for anc in ancs]

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
def get_ancestral_seqs(taf_file, target, ancs, size):

    taf_index = TafIndex(taf_file + '.tai', is_maf=False)

    #this is a temporary fix to deal with the taffy view error where the first chromosome isn't found in the index
    #change this to potentially read from the .tai actually  
    #this won't be needed when pull request goes through
    start=0
    if target == 'hg38.chr1':
        start = 10000
    
    with AlignmentReader(taf_file, taf_index=taf_index, sequence_name=target, start=start, length=size) as mp:
        print('iterating', flush=True)
        ancestral_seqs = []
        total_blocks = 0
        valid_blocks = 0
        for block in mp:
            
            #The first row of each alignment block consists of a sequence from the reference genome.
            #The first row is on the forward (+) strand
            #this is in order of actual column sequence
            total_blocks += 1
            col_names = block.get_column_sequences()
            row = block.first_row()
            ancestor_names, valid = get_ancs_in_block(ancs, block, row)

            if valid:
                valid_blocks+=1

                anc_to_indices = {anc:[] for anc in ancs}
                for i, name in enumerate(col_names):
                    prefix = remove_last_numbers(name)
                    if prefix in ancs:
                        anc_to_indices[prefix] = i

                #print(col_names)
                #print(anc_to_indices)
                anc_alleles = []
                for i in range(block.column_number()):
                    anc_alleles.append(get_anc_allele(block.get_column(i), ancs, anc_to_indices))
                    #anc_alleles.append(get_anc_allele(block.get_column(i), ancs, anc_to_indices, row))
                #print(anc_alleles)
                anc_seq = AncSeq(seq=build_seq(anc_alleles, row), start=row.start(), end=row.start()+row.length(), strand=row.strand())
                ancestral_seqs.append(anc_seq)

    print(f'valid blocks: {valid_blocks}')
    print(f'total blocks: {total_blocks}')
    return ancestral_seqs

#ancestral seqs is ordered
#return a generator 
def get_fasta(anc_seqs, size, line_length=80):
    #zero based
    total_chars = 0
    gap_chars_written = 0
    last_end = 0
    for seq in anc_seqs:
        gap_str = '.' * (seq.start - last_end)
        assert (seq.start - last_end) >= 0

        last_end = seq.end

        total_chars += len(gap_str)
        total_chars += len(seq.seq)
        gap_chars_written += len(gap_str)

        if len(gap_str) > 0:
            for i in range(0, len(gap_str), line_length):
                yield gap_str[i:i + line_length]
        for i in range(0, len(seq.seq), line_length):
            yield seq.seq[i:i + line_length]

    print(f'gap chars: {gap_chars_written}')
    print(f'total chars: {total_chars}')
    #print(f'last-last-end: {last_end}')
    print(f'size {size}')
    #fill end of contig
    if last_end < size:
        gap_str = '.' * (size - last_end)
        total_chars += len(gap_str)
        for i in range(0, len(gap_str), line_length):
            yield gap_str[i:i + line_length]
    
    print(f'total chars after blocks {total_chars}')


def write_fasta(anc_seqs, target, size, out, line_length=80):
    
    total_chars = 0
    with open(out, 'w') as f:
        f.write(f'>{target}\n')
        buffer = ''

        for base in get_fasta(anc_seqs, size):
            buffer += base
            total_chars += len(base)

            if len(buffer) >= line_length:
                f.write(buffer[:line_length] + '\n')
                buffer = buffer[line_length:]

        if buffer:
            f.write(buffer + '\n')

    return total_chars


def main():
    args = get_args()
    #iterate through each block in sorted order- fill gaps with . to represent no alignment data
    size = get_chr_size(args.target, args.contigs_info)
    if size is None:
        print('error could not get size of target chr')
        sys.exit()

    print(f'{args.target} size: {size}')
    
    start_time = time.time()
    anc_seqs = get_ancestral_seqs(args.taf, args.target, args.ancs, size)
    end_time = time.time()
    
    print(f'time: {end_time - start_time}')

    #check that it's sorted
    assert(sorted(anc_seqs, key = lambda x:x.start) == anc_seqs)
    #check for only positive strands in reference
    assert all(s.strand for s in anc_seqs)

    print(f'max end coord in taf {max(s.end for s in anc_seqs)}')

    total_chars = write_fasta(anc_seqs, args.target, size, args.o)
    print(f'{total_chars} sequence characters written to fasta')

if __name__ == '__main__':
    main()
