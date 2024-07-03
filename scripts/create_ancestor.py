import argparse
import re
import time
import sys
from taffy.lib import AlignmentReader, TafIndex


#ancseq of snps for reference
class AncSeq:
    def __init__(self, seq, start, end, strand):
        self.seq = seq 
        self.start = start
        self.end = end 
        self.strand = strand
        #may want to add strand direction

#variants class for writing out vcf 
class Variant:
    def __init__(self, chrom):
        self.chrom = chrom
        self.ref = ''
        self.alt = '' 
        self.start = None

        
def get_args():
    parser = argparse.ArgumentParser(description='create fasta containing ancestral alleles')
    #parser.add_argument('--ref', type=str, required=True, help='reference genome fasta (chm13, hg38)')
    parser.add_argument('--taf', type=str, required=True, help='taf (hg38/CHM13)')
    parser.add_argument('--target', type=str, required=True, help='ex. hg38.chr20')
    parser.add_argument('--ancs', type=str, required=True, nargs='+', help='Anc4.Anc4.refChr, Anc3.Anc3.refChr etc.')
    #parser.add_argument('--chr', type=int, required=True, help='chr corresponding to taf')
    parser.add_argument('--contigs-info', type=str, required=True, help='tsv containing contig, size, chr corresponding to reference (hg38/CHM13)')
    parser.add_argument('--o', type=str, required=True, help='indel vcf, fasta to output for contig')

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

In our case let (a) be the human-chimp-bonobo, (b) chimp-bonobo (c) human-chimp-bonobo-gorilla

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
def get_consensus(a, b, c):
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

#this breaks if MAF/TAF isn't properly filtered to contain single copies of contigs 
def get_anc_allele(col, ancs, anc_to_indices):
    #ancs is in order of a, b, c 
    a, b, c = [col[anc_to_indices[anc]] for anc in ancs]

    return get_consensus(a, b, c)


#do we need to realign seqs for context like Ensembl?
#since we need an hg38 based fasta, we need only original aligned bases
def build_seq(anc_alleles, row):
    
    #get indices or original bases
    indices = []
    bases = row.bases()
    start = row.start()

    #these are original hg38/chm13 bases
    for i, base in enumerate(bases):
        if base != '-':
            indices.append(i)
    
    seq = ''.join([anc_alleles[idx] for idx in indices])
    #row.length() is the number of bases in the row of our reference
    assert len(seq) == row.length()
    return seq


#remove numbers from ancestor contig 
def remove_last_numbers(input_string):

    return re.sub(r'\d+$', '', input_string)


def get_alleles(chrom, row, block, ancs, anc_to_indices):
    #ref will always be first in col
    curr_var = Variant(chrom=chrom)
    ref_bases = row.bases()
    indels = []
    anc_alleles = []
    
    #init ref_position to first valid_base in block
    for idx, base in enumerate(ref_bases):
        if base != '-':
            ref_position = row.start() + idx + 1
            break

    assert len(ref_bases) == block.column_number()
    for i in range(0, block.column_number()):
        col = block.get_column(i)
        a, b, c = [col[anc_to_indices[anc]] for anc in ancs]
        


        anc_allele = get_consensus(a, b, c)
        anc_alleles.append(anc_allele)

        #now classify as indel or snp relative to ref 

        #anc allele is insertion relative to ref - should we go beyond 1bp (why not)?
        #if ref_bases[i] == '-' and anc_allele != 'N' and anc_allele != '-':
        #filter insertions and deletions later
        if ref_bases[i] == '-':
            #build insertion
            #print('insert')
            #print(ref_bases[i-3:i+1])
            #print(anc_allele)
            if curr_var.ref == '':
                curr_var.start = ref_position - 1#check that is based correctly
                curr_var.ref = ref_bases[i-1]
                curr_var.alt += ref_bases[i-1] + anc_allele
            #longer than 1bp 
            else:
                curr_var.alt += anc_allele

        #anc allele is deletion relative to ref
        elif ref_bases[i] != '-' and anc_allele == '-':
            #print(ref_bases[i-3:i+1])
            if curr_var.ref == '':
                #make 1-based
                curr_var.start = row.start() + i - 1
                curr_var.ref = ref_bases[i-1] + ref_bases[i]
                curr_var.alt = ref_bases[i-1]
                #longer than 1bp
                #print(ref_bases[i-3:i+3])
                #print(curr_var.start, curr_var.ref, curr_var.alt)
                #sys.exit()
            else:
                curr_var.ref += ref_bases[i]

        #indel ends if there aren't any more gaps
        elif curr_var.ref != '':
            indels.append(curr_var)
            #print(curr_var.ref, curr_var.alt, curr_var.start)
            #sys.exit()
            curr_var = Variant(chrom=chrom)
        
        if ref_bases[i] != '-':
            ref_position += 1

    return indels, anc_alleles
            

#iterate through blocks, getting the inferred ancestral sequence and the interval of the sequence
#we may want to refactor this to just use the block obj later.
def get_ancestral_seqs(taf_file, target, ancs, size):

    taf_index = TafIndex(taf_file + '.tai', is_maf=False)

    start = 0
    #chrom = target.split('.')[1]
    chrom = target.split('.')[1]
    #chrom = re.search(r'\d+$', target.split('.')[1]).group(0)
    print(f'chr: {chrom}')

    #this is a temporary fix to deal with the taffy view error where the first chromosome isn't found in the index
    #this is fixed in newest taffy update 
    if target == 'hg38.chr1':
        start = 10000
    
    with AlignmentReader(taf_file, taf_index=taf_index, sequence_name=target, start=start, length=size) as mp:
        #print('iterating', flush=True)
        ancestral_seqs = []
        total_blocks = 0
        valid_blocks = 0

        indels = []

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


                new_indels, anc_alleles = get_alleles(chrom, row, block, ancs, anc_to_indices)
                if new_indels:
                    indels.extend(new_indels)
                
                anc_seq = AncSeq(seq=build_seq(anc_alleles, row), start=row.start(), end=row.start()+row.length(), strand=row.strand())
                #we can change this to yield 
                ancestral_seqs.append(anc_seq)

    print(f'valid blocks: {valid_blocks}')
    print(f'total blocks: {total_blocks}')
    return indels, ancestral_seqs

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

#write out indels in vcf format
def write_vcf(indels, out):

    #write header
    with open(out+'.vcf', 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        
        #id will always be blank until we get rs numbers from dbSNP
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for indel in indels:
            if '-' not in indel.alt and 'N' not in indel.alt:
                #check for 1bp indels
                if (len(indel.ref) == 2 and len(indel.alt) == 1) or (len(indel.ref) == 1 and len(indel.alt) == 2):
                    f.write(f'{indel.chrom}\t{indel.start}\t.\t{indel.ref}\t{indel.alt}\n')


def write_fasta(anc_seqs, target, size, out, line_length=80):
    
    total_chars = 0
    with open(out+'.fa', 'w') as f:
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
    indels, anc_seqs = get_ancestral_seqs(args.taf, args.target, args.ancs, size)
    end_time = time.time()
    
    print(f'time: {end_time - start_time}')
    
    print(indels[0].ref, indels[0].alt, indels[0].start)
    
    #check that it's sorted
    assert(sorted(anc_seqs, key = lambda x:x.start) == anc_seqs)
    #check for only positive strands in reference
    assert all(s.strand for s in anc_seqs)

    print(f'max end coord in taf {max(s.end for s in anc_seqs)}')
    

    #should be no difference from our other fasta 
    total_chars = write_fasta(anc_seqs, args.target, size, args.o)
    print(f'{total_chars} sequence characters written to fasta')

    #write out vcf with indels
    write_vcf(indels, args.o)

if __name__ == '__main__':
    main()
