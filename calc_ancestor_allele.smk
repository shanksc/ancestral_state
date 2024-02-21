import sys
import os

cpu_count = int(os.environ['SLURM_CPUS_PER_TASK'])
assert cpu_count >= 22

print(config['taf'])
print(config['ref'])

#pass config file
rule all:
    input:
        f"{config['working']}/{config['ref']}.done.txt"
        #f"{config['taf']}.filtered"

#just do autosomes
#we are always splitting by the reference (hg38 or chm13).
rule get_contigs:
    input:
        f"{config['taf']}.filtered",
        config['working']
    output:
        f"{config['working']}/{config['ref']}.done.txt"
    shell:
        #split by contigs with taffy
        "bash scripts/split_by_contig.sh {config[ref]} {config[taf]} {config[working]}"

#def get_contigs()
#for each base in GRCh38 or CHM13, compute the ancestral state for each contig
#rule calc_anc:
#    input:
#        f"{config['working']}/{config['ref']}.done.txt"
        #check that all autosomes are ready, otherwis throw error 
