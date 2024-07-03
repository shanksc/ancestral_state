import sys
import os

#change this in future

#cpu_count = int(os.environ['SLURM_CPUS_PER_TASK'])
#assert cpu_count >= 22
#just ensuring that our cpu count is utilized by the 

print(config['taf'])
print(config['ref'])
print(config['ancestors'])
print(config['contig-info'])
#print(config['compare'])


chroms=list(range(1,23)) + ['X', 'Y']
#pass config file
rule all:
    input:
        #f"{config['working']}/{config['ref']}.ready.txt",
        #expand('{prefix}_{autosome}.txt', autosome=range(1,23), prefix=f"{config['working']}/{config['ref']}"),
        expand('{prefix}_chr{autosome}.fa', autosome=chroms, prefix=f"{config['working']}/{config['ref']}"),
        #add rule to produce the felsing version, expand the config to support this
        #get stats for each of these as well 
	#expand('{prefix}_chr{autosome}_stats.txt', autosome=chroms, prefix=f"{config['working']}/{config['ref']}")

#computing only autosomes
#create empty stats files for each chr 
rule create_stats:
    input:
        config['taf']
    output:
        expand('{prefix}_{autosome}.txt', autosome=range(1,23), prefix=f"{config['working']}/{config['ref']}")
        #f"{config['working']}/{config['ref']}.ready.txt"
    run:
        # Create empty files for autosomes
        #open(f"{config['working']}/{config['ref']}.ready.txt", 'w').close()
        open(output[0], 'w').close()
        for autosome in range(1, 23):
            autosome_file = os.path.join(config['working'], f"{config['ref']}_{autosome}.txt")
            open(autosome_file, 'w').close()
        open(output[0], 'w').close()

#create hg38 based fastas of SNPs 
rule calc_anc_snps:
    input:  
        #"{prefix}_{autosome}.txt",
        config['contig-info']
    output: 
        "{prefix}_chr{autosome}.fa"
    	#"{prefix_chr{autosome}_indels.vcf"
    shell:
        """
        python3 scripts/create_ancestor.py --taf {config[taf]} --target {config[ref]}.chr{wildcards.autosome}\
        --ancs {config[ancestors][anc3]} {config[ancestors][anc2]} {config[ancestors][anc1]}\
        --contigs-info {config[contig-info]} --o {wildcards.prefix}_chr{wildcards.autosome}
        """

rule compute_stats:
    input:
        "{prefix}_{autosome}.fa"
    output:
        "{prefix}_{autosome}_stats.txt"
    shell:
        "bash scripts/get_base_counts.sh {input} > {output}"

