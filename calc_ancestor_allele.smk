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

#pass config file
rule all:
    input:
        #f"{config['working']}/{config['ref']}.ready.txt",
        #expand("{config['working']}/{config['ref']}_{autosome}_processed.txt", autosome=range(1, 23))
        #expand("{config['working']}/{config['ref']}_{autosome}.fa", autosome=range(1, 23))
        expand('{prefix}_chr{autosome}.fa', autosome=range(1,23), prefix=f"{config['working']}/{config['ref']}"),
        expand('{prefix}_chr{autosome}_stats.txt', autosome=range(1,23), prefix=f"{config['working']}/{config['ref']}")
        #f"{config['taf']}.filtered"

#computing only autosomes
#create empty stats files for each chr 
checkpoint create_stats:
    input:
        config['taf']
    output:
        f"{config['working']}/{config['ref']}.ready.txt"
    run:
        # Create empty files for autosomes
        #open(f"{config['working']}/{config['ref']}.ready.txt", 'w').close()
        open(output[0], 'w').close()
        for autosome in range(1, 23):
            autosome_file = os.path.join(config['working'], f"{config['ref']}_{autosome}.txt")
            open(autosome_file, 'w').close()
        open(output[0], 'w').close()


rule calc_anc:
    input:  
        "{prefix}_{autosome}.txt",
        config['contig-info']
    output: 
        "{prefix}_chr{autosome}.fa"
    shell:
        """
        python3 scripts/create_ancestor.py --taf {config[taf]} --target {config[ref]}.chr{wildcards.autosome}\
        --ancs {config[ancestors][anc3]} {config[ancestors][anc2]} {config[ancestors][anc1]}\
        --contigs-info {config[contig-info]} --o {output}
        """

rule compute_stats:
    input:
        "{prefix}_{autosome}.fa"
    output:
        "{prefix}_{autosome}_stats.txt"
    shell:
        "bash scripts/get_base_counts.sh {input} > {output}"

#rule 